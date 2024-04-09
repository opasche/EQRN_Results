# This implementation is based on the approach in "Extreme Quantile Estimation for Autoregressive Models"
# by Deyuan Li and Huixia Judy Wang, as well as on some direct implementation insights from the authors.
# It generalizes it to multi-process covariates, and wraps the EQRN-style UI.
# Please ask the authors rather than using this implementation if you want to apply their method.
# Olivier PASCHE, 2023

library(quantreg)


EXQAR_fit <- function(y, X=NULL, lag, delta1="default", delta2="default") {
  
  if(is.null(X)){
    X_feat <- lagged_features(X=matrix(y, ncol=1), max_lag=lag, drop_present=TRUE)
  }else{
    X_feat <- lagged_features(X=cbind(y,X), max_lag=lag, drop_present=TRUE)
  }
  n <- nrow(X_feat)
  y_resp <- y[lag+(1:n)]
  stopifnot("Dimension error."=(n==(length(y)-lag)))
  
  if(delta1=="default"){delta1 <- n^(-2/3)*4.5}
  if(delta2=="default"){delta2 <- n^(-0.9)}
  
  k1 <- as.integer(n*delta1)
  k2 <- as.integer(n*delta2)
  tau_grid <- seq(1-delta1, 1-delta2, length=k1-k2)
  
  # Estimate the cond. quantile coefficients of y given ypast and xpast at levels tau_grid
  beta_coeffs <- quantreg::rq(y_resp ~ X_feat, tau_grid)$coef
  
  EXQAR_obj <- list(beta_coeffs=beta_coeffs, lag=lag, delta1=delta1, delta2=delta2, k1=k1, k2=k2, tau_grid=tau_grid)
  class(EXQAR_obj) <- c("EXQAR")
  return(EXQAR_obj)
}


EXQAR_predict <- function(EXQAR_obj, y, X=NULL, prob_lvls_predict, tol=1e-4, min_prop=0.3, return_infos=FALSE){
  
  if(!is.null(y)){
    if(is.null(X)){
      X_feat <- lagged_features(X=matrix(y, ncol=1), max_lag=EXQAR_obj$lag, drop_present=TRUE)
    }else{
      X_feat <- lagged_features(X=cbind(y,X), max_lag=EXQAR_obj$lag, drop_present=TRUE)
    }
  }else{
    if(is.null(X)){
      stop("Must provide at least one of 'y' or 'X' in 'EXQAR_predict'.\n")
    }else{
      X_feat <- matrix(X)
    }
  }
  n_test <- nrow(X_feat)
  nb_prob_lvls_predict <- length(prob_lvls_predict)
  
  # Estimate the cond. quantiles of y given ypast_test and xpast_test at levels EXQAR_obj$tau_grid
  Q <- cbind(1, X_feat) %*% as.matrix(EXQAR_obj$beta_coeffs)
  
  # Obtain the time-dependent shape (EVI) and extreme quantile estimates
  Q_ex <- matrix(as.double(NA), nrow=n_test, ncol=nb_prob_lvls_predict)
  shape <- rep(as.double(NA), n_test)
  am <- rep(as.double(NA), n_test)
  Ut <- rep(as.double(NA), n_test)
  for(j in 1:n_test){
    tail_est <- tail_moment_estim(q=Q[j,], interm_lvls=EXQAR_obj$tau_grid, prob_lvls_predict=prob_lvls_predict, tol=tol, min_prop=min_prop)
    Q_ex[j,] <- tail_est$Q_ex
    shape[j] <- tail_est$shape
    am[j] <- tail_est$am
    Ut[j] <- tail_est$Ut
  }
  
  if(return_infos){
    out <- list(prediction=Q_ex, shape=shape, am=am, Ut=Ut)
    return(out)
  }else{
    return(Q_ex)
  }
}


EXQAR_excess_probability <- function(EXQAR_obj, val, y, X=NULL, body_proba="default", proba_type=c("excess","cdf"),
                                     tol=1e-4, min_prop=0.3, return_infos=FALSE){
  proba_type <- match.arg(proba_type)
  preds <- EXQAR_predict(EXQAR_obj, y=y, X=X, prob_lvls_predict=NULL, tol=tol, min_prop=min_prop, return_infos=TRUE)
  
  Probs <- GPD_excess_probability(val, sigma=preds$am, xi=preds$shape, interm_threshold=preds$Ut,
                                  threshold_p=1-EXQAR_obj$delta1, body_proba=body_proba, proba_type=proba_type)
  if(return_infos){
    return(list(Probs=c(Probs), pred_object=preds))
  }else{
    return(c(Probs))
  }
}


predict.EXQAR <- function(EXQAR_obj, ...){
  # The 'predict' method for class "EXQAR".
  # See 'EXQAR_predict' for details.
  return(EXQAR_predict(EXQAR_obj, ...))
}

excess_probability.EXQAR <- function(EXQAR_obj, ...){
  # The 'excess_probability' prediction method for class "EXQAR".
  # See 'EXQAR_excess_probability' for details.
  return(EXQAR_excess_probability(EXQAR_obj, ...))
}


rearrange_list <- function(f, xmin, xmax){
  # Modification of `quantreg::rearrange` to return a list instead of a `stepfun`.
  # Also fixing the TRUE() | T() argument check
  if (is.list(f)) 
    lapply(f, rearrange_list)
  else {
    if (!is.stepfun(f)) 
      stop("Only stepfuns can be rearranged.\n")
    call <- attributes(f)$call
    right <- call[match("right", names(call))] == "TRUE()" | 
              call[match("right", names(call))] == "T()"
    x <- knots(f)
    n <- length(x)
    if (missing(xmin)) 
      xmin <- x[1]
    if (missing(xmax)) 
      xmax <- x[n]
    x <- x[(x >= xmin) & (x <= xmax)]
    x <- c(xmin, x, xmax)
    n <- length(x)
    y <- f(x)
    o <- ifelse(rep(right, n - 1), order(y[-1]) + 1, order(y[-n]))
    x <- cumsum(c(x[1], diff(x)[o - right]))
    y <- y[o]
    y <- c(y[1], y, max(y))
    return(list(x=x, y=y))
  }
}


tail_moment_estim <- function(q, interm_lvls, prob_lvls_predict=NULL, tol=1e-4, min_prop=0.3){
  diff_prop <- mean(diff(q) > tol, na.rm=T)
  q_posi <- which(!is.na(q) & q > 0)
  if(!is.na(diff_prop) & (diff_prop > min_prop) & (length(q_posi) > 2)){
    q_pos <- q[q_posi]
    interm_lvls2 <- interm_lvls[q_posi]
    q_mono <- rearrange_list(stepfun(interm_lvls2[-length(interm_lvls2)], q_pos), xmax=1)$y
    Ut <- q_mono[1]
    Mn1 <- mean((log(q_mono/Ut)))
    Mn2 <- mean((log(q_mono/Ut))^2)
    Mn_ratio <- 0.5/(1-Mn1^2/Mn2)
    shape <- Mn1 + 1 - Mn_ratio
    shape <- pmin(pmax(shape, -1),1) # bound shape estimate for stability
    am <- Ut * Mn1 * Mn_ratio
    if(!is.null(prob_lvls_predict)){
      Q_ex <- Ut + am*(((1-interm_lvls[1])/(1-prob_lvls_predict))^shape-1)/shape
    }else{
      Q_ex <- NULL
    }
  }else{
    shape <- Q_ex <- am <- Ut <- NA
  } 
  return(list(shape=shape, Q_ex=Q_ex, am=am, Ut=Ut))
}



QAR_fit <- function(y, X=NULL, lag){
  
  if(is.null(X)){
    X_feat <- lagged_features(X=matrix(y, ncol=1), max_lag=lag, drop_present=TRUE)
  }else{
    X_feat <- lagged_features(X=cbind(y,X), max_lag=lag, drop_present=TRUE)
  }
  n <- nrow(X_feat)
  y_resp <- y[lag+(1:n)]
  stopifnot("Dimension error."=(n==(length(y)-lag)))
  
  #QAR estimation
  coef <- quantreg::rq(y_resp ~ X_feat, prob_lvls_predict)$coef
  
  QAR_obj <- list(coef=coef, lag=lag)
  class(QAR_obj) <- c("QAR")
  return(QAR_obj)
}


QAR_predict <- function(QAR_obj, y, X=NULL, prob_lvls_predict, return_infos=FALSE) {
  
  if(!is.null(y)){
    if(is.null(X)){
      X_feat <- lagged_features(X=matrix(y, ncol=1), max_lag=QAR_obj$lag, drop_present=TRUE)
    }else{
      X_feat <- lagged_features(X=cbind(y,X), max_lag=QAR_obj$lag, drop_present=TRUE)
    }
  }else{
    if(is.null(X)){
      stop("Must provide at least one of 'y' or 'X' in 'EXQAR_predict'.\n")
    }else{
      X_feat <- matrix(X)
    }
  }
  #n_test <- nrow(X_feat)
  #nb_prob_lvls_predict <- length(prob_lvls_predict)
  
  QAR_pred <- cbind(1, X_feat) %*% as.matrix(QAR_obj$coef)
  
  return(QAR_pred)
}


predict.QAR <- function(QAR_obj, ...){
  # The 'predict' method for class "QAR".
  # See 'QAR_predict' for details.
  return(QAR_predict(QAR_obj, ...))
}

