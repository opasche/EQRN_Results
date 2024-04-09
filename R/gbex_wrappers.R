# Wrappers around the gbex EQR model for EQRN-like UI.
# They extend the gbex approach to use the conditional intermediate quantiles
# as an accuracy-improving covariate.
# Olivier PASCHE, 2022

#devtools::install_github("JVelthoen/gbex")
library(tidyverse)
library(gbex)
library(future)
library(doFuture)

gbex_fit <- function(X, y, intermediate_quantiles, interm_lvl, intermediate_q_feature=FALSE, scale_features=FALSE, ...){
  data_excesses <- get_excesses(X=X, y=y, quantiles=intermediate_quantiles,
                                intermediate_q_feature=intermediate_q_feature,
                                scale_features=scale_features, X_scaling=NULL)
  Y_excesses <- data_excesses$Y_excesses
  X_feats_excesses <- data_excesses$X_excesses
  X_scaling <- data_excesses$X_scaling
  fit_gbex <- gbex::gbex(y=Y_excesses, X=data.frame(X_feats_excesses), ...)
  
  gbex_obj <- list(fit_gbex=fit_gbex, intermediate_q_feature=intermediate_q_feature, 
                   interm_lvl=interm_lvl, scale_features=scale_features, X_scaling=X_scaling)
  class(gbex_obj) <- "gbex_eqr"
  return(gbex_obj)
}


gbex_predict <- function(fitted_gbex, X_test, to_predict, intermediate_quantiles, interm_lvl=fitted_gbex$interm_lvl){
  
  if(interm_lvl!=fitted_gbex$interm_lvl){stop("gbex intermediate quantiles interm_lvl does not match in train and predict.")}
  
  X_feats <- X_test
  if(fitted_gbex$intermediate_q_feature){
    if(is.null(intermediate_quantiles)){stop("intermediate_quantiles needed in gbex_predict for intermediate_q_feature.")}
    X_feats <- cbind(X_feats, intermediate_quantiles)
  }
  X_feats <- perform_scaling(X=X_feats, X_scaling=fitted_gbex$X_scaling, scale_features=fitted_gbex$scale_features)$X_scaled
  
  pred_params <- predict(fitted_gbex$fit_gbex, newdata=data.frame(X_feats), what="par")
  
  if(any(to_predict=="par")){
    return(pred_params)
  }
  
  if(length(dim(to_predict))>1){
    stop("Please provide a single value or 1D vector as to_predict in gbex_predict")
  }
  
  if(length(to_predict)==1){
    return(GPD_quantiles(to_predict, interm_lvl, intermediate_quantiles, pred_params$s, pred_params$g))
  } else if(length(to_predict)>1){
    nb_prob_lvls_predict <- length(to_predict)
    gbex_quantiles <- matrix(as.double(NA), nrow=nrow(X_test), ncol=nb_prob_lvls_predict)
    for(i in 1:nb_prob_lvls_predict){
      gbex_quantiles[,i] <- GPD_quantiles(to_predict[i], interm_lvl, intermediate_quantiles, pred_params$s, pred_params$g)
    }
    return(gbex_quantiles)
  } else {
    stop("Please provide a single value or 1D vector as to_predict in gbex_predict")
  }
}

predict.gbex_eqr <- function(fitted_gbex, ...){
  # The 'predict' method for class "gbex_eqr".
  # See 'gbex_predict' for details.
  return(gbex_predict(fitted_gbex, ...))
}


gbex_CV <- function(X, y, intermediate_quantiles, interm_lvl, intermediate_q_feature=FALSE, scale_features=FALSE,
                    num_folds=5, Bmax=500, grid_lambda_ratio=c(5,6,7,8,9,10), grid_depth=list(c(1,0),c(1,1),c(2,1),c(2,2),c(3,1),c(3,2),c(3,3)),
                    stratified=T, lambda_scale=0.01, min_leaf_size=NULL, sf=0.75,
                    parallel_strat=c("sequential", "multisession", "multicore"), n_workers=(availableCores() - 1),
                    seed=NULL, err_handling=c("stop", "remove", "pass"), ...){
  parallel_strat <- match.arg(parallel_strat)
  err_handling <- match.arg(err_handling)
  
  data_excesses <- get_excesses(X=X, y=y, quantiles=intermediate_quantiles,
                                intermediate_q_feature=intermediate_q_feature,
                                scale_features=scale_features, X_scaling=NULL)
  Y_excesses <- data_excesses$Y_excesses
  X_feats_excesses <- data_excesses$X_excesses
  X_scaling <- data_excesses$X_scaling
  #map cross params
  if(is.null(min_leaf_size)){min_leaf_size <- rep(max(10, length(Y_excesses)/100), 2)}
  
  n_workers <- min(n_workers,length(grid_depth))
  `%fun%` <- set_doFuture_strategy(parallel_strat, n_workers=n_workers)
  
  results_grid_fit <- foreach(depth=grid_depth, .errorhandling=err_handling, .combine=rbind) %fun% {
    if(!is.null(seed)){set.seed(seed)}
    CV_lambda_ratio = gbex::CV_gbex(y=Y_excesses, X=data.frame(X_feats_excesses), num_folds=num_folds, Bmax=Bmax, par_name="lambda_ratio", par_grid=grid_lambda_ratio,
                                    stratified=T, lambda_scale=lambda_scale, depth=depth, min_leaf_size=min_leaf_size, sf=sf, silent=T,ncores=1, ...)
    
    lrr_opt <- CV_lambda_ratio$par_CV
    B_opt <- CV_lambda_ratio$B_opt
    dev_best <- CV_lambda_ratio$dev_all[B_opt, which(lrr_opt==CV_lambda_ratio$par_grid)]
    
    output <- c(depth_scale=depth[1], depth_shape=depth[2], lrr_opt=lrr_opt, B_opt=B_opt, dev_best=dev_best)
    output
  }
  results_tibble <- tibble::as_tibble(results_grid_fit)
  
  end_doFuture_strategy()
  
  best_params <- results_tibble %>% slice_min(dev_best) %>%
    slice_min(depth_shape) %>% slice_min(depth_scale) %>% slice_max(lrr_opt) %>% as.list()
  
  best_params_renamed <- list(depth=c(best_params$depth_scale, best_params$depth_shape), lambda_ratio=best_params$lrr_opt, B=best_params$B_opt)
  
  return(list(best_params=best_params_renamed, results_tibble=results_tibble,
              intermediate_q_feature=intermediate_q_feature, interm_lvl=interm_lvl, scale_features=scale_features, X_scaling=X_scaling))
}


gbex_excess_probability <- function(fitted_gbex, val, X_test, intermediate_quantiles, interm_lvl=fitted_gbex$interm_lvl,
                                    body_proba="default", proba_type=c("excess","cdf")){
  
  proba_type <- match.arg(proba_type)
  if(interm_lvl!=fitted_gbex$interm_lvl){stop("gbex intermediate quantiles interm_lvl does not match in train and predict.")}
  
  X_feats <- X_test
  if(fitted_gbex$intermediate_q_feature){
    if(is.null(intermediate_quantiles)){stop("intermediate_quantiles needed in gbex_predict for intermediate_q_feature.")}
    X_feats <- cbind(X_feats, intermediate_quantiles)
  }
  X_feats <- perform_scaling(X=X_feats, X_scaling=fitted_gbex$X_scaling, scale_features=fitted_gbex$scale_features)$X_scaled
  
  pred_params <- predict(fitted_gbex$fit_gbex, newdata=data.frame(X_feats), what="par")
  
  Probs <- GPD_excess_probability(val, sigma=pred_params$s, xi=pred_params$g, interm_threshold=intermediate_quantiles,
                                  threshold_p=interm_lvl, body_proba=body_proba, proba_type=proba_type)
  return(c(Probs))
}

excess_probability.gbex_eqr <- function(fitted_gbex, ...){
  # The 'excess_probability' prediction method for class "gbex_eqr".
  # See 'gbex_excess_probability' for details.
  return(gbex_excess_probability(fitted_gbex, ...))
}

