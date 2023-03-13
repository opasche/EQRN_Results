# Dependencies
library(tidyverse)
library(evgam)
library(sfsmisc)
#source("R/EQRN_loader.R")

# Functions
fit_gpd_gam <- function(X, y, intermediate_quantiles, interm_lvl=0.8, model_shape=FALSE,
                        intermediate_q_feature=FALSE, scale_features=FALSE, ...) {
  ## fit an extreme GAM
  
  data_excesses <- get_excesses(X=X, y=y, quantiles=intermediate_quantiles,
                                intermediate_q_feature=intermediate_q_feature,
                                scale_features=scale_features, X_scaling=NULL)
  Y_excesses <- data_excesses$Y_excesses
  X_feats_excesses <- data_excesses$X_excesses
  X_scaling <- data_excesses$X_scaling
  
  # reformat data
  p <- ncol(X_feats_excesses)
  
  data <- bind_cols(
    as_tibble(X_feats_excesses, .name_repair = ~ paste0("X", 1:p)),
    Y_excesses = Y_excesses
  )
  
  
  # Prepare formula for evgam
  
  fmla_gpd <- list(
    sfsmisc::wrapFormula(
      Y_excesses ~ .,
      data %>% select(-Y_excesses)
    ),
    if (model_shape) {
      sfsmisc::wrapFormula(
        Y_excesses ~ .,
        data %>% select(-Y_excesses)
      )
    } else {
      Y_excesses ~ 1
    }
  )
  
  # return fitted gpd_gam
  structure(list(
    "evgam" = evgam::evgam(fmla_gpd, data, family = "gpd", ...),
    "interm_lvl" = interm_lvl,
    intermediate_q_feature=intermediate_q_feature,
    scale_features=scale_features,
    X_scaling=X_scaling
  ),
  class = "gpd_gam"
  )
}

predict_gpd_gam <- function(fitted_gpd_gam, X_test, to_predict=c(0.95, 0.99),
                            intermediate_quantiles, interm_lvl=fitted_gpd_gam$interm_lvl) {
  ## predicts extreme quantiles using evgam
  if(interm_lvl!=fitted_gpd_gam$interm_lvl){stop("gpd_gam intermediate quantiles interm_lvl does not match in train and predict.")}
  
  X_feats <- X_test
  if(fitted_gpd_gam$intermediate_q_feature){
    if(is.null(intermediate_quantiles)){stop("intermediate_quantiles needed in predict_gpd_gam for intermediate_q_feature.")}
    X_feats <- cbind(X_feats, intermediate_quantiles)
  }
  X_feats <- perform_scaling(X=X_feats, X_scaling=fitted_gpd_gam$X_scaling, scale_features=fitted_gpd_gam$scale_features)$X_scaled
  
  # validate and reformat newdata
  p <- ncol(X_feats)
  data_test <- as_tibble(
    X_feats,
    .name_repair = ~ paste0("X", 1:p)
  )
  
  # compute optimal GPD parameters
  pred_params <- predict(fitted_gpd_gam$evgam, data_test)
  pred_params[, 1] <- exp(pred_params[, 1])
  
  # write is as tibble
  pred_params <- tibble::tibble(
    "sigma" = pred_params[, 1],
    "xi" = pred_params[, 2]
  )
  
  if(any(to_predict=="par")){
    return(pred_params)
  }
  
  if(length(dim(to_predict))>1){
    stop("Please provide a single value or 1D vector as to_predict in predict_gpd_gam")
  }
  
  if(length(to_predict)==1){
    return(GPD_quantiles(to_predict, interm_lvl, intermediate_quantiles, pred_params$sigma, pred_params$xi))
  } else if(length(to_predict)>1){
    nb_quantiles_predict <- length(to_predict)
    egam_quantiles <- matrix(as.double(NA), nrow=nrow(X_test), ncol=nb_quantiles_predict)
    for(i in 1:nb_quantiles_predict){
      egam_quantiles[,i] <- GPD_quantiles(to_predict[i], interm_lvl, intermediate_quantiles, pred_params$sigma, pred_params$xi)
    }
    return(egam_quantiles)
  } else {
    stop("Please provide a single value or 1D vector as to_predict in predict_gpd_gam")
  }
}

excess_probability_gpd_gam <- function(fitted_gpd_gam, val, X_test, intermediate_quantiles, interm_lvl=fitted_gpd_gam$interm_lvl,
                                       body_proba="default", proba_type=c("excess","cdf")){
  ## predicts extreme quantiles using evgam
  proba_type <- match.arg(proba_type)
  if(interm_lvl!=fitted_gpd_gam$interm_lvl){stop("gpd_gam intermediate quantiles interm_lvl does not match in train and predict.")}
  
  X_feats <- X_test
  if(fitted_gpd_gam$intermediate_q_feature){
    if(is.null(intermediate_quantiles)){stop("intermediate_quantiles needed in predict_gpd_gam for intermediate_q_feature.")}
    X_feats <- cbind(X_feats, intermediate_quantiles)
  }
  X_feats <- perform_scaling(X=X_feats, X_scaling=fitted_gpd_gam$X_scaling, scale_features=fitted_gpd_gam$scale_features)$X_scaled
  
  # validate and reformat newdata
  p <- ncol(X_feats)
  data_test <- as_tibble(
    X_feats,
    .name_repair = ~ paste0("X", 1:p)
  )
  
  # compute optimal GPD parameters
  pred_params <- predict(fitted_gpd_gam$evgam, data_test)
  pred_params[, 1] <- exp(pred_params[, 1])
  
  # write as tibble
  pred_params <- tibble::tibble(
    "sigma" = pred_params[, 1],
    "xi" = pred_params[, 2]
  )
  
  Probs <- GPD_excess_probability(val, sigma=pred_params$sigma, xi=pred_params$xi, interm_threshold=intermediate_quantiles,
                                  threshold_p=interm_lvl, body_proba=body_proba, proba_type=proba_type)
  return(c(Probs))
}

