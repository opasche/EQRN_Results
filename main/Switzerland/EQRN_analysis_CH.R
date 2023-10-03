library(tidyverse)
library(here)
library(ggpubr)
library(grf)
source("R/gbex_wrappers.R")
source("R/EGAM_wrappers.R")
source("R/EXQAR/EXQAR.R")
source("R/EQRN_loader.R")

## =============================== PARAMETERS ===============================
save_path <- "Results/Switzerland/EQRN_analysis_results/"
save_plots <- TRUE
force_refit <- FALSE
do_competitors <- TRUE

seedR <- 0
seedGRF <- 1
seedT <- seedR

# PARAM: Data
Y_name <- c("station_62")
X_disch_names <- c("station_43")
precip_names <- c("LTB","BRZ","MER","THU","GHS","BEP")
use_stl_remainder <- FALSE
use_mean_precip <- FALSE
prop_test <- 2/3
prop_valid <- 1/4

# PARAM: General
intermediate_method <- "qrn"
interm_lvl = 0.8
quantiles_predict = c(0.995,0.999,0.9995)
years_return <- 100
q_levels_nex <- 1-(1e-1*(10^seq(0,-3,-0.5)))

# Params: QRNN
interm_path <- "data/Switzerland/qrn_intermediate_quantile_best/"
par_qrn <- list(
  nb_fits=3,
  rnn_type="lstm",
  num_layers=2,
  hidden_size=256,
  p_drop=0,
  L2_pen=1e-6,
  seq_len=10,
  learning_rate=1e-3,
  n_epochs=1e3,
  batch_size=256,
  lr_decay=0.4,
  patience_decay=5,
  min_lr=1e-5,
  patience_stop=20,
  scale_features=TRUE,
  tol=1e-4
)

# ============== Params: EQRN ==============
path_eqrn <- "Results/Switzerland/EQRN_gridsearch_sf/"
par_eqrn <- list(
  nb_fits = 3,
  shape_fixed=TRUE,
  rnn_type="lstm",
  num_layers=2,
  hidden_size=16,
  p_drop=0,
  intermediate_q_feature=TRUE,
  L2_pen=1e-6,
  shape_penalty=0,
  seq_len=par_qrn$seq_len,
  learning_rate=1e-3,
  scale_features=TRUE,
  orthogonal_gpd=TRUE,
  n_epochs=1000,
  batch_size=256,
  lr_decay=0.4,
  patience_decay=5,
  min_lr=1e-5,
  patience_stop=20,
  patience_lag=5,
  tol=1e-5
)

# Params: GRF
grf_pars <- list(
  timedep = par_qrn$seq_len,
  num.trees = 5e3,
  quantiles_fit = c(0.1, 0.5, 0.9),
  sample.fraction = 0.5,
  #mtry = min(ceiling(sqrt(p) + 20), p),
  min.node.size = 5,#5
  honesty = TRUE,
  honesty.fraction = 0.5,
  honesty.prune.leaves = TRUE,
  alpha = 0.05,
  imbalance.penalty = 0
)

# Params: gbex
interm_method_competitors <- intermediate_method
gbex_params <- list(
  B=1293,
  lambda=NULL,
  lambda_ratio=5,
  lambda_scale=0.01,
  depth=c(2, 1),
  sf=0.75,
  intermediate_q_feature=TRUE,
  scale_features=FALSE)

# Params: egam
egam_params <- list(
  model_shape=FALSE,
  intermediate_q_feature=gbex_params$intermediate_q_feature,
  scale_features=gbex_params$scale_features)

# Params: EXQAR
exqar_params <- list(
  delta1 = 1-interm_lvl,
  delta2 = "default"
)

# Events and dates
ids_2005 <- (17227:17333)
date_2005_first <- as.Date("2005-08-22")
date_2005_flood <- as.Date("2005-08-23")
ids_datplt <- (1340:3165)

# Plots
mp_lwdt <- 150
mp_lhgt <- 30
mp_wdth <- 60
ap_wdth <- 210
ap_hght <- 82

## =============================== END PARAMETERS ===============================

set.seed(seedR)
check_directory(save_path, recursive=TRUE)
ts_datp_path <- paste0(save_path, "data_plots/")
check_directory(ts_datp_path, recursive=TRUE)
ex_proba_path <- paste0(save_path, "exceedence_probabilities/")
check_directory(ex_proba_path, recursive=TRUE)
if(do_competitors){check_directory(paste0(ex_proba_path,"competitors/"), recursive=TRUE)}

qrn_params_string <- paste0("qrrnn_", par_qrn$rnn_type, "_", par_qrn$num_layers, "x", par_qrn$hidden_size, "_s", par_qrn$seq_len, "_do", par_qrn$p_drop*100,
                            "_L2", str_replace(toString(par_qrn$L2_pen),"([.])","p"), "_lr", str_replace(toString(par_qrn$learning_rate),"([.])","p"))

params_string <- paste0("reqrnn_", par_eqrn$rnn_type, "_", par_eqrn$num_layers, "x", par_eqrn$hidden_size, "_",
                        "n"[!par_eqrn$intermediate_q_feature], "u_s", par_eqrn$seq_len, "_do", par_eqrn$p_drop*100,
                        "_L2", str_replace(toString(par_eqrn$L2_pen),"([.])","p"), "_lr", str_replace(toString(par_eqrn$learning_rate),"([.])","p"))

# For English month names in plots:
Sys.setlocale("LC_TIME", "English_UK")

start_time <- Sys.time()

# Data
if(!file.exists(paste0(interm_path,"Data_backup.rds"))){
  #warning("File not found, new data is created")
  stop("Data file not found")
}else{
  data_save <- readRDS(paste0(interm_path,"Data_backup.rds"))
  X <- data_save$X
  Y <- data_save$Y
  Dates <- data_save$Dates
  n_train <- data_save$n_train
  n_valid <- data_save$n_valid
  n_test <- data_save$n_test
  seq_len <- data_save$seq_len
  n_train_all <- n_train+n_valid
  n_tot <- n_train_all+n_test
}
rm(data_save)

if(par_qrn$seq_len!=seq_len){stop("different seq_len.")}

X_train <- X[1:n_train, , drop=F]
y_train <- Y[1:n_train]
X_valid <- X[(n_train+1-seq_len):(n_train+n_valid), , drop=F]
y_valid <- Y[(n_train+1-seq_len):(n_train+n_valid)]
X_test <- X[(n_train+n_valid+1-seq_len):(n_train+n_valid+n_test), , drop=F]
y_test <- Y[(n_train+n_valid+1-seq_len):(n_train+n_valid+n_test)]

X_train_all <- X[1:(n_train+n_valid), , drop=F]
y_train_all <- Y[1:(n_train+n_valid)]


for(ido in list(list(id=ids_datplt, nm="typical", event=NULL, db="12 months", dmb="1 month"),
                list(id=ids_2005, nm="2005_first", event=date_2005_first, db="1 month", dmb="1 day"))){
  for(pri in seq_along(precip_names)){
    plt_dat <- plot_series_comparison3(y_test[seq_len+ido$id], X_test[seq_len+ido$id,1], X_test[seq_len+ido$id,(pri+1)],
                                       Date_index=Dates[n_train_all+(ido$id)],
                                       var_names=c(latex2exp::TeX(r'(Bern $\[m^3s^{-1}\]$)'), latex2exp::TeX(r'(Gsteig $\[m^3s^{-1}\]$)'),
                                                   latex2exp::TeX(paste0(precip_names[pri], r'( $\[m^3s^{-1}\]$)'))),
                                       date_breaks=ido$db, date_minor_breaks=ido$dmb, event_dates=ido$event, show_top=0)
    # plot(plt_dat)
    if(!is.null(save_path)){
      save_myplot(plt=plt_dat, plt_nm=paste0(ts_datp_path, "ts_data_", ido$nm, "_62_43_",precip_names[pri],".pdf"),
                  width = mp_lwdt, height = mp_lhgt, cairo=FALSE)
    }
  }
}

## ======= INTERMEDIATE QUANTILES =======

if(intermediate_method=="qrn"){
  # QRNN fit
  foldwise_obj <- readRDS(paste0(interm_path, "Results_", qrn_params_string, ".rds"))
  foldwise_obj$fits <- list()
  for (k in 1:(foldwise_obj$n_folds)) {
    foldwise_obj$fits[[k]] <- EQRN_load(paste0(interm_path, "networks/"), paste0(qrn_params_string, "_f", k))
  }
  #Intermediate quantiles on all data
  interm_quantiles_all <- matrix(foldwise_obj$predictions, ncol=1)
  thresh_quant_all <- rbind(interm_quantiles_all, foldwise_obj$test_predictions[,(foldwise_obj$n_folds+1), drop=F])
}else if(intermediate_method=="oracle"){
  stop("No ground truth for real data.")
}
intermediate_quantiles <- thresh_quant_all[1:n_train, , drop=F]
valid_quantiles <- thresh_quant_all[(n_train+1-par_qrn$seq_len):(n_train+n_valid), , drop=F]
interm_quantiles_all <- thresh_quant_all[1:(n_train+n_valid), , drop=F]
#Intermediate quantile predictions on X_test
pred_interm <- thresh_quant_all[(n_train+n_valid+1-par_qrn$seq_len):(n_train+n_valid+n_test), , drop=F]


## ======== EQRN FIT ========

#Fit EQRN with intermediate (or oracle) quantiles
if((!dir.exists(paste0(path_eqrn, "networks/", params_string,"/"))) | force_refit){
  torch_manual_seed(seedT)
  fit_eqrn <- EQRN_fit_restart(X_train, y_train, intermediate_quantiles=intermediate_quantiles, interm_lvl=interm_lvl, number_fits=par_eqrn$nb_fits,
                                shape_fixed=par_eqrn$shape_fixed, hidden_size=par_eqrn$hidden_size, num_layers=par_eqrn$num_layers,
                                rnn_type=par_eqrn$rnn_type, p_drop=par_eqrn$p_drop, intermediate_q_feature=par_eqrn$intermediate_q_feature,
                                learning_rate=par_eqrn$learning_rate, L2_pen=par_eqrn$L2_pen, seq_len=par_eqrn$seq_len,
                                shape_penalty=par_eqrn$shape_penalty, scale_features=par_eqrn$scale_features, n_epochs=par_eqrn$n_epochs,
                                batch_size=par_eqrn$batch_size, X_valid=X_valid, y_valid=y_valid, quant_valid=valid_quantiles,
                                lr_decay=par_eqrn$lr_decay, patience_decay=par_eqrn$patience_decay, min_lr=par_eqrn$min_lr,
                                patience_stop=par_eqrn$patience_stop, tol=par_eqrn$tol, orthogonal_gpd=par_eqrn$orthogonal_gpd,
                                patience_lag=par_eqrn$patience_lag, data_type="seq")
  #EQRN_save(fit_eqrn, paste0(path_eqrn, "networks/"), params_string)
}else{
  fit_eqrn <- EQRN_load(paste0(path_eqrn, "networks/"), params_string)
  #if(file.exists(paste0(path_eqrn, "Results_",params_string,".rds"))){
  RESULTS_gs <- readRDS(paste0(path_eqrn, "Results_",params_string,".rds"))
  if(any(unlist(par_eqrn) != RESULTS_gs$output[1:length(par_eqrn)])){
    stop("Different parameters between saved and desired EQRN model.")
  }
  #}
}


cat("\nElapsed time (EQRN fit):\n")
print(Sys.time() - start_time)


## ======== TESTING (at levels 'quantiles_predict') ========

#Final EQRN quantile predictions on X_test (and X_valid)
pred_eqrn_val <- EQRN_predict_seq(fit_eqrn, X_valid, y_valid, quantiles_predict, valid_quantiles, interm_lvl, crop_predictions=TRUE)
pred_eqrn <- EQRN_predict_seq(fit_eqrn, X_test, y_test, quantiles_predict, pred_interm, interm_lvl, crop_predictions=TRUE)#[params$seq_len+(1:n_test)]
pred_eqrn_params <- EQRN_predict_params_seq(fit_eqrn, X_test, y_test, intermediate_quantiles=pred_interm,
                                            return_parametrization="classical")


# UNCONDITIONAL predicted quantile(s) (Y quantile on X_train)
pred_unc <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = quantiles_predict, Y = y_train_all, ntest = n_test)

# SEMI-CONDITIONAL predicted quantiles
pred_semicond <- predict_GPD_semiconditional(Y=y_train_all[(seq_len+1):n_train_all], interm_lvl=interm_lvl, thresh_quantiles=interm_quantiles_all[(seq_len+1):n_train_all],
                                             interm_quantiles_test=pred_interm[seq_len+(1:n_test)], quantiles_predict=quantiles_predict)

# Unconditional parameters and losses (for comparison with EQRN fit losses)
uncond_losses_fixed <- unconditional_train_valid_GPD_loss(Y_train=y_train[(seq_len+1):n_train], interm_lvl=interm_lvl, Y_valid=y_valid[seq_len+(1:n_valid)])
uncond_losses_interm <- semiconditional_train_valid_GPD_loss(Y_train=y_train[(seq_len+1):n_train], Y_valid=y_valid[seq_len+(1:n_valid)],
                                                             interm_quant_train=intermediate_quantiles[(seq_len+1):n_train],
                                                             interm_quant_valid=valid_quantiles[seq_len+(1:n_valid)])

## ======== Block maxima PoT and empirical fits for 100y event ========

p_daily_100y <- (1 - 1/(years_return*365.25)) # corresponding to 100y event

# EQRN prediciton at level 'p_daily_100y'
pred_eqrn_rl <- EQRN_predict_seq(fit_eqrn, X_test, y_test, p_daily_100y, pred_interm, interm_lvl, crop_predictions=TRUE)

#Select observations and date index for block maxima (and PoT, empirical) approaches
Y_rl <- y_train_all
Dates_rl <- Dates[1:n_train_all]

Ymax_table <- tibble(Year=lubridate::year(Dates_rl), Y=Y_rl) %>% group_by(Year) %>% summarise(ymax=max(Y))
Ymax <- Ymax_table$ymax[1:length(Ymax_table$ymax)-1]
gev_fit <- ismev::gev.fit(Ymax, method="Nelder-Mead", maxit=1e6, show=FALSE)
return_lvl_100y <- evd::qgev(1-1/years_return, loc=gev_fit$mle[1], scale=gev_fit$mle[2], shape=gev_fit$mle[3])
PoT_100yrl <- c(predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = p_daily_100y, Y = Y_rl, ntest = 1)$predictions)
empquant_100yrl <- quantile(Y_rl,p_daily_100y)


## ========  EQRN 100y excess probability predictions on X_test (and X_valid) for 100y event ========
excess_prob_100y_test <- EQRN_excess_probability_seq(val=return_lvl_100y, fit_eqrn=fit_eqrn, X=X_test, Y=y_test,
                                                     intermediate_quantiles=pred_interm, interm_lvl=interm_lvl,
                                                     crop_predictions=TRUE, body_proba="default", proba_type="excess")
excess_prob_100y_valid <- EQRN_excess_probability_seq(val=return_lvl_100y, fit_eqrn=fit_eqrn, X=X_valid, Y=y_valid,
                                                      intermediate_quantiles=valid_quantiles, interm_lvl=interm_lvl,
                                                      crop_predictions=TRUE, body_proba="default", proba_type="excess")



## ======== Competitor FITS =========
if(do_competitors){
  lagged_train_all <- lagged_features(X=cbind(y_train_all,X_train_all), max_lag=grf_pars$timedep, drop_present=TRUE)
  lagged_test <- lagged_features(X=cbind(y_test,X_test), max_lag=grf_pars$timedep, drop_present=TRUE)
  
  if(file.exists(paste0(save_path, "competitor_fits/", "competitor_preds.rds"))){
    competitor_preds <- readRDS(paste0(save_path, "competitor_fits/", "competitor_preds.rds"))
    for(nm in names(competitor_preds)){assign(nm, competitor_preds[[nm]])}
    rm(competitor_preds)
  }else{
    
    fit_grf <- quantile_forest(lagged_train_all, y_train_all[(grf_pars$timedep+1):length(y_train_all)],
                               num.trees=grf_pars$num.trees, quantiles=grf_pars$quantiles_fit, sample.fraction=grf_pars$sample.fraction,# mtry=grf_pars$mtry,
                               min.node.size=grf_pars$min.node.size, honesty=grf_pars$honesty,
                               honesty.fraction=grf_pars$honesty.fraction, honesty.prune.leaves=grf_pars$honesty.prune.leaves,
                               alpha=grf_pars$alpha, imbalance.penalty=grf_pars$imbalance.penalty, seed=seedGRF)#,compute.oob.predictions = FALSE
    
    if(interm_method_competitors=="qrn"){
      # QRNN fit
      foldwise_obj <- readRDS(paste0(interm_path, "Results_", qrn_params_string, ".rds"))
      foldwise_obj$fits <- list()
      for (k in 1:(foldwise_obj$n_folds)) {
        foldwise_obj$fits[[k]] <- EQRN_load(paste0(interm_path, "networks/"), paste0(qrn_params_string, "_f", k))
      }
      #Intermediate quantiles on all data
      thresh_quant_all_c <- rbind(matrix(foldwise_obj$predictions, ncol=1),
                                  foldwise_obj$test_predictions[,(foldwise_obj$n_folds+1), drop=F])
    }else if(interm_method_competitors=="grf"){
      #Intermediate quantiles all data
      pred_grf_interm_train_all <- predict(fit_grf, newdata=NULL, quantiles = c(interm_lvl))$predictions
      pred_grf_interm_test <- predict(fit_grf, newdata=lagged_test, quantiles = c(interm_lvl))$predictions
      thresh_quant_all_c <- rbind(matrix(as.double(NA), nrow=grf_pars$timedep, ncol=1), pred_grf_interm_train_all, pred_grf_interm_test)
    }else if(interm_method_competitors=="oracle"){
      stop("No ground truth for real data.")
    }
    interm_quantiles_all_c <- thresh_quant_all_c[1:(n_train+n_valid), , drop=F]
    pred_interm_c <- thresh_quant_all_c[(n_train+n_valid+1-par_qrn$seq_len):(n_train+n_valid+n_test), , drop=F]
    
    lagged_interm_q_trall <- interm_quantiles_all_c[(grf_pars$timedep+1):length(interm_quantiles_all_c), , drop=F]
    laged_interm_q_test <- pred_interm_c[(grf_pars$timedep+1):length(pred_interm_c), , drop=F]
    #High quantile prediction with GRF
    pred_grf_test <- predict(fit_grf, newdata=lagged_test, quantiles = quantiles_predict)$predictions
    pred_grf_rl <- predict(fit_grf, newdata=lagged_test, quantiles = p_daily_100y)$predictions
    pred_grf_nex <- predict(fit_grf, newdata=lagged_test, quantiles = q_levels_nex)$predictions
    #pred_grf_train_all <- predict(fit_grf, newdata=NULL, quantiles = quantiles_predict)$predictions
    
    fit_gbex <- gbex_fit(X=lagged_train_all, y=y_train_all[(grf_pars$timedep+1):length(y_train_all)], intermediate_quantiles=lagged_interm_q_trall,
                         interm_lvl=interm_lvl, intermediate_q_feature=gbex_params$intermediate_q_feature, scale_features=gbex_params$scale_features,
                         B=gbex_params$B, lambda=gbex_params$lambda, lambda_ratio=gbex_params$lambda_ratio,
                         lambda_scale=gbex_params$lambda_scale, depth=gbex_params$depth, sf=gbex_params$sf)
    fit_egam <- fit_gpd_gam(X=lagged_train_all, y=y_train_all[(grf_pars$timedep+1):length(y_train_all)],
                            intermediate_quantiles=lagged_interm_q_trall, interm_lvl=interm_lvl, model_shape=egam_params$model_shape,
                            intermediate_q_feature=egam_params$intermediate_q_feature, scale_features=egam_params$scale_features)
    fit_exqar <- EXQAR_fit(y_train_all, X=X_train_all, lag=grf_pars$timedep,
                           delta1=exqar_params$delta1, delta2=exqar_params$delta2)
    
    pred_gbex <- gbex_predict(fit_gbex, lagged_test, to_predict=quantiles_predict, intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
    pred_gbex_rl <- gbex_predict(fit_gbex, lagged_test, to_predict=p_daily_100y, intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
    pred_egam <- predict_gpd_gam(fit_egam, lagged_test, to_predict=quantiles_predict,
                                 intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
    pred_egam_rl <- predict_gpd_gam(fit_egam, lagged_test, to_predict=p_daily_100y,
                                    intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
    pred_exqar <- EXQAR_predict(fit_exqar, y_test, X=X_test, quantiles_predict, tol=1e-4, min_prop=0.3, return_infos=FALSE)
    pred_exqar_rl <- EXQAR_predict(fit_exqar, y_test, X=X_test, p_daily_100y, tol=1e-4, min_prop=0.3, return_infos=FALSE)
    pred_gbex_nex <- gbex_predict(fit_gbex, lagged_test, to_predict=q_levels_nex, intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
    pred_egam_nex <- predict_gpd_gam(fit_egam, lagged_test, to_predict=q_levels_nex,
                                     intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
    pred_exqar_nex <- EXQAR_predict(fit_exqar, y_test, X=X_test, q_levels_nex, tol=1e-4, min_prop=0.3, return_infos=FALSE)
    
    # ==== excess probs competitors ====
    
    excess_prob_100y_gbex <- gbex_excess_probability(fit_gbex, val=return_lvl_100y, X_test=lagged_test, intermediate_quantiles=laged_interm_q_test,
                                                     interm_lvl=interm_lvl, body_proba="default", proba_type="excess")
    excess_prob_100y_egam <- excess_probability_gpd_gam(fit_egam, val=return_lvl_100y, X_test=lagged_test, intermediate_quantiles=laged_interm_q_test,
                                                        interm_lvl=interm_lvl, body_proba="default", proba_type="excess")
    excess_prob_100y_exqar <- EXQAR_excess_probability(fit_exqar, val=return_lvl_100y, y=y_test, X=X_test, body_proba=1-interm_lvl,#as.double(NA),
                                                       proba_type="excess", tol=1e-4, min_prop=0.3, return_infos=FALSE)
    
    competitor_preds <- list(pred_gbex=pred_gbex, pred_gbex_rl=pred_gbex_rl, pred_gbex_nex=pred_gbex_nex,
                             pred_egam=pred_egam, pred_egam_rl=pred_egam_rl, pred_egam_nex=pred_egam_nex,
                             pred_exqar=pred_exqar, pred_exqar_rl=pred_exqar_rl, pred_exqar_nex=pred_exqar_nex,
                             pred_grf_test=pred_grf_test, pred_grf_rl=pred_grf_rl, pred_grf_nex=pred_grf_nex,
                             excess_prob_100y_gbex=excess_prob_100y_gbex, excess_prob_100y_egam=excess_prob_100y_egam,
                             excess_prob_100y_exqar=excess_prob_100y_exqar)
    safe_save_rds(competitor_preds, paste0(save_path, "competitor_fits/", "competitor_preds.rds"))
    safe_save_rds(fit_exqar, paste0(save_path, "competitor_fits/", "exqar.rds"))
    safe_save_rds(fit_egam, paste0(save_path, "competitor_fits/", "egam.rds"))
    # safe_save_rds(fit_gbex, paste0(save_path, "competitor_fits/", "gbex.rds"))
    rm(list=c("fit_exqar","fit_egam","fit_gbex","competitor_preds"))
  }
}

excess_prob_100y_scond <- GPD_excess_probability(val=return_lvl_100y, sigma=pred_semicond$pars[[1]], xi=pred_semicond$pars[[2]],
                                                 interm_threshold=pred_interm[seq_len+(1:n_test)],
                                                 threshold_p=interm_lvl, body_proba="default", proba_type="excess")
excess_prob_100y_uncond <- GPD_excess_probability(val=return_lvl_100y, sigma=pred_unc$pars[,1], xi=pred_unc$pars[,2],
                                                  interm_threshold=rep(quantile(y_train_all, interm_lvl), n_test),
                                                  threshold_p=interm_lvl, body_proba="default", proba_type="excess")

pred_unc_rl <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = p_daily_100y, Y = y_train_all, ntest = n_test)
pred_unc_nex <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = q_levels_nex, Y = y_train_all, ntest = n_test)
pred_semicond_rl <- predict_GPD_semiconditional(Y=y_train_all[(seq_len+1):n_train_all], interm_lvl=interm_lvl, thresh_quantiles=interm_quantiles_all[(seq_len+1):n_train_all],
                                                interm_quantiles_test=pred_interm[seq_len+(1:n_test)], quantiles_predict=p_daily_100y)
pred_semicond_nex <- predict_GPD_semiconditional(Y=y_train_all[(seq_len+1):n_train_all], interm_lvl=interm_lvl, thresh_quantiles=interm_quantiles_all[(seq_len+1):n_train_all],
                                                interm_quantiles_test=pred_interm[seq_len+(1:n_test)], quantiles_predict=q_levels_nex)


## ============= Quantile and exceedance proba prediction plots ==============

zoom_windows <- list(list(name="2005_first", ind=ids_2005, event=date_2005_first, db="1 month"),
                     list(name="2005_flood_day", ind=ids_2005, event=date_2005_flood, db="1 month"))

if(do_competitors){
  comp_excesses_itterator <- list(list(name=params_string, p=excess_prob_100y_test, q=c(pred_eqrn_rl), f="", c="EQRN"),
                                  list(name="gbex", p=excess_prob_100y_gbex, q=c(pred_gbex_rl), f="competitors/", c="GBEX"),
                                  list(name="egam", p=excess_prob_100y_egam, q=c(pred_egam_rl), f="competitors/", c="EGAM"),
                                  list(name="exqar", p=excess_prob_100y_exqar, q=c(pred_exqar_rl), f="competitors/", c="EXQAR"),
                                  list(name="semicond", p=excess_prob_100y_scond, q=c(pred_semicond_rl$predictions), f="competitors/", c="Semicond"),
                                  list(name="uncond", p=excess_prob_100y_uncond, q=c(pred_unc_rl$predictions), f="competitors/", c="Uncond"))
}else{
  comp_excesses_itterator <- list(list(name=params_string, p=excess_prob_100y_test, q=c(pred_eqrn_rl), f=""))
}
for(j in seq_along(comp_excesses_itterator)){
  cei <- comp_excesses_itterator[[j]]
  for(i in seq_along(zoom_windows)){
    zw <- zoom_windows[[i]]
    p_pvrl <- plot_pred_vs_return_lvl_ts(cei$q[zw$ind], return_lvl_100y, PoT_100yrl, NaN,#empquant_100yrl,
                                         latex2exp::TeX(r'(Discharge $\[m^3s^{-1}\]$)'), y_obs=y_test[seq_len+(zw$ind)],
                                         Date_index=Dates[n_train_all+(zw$ind)], date_breaks=zw$db,
                                         event_dates=zw$event, legend.position="none", color_met=my_palette_methods[[cei$c]])
    plot(p_pvrl)
    p_evpts <- plot_exceedence_proba_ts(cei$p[zw$ind], p_daily_100y, Date_index=Dates[n_train_all+(zw$ind)],
                                        type_probas="probability", proba_ratio=TRUE, date_breaks=zw$db, event_dates=zw$event,
                                        color_met=my_palette_methods[[cei$c]])
    plot(p_evpts)
    p_predrisk <- plot_pred_quant_risk_ts(cei$q[zw$ind], cei$p[zw$ind], p_daily_100y, return_lvl_100y, PoT_100yrl, NaN,#empquant_100yrl,
                                          latex2exp::TeX(r'(Discharge $\[m^3s^{-1}\]$)'), y_obs=y_test[seq_len+(zw$ind)],
                                          Date_index=Dates[n_train_all+(zw$ind)], type_probas="probability",
                                          proba_ratio=TRUE, date_breaks=zw$db, event_dates=zw$event, legend.position="none",
                                          p_lab="'Probability ratio'", color_met=my_palette_methods[[cei$c]])
    plot(p_predrisk)
    if(save_plots){
      save_myplot(plt=p_pvrl, plt_nm=paste0(ex_proba_path, cei$f, cei$name, "_",years_return,"y_pred_vs_return_lvl_", zw$name,".pdf"),
                  width = mp_lwdt, height = 50, cairo=FALSE)
      save_myplot(plt=p_evpts, plt_nm=paste0(ex_proba_path, cei$f, cei$name, "_",years_return,"y_exceedence_proba_ts_", zw$name,".pdf"),
                  width = mp_lwdt, height = 50, cairo=FALSE)
      save_myplot(plt=p_predrisk, plt_nm=paste0(ex_proba_path, cei$f, cei$name, "_",years_return,"y_pred_both_", zw$name,"_pratio.pdf"),
                  width = mp_lwdt, height = mp_lhgt, cairo=FALSE)
      save_myplot(plt=p_predrisk, plt_nm=paste0(ex_proba_path, cei$f, cei$name, "_",years_return,"y_pred_both_", zw$name,"_pratio_w.pdf"),
                  width = mp_lwdt, height = 40, cairo=FALSE)
    }
  }
}

if(do_competitors){
  fit_exqar <- readRDS(paste0(save_path, "competitor_fits/exqar.rds"))
  pred_exqar_rldet <- EXQAR_predict(fit_exqar, y_test, X=X_test, p_daily_100y, tol=1e-4, min_prop=0.3, return_infos=TRUE)
  p_exqar_expl <- plot_series_comparison4(c(pred_exqar_rl)[ids_2005], pred_exqar_rldet$EVI[ids_2005], pred_exqar_rldet$am[ids_2005],
                                          pred_exqar_rldet$Ut[ids_2005], Date_index=Dates[n_train_all+(ids_2005)],
                                          var_names=c("'EXQAR EQ-pred'", "'EXQAR EVI pred'", "'EXQAR scale pred'", "'EXQAR Q0-pred'"),
                                          date_breaks="1 month", date_minor_breaks="1 day", event_dates=date_2005_first, show_top=0)
  save_myplot(plt=p_exqar_expl, plt_nm=paste0(ex_proba_path, "competitors/exqar_", years_return,"y_expl_2005_first.pdf"),
              width = mp_lwdt, height = mp_lhgt, cairo=FALSE)
}

Dates_pred_test <- Dates[n_train_all+(1:n_test)]
zw_featspred <- 17246:17287
for(pri in seq_along(precip_names)){
  plt_feat_pred <- plot_features_for_pred_ts(y_test[seq_len+zw_featspred], X_test[seq_len+zw_featspred,1], X_test[seq_len+zw_featspred,(pri+1)],
                                             pred_eqrn_rl[zw_featspred], return_lvl_100y, event_ind=date_2005_first,
                                             Date_index=Dates[n_train_all+zw_featspred], seq_len=seq_len,
                                             var_names=c(latex2exp::TeX(r'(Bern $\[m^3s^{-1}\]$)'), latex2exp::TeX(r'(Gsteig $\[m^3s^{-1}\]$)'),
                                                         latex2exp::TeX(paste0(precip_names[pri], r'( $\[m^3s^{-1}\]$)'))),
                                             date_breaks="1 month", date_minor_breaks="1 day", legend.position="none")
  # plot(plt_dat)
  if(!is.null(save_path)){
    save_myplot(plt=plt_feat_pred, plt_nm=paste0(ts_datp_path,"pred_feat_62_43_",precip_names[pri],".pdf"),
                width = mp_lwdt, height = mp_lhgt, cairo=FALSE)
  }
}


## ============= Quantile estimation windows for non-stationarity ==============

Ymax_table_all <- tibble(Year=lubridate::year(Dates), Y=Y) %>% group_by(Year) %>% summarise(ymax=max(Y))

quantiles_ns_window <- foreach(yr=(1959:2014), .errorhandling="stop", .combine=bind_rows) %do% {
  inds_yr <- which(lubridate::year(Dates)==yr)
  start_yr <- inds_yr[1]
  end_yr <- last_elem(inds_yr)
  Ymax_l <- filter(Ymax_table_all, Year<yr)$ymax
  gev_fit_l <- ismev::gev.fit(Ymax_l, method="Nelder-Mead", maxit=1e6, show=FALSE)
  return_lvl_100y_l <- evd::qgev(1-1/years_return, loc=gev_fit_l$mle[1], scale=gev_fit_l$mle[2], shape=gev_fit_l$mle[3])
  PoT_100yrl_l80 <- c(predict_unconditional_quantiles(interm_lvl = 0.8, quantiles = p_daily_100y, Y = Y[1:(start_yr-1)], ntest = 1)$predictions)
  PoT_100yrl_l90 <- c(predict_unconditional_quantiles(interm_lvl = 0.9, quantiles = p_daily_100y, Y = Y[1:(start_yr-1)], ntest = 1)$predictions)
  PoT_100yrl_l95 <- c(predict_unconditional_quantiles(interm_lvl = 0.95, quantiles = p_daily_100y, Y = Y[1:(start_yr-1)], ntest = 1)$predictions)
  empquant_100yrl_l <- quantile(Y[1:(start_yr-1)],p_daily_100y)
  
  w_year <- yr - 29
  w_inds_yr <- which(lubridate::year(Dates)==w_year)
  start_window <- w_inds_yr[1]
  Ymax_w <- filter(Ymax_table_all, Year<yr & Year>=w_year)$ymax
  gev_fit_w <- ismev::gev.fit(Ymax_w, method="Nelder-Mead", maxit=1e6, show=FALSE)
  return_lvl_100y_w <- evd::qgev(1-1/years_return, loc=gev_fit_w$mle[1], scale=gev_fit_w$mle[2], shape=gev_fit_w$mle[3])
  PoT_100yrl_w80 <- c(predict_unconditional_quantiles(interm_lvl = 0.8, quantiles = p_daily_100y, Y = Y[start_window:(start_yr-1)], ntest = 1)$predictions)
  PoT_100yrl_w90 <- c(predict_unconditional_quantiles(interm_lvl = 0.9, quantiles = p_daily_100y, Y = Y[start_window:(start_yr-1)], ntest = 1)$predictions)
  PoT_100yrl_w95 <- c(predict_unconditional_quantiles(interm_lvl = 0.95, quantiles = p_daily_100y, Y = Y[start_window:(start_yr-1)], ntest = 1)$predictions)
  empquant_100yrl_w <- quantile(Y[start_window:(start_yr-1)],p_daily_100y)
  
  pred_eqrn_rl_l <- EQRN_predict_seq(fit_eqrn, X[(start_yr-seq_len):end_yr, , drop=F], Y[(start_yr-seq_len):end_yr], p_daily_100y,
                                     thresh_quant_all[(start_yr-seq_len):end_yr, , drop=F], interm_lvl, crop_predictions=TRUE)
  pred_eqrn_rl_lav <- mean(pred_eqrn_rl_l)
  pred_eqrn_rl_lmax <- max(pred_eqrn_rl_l)
  pred_eqrn_rl_lmin <- min(pred_eqrn_rl_l)
  pred_eqrn_rl_lmed <- median(pred_eqrn_rl_l)
  pred_eqrn_rl_lm5 <- mean(tail(sort(pred_eqrn_rl_l),5))
  Ym1 <- mean(tail(sort(Y[start_yr:end_yr]),1))
  
  tibble(Year=yr, EQRN_av=pred_eqrn_rl_lav, EQRN_max=pred_eqrn_rl_lmax, EQRN_min=pred_eqrn_rl_lmin, EQRN_med=pred_eqrn_rl_lmed,
         EQRN_m5 = pred_eqrn_rl_lm5, return_lvl=return_lvl_100y_l, PoT80=PoT_100yrl_l80,
         PoT90=PoT_100yrl_l90, PoT95=PoT_100yrl_l95, empirical=empquant_100yrl_l, rl0=return_lvl_100y, Ym1=Ym1,
         return_lvl_w=return_lvl_100y_w, PoT80w=PoT_100yrl_w80, PoT90w=PoT_100yrl_w90, PoT95w=PoT_100yrl_w95, empiricalw=empquant_100yrl_w)
}
quantiles_ns_window <- quantiles_ns_window %>% full_join(Ymax_table_all,by=c("Year"="Year","Ym1"="ymax")) %>% arrange(Year)
quantiles_ns_window$rl0 <- return_lvl_100y
qevol_plt <- quantiles_ns_window %>% select(Year, return_lvl, rl0, Ym1) %>% #, empirical
  rename(`GEV evol`=return_lvl, `GEV train`=rl0) %>% #, `empirical`=empirical
  gather(key="Method", value="Quantile", -Year, -Ym1, factor_key=TRUE) %>%
  ggplot(aes(x=Year, y=Quantile, group=Method, color=Method, linetype=Method)) +
  geom_line(size=1) + scale_linetype_manual(values=c("solid","dashed")) +
  scale_color_manual(values=c(my_palette$green,my_palette$red)) +
  geom_point(aes(y=Ym1, fill="Ymax"), color="darkslategrey", alpha=0.8) + scale_fill_manual("", breaks="Ymax", values=c(Ymax="darkslategrey")) +
  labs(x="Year", y=latex2exp::TeX(r'(Discharge $\[m^3s^{-1}\]$)'), color=NULL, linetype=NULL) +
  scale_x_continuous(expand=c(0.01,0)) + scale_y_continuous(expand=c(0.02,0))
plot(qevol_plt)

qwindow_plt <- quantiles_ns_window %>% select(Year, return_lvl_w, rl0, Ym1) %>% #, empiricalw
  rename(`GEV evol`=return_lvl_w, `GEV train`=rl0) %>% #, `empirical`=empiricalw
  gather(key="Method", value="Quantile", -Year, -Ym1, factor_key=TRUE) %>%
  ggplot(aes(x=Year, y=Quantile, group=Method, color=Method, linetype=Method)) +
  geom_line(size=1) + scale_linetype_manual(values=c("solid","dashed")) +
  scale_color_manual(values=c(my_palette$green,my_palette$red)) +
  geom_point(aes(y=Ym1, fill="Ymax"), color="darkslategrey", alpha=0.8) + scale_fill_manual("", breaks="Ymax", values=c(Ymax="darkslategrey")) +
  labs(x="Year", y=latex2exp::TeX(r'(Discharge $\[m^3s^{-1}\]$)'), color=NULL, linetype=NULL) +
  scale_x_continuous(expand=c(0.01,0)) + scale_y_continuous(expand=c(0.02,0))
plot(qwindow_plt)
if(save_plots){
  save_myplot(plt=qevol_plt, plt_nm=paste0(save_path,"evolution_quantiles.pdf"),
              width = mp_lwdt, height = 75, cairo=FALSE)
  save_myplot(plt=qevol_plt, plt_nm=paste0(save_path,"evolution_quantiles_square.pdf"),
              width = mp_wdth, height = mp_wdth, cairo=FALSE)
  save_myplot(plt=qwindow_plt, plt_nm=paste0(save_path,"running_window_quantiles.pdf"),
              width = mp_lwdt, height = 75, cairo=FALSE)
}


## ============= Quantile calibration plot for EQRN ==============

pred_eqrn_nex <- EQRN_predict_seq(fit_eqrn, X_test, y_test, q_levels_nex, pred_interm, interm_lvl, crop_predictions=TRUE)
nex_plt <- plot_exceedences_quantile(pred_eqrn_nex, y_test[seq_len+(1:n_test)], quantile_levels=roundm(q_levels_nex,4),
                                     show_legend=FALSE)
plot(nex_plt)
if(!is.null(save_path)){
  save_myplot(plt=nex_plt, plt_nm=paste0(save_path,"number_exceedences_quant.pdf"),
              width = mp_wdth, height = mp_wdth, cairo=FALSE)
}

if(do_competitors){
  nex_plt_c <- plot_exceedences_quantile_comp(pred_eqrn_nex, pred_grf_nex, pred_unc_nex$predictions, pred_semicond_nex$predictions, 
                                              pred_exqar_nex, pred_gbex_nex, pred_egam_nex,
                                              y=y_test[seq_len+(1:n_test)], quantile_levels=roundm(q_levels_nex,4), 
                                              legend.position="bottom")
  nex_plt_d <- plot_exceedences_quantile_diff(pred_eqrn_nex, pred_grf_nex, pred_unc_nex$predictions, pred_semicond_nex$predictions, 
                                              pred_exqar_nex, pred_gbex_nex, pred_egam_nex,
                                              y=y_test[seq_len+(1:n_test)], quantile_levels=roundm(q_levels_nex,4), type_diff="diff", 
                                              legend.position="bottom")
  plot(nex_plt_c)
  if(!is.null(save_path)){
    save_myplot(plt=nex_plt_c, plt_nm=paste0(save_path,"number_exceedences_quant_c.pdf"),
                width = mp_wdth, height = mp_wdth, cairo=FALSE)
    save_myplot(plt=nex_plt_d, plt_nm=paste0(save_path,"number_exceedences_quant_d.pdf"),
                width = mp_wdth, height = mp_wdth, cairo=FALSE)
  }
}

## ============= Validation curve during EQRN training ==============

valid_plot <- validation_plot_eqrn(fit_eqrn, uncond_losses_interm, show_legend=FALSE)
plot(valid_plot)
if(save_plots){
  save_myplot(plt=valid_plot, plt_nm=paste0(save_path, params_string, "_validation.pdf"),
              width = mp_wdth, height = mp_wdth, cairo=FALSE)
}


## ============= Identifying other clusters of 100y exceedances ==============


cat("\nExceedance cluster identification with EQRN:\n")

dft <- data.frame(Date=Dates_pred_test, p=excess_prob_100y_test, q=c(pred_eqrn_rl), y=y_test[seq_len+(1:n_test)])

p_ratio100 <- dft$p / (1-p_daily_100y)


pratio_thresh <- 100
which(p_ratio100>pratio_thresh)
inds_clust <- which(diff(p_ratio100>pratio_thresh)>0)+1
size_clusters <- (which(diff(p_ratio100>pratio_thresh)<0)+1) - inds_clust
names(size_clusters) <- inds_clust
cat("\n Test set: on average", length(size_clusters)/n_test * 365.25,
    "early warnings are triggered per year with a ratio above", pratio_thresh,"\n")


first_cluster_exceedences_ind <- which(diff(dft$y > return_lvl_100y)>0)+1 # 14986 16834 17279 17996
cat("\n Test set:", sum((p_ratio100>pratio_thresh)[first_cluster_exceedences_ind]), "of the",
    length(first_cluster_exceedences_ind), "clusters of exceedences were forecast early with a ratio above", pratio_thresh,"\n")
print((p_ratio100>pratio_thresh)[first_cluster_exceedences_ind])



end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)



