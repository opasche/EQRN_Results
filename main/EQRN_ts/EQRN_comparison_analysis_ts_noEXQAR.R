# This is a version of the script 'EQRN_comparison_analysis_ts.R' without the competitor EXQAR,
# for reproducibility until we get permission to share our implementation of the method.

library(tidyverse)
library(here)
library(ggpubr)
library(grf)
source("R/gbex_wrappers.R")
source("R/EGAM_wrappers.R")
source("R/EQRN_loader.R")

## =============================== PARAMETERS ===============================
save_path <- "Results/ts/EQRN_comp_pen_nopen_comprnn_noEXQAR/"
check_directory(save_path, recursive=TRUE)
save_plots <- TRUE
force_refit <- FALSE

seedR <- 0
seedGRF <- 1
seedT <- seedR
set.seed(seedR)

# PARAM: Data
n <- 5e3
ntest <- 5e3
n_valid <- 2e3
X_distr="foldnormal"
Y_distr="foldnormal"
df <- 4
alphas=c(.2,.1,.1,.1,.1)
betas=c(0)
sX=c(.3,.2,.1,.1,.1)
S0=1
AR=c(0)
MA=c(0)
muX=c(0)
mu0=0
ARX=c(0.4)
seasonal_hetero=0

# PARAM: General
intermediate_method <- "qrn"#"qrn""oracle"
interm_lvl = 0.8
quantiles_predict = c(0.995,0.999,0.9995)

# Params: QRNN
interm_path <- "data/simulations/ts_intermediate_quantile_best/"
par_qrn <- list(
  nb_fits=3,
  rnn_type="lstm",
  num_layers=1,
  hidden_size=128,
  p_drop=0,
  L2_pen=0,
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
path_eqrn <- "Results/ts/EQRN_gridsearch_ts/"#save_path
par_eqrn_p <- list(
  nb_fits = 3,
  shape_fixed=TRUE,
  rnn_type="lstm",
  num_layers=1,
  hidden_size=128,
  p_drop=0,
  intermediate_q_feature=TRUE,
  L2_pen=1e-4,
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
  tol=1e-5
)
par_eqrn_np <- list(
  nb_fits = 3,
  shape_fixed=TRUE,
  rnn_type="lstm",
  num_layers=2,
  hidden_size=128,
  p_drop=0,
  intermediate_q_feature=TRUE,
  L2_pen=0,
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
  tol=1e-5
)

hid_f_str <- "tanh"
lr_str <- "lrd4"

# Params: GRF
grf_pars <- list(
  timedep = par_qrn$seq_len,
  num.trees = 5e3,
  quantiles_fit = c(0.1, 0.5, 0.9),
  sample.fraction = 0.5,
  min.node.size = 5,
  honesty = TRUE,
  honesty.fraction = 0.5,
  honesty.prune.leaves = TRUE,
  alpha = 0.05,
  imbalance.penalty = 0
)

# Params: GBEX
interm_method_competitors <- intermediate_method
gbex_params <- list(
  B=556,
  lambda=NULL,
  lambda_ratio=6,
  lambda_scale=0.01,
  depth=c(1, 0),
  sf=0.75,
  intermediate_q_feature=TRUE,
  scale_features=FALSE)

# Params: EGAM
egam_params <- list(
  model_shape=FALSE,
  intermediate_q_feature=gbex_params$intermediate_q_feature,
  scale_features=gbex_params$scale_features)

# Plots
mp_lwdt <- 150
mp_lhgt <- 30
mp_wdth <- 60
ap_wdth <- 210
ap_hght <- 82


## =============================== END PARAMETERS ===============================

qrn_params_string <- paste0("qrrnn_", par_qrn$rnn_type, "_", par_qrn$num_layers, "x", par_qrn$hidden_size, "_s", par_qrn$seq_len, "_do", par_qrn$p_drop*100,
                            "_L2", str_replace(toString(par_qrn$L2_pen),"([.])","p"), "_lr", str_replace(toString(par_qrn$learning_rate),"([.])","p"))

params_string_p <- paste0("reqrnn_", par_eqrn_p$rnn_type, "_", par_eqrn_p$num_layers, "x", par_eqrn_p$hidden_size, "_",
                          "n"[!par_eqrn_p$intermediate_q_feature], "u_s", par_eqrn_p$seq_len, "_do", par_eqrn_p$p_drop*100,
                          "_L2", str_replace(toString(par_eqrn_p$L2_pen),"([.])","p"), "_lr", str_replace(toString(par_eqrn_p$learning_rate),"([.])","p"))
params_string_np <- paste0("reqrnn_", par_eqrn_np$rnn_type, "_", par_eqrn_np$num_layers, "x", par_eqrn_np$hidden_size, "_",
                           "n"[!par_eqrn_np$intermediate_q_feature], "u_s", par_eqrn_np$seq_len, "_do", par_eqrn_np$p_drop*100,
                           "_L2", str_replace(toString(par_eqrn_np$L2_pen),"([.])","p"), "_lr", str_replace(toString(par_eqrn_np$learning_rate),"([.])","p"))


start_time <- Sys.time()

# Data
if(!file.exists(paste0(interm_path,"Data_backup.rds"))){
  warning("File not found, new data is generated.")
  dat <- generate_series_model(n=n+n_valid+ntest, df=df, AR=AR, MA=MA, muX=muX, mu0=mu0,
                               alphas=alphas, betas=betas, sX=sX, S0=S0, ARX=ARX,
                               X_distr=X_distr, Y_distr=Y_distr, seasonal_hetero=seasonal_hetero)
  
  data_save <- list(dat=dat, n=n, n_valid=n_valid, ntest=ntest, df=df, AR=AR, MA=MA, muX=muX, mu0=mu0,
                    alphas=alphas, betas=betas, sX=sX, S0=S0, ARX=ARX, X_distr=X_distr, Y_distr=Y_distr, seasonal_hetero=seasonal_hetero)
  #safe_save_rds(data_save, paste0(interm_path,"Data_backup.rds"))
}else{
  data_save <- readRDS(paste0(interm_path,"Data_backup.rds"))
  dat <- data_save$dat
  verif_dat <- c(n==data_save$n, n_valid==data_save$n_valid, ntest==data_save$ntest, df==data_save$df,
                 AR==data_save$AR, MA==data_save$MA, muX==data_save$muX, mu0==data_save$mu0,
                 alphas==data_save$alphas, betas==data_save$betas, sX==data_save$sX, S0==data_save$S0, ARX==data_save$ARX,
                 X_distr==data_save$X_distr, Y_distr==data_save$Y_distr)
  if(any(!verif_dat)){stop("Issue with data generating parameters")}
}
rm(data_save)

n_train_all <- n+n_valid
true_quantiles <- series_theoretical_quantiles(quantiles_predict, dat, Y_distr=Y_distr)
seq_len <- par_qrn$seq_len
X_train <- dat$X[1:n, , drop=F]
y_train <- dat$Y[1:n]
true_quantiles_train <- true_quantiles[1:n, , drop=F]
X_valid <- dat$X[(n+1-seq_len):(n_train_all), , drop=F]
y_valid <- dat$Y[(n+1-seq_len):(n_train_all)]
true_quantiles_valid <- true_quantiles[(n+1):(n_train_all), , drop=F]
X_test <- dat$X[(n_train_all+1-seq_len):(n_train_all+ntest), , drop=F]
y_test <- dat$Y[(n_train_all+1-seq_len):(n_train_all+ntest)]
true_quantiles_test <- true_quantiles[(n_train_all+1):(n_train_all+ntest), , drop=F]

X_train_all <- dat$X[1:(n_train_all), , drop=F]
y_train_all <- dat$Y[1:(n_train_all)]
true_quantiles_train_all <- true_quantiles[1:(n_train_all), , drop=F]

plt_dat <- plot_series_comparison2(y_train[1:1000], X_train[1:1000, ], var_names=c("Y","X"))
plot(plt_dat)
if(!is.null(save_path)){
  save_myplot(plt=plt_dat, plt_nm=paste0(save_path, "simulated_ts_data.pdf"),
              width = mp_lwdt, height = mp_lhgt, cairo=FALSE)
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
  #Intermediate quantiles all data
  thresh_quant_all <- series_theoretical_quantiles(interm_lvl, dat, Y_distr=Y_distr)
}
intermediate_quantiles <- thresh_quant_all[1:n, , drop=F]
valid_quantiles <- thresh_quant_all[(n+1-par_qrn$seq_len):(n_train_all), , drop=F]
interm_quantiles_all <- thresh_quant_all[1:(n_train_all), , drop=F]
#Predict intermediate quantiles on X_test
pred_interm <- thresh_quant_all[(n_train_all+1-par_qrn$seq_len):(n_train_all+ntest), , drop=F]


## ======== EQRN FIT ========

#Fit EQRN with intermediate (or oracle) quantiles
if((!dir.exists(paste0(path_eqrn, "networks/", params_string_p,"/"))) | (!dir.exists(paste0(path_eqrn, "networks/", params_string_np,"/"))) | force_refit){
  stop("no refit in 'EQRN_comparison_analysis.R'")
}else{
  fit_eqrn_p <- EQRN_load(paste0(path_eqrn, "networks/"), params_string_p)
  fit_eqrn_np <- EQRN_load(paste0(path_eqrn, "networks/"), params_string_np)
}


cat("\nElapsed time (EQRN fit):\n")
print(Sys.time() - start_time)


## ======== TESTING ========

#Final EQRN predictions on X_test
pred_eqrn_p_val <- EQRN_predict_seq(fit_eqrn_p, X_valid, y_valid, quantiles_predict, valid_quantiles, interm_lvl, crop_predictions=TRUE)
pred_eqrn_p <- EQRN_predict_seq(fit_eqrn_p, X_test, y_test, quantiles_predict, pred_interm, interm_lvl, crop_predictions=TRUE)
pred_eqrn_np_val <- EQRN_predict_seq(fit_eqrn_np, X_valid, y_valid, quantiles_predict, valid_quantiles, interm_lvl, crop_predictions=TRUE)
pred_eqrn_np <- EQRN_predict_seq(fit_eqrn_np, X_test, y_test, quantiles_predict, pred_interm, interm_lvl, crop_predictions=TRUE)


# UNCONDITIONAL predicted quantile(s) (Y quantile on X_train)
pred_unc <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = quantiles_predict, Y = y_train_all, ntest = ntest)

# SEMI-CONDITIONAL predicted quantiles
pred_semicond <- predict_GPD_semiconditional(Y=y_train[(seq_len+1):n], interm_lvl=interm_lvl, thresh_quantiles=intermediate_quantiles[(seq_len+1):n],
                                             interm_quantiles_test=pred_interm[seq_len+(1:ntest)], quantiles_predict=quantiles_predict)

# GROUND-TRUTH (y_test)
pred_true <- true_quantiles_test

# Unconditional parameters and losses (for comparison with EQRN fit)
uncond_losses_fixed <- unconditional_train_valid_GPD_loss(Y_train=y_train[(seq_len+1):n], interm_lvl=interm_lvl, Y_valid=y_valid[seq_len+(1:n_valid)])
uncond_losses_interm <- semiconditional_train_valid_GPD_loss(Y_train=y_train[(seq_len+1):n], Y_valid=y_valid[seq_len+(1:n_valid)],
                                                             interm_quant_train=intermediate_quantiles[(seq_len+1):n],
                                                             interm_quant_valid=valid_quantiles[seq_len+(1:n_valid)])

## ======== GRF GBEX and EGAM FITS =========

lagged_train_all <- lagged_features(X=cbind(y_train_all,X_train_all), max_lag=grf_pars$timedep, drop_present=TRUE)
lagged_test <- lagged_features(X=cbind(y_test,X_test), max_lag=grf_pars$timedep, drop_present=TRUE)

fit_grf <- quantile_forest(lagged_train_all, y_train_all[(grf_pars$timedep+1):length(y_train_all)],
                           num.trees=grf_pars$num.trees, quantiles=grf_pars$quantiles_fit, sample.fraction=grf_pars$sample.fraction,
                           min.node.size=grf_pars$min.node.size, honesty=grf_pars$honesty,
                           honesty.fraction=grf_pars$honesty.fraction, honesty.prune.leaves=grf_pars$honesty.prune.leaves,
                           alpha=grf_pars$alpha, imbalance.penalty=grf_pars$imbalance.penalty, seed=seedGRF)

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
  #Intermediate quantiles all data
  thresh_quant_all_c <- series_theoretical_quantiles(interm_lvl, dat, Y_distr=Y_distr)
}
interm_quantiles_all_c <- thresh_quant_all_c[1:(n_train_all), , drop=F]
pred_interm_c <- thresh_quant_all_c[(n_train_all+1-par_qrn$seq_len):(n_train_all+ntest), , drop=F]

lagged_interm_q_trall <- interm_quantiles_all_c[(grf_pars$timedep+1):length(interm_quantiles_all_c), , drop=F]
laged_interm_q_test <- pred_interm_c[(grf_pars$timedep+1):length(pred_interm_c), , drop=F]
#High quantile prediction with GRF
pred_grf_test <- predict(fit_grf, newdata=lagged_test, quantiles = quantiles_predict)$predictions

# GBEX and EGAM fits and prediction
fit_gbex <- gbex_fit(X=lagged_train_all, y=y_train_all[(grf_pars$timedep+1):length(y_train_all)], intermediate_quantiles=lagged_interm_q_trall,
                     interm_lvl=interm_lvl, intermediate_q_feature=gbex_params$intermediate_q_feature, scale_features=gbex_params$scale_features,
                     B=gbex_params$B, lambda=gbex_params$lambda, lambda_ratio=gbex_params$lambda_ratio,
                     lambda_scale=gbex_params$lambda_scale, depth=gbex_params$depth, sf=gbex_params$sf)
pred_gbex <- gbex_predict(fit_gbex, lagged_test, to_predict=quantiles_predict, intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
fit_egam <- fit_gpd_gam(X=lagged_train_all, y=y_train_all[(grf_pars$timedep+1):length(y_train_all)],
                        intermediate_quantiles=lagged_interm_q_trall, interm_lvl=interm_lvl, model_shape=egam_params$model_shape,
                        intermediate_q_feature=egam_params$intermediate_q_feature, scale_features=egam_params$scale_features)
pred_egam <- predict_gpd_gam(fit_egam, lagged_test, to_predict=quantiles_predict,
                             intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
## ===========================


train_plot_p <- training_plot_eqrn(fit_eqrn_p, uncond_losses_interm, show_legend=TRUE)
train_plot_np <- training_plot_eqrn(fit_eqrn_np, uncond_losses_interm, show_legend=TRUE)
valid_plot_p <- validation_plot_eqrn(fit_eqrn_p, uncond_losses_interm, show_legend=FALSE)
valid_plot_np <- validation_plot_eqrn(fit_eqrn_np, uncond_losses_interm, show_legend=FALSE)
plot(train_plot_p)
plot(train_plot_np)
plot(valid_plot_p)
plot(valid_plot_np)
if(save_plots){
  save_myplot(plt=train_plot_p, plt_nm=paste0(save_path, params_string_p, "_training_p.pdf"),
              width = mp_wdth, height = mp_wdth, cairo=FALSE)
  save_myplot(plt=train_plot_np, plt_nm=paste0(save_path, params_string_np, "_training_np.pdf"),
              width = mp_wdth, height = mp_wdth, cairo=FALSE)
  save_myplot(plt=valid_plot_p, plt_nm=paste0(save_path, params_string_p, "_validation_p.pdf"),
              width = mp_wdth, height = mp_wdth, cairo=FALSE)
  save_myplot(plt=valid_plot_np, plt_nm=paste0(save_path, params_string_np, "_validation_np.pdf"),
              width = mp_wdth, height = mp_wdth, cairo=FALSE)
}

# Compute losses for desired predicted quantiles
nb_quantiles_predict <- length(quantiles_predict)

quant_pred_plot <- c(0.8, 0.85, 0.9, 0.95, 0.99, 0.995, 0.999, 0.9995)
quant_pred_plot_smooth <- c(0.8,1-10^c(-1,-1.5,-2,-2.5,-3,-3.5,-4))
for (Q_pred_levels in list(list(qs=quant_pred_plot, nm="dq"),
                           list(qs=quant_pred_plot_smooth, nm="sq"))){
  
  true_q_test_plot <- series_theoretical_quantiles(Q_pred_levels$qs, dat, Y_distr=Y_distr)[(n_train_all+1):(n_train_all+ntest), , drop=F]
  err_quant_plt_p <- plot_rmse_quantile_ts(fit_eqrn_p, fit_grf, fit_gbex, fit_egam, y_train, X_test, y_test, Q_pred_levels$qs, pred_interm,
                                           interm_lvl, intermediate_quantiles, true_q_test_plot, factorize=TRUE, test_data="other",
                                           pred_interm_c=pred_interm_c, bias_variance=TRUE, crop_Rbv=c(0,0,0))
  err_quant_plt_np <- plot_rmse_quantile_ts(fit_eqrn_np, fit_grf, fit_gbex, fit_egam, y_train, X_test, y_test, Q_pred_levels$qs, pred_interm,
                                            interm_lvl, intermediate_quantiles, true_q_test_plot, factorize=TRUE, test_data="other",
                                            pred_interm_c=pred_interm_c, bias_variance=TRUE, crop_Rbv=c(0,0,0))
  err_quant_comp_plt <- plot_rmse_quantile_comp_ts_old(fit_eqrn_p, fit_eqrn_np, fit_grf, fit_gbex, fit_egam, y_train, X_test, y_test, Q_pred_levels$qs, pred_interm,
                                                       interm_lvl, intermediate_quantiles, true_q_test_plot, factorize=TRUE, test_data="other",
                                                       pred_interm_c=pred_interm_c, names_EQRN=c("EQRN best","EQRN unpen."), bias_variance=TRUE,
                                                       crop_Rbv=c(0,0,0))
  plot(err_quant_plt_p$rmse_plot)
  plot(err_quant_plt_np$rmse_plot)
  plot(err_quant_comp_plt$rmse_plot)
  if(save_plots){
    save_myplot(plt=err_quant_plt_p$rmse_plot, plt_nm=paste0(save_path, params_string_p, "_errquant_p_",Q_pred_levels$nm,".pdf"),
                width = mp_wdth, height = mp_wdth, cairo=FALSE)
    save_myplot(plt=err_quant_plt_np$rmse_plot, plt_nm=paste0(save_path, params_string_np, "_errquant_np_",Q_pred_levels$nm,".pdf"),
                width = mp_wdth, height = mp_wdth, cairo=FALSE)
    save_myplot(plt=err_quant_comp_plt$rmse_plot, plt_nm=paste0(save_path, "compar_errquant_",Q_pred_levels$nm,".pdf"),
                width = mp_wdth, height = mp_wdth, cairo=FALSE)
    save_myplot(plt=err_quant_plt_p$Rbv_plot, plt_nm=paste0(save_path, params_string_p, "_Rbvquant_p_",Q_pred_levels$nm,".pdf"),
                width = ap_wdth, height = ap_hght, cairo=FALSE)
    save_myplot(plt=err_quant_plt_np$Rbv_plot, plt_nm=paste0(save_path, params_string_np, "_Rbvquant_np_",Q_pred_levels$nm,".pdf"),
                width = ap_wdth, height = ap_hght, cairo=FALSE)
    save_myplot(plt=err_quant_comp_plt$Rbv_plot, plt_nm=paste0(save_path, "compar_Rbvquant_",Q_pred_levels$nm,".pdf"),
                width = ap_wdth, height = ap_hght, cairo=FALSE)
  }
}

for(i in 1:nb_quantiles_predict){
  pred_vs_true_comp_sc <- plot_predictions_vs_truth_comp(pred_eqrn_np[,i], pred_eqrn_p[,i], pred_semicond$predictions[,i],
                                                         pred_true[,i], names_EQRN=c("EQRN unpen.","EQRN best"),
                                                         name_other="Semi_cond", legend.position="none")#"none""bottom"
  pred_vs_true_comp_solo <- plot_predictions_vs_truth_comp(pred_eqrn_np[,i], pred_eqrn_p[,i], NaN,
                                                           pred_true[,i], names_EQRN=c("EQRN unpen.","EQRN best"),
                                                           name_other="null", legend.position="none")#"none""bottom"
  plot(pred_vs_true_comp_sc)
  plot(pred_vs_true_comp_solo)
  if(save_plots){
    save_myplot(plt=pred_vs_true_comp_sc, plt_nm=paste0(save_path, "predVStruth_comp_sc_q",quantiles_predict[i]*10000,".pdf"),
                width = mp_wdth, height = mp_wdth, cairo=FALSE)
    save_myplot(plt=pred_vs_true_comp_solo, plt_nm=paste0(save_path, "predVStruth_comp_solo_q",quantiles_predict[i]*10000,".pdf"),
                width = mp_wdth, height = mp_wdth, cairo=FALSE)
  }
}


end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)

