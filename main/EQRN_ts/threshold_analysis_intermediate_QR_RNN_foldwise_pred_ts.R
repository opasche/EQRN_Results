library(tidyverse)
library(here)
library(ggpubr)
library(grf)
source("R/EQRN_loader.R")

## =============================== PARAMETERS ===============================
save_path <- "data/simulations/ts_intermediate_quantile_threshold_analysis/"
results_path <- "Results/ts/intermediate_quantiles_threshold_analysis/"
check_directory(save_path, recursive=TRUE)
check_directory(results_path, recursive=TRUE)
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
interm_lvl = 0.8 #0.7 0.75 0.8 0.85 0.9 0.95
n_folds <- 4

# Params: QRNN
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

## =============================== END PARAMETERS ===============================

params_string <- paste0("qrrnn_", par_qrn$rnn_type, "_", par_qrn$num_layers, "x", par_qrn$hidden_size, "_s", par_qrn$seq_len, "_do", par_qrn$p_drop*100,
                        "_L2", str_replace(toString(par_qrn$L2_pen),"([.])","p"), "_lr", str_replace(toString(par_qrn$learning_rate),"([.])","p"),
                        "_ilv",interm_lvl*100)

start_time <- Sys.time()

# Data

if(!file.exists(paste0(save_path,"Data_backup.rds"))){
  stop("Data file 'Data_backup.rds' not found.")
  warning("Data file 'Data_backup.rds' not found, new data is generated.")
  dat <- generate_series_model(n=n+n_valid+ntest, df=df, AR=AR, MA=MA, muX=muX, mu0=mu0,
                               alphas=alphas, betas=betas, sX=sX, S0=S0, ARX=ARX,
                               X_distr=X_distr, Y_distr=Y_distr, seasonal_hetero=seasonal_hetero)
  
  data_save <- list(dat=dat, n=n, n_valid=n_valid, ntest=ntest, df=df, AR=AR, MA=MA, muX=muX, mu0=mu0,
                    alphas=alphas, betas=betas, sX=sX, S0=S0, ARX=ARX, X_distr=X_distr, Y_distr=Y_distr, seasonal_hetero=seasonal_hetero)
  safe_save_rds(data_save, paste0(save_path, "Data_","backup",".rds"))
}else{
  data_save <- readRDS(paste0(save_path,"Data_backup.rds"))
  dat <- data_save$dat
  verif_dat <- c(n==data_save$n, n_valid==data_save$n_valid, ntest==data_save$ntest, df==data_save$df,
                 AR==data_save$AR, MA==data_save$MA, muX==data_save$muX, mu0==data_save$mu0,
                 alphas==data_save$alphas, betas==data_save$betas, sX==data_save$sX, S0==data_save$S0, ARX==data_save$ARX,
                 X_distr==data_save$X_distr, Y_distr==data_save$Y_distr, seasonal_hetero==data_save$seasonal_hetero)
  if(any(!verif_dat)){stop("Issue with data generating parameters")}
}
rm(data_save)

true_quantiles <- series_theoretical_quantiles(interm_lvl, dat, Y_distr=Y_distr)
seq_len <- par_qrn$seq_len
X_train <- dat$X[1:n, , drop=F]
y_train <- dat$Y[1:n]
true_quantiles_train <- true_quantiles[1:n]
X_valid <- dat$X[(n+1-seq_len):(n+n_valid), , drop=F]
y_valid <- dat$Y[(n+1-seq_len):(n+n_valid)]
true_quantiles_valid <- true_quantiles[(n+1):(n+n_valid)]
X_test <- dat$X[(n+n_valid+1-seq_len):(n+n_valid+ntest), , drop=F]
y_test <- dat$Y[(n+n_valid+1-seq_len):(n+n_valid+ntest)]
true_quantiles_test <- true_quantiles[(n+n_valid+1):(n+n_valid+ntest)]

#plot(dat$Y)
#lines(true_quantiles)
X_train_all <- dat$X[1:(n+n_valid), , drop=F]
y_train_all <- dat$Y[1:(n+n_valid)]
true_quantiles_train_all <- true_quantiles[1:(n+n_valid)]


if((!file.exists(paste0(save_path, "Results_", params_string, ".rds"))) | force_refit){
  cat("Save files not found, new networks are fitted.\n")
  # QRNN fit
  torch_manual_seed(seedT)
  foldwise_obj <- QRN_seq_predict_foldwise(X_train_all, y_train_all, q_level=interm_lvl, n_folds=n_folds, number_fits=par_qrn$nb_fits, seq_len=par_qrn$seq_len,
                                           hidden_size=par_qrn$hidden_size, num_layers=par_qrn$num_layers,
                                           rnn_type=par_qrn$rnn_type, p_drop=par_qrn$p_drop, learning_rate=par_qrn$learning_rate, L2_pen=par_qrn$L2_pen,
                                           scale_features=par_qrn$scale_features, n_epochs=par_qrn$n_epochs,
                                           batch_size=par_qrn$batch_size, lr_decay=par_qrn$lr_decay,
                                           patience_decay=par_qrn$patience_decay, min_lr=par_qrn$min_lr, patience_stop=par_qrn$patience_stop, tol=par_qrn$tol)
  
  cat("\nElapsed time (QRN foldwise fit predictions):\n")
  print(Sys.time() - start_time)
  
  for (k in 1:n_folds) {
    EQRN_save(foldwise_obj$fits[[k]], paste0(save_path, "networks/"), paste0(params_string, "_f", k))
  }
  safe_save_rds(foldwise_obj[names(foldwise_obj)!="fits"], paste0(save_path, "Results_", params_string, ".rds"))#checkpoint
  
  # Diagnostics
  
  i_first_pred <- seq_len+1
  n_trall <- n + n_valid
  SEs_valid_all <- mean_squared_error(true_quantiles_train_all[i_first_pred:n_trall],
                                      foldwise_obj$predictions[i_first_pred:n_trall], return_agg = "vector")
  MSE_valid_all <- mean(SEs_valid_all)
  AEs_valid_all <- mean_absolute_error(true_quantiles_train_all[i_first_pred:n_trall],
                                       foldwise_obj$predictions[i_first_pred:n_trall], return_agg = "vector")
  MAE_valid_all <- mean(AEs_valid_all)
  qlosses_valid_all <- quantile_loss(y=y_train_all[i_first_pred:n_trall],
                                     y_hat=foldwise_obj$predictions[i_first_pred:n_trall], q=interm_lvl, return_agg = "vector")
  qloss_valid_all <- mean(qlosses_valid_all)
  
  qrn_pred_test <- matrix(as.double(NA), nrow=ntest, ncol=n_folds)
  fold_MSE_valid <- rep(as.double(NA), n_folds)
  fold_MAE_valid <- rep(as.double(NA), n_folds)
  fold_qloss_valid <- rep(as.double(NA), n_folds)
  test_MSE <- rep(as.double(NA), n_folds)
  test_MAE <- rep(as.double(NA), n_folds)
  test_qloss <- rep(as.double(NA), n_folds)
  for (k in 1:n_folds) {
    # Train and validation losses
    fold_MSE_valid[k] <- mean_squared_error(true_quantiles_train_all[foldwise_obj$cuts==k],
                                            foldwise_obj$predictions[foldwise_obj$cuts==k],
                                            return_agg = "mean", na.rm=(k==1))
    fold_MAE_valid[k] <- mean_absolute_error(true_quantiles_train_all[foldwise_obj$cuts==k],
                                             foldwise_obj$predictions[foldwise_obj$cuts==k],
                                             return_agg = "mean", na.rm=(k==1))
    fold_qloss_valid[k] <- quantile_loss(y=y_train_all[foldwise_obj$cuts==k],
                                         y_hat=foldwise_obj$predictions[foldwise_obj$cuts==k], q=interm_lvl,
                                         return_agg = "mean", na.rm=(k==1))
    # Test predictions
    qrn_pred_test[,k] <- QRN_seq_predict(foldwise_obj$fits[[k]], X_test, y_test, q_level=interm_lvl, crop_predictions=TRUE)
    test_MSE[k] <- mean_squared_error(true_quantiles_test, qrn_pred_test[,k])
    test_MAE[k] <- mean_absolute_error(true_quantiles_test, qrn_pred_test[,k])
    test_qloss[k] <- quantile_loss(y=y_test[seq_len+(1:ntest)], y_hat=qrn_pred_test[,k], q=interm_lvl)
  }
  qrn_mean_pred_test <- rowMeans(qrn_pred_test)
  qrn_pred_test <- cbind(qrn_pred_test,qrn_mean_pred_test)
  mean_test_MSE <- mean_squared_error(true_quantiles_test, qrn_mean_pred_test)
  mean_test_MAE <- mean_absolute_error(true_quantiles_test, qrn_mean_pred_test)
  mean_test_qloss <- quantile_loss(y=y_test[seq_len+(1:ntest)], y_hat=qrn_mean_pred_test, q=interm_lvl)
  
  
  results_tibble <- tibble(fold=paste0("f_",1:n_folds), train_loss=foldwise_obj$train_losses, valid_loss=foldwise_obj$valid_losses,
                           min_valid_loss=foldwise_obj$min_valid_losses, fold_MSE_valid=fold_MSE_valid, fold_MAE_valid=fold_MAE_valid,
                           test_MSE=test_MSE, test_MAE=test_MAE, test_qloss=test_qloss,
                           fold_qloss_valid=fold_qloss_valid, min_valid_e=foldwise_obj$min_valid_e)
  results_tibble <- results_tibble %>% add_row(fold="aggr", train_loss=NA, valid_loss=qloss_valid_all,
                                               min_valid_loss=NA, fold_MSE_valid=MSE_valid_all, fold_MAE_valid=MAE_valid_all,
                                               test_MSE=mean_test_MSE, test_MAE=mean_test_MAE, test_qloss=mean_test_qloss,
                                               fold_qloss_valid=qloss_valid_all, min_valid_e=NA)
  
  filename <- paste0(save_path,"Disnostics_QR_RNN_foldwise_predictions_",params_string,".csv")
  write_csv(results_tibble, file=filename)
  
  foldwise_obj$test_predictions <- qrn_pred_test
  foldwise_obj$results_tibble <- results_tibble
  foldwise_obj$qrn_parameters <- par_qrn
  safe_save_rds(foldwise_obj[names(foldwise_obj)!="fits"], paste0(save_path, "Results_", params_string, ".rds"))
  
} else {
  
  foldwise_obj <- readRDS(paste0(save_path, "Results_", params_string, ".rds"))
  foldwise_obj$fits <- list()
  for (k in 1:n_folds) {
    foldwise_obj$fits[[k]] <- EQRN_load(paste0(save_path, "networks/"), paste0(params_string, "_f", k))
  }
}

## DIAGNOSTIC PLOTS

for (k in 1:n_folds) {
  uncond_losses_interm <- list(train_loss=quantile_loss(y_train_all[foldwise_obj$cuts!=k], quantile(y_train_all[foldwise_obj$cuts!=k], interm_lvl), interm_lvl),
                               valid_loss=quantile_loss(y_train_all[foldwise_obj$cuts==k], quantile(y_train_all[foldwise_obj$cuts!=k], interm_lvl), interm_lvl))
  train_plot <- training_plot_eqrn(foldwise_obj$fits[[k]], uncond_losses_interm, title="QR RNN training", y_lab="Quantile loss")
  plot(train_plot)
  if(save_plots){
    ggsave(paste0(params_string, "_training_f_",k,".pdf"), plot=train_plot, device="pdf",
           path=results_path, width=160, height=120, units="mm", dpi=300)
  }
}

## Train foldwise intermediate quantiles
foldwise_obj$predictions

predplot <- plot_predictions_ts(foldwise_obj$predictions, NaN, true_quantiles_train_all, NaN, interm_lvl)
plot(predplot)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_preds_interm_folds.pdf"), plot=predplot, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}

predplotdiff <- plot_predictions_diff_ts(foldwise_obj$predictions, NaN, true_quantiles_train_all, interm_lvl)
plot(predplotdiff)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_preds_interm_folds_diff.pdf"), plot=predplotdiff, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}

pred_vs_true_plot_interm <- plot_predictions_vs_truth(matrix(foldwise_obj$predictions), matrix(NaN, nrow=(n+n_valid)),
                                                      matrix(rep(quantile(y_train_all, interm_lvl)[[1]], (n+n_valid))),
                                                      matrix(true_quantiles_train_all), interm_lvl, name_other="NULL")
plot(pred_vs_true_plot_interm)
if(save_plots){
  ggsave(paste0(params_string, "_predVStruth_interm_folds.pdf"), plot=pred_vs_true_plot_interm, device="pdf",
         path=results_path, width=100, height=100, units="mm", dpi=300)
}

resid_box2 <- residuals_boxplot_eqrn2(foldwise_obj$predictions, NaN, NaN, true_quantiles_train_all, interm_lvl)
plot(resid_box2)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_residbox_interm_folds.pdf"), plot=resid_box2, device="pdf",
         path=results_path, width=200, height=150, units="mm", dpi=300)
}


## Aggregated test predictions

predplot <- plot_predictions_ts(foldwise_obj$test_predictions[,foldwise_obj$n_folds+1], NaN, true_quantiles_test, NaN, interm_lvl)
plot(predplot)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_preds_interm_test.pdf"), plot=predplot, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}

predplotdiff <- plot_predictions_diff_ts(foldwise_obj$test_predictions[,foldwise_obj$n_folds+1], NaN, true_quantiles_test, interm_lvl)
plot(predplotdiff)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_preds_interm_test_diff.pdf"), plot=predplotdiff, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}

pred_vs_true_plot_interm <- plot_predictions_vs_truth(foldwise_obj$test_predictions[,foldwise_obj$n_folds+1,drop=F], matrix(NaN, nrow=ntest),
                                                      matrix(rep(quantile(y_train_all, interm_lvl)[[1]], ntest)),
                                                      matrix(true_quantiles_test), interm_lvl, name_other="NULL")
plot(pred_vs_true_plot_interm)
if(save_plots){
  ggsave(paste0(params_string, "_predVStruth_interm_test.pdf"), plot=pred_vs_true_plot_interm, device="pdf",
         path=results_path, width=100, height=100, units="mm", dpi=300)
}

resid_box2 <- residuals_boxplot_eqrn2(foldwise_obj$test_predictions[,foldwise_obj$n_folds+1], NaN, NaN, true_quantiles_test, interm_lvl)
plot(resid_box2)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_residbox_interm_test.pdf"), plot=resid_box2, device="pdf",
         path=results_path, width=200, height=150, units="mm", dpi=300)
}


end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)


