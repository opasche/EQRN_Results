library(tidyverse)
library(evd)
library(ismev)
library(here)
library(tsibble)
library(feasts)
library(fable)
library(ggpubr)
library(latex2exp)
library(stlplus)
library(grf)
source("R/EQRN_loader.R")


## =============================== PARAMETERS ===============================
save_path <- "data/Switzerland/qrn_intermediate_quantile_2_256_1e6_new/"
results_path <- "Results/Switzerland/intermediate_quantiles_foldwise_analysis/QR_RNN_foldwise_2_256_1e6/"
check_directory(save_path, recursive=TRUE)
check_directory(results_path, recursive=TRUE)
save_plots <- TRUE
force_refit <- FALSE

seedR <- 0
seedGRF <- 1
seedT <- seedR
set.seed(seedR)

# PARAM: Data
Y_name <- c("station_62")
X_disch_names <- c("station_43")
precip_names <- c("LTB","BRZ","MER","THU","GHS","BEP")
use_stl_remainder <- FALSE
use_mean_precip <- FALSE
prop_test <- 2/3
prop_valid <- 1/4

# PARAM: General
interm_lvl = 0.8
n_folds <- 4

# Params: QRN
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

## =============================== END PARAMETERS ===============================


start_time <- Sys.time()

disch_names <- c(Y_name, X_disch_names)
sgroup <- c(disch_names, precip_names)

## =============================== Data Processing ===============================
if(!file.exists(paste0(save_path,"Data_backup.rds"))){
  warning("File not found, new data is generated.")
  
  #Load Discharge data
  Disch <- read_csv("data/data_wrangled/Discharges_all.csv",guess_max = 10000)
  Disch <- Disch %>% arrange(Date)
  inds <- Disch[["Date"]]
  #Load precipitation data
  Precipitations <- read_csv("data/Data/precip.csv",guess_max = 10000)
  names(Precipitations)[1] <- "Date"
  
  if(use_stl_remainder){
    #Create stl remainder Discharge data (seasonality in trend removal)
    Disch_ts <- ts(select(Disch, -Date), start = c(as.numeric(format(inds[1], "%Y")), as.numeric(format(inds[1], "%j"))), frequency = 365.25)
    Disch_stl_remainder <- Disch["Date"]
    for (stati in names(select(Disch, -Date))){
      Disch_ts_s <- Disch_ts[,stati]
      stl_s <- stlplus(Disch_ts_s, s.window="periodic", t.window=1e3*nrow(Disch))
      rem <- tibble(Disch_stl_remainder["Date"], stl_s$data$remainder)
      names(rem)[2] <- stati
      Disch_stl_remainder <- full_join(Disch_stl_remainder, rem, by="Date")
    }
    #Merge stl-remainder Discharge and Precipitation data
    Disch_Precip <- full_join(Disch_stl_remainder, Precipitations, by="Date") %>% arrange(Date)#Disch_stl_Precip
  }else{
    #Merge Discharge and Precipitation data
    Disch_Precip <- full_join(Disch, Precipitations, by="Date") %>% arrange(Date)
  }
  
  ## This is the meta-info for all stations, some scraped from the admin DAFU website
  #info_disch <- read_csv("data/data_wrangled/Discharges_info_merged.csv",guess_max = 10000)
  
  incomplete_stations_all <- c()
  for(stati in colnames(Precipitations)[-1]){
    droped_stati <- Disch_Precip[c("Date",stati)] %>% drop_na()
    daterange_stati <- seq(droped_stati[[1,"Date"]], (droped_stati %>% tail(1))[["Date"]], by = "day")
    if(nrow(droped_stati)!=length(daterange_stati)){
      incomplete_stations_all <- c(incomplete_stations_all,stati)
    }
  }
  #print(incomplete_stations_all)
  
  if(use_mean_precip){
    if(all(is.element(incomplete_stations_all, sgroup))){
      stop("All incomplete precipitation stations in selection.")
    }
    dat_tb <- Disch_Precip %>%
      mutate(precip_mean=rowMeans(select(., all_of(precip_names)), na.rm=TRUE)) %>%
      select(Date, all_of(disch_names), precip_mean) %>% drop_na()
  }else{
    if(any(is.element(incomplete_stations_all, sgroup))){
      stop("Incomplete precipitation stations in selection.")
    }
    dat_tb <- Disch_Precip[c("Date",setdiff(sgroup,incomplete_stations_all))] %>% drop_na()
  }
  dat <- as.matrix(dat_tb %>% select(-Date))
  X <- dat[, -1, drop=F]
  Y <- dat[, 1]
  Dates <- dat_tb[["Date"]]
  
  n_tot <- nrow(X)
  n_test <- floor(prop_test*n_tot)
  n_train_all <- n_tot-n_test
  n_valid <- floor(prop_valid*n_train_all)
  n_train <- n_train_all-n_valid
  
  seq_len <- par_qrn$seq_len
  
  data_save <- list(X=X, Y=Y, Dates=Dates, n_train=n_train, n_valid=n_valid, n_test=n_test, seq_len=seq_len)
  safe_save_rds(data_save, paste0(save_path, "Data_","backup",".rds"))
  
  suppressWarnings(rm(Disch, Precipitations, Disch_stl_remainder, Disch_Precip, dat_tb, dat))
}else{
  data_save <- readRDS(paste0(save_path,"Data_backup.rds"))
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

## =============================== QRN ===============================

params_string <- paste0("qrrnn_", par_qrn$rnn_type, "_", par_qrn$num_layers, "x", par_qrn$hidden_size, "_s", par_qrn$seq_len, "_do", par_qrn$p_drop*100,
                        "_L2", str_replace(toString(par_qrn$L2_pen),"([.])","p"), "_lr", str_replace(toString(par_qrn$learning_rate),"([.])","p"))

# Unconditional parameters and losses (for comparison with EQRNN fit)
uncond_losses_interm <- list(train_loss=quantile_loss(y_train, quantile(y_train, interm_lvl), interm_lvl),
                             valid_loss=quantile_loss(y_valid, quantile(y_train, interm_lvl), interm_lvl))


if((!file.exists(paste0(save_path, "Results_", params_string, ".rds"))) | force_refit){
  # QRN fit
  cat("Fitting foldwise model")
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
  safe_save_rds(foldwise_obj[names(foldwise_obj)!="fits"], paste0(save_path, "Results_", params_string, ".rds"))
  
  # Diagnostics
  
  i_first_pred <- seq_len+1
  qlosses_valid_all <- quantile_loss(y=y_train_all[i_first_pred:n_train_all],
                                     y_hat=foldwise_obj$predictions[i_first_pred:n_train_all], q=interm_lvl, return_agg = "vector")
  qloss_valid_all <- mean(qlosses_valid_all)
  
  qrn_pred_test <- matrix(as.double(NA), nrow=n_test, ncol=n_folds)
  fold_qloss_valid <- rep(as.double(NA), n_folds)
  test_qloss <- rep(as.double(NA), n_folds)
  for (k in 1:n_folds) {
    # Train and validation losses
    fold_qloss_valid[k] <- quantile_loss(y=y_train_all[foldwise_obj$cuts==k],
                                         y_hat=foldwise_obj$predictions[foldwise_obj$cuts==k], q=interm_lvl,
                                         return_agg = "mean", na.rm=(k==1))
    # Test predictions
    qrn_pred_test[,k] <- QRN_seq_predict(foldwise_obj$fits[[k]], X_test, y_test, q_level=interm_lvl, crop_predictions=TRUE)
    test_qloss[k] <- quantile_loss(y=y_test[seq_len+(1:n_test)], y_hat=qrn_pred_test[,k], q=interm_lvl)
  }
  qrn_mean_pred_test <- rowMeans(qrn_pred_test)
  qrn_pred_test <- cbind(qrn_pred_test,qrn_mean_pred_test)
  mean_test_qloss <- quantile_loss(y=y_test[seq_len+(1:n_test)], y_hat=qrn_mean_pred_test, q=interm_lvl)
  
  
  results_tibble <- tibble(fold=paste0("f_",1:n_folds), train_loss=foldwise_obj$train_losses, valid_loss=foldwise_obj$valid_losses,
                           min_valid_loss=foldwise_obj$min_valid_losses,
                           test_qloss=test_qloss,
                           fold_qloss_valid=fold_qloss_valid, min_valid_e=foldwise_obj$min_valid_e)
  results_tibble <- results_tibble %>% add_row(fold="aggr", train_loss=NA, valid_loss=qloss_valid_all,
                                               min_valid_loss=NA,
                                               test_qloss=mean_test_qloss,
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

# Other models
lagged_train_all <- lagged_features(X=cbind(y_train_all,X_train_all), max_lag=grf_pars$timedep, drop_present=TRUE)
lagged_test <- lagged_features(X=cbind(y_test,X_test), max_lag=grf_pars$timedep, drop_present=TRUE)
fit_grf <- quantile_forest(lagged_train_all, y_train_all[(grf_pars$timedep+1):length(y_train_all)],
                           num.trees=grf_pars$num.trees, quantiles=grf_pars$quantiles_fit, sample.fraction=grf_pars$sample.fraction,# mtry=grf_pars$mtry,
                           min.node.size=grf_pars$min.node.size, honesty=grf_pars$honesty,
                           honesty.fraction=grf_pars$honesty.fraction, honesty.prune.leaves=grf_pars$honesty.prune.leaves,
                           alpha=grf_pars$alpha, imbalance.penalty=grf_pars$imbalance.penalty, seed=seedGRF)#,compute.oob.predictions = FALSE
pred_grf_train_all <- predict(fit_grf, newdata=NULL, quantiles = c(interm_lvl))$predictions
pred_grf_test <- predict(fit_grf, newdata=lagged_test, quantiles = c(interm_lvl))$predictions

# plots test
qrn_pred_test <- foldwise_obj$test_predictions[, (n_folds+1)]
predplot <- plot_predictions_ts(qrn_pred_test, pred_grf_test, rep(NaN,n_test), y_test[par_qrn$seq_len+(1:n_test)], interm_lvl)
plot(predplot)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_preds_interm.pdf"), plot=predplot, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}

predplotdiff <- plot_predictions_diff_ts(qrn_pred_test, pred_grf_test, y_test[par_qrn$seq_len+(1:n_test)], interm_lvl)
plot(predplotdiff)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_preds_interm_diff_obs.pdf"), plot=predplotdiff, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}

predplotzoom <- plot_predictions_ts(qrn_pred_test[1:5000], pred_grf_test[1:5000], rep(NaN,5000), y_test[par_qrn$seq_len+(1:5000)], interm_lvl)
plot(predplotzoom)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_preds_interm_zoom.pdf"), plot=predplotzoom, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}

predplotdiffzoom <- plot_predictions_diff_ts(qrn_pred_test[1:5000], pred_grf_test[1:5000], y_test[par_qrn$seq_len+(1:5000)], interm_lvl)
plot(predplotdiffzoom)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_preds_interm_diff_obs_zoom.pdf"), plot=predplotdiffzoom, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}

#plots folds
fpreds <- foldwise_obj$predictions[(par_qrn$seq_len+1):n_train_all]
fpredplot <- plot_predictions_ts(fpreds, pred_grf_train_all, rep(NaN,n_train_all-par_qrn$seq_len), y_train_all[(par_qrn$seq_len+1):n_train_all], interm_lvl)
plot(fpredplot)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_fpreds_interm.pdf"), plot=fpredplot, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}

fpredplotdiff <- plot_predictions_diff_ts(fpreds, pred_grf_train_all, y_train_all[(par_qrn$seq_len+1):n_train_all], interm_lvl)
plot(fpredplotdiff)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_fpreds_interm_diff_obs.pdf"), plot=fpredplotdiff, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}

fpredplotzoom <- plot_predictions_ts(fpreds[1:5000], pred_grf_train_all[1:5000], rep(NaN,5000), y_train_all[par_qrn$seq_len+(1:5000)], interm_lvl)
plot(fpredplotzoom)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_fpreds_interm_zoom.pdf"), plot=fpredplotzoom, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}

fpredplotdiffzoom <- plot_predictions_diff_ts(fpreds[1:5000], pred_grf_train_all[1:5000], y_train_all[par_qrn$seq_len+(1:5000)], interm_lvl)
plot(fpredplotdiffzoom)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_fpreds_interm_diff_obs_zoom.pdf"), plot=fpredplotdiffzoom, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}

f_ex_ind <- which(y_train_all[(par_qrn$seq_len+1):n_train_all]>fpreds)
fexcesses <- pmax(0,y_train_all[(par_qrn$seq_len+1):n_train_all]-fpreds)[f_ex_ind]

t_ex_ind <- which(y_test[par_qrn$seq_len+(1:n_test)]>qrn_pred_test)
texcesses <- pmax(0,y_test[par_qrn$seq_len+(1:n_test)]-qrn_pred_test)[t_ex_ind]

fex_folds <- data.frame(i=f_ex_ind, excesses=fexcesses) %>% ggplot(aes(x=i, y=excesses)) +
  geom_point() + labs(y="Conditional excesses", x="Observation in the training set", title=NULL)
plot(fex_folds)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_fpreds_interm_excesses.pdf"), plot=fex_folds, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}
fex_folds_zoom <- data.frame(i=f_ex_ind[1:1000], excesses=fexcesses[1:1000]) %>% ggplot(aes(x=i, y=excesses)) +
  geom_point() + labs(y="Conditional excesses", x="Observation in the training set", title=NULL)
plot(fex_folds_zoom)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_fpreds_interm_excesses_zoom.pdf"), plot=fex_folds_zoom, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}

fex_folds <- data.frame(i=t_ex_ind, excesses=texcesses) %>% ggplot(aes(x=i, y=excesses)) +
  geom_point() + labs(y="Conditional excesses", x="Observation in the test set", title=NULL)
plot(fex_folds)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_preds_interm_excesses.pdf"), plot=fex_folds, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}
fex_folds_zoom <- data.frame(i=t_ex_ind[1:1000], excesses=texcesses[1:1000]) %>% ggplot(aes(x=i, y=excesses)) +
  geom_point() + labs(y="Conditional excesses", x="Observation in the test set", title=NULL)
plot(fex_folds_zoom)
if(save_plots){
  ggsave(paste0(params_string, "_q",interm_lvl*10000,"_preds_interm_excesses_zoom.pdf"), plot=fex_folds_zoom, device="pdf",
         path=results_path, width=300, height=200, units="mm", dpi=300)
}

sink(paste0(results_path,"model_errors_",params_string,".txt"))
cat("Train set (out of sample):\n")
cat("Quantile loss QRN: ", quantile_loss(y=y_train_all[(par_qrn$seq_len+1):n_train_all], y_hat=fpreds, q=interm_lvl),"\n")
cat("Quantile loss GRF: ", quantile_loss(y=y_train_all[(par_qrn$seq_len+1):n_train_all], y_hat=pred_grf_train_all, q=interm_lvl),"\n")
cat("Obs Below QRN: ", sum(y_train_all[(par_qrn$seq_len+1):n_train_all]<=fpreds)/n_train_all,"\n")
cat("Obs Below GRF: ", sum(y_train_all[(par_qrn$seq_len+1):n_train_all]<=pred_grf_train_all)/n_train_all,"\n")
cat("Q Pred Error QRN: ", quantile_prediction_error(y_train_all[(par_qrn$seq_len+1):n_train_all], fpreds, interm_lvl),"\n")
cat("Q Pred Error GRF: ", quantile_prediction_error(y_train_all[(par_qrn$seq_len+1):n_train_all], pred_grf_train_all, interm_lvl),"\n")

cat("Test set:\n")
cat("Quantile loss QRN: ", quantile_loss(y=y_test[seq_len+(1:n_test)], y_hat=qrn_pred_test, q=interm_lvl),"\n")
cat("Quantile loss GRF: ", quantile_loss(y=y_test[seq_len+(1:n_test)], y_hat=pred_grf_test, q=interm_lvl),"\n")
cat("Obs Below QRN: ", sum(y_test[seq_len+(1:n_test)]<=qrn_pred_test)/n_test,"\n")
cat("Obs Below GRF: ", sum(y_test[seq_len+(1:n_test)]<=pred_grf_test)/n_test,"\n")
cat("Q Pred Error QRN: ", quantile_prediction_error(y_test[seq_len+(1:n_test)], qrn_pred_test, interm_lvl),"\n")
cat("Q Pred Error GRF: ", quantile_prediction_error(y_test[seq_len+(1:n_test)], pred_grf_test, interm_lvl),"\n")
sink()


## Train foldwise intermediate quantiles
foldwise_obj$predictions


end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)




