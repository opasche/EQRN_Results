library(tidyverse)
library(here)
library(ggpubr)
library(grf)
library(evd)
library(ismev)
library(future)
library(doFuture)
source("R/gbex_wrappers.R")
source("R/EQRN_loader.R")

## =============================== PARAMETERS ===============================
save_path <- "Results/competitors/gbex/Switzerland/"
check_directory(save_path, recursive=TRUE)
parallel_strat <- "multisession"#"sequential", "multisession", "multicore"
n_workers <- 12#(availableCores() - 1)#PARAM
err_handling <- "remove"#"stop", "remove", "pass"

seedR <- 0
seedGRF <- 1
seedgbex <- 1
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
intermediate_method <- "qrn"#"qrn""grf""oracle"
interm_lvl = 0.8
prob_lvls_predict = c(0.995,0.999,0.9995)#c(interm_lvl,0.995,0.999,0.9995)

# Params: QRN
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

# Params: gbex
num_folds=5
gbex_params <- list(
  intermediate_q_feature=TRUE,
  scale_features=FALSE,
  Bmax=3000,
  grid_lambda_ratio=c(5,6,7,8,9,10),
  grid_depth=list(c(1,0),c(1,1),c(2,1),c(2,2),c(3,1),c(3,2),c(3,3)),
  stratified=TRUE,
  lambda_scale=0.01,
  min_leaf_size=NULL,
  sf=0.75)

## =============================== END PARAMETERS ===============================

qrn_params_string <- paste0("qrrnn_", par_qrn$rnn_type, "_", par_qrn$num_layers, "x", par_qrn$hidden_size, "_s", par_qrn$seq_len, "_do", par_qrn$p_drop*100,
                            "_L2", str_replace(toString(par_qrn$L2_pen),"([.])","p"), "_lr", str_replace(toString(par_qrn$learning_rate),"([.])","p"))


start_time <- Sys.time()

# Data
if(!file.exists(paste0(interm_path,"Data_backup.rds"))){
  warning("File not found, new data is generated.")
  stop("insert data creation here")
}else{
  data_save <- readRDS(paste0(interm_path,"Data_backup.rds"))
  X <- data_save$X
  Y <- data_save$Y
  Dates <- data_save$Dates
  n_train <- data_save$n_train
  n_valid <- data_save$n_valid
  n_test <- data_save$n_test
  seq_len <- data_save$seq_len
}
rm(data_save)

seq_len <- par_qrn$seq_len
X_train <- X[1:n_train, , drop=F]
y_train <- Y[1:n_train]
X_valid <- X[(n_train+1-seq_len):(n_train+n_valid), , drop=F]
y_valid <- Y[(n_train+1-seq_len):(n_train+n_valid)]
X_test <- X[(n_train+n_valid+1-seq_len):(n_train+n_valid+n_test), , drop=F]
y_test <- Y[(n_train+n_valid+1-seq_len):(n_train+n_valid+n_test)]

X_train_all <- X[1:(n_train+n_valid), , drop=F]
y_train_all <- Y[1:(n_train+n_valid)]


lagged_train_all <- lagged_features(X=cbind(y_train_all,X_train_all), max_lag=grf_pars$timedep, drop_present=TRUE)
lagged_test <- lagged_features(X=cbind(y_test,X_test), max_lag=grf_pars$timedep, drop_present=TRUE)


## ======= INTERMEDIATE QUANTILES =======

if(intermediate_method=="qrn"){
  # QRN fit
  foldwise_obj <- readRDS(paste0(interm_path, "Results_", qrn_params_string, ".rds"))
  foldwise_obj$fits <- list()
  for (k in 1:(foldwise_obj$n_folds)) {
    foldwise_obj$fits[[k]] <- EQRN_load(paste0(interm_path, "networks/"), paste0(qrn_params_string, "_f", k))
  }
  #Intermediate quantiles on all data
  interm_quantiles_all <- matrix(foldwise_obj$predictions, ncol=1)
  thresh_quant_all <- rbind(interm_quantiles_all, foldwise_obj$test_predictions[,(foldwise_obj$n_folds+1), drop=F])
  
}else if(intermediate_method=="grf"){
  #Intermediate quantiles all data
  fit_grf <- quantile_forest(lagged_train_all, y_train_all[(grf_pars$timedep+1):length(y_train_all)],
                             num.trees=grf_pars$num.trees, quantiles=grf_pars$quantiles_fit, sample.fraction=grf_pars$sample.fraction,# mtry=grf_pars$mtry,
                             min.node.size=grf_pars$min.node.size, honesty=grf_pars$honesty,
                             honesty.fraction=grf_pars$honesty.fraction, honesty.prune.leaves=grf_pars$honesty.prune.leaves,
                             alpha=grf_pars$alpha, imbalance.penalty=grf_pars$imbalance.penalty, seed=seedGRF)#,compute.oob.predictions = FALSE
  pred_grf_train_all <- predict(fit_grf, newdata=NULL, quantiles = c(interm_lvl))$predictions
  pred_grf_test <- predict(fit_grf, newdata=lagged_test, quantiles = c(interm_lvl))$predictions
  thresh_quant_all <- rbind(matrix(as.double(NA), nrow=grf_pars$timedep, ncol=1), pred_grf_train_all, pred_grf_test)
}else if(intermediate_method=="oracle"){
  stop("No ground truth for real data.")
}
intermediate_quantiles <- thresh_quant_all[1:n_train, , drop=F]
valid_quantiles <- thresh_quant_all[(n_train+1-par_qrn$seq_len):(n_train+n_valid), , drop=F]
interm_quantiles_all <- thresh_quant_all[1:(n_train+n_valid), , drop=F]
#Predict intermediate quantiles on X_test
pred_interm <- thresh_quant_all[(n_train+n_valid+1-par_qrn$seq_len):(n_train+n_valid+n_test), , drop=F]


lagged_interm_q_trall <- interm_quantiles_all[(grf_pars$timedep+1):length(interm_quantiles_all), , drop=F]
laged_interm_q_test <- pred_interm[(grf_pars$timedep+1):length(pred_interm), , drop=F]



## ======== GBEX CV =========

CV_results <- gbex_CV(X=lagged_train_all, y=y_train_all[(grf_pars$timedep+1):length(y_train_all)], intermediate_quantiles=lagged_interm_q_trall,
                      interm_lvl=interm_lvl, intermediate_q_feature=gbex_params$intermediate_q_feature, scale_features=gbex_params$scale_features,
                      num_folds=num_folds, Bmax=gbex_params$Bmax, grid_lambda_ratio=gbex_params$grid_lambda_ratio, grid_depth=gbex_params$grid_depth,
                      stratified=gbex_params$stratified, lambda_scale=gbex_params$lambda_scale, min_leaf_size=gbex_params$min_leaf_size, sf=gbex_params$sf,
                      parallel_strat=parallel_strat, n_workers=n_workers, seed=seedgbex, err_handling=err_handling)

filename <- paste0(save_path,"results_gbex_ts_CV_",format(Sys.time(),'%Y%m%d_%H%M%S'),".csv")

write_csv(CV_results$results_tibble, file=filename)

cat(paste0("\n Best params:\n",CV_results$best_params,"\n\n"))

end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)


