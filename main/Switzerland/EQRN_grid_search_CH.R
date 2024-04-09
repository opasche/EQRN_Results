library(tidyverse)
library(here)
library(ggpubr)
library(future)
library(doFuture)
source("R/EQRN_loader.R")

## =============================== PARAMETERS ===============================
save_path <- "Results/Switzerland/EQRN_gridsearch/"
check_directory(save_path, recursive=TRUE)
parallel_strat <- "sequential"#"sequential", "multisession", "multicore"
n_workers <- 8#(availableCores() - 1)#PARAM
save_plots <- TRUE

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
intermediate_method <- "qrn"
interm_lvl = 0.8
prob_lvls_predict = c(0.995,0.999,0.9995)#c(interm_lvl,0.995,0.999,0.9995)

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

# Params: EQRN
params_list <- list(
  nb_fits = 3,
  shape_fixed=TRUE,
  rnn_type=list("lstm"),#"gru","lstm"
  num_layers=list(1,2,3),
  hidden_size=list(8,16,32,64,128,256),
  p_drop=0,
  intermediate_q_feature=TRUE,
  L2_pen=list(0,1e-6,1e-5,1e-4),
  shape_penalty=0,
  seq_len=par_qrn$seq_len,
  learning_rate=list(1e-3),
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


## =============================== END PARAMETERS ===============================

qrn_params_string <- paste0("qrrnn_", par_qrn$rnn_type, "_", par_qrn$num_layers, "x", par_qrn$hidden_size, "_s", par_qrn$seq_len, "_do", par_qrn$p_drop*100,
                            "_L2", str_replace(toString(par_qrn$L2_pen),"([.])","p"), "_lr", str_replace(toString(par_qrn$learning_rate),"([.])","p"))


start_time <- Sys.time()

disch_names <- c(Y_name, X_disch_names)
sgroup <- c(disch_names, precip_names)

# Data
if(!file.exists(paste0(interm_path,"Data_backup.rds"))){
  #warning("File not found, new data is generated.")
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
  # n_train_all <- n_train+n_valid
  # n_tot <- n_train_all+n_test
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
#Predict intermediate quantiles on X_test
pred_interm <- thresh_quant_all[(n_train+n_valid+1-par_qrn$seq_len):(n_train+n_valid+n_test), , drop=F]


## ======== For TESTING ========

# UNCONDITIONAL predicted quantile(s) (Y quantile on X_train)
pred_unc <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = prob_lvls_predict, Y = y_train_all, ntest = n_test)

# SEMI-CONDITIONAL predicted quantiles
pred_semicond <- predict_GPD_semiconditional(Y=y_train[(seq_len+1):n_train], interm_lvl=interm_lvl, thresh_quantiles=intermediate_quantiles[(seq_len+1):n_train],
                                             interm_quantiles_test=pred_interm[seq_len+(1:n_test)], prob_lvls_predict=prob_lvls_predict)

# Unconditional parameters and losses (for comparison with EQRN fit)
uncond_losses_fixed <- unconditional_train_valid_GPD_loss(Y_train=y_train[(seq_len+1):n_train], interm_lvl=interm_lvl, Y_valid=y_valid[seq_len+(1:n_valid)])
uncond_losses_interm <- semiconditional_train_valid_GPD_loss(Y_train=y_train[(seq_len+1):n_train], Y_valid=y_valid[seq_len+(1:n_valid)],
                                                                             interm_quant_train=intermediate_quantiles[(seq_len+1):n_train],
                                                                             interm_quant_valid=valid_quantiles[seq_len+(1:n_valid)])


## ======== EQRN FIT ========

grid_l <- purrr::cross(params_list, .filter = NULL)

`%fun%` <- set_doFuture_strategy(parallel_strat, n_workers=n_workers)

results_grid_fit <- foreach(params=grid_l, .errorhandling="remove", .combine=rbind) %fun% {
  #source("R/EQRN_loader.R")
  
  params_string <- paste0("reqrnn_", params$rnn_type, "_", params$num_layers, "x", params$hidden_size, "_",
                          "n"[!params$intermediate_q_feature], "u_s", params$seq_len, "_do", params$p_drop*100,
                          "_L2", str_replace(toString(params$L2_pen),"([.])","p"), "_lr", str_replace(toString(params$learning_rate),"([.])","p"))
  cat("======== Start: ", params_string, " ========\n")
  #Fit EQRN with intermediate OOB (or oracle) quantiles
  torch_manual_seed(seedT)
  fit_eqrn <- EQRN_fit_restart(X_train, y_train, intermediate_quantiles=intermediate_quantiles, interm_lvl=interm_lvl, number_fits=params$nb_fits,
                                shape_fixed=params$shape_fixed, hidden_size=params$hidden_size, num_layers=params$num_layers,
                                rnn_type=params$rnn_type, p_drop=params$p_drop, intermediate_q_feature=params$intermediate_q_feature,
                                learning_rate=params$learning_rate, L2_pen=params$L2_pen, seq_len=params$seq_len,
                                shape_penalty=params$shape_penalty, scale_features=params$scale_features, n_epochs=params$n_epochs,
                                batch_size=params$batch_size, X_valid=X_valid, y_valid=y_valid, quant_valid=valid_quantiles,
                                lr_decay=params$lr_decay, patience_decay=params$patience_decay, min_lr=params$min_lr,
                                patience_stop=params$patience_stop, tol=params$tol, orthogonal_gpd=params$orthogonal_gpd,
                                patience_lag=params$patience_lag, data_type="seq")
  
  try(EQRN_save(fit_eqrn, paste0(save_path, "networks/"), params_string))
  
  cat("\nElapsed time (EQRN fit):\n")
  print(Sys.time() - start_time)
  
  
  ## ======== TESTING ========
  
  #Final EQRN predictions on X_valid and X_test
  pred_eqrn_val <- EQRN_predict_seq(fit_eqrn, X_valid, y_valid, prob_lvls_predict, valid_quantiles, interm_lvl, crop_predictions=TRUE)
  pred_eqrn <- EQRN_predict_seq(fit_eqrn, X_test, y_test, prob_lvls_predict, pred_interm, interm_lvl, crop_predictions=TRUE)#[params$seq_len+(1:n_test)]
  
  
  train_plot <- training_plot_eqrn(fit_eqrn, uncond_losses_interm)
  #plot(train_plot)
  if(save_plots){
    ggsave(paste0(params_string, "_training.pdf"), plot=train_plot, device="pdf",#substr(sims_codes[1],1,str_length(sims_codes[1])-6)
           path=save_path, width=160, height=120, units="mm", dpi=300)
  }
  
  # Compute losses for desired predicted quantiles
  nb_prob_lvls_predict <- length(prob_lvls_predict)
  qPredErrs_v <- multilevel_q_pred_error(y=y_valid[seq_len+(1:n_valid)], Pred_Q=pred_eqrn_val,
                                         proba_levels=prob_lvls_predict, prefix="valid_", na.rm=FALSE)#rep(as.double(NA), nb_prob_lvls_predict)
  qPredErrs_t <- multilevel_q_pred_error(y=y_test[seq_len+(1:n_test)], Pred_Q=pred_eqrn,
                                         proba_levels=prob_lvls_predict, prefix="test_", na.rm=FALSE)
  qlosses_v <- multilevel_q_loss(y=y_valid[seq_len+(1:n_valid)], Pred_Q=pred_eqrn_val,
                                 proba_levels=prob_lvls_predict, prefix="valid_", na.rm=FALSE)
  qlosses_t <- multilevel_q_loss(y=y_test[seq_len+(1:n_test)], Pred_Q=pred_eqrn,
                                 proba_levels=prob_lvls_predict, prefix="test_", na.rm=FALSE)
  propb_v <- multilevel_prop_below(y=y_valid[seq_len+(1:n_valid)], Pred_Q=pred_eqrn_val,
                                   proba_levels=prob_lvls_predict, prefix="valid_", na.rm=FALSE)
  propb_t <- multilevel_prop_below(y=y_test[seq_len+(1:n_test)], Pred_Q=pred_eqrn,
                                   proba_levels=prob_lvls_predict, prefix="test_", na.rm=FALSE)
  
  test_loss <- compute_EQRN_seq_GPDLoss(fit_eqrn, X_test, y_test, intermediate_quantiles=pred_interm, interm_lvl=interm_lvl)
  output <- c(unlist(params),
              Train_loss=fit_eqrn$train_loss[length(fit_eqrn$train_loss)], Valid_loss=fit_eqrn$valid_loss[length(fit_eqrn$valid_loss)],
              min_valid_loss=min(fit_eqrn$valid_loss, na.rm=TRUE), min_valid_e=which.min(fit_eqrn$valid_loss),
              test_loss=test_loss, qPredErrs_v, qlosses_v, propb_v, qPredErrs_t, qlosses_t, propb_t)
  RESULTS <- list(params_string=params_string, output=output)
  try(safe_save_rds(RESULTS, paste0(save_path, "Results_",params_string,".rds")))
  cat("==== End: ", params_string, " ====\n")
  output
}
results_tibble <- tibble::as_tibble(results_grid_fit)
results_tibble <- bind_cols(use_stl_remainder=use_stl_remainder, use_mean_precip=use_mean_precip,
                            n_train=n_train, p=ncol(X_train), interm_lvl=interm_lvl, results_tibble)

filename <- paste0(save_path,"results_EQRN_ts_grid_fit_",format(Sys.time(),'%Y%m%d_%H%M%S'),".csv")

write_csv(results_tibble, file=filename)

end_doFuture_strategy()


end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)


#Sys.sleep(300)
#system('shutdown -t 30 -s')



