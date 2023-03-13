library(tidyverse)
library(here)
library(ggpubr)
library(future)
library(doFuture)
source("R/EQRN_loader.R")

## =============================== PARAMETERS ===============================
save_path <- "Results/ts/EQRN_gridsearch_ts/"
check_directory(save_path, recursive=TRUE)
parallel_strat <- "multisession"#"sequential", "multisession", "multicore"
n_workers <- 12#(availableCores() - 1)#PARAM
save_plots <- TRUE

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
quantiles_predict = c(0.995,0.999,0.9995)#c(interm_lvl,0.995,0.999,0.9995)

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

# Params: EQRN
params_list <- list(
  nb_fits = 3,
  shape_fixed=TRUE,
  rnn_type=list("lstm"),#"gru","lstm"
  num_layers=list(1,2,3),
  hidden_size=list(32,64,128,256),
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
  #mtry = min(ceiling(sqrt(p) + 20), p),
  min.node.size = 5,#5
  honesty = TRUE,
  honesty.fraction = 0.5,
  honesty.prune.leaves = TRUE,
  alpha = 0.05,
  imbalance.penalty = 0
)

## =============================== END PARAMETERS ===============================

qrn_params_string <- paste0("qrrnn_", par_qrn$rnn_type, "_", par_qrn$num_layers, "x", par_qrn$hidden_size, "_s", par_qrn$seq_len, "_do", par_qrn$p_drop*100,
                            "_L2", str_replace(toString(par_qrn$L2_pen),"([.])","p"), "_lr", str_replace(toString(par_qrn$learning_rate),"([.])","p"))


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

true_quantiles <- series_theoretical_quantiles(quantiles_predict, dat, Y_distr=Y_distr)
seq_len <- par_qrn$seq_len
X_train <- dat$X[1:n, , drop=F]
y_train <- dat$Y[1:n]
true_quantiles_train <- true_quantiles[1:n, , drop=F]
X_valid <- dat$X[(n+1-seq_len):(n+n_valid), , drop=F]
y_valid <- dat$Y[(n+1-seq_len):(n+n_valid)]
true_quantiles_valid <- true_quantiles[(n+1):(n+n_valid), , drop=F]
X_test <- dat$X[(n+n_valid+1-seq_len):(n+n_valid+ntest), , drop=F]
y_test <- dat$Y[(n+n_valid+1-seq_len):(n+n_valid+ntest)]
true_quantiles_test <- true_quantiles[(n+n_valid+1):(n+n_valid+ntest), , drop=F]

X_train_all <- dat$X[1:(n+n_valid), , drop=F]
y_train_all <- dat$Y[1:(n+n_valid)]
true_quantiles_train_all <- true_quantiles[1:(n+n_valid), , drop=F]


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
valid_quantiles <- thresh_quant_all[(n+1-par_qrn$seq_len):(n+n_valid), , drop=F]
interm_quantiles_all <- thresh_quant_all[1:(n+n_valid), , drop=F]
#Predict intermediate quantiles on X_test
pred_interm <- thresh_quant_all[(n+n_valid+1-par_qrn$seq_len):(n+n_valid+ntest), , drop=F]


## ======== For TESTING ========

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
                                  patience_stop=params$patience_stop, tol=params$tol, orthogonal_gpd=params$orthogonal_gpd, data_type="seq")
  
  try(EQRN_save(fit_eqrn, paste0(save_path, "networks/"), params_string))
  
  cat("\nElapsed time (EQRN fit):\n")
  print(Sys.time() - start_time)
  
  
  ## ======== TESTING ========
  
  #Final EQRN predictions on X_test
  pred_eqrn <- EQRN_predict_seq(fit_eqrn, X_test, y_test, quantiles_predict, pred_interm, interm_lvl, crop_predictions=TRUE)#[params$seq_len+(1:ntest)]
  
  
  train_plot <- training_plot_eqrn(fit_eqrn, uncond_losses_interm)
  plot(train_plot)
  if(save_plots){
    ggsave(paste0(params_string, "_training.pdf"), plot=train_plot, device="pdf",
           path=save_path, width=160, height=120, units="mm", dpi=300)
  }
  
  # Compute losses for desired predicted quantiles
  MSE_losses <- multilevel_MSE(pred_true, pred_eqrn, quantiles_predict, prefix="test_", give_names=TRUE)
  MAE_losses <- multilevel_MAE(pred_true, pred_eqrn, quantiles_predict, prefix="test_", give_names=TRUE)
  
  #test_loss <- compute_EQRN_seq_GPDLoss(fit_eqrn, X_test, y_test, intermediate_quantiles=pred_interm, interm_lvl=interm_lvl)
  output <- c(unlist(params),
              Train_loss=fit_eqrn$train_loss[length(fit_eqrn$train_loss)], Valid_loss=fit_eqrn$valid_loss[length(fit_eqrn$valid_loss)],
              min_valid_loss=min(fit_eqrn$valid_loss, na.rm=TRUE), min_valid_e=which.min(fit_eqrn$valid_loss),
              #test_loss=test_loss,
              MAE_losses, MSE_losses)
  RESULTS <- list(params_string=params_string, output=output)
  try(safe_save_rds(RESULTS, paste0(save_path, "Results_",params_string,".rds")))
  cat("==== End: ", params_string, " ====\n")
  output
}
results_tibble <- tibble::as_tibble(results_grid_fit)
results_tibble <- bind_cols(Y_distr=Y_distr, n=n, df=df, interm_lvl=interm_lvl, results_tibble)

filename <- paste0(save_path,"results_EQRN_ts_grid_fit_",format(Sys.time(),'%Y%m%d_%H%M%S'),".csv")

write_csv(results_tibble, file=filename)

end_doFuture_strategy()


end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)


#Sys.sleep(300)
#system('shutdown -t 30 -s')


