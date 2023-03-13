library(tidyverse)
library(here)
library(ggpubr)
source("R/EQRN_loader.R")

## =============================== PARAMETERS ===============================
save_path <- "Results/ts/intermediate_QR_RNN_gridsearch/"
check_directory(save_path, recursive=TRUE)
parallel_strat <- "multicore"#"sequential", "multisession", "multicore"
n_workers <- 12#(availableCores() - 1)#PARAM
save_plots <- TRUE

seedR <- 0
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
interm_lvl = 0.8

# Params: QRN
par_qrn <- list(
  rnn_type=list("lstm"),#"gru","lstm"
  num_layers=list(1,2,3),
  hidden_size=list(32,64,128,256),
  p_drop=0,
  L2_pen=list(0,1e-6,1e-5,1e-4),
  seq_len=10,
  learning_rate=list(1e-3),
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


start_time <- Sys.time()

# Data
dat <- generate_series_model(n=n+n_valid+ntest, df=df, AR=AR, MA=MA, muX=muX, mu0=mu0,
                             alphas=alphas, betas=betas, sX=sX, S0=S0, ARX=ARX,
                             X_distr=X_distr, Y_distr=Y_distr, seasonal_hetero=seasonal_hetero)
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

# Unconditional parameters and losses (for comparison with EQRNN fit)
uncond_losses_interm <- list(train_loss=quantile_loss(y_train, quantile(y_train, interm_lvl), interm_lvl),
                             valid_loss=quantile_loss(y_valid, quantile(y_train, interm_lvl), interm_lvl))

#plot(dat$Y)
#lines(true_quantiles)
X_train_all <- dat$X[1:(n+n_valid), , drop=F]
y_train_all <- dat$Y[1:(n+n_valid)]
true_quantiles_train_all <- true_quantiles[1:(n+n_valid)]

data_save <- list(dat=dat, n=n, n_valid=n_valid, ntest=ntest, df=df, AR=AR, MA=MA, muX=muX, mu0=mu0,
                  alphas=alphas, betas=betas, sX=sX, S0=S0, ARX=ARX, X_distr=X_distr, Y_distr=Y_distr)
safe_save_rds(data_save, paste0(save_path, "Data_","backup",".rds"))

grid_l <- purrr::cross(par_qrn, .filter = NULL)

`%fun%` <- set_doFuture_strategy(parallel_strat, n_workers=n_workers)

results_grid_fit <- foreach(params=grid_l, .errorhandling="remove", .combine=rbind) %fun% {
  
  params_string <- paste0("qrrnn_", params$rnn_type, "_", params$num_layers, "x", params$hidden_size, "_s", params$seq_len, "_do", params$p_drop*100,
                          "_L2", str_replace(toString(params$L2_pen),"([.])","p"), "_lr", str_replace(toString(params$learning_rate),"([.])","p"))
  cat("======== Start: ", params_string, " ========\n")
  
  # QRNN fit
  fit_qrn <- QRN_seq_fit(X=X_train, Y=y_train, q_level=interm_lvl, hidden_size=params$hidden_size, num_layers=params$num_layers,
                         rnn_type=params$rnn_type, p_drop=params$p_drop, learning_rate=params$learning_rate, L2_pen=params$L2_pen,
                         seq_len=params$seq_len, scale_features=params$scale_features, n_epochs=params$n_epochs,
                         batch_size=params$batch_size, X_valid=X_valid, Y_valid=y_valid, lr_decay=params$lr_decay,
                         patience_decay=params$patience_decay, min_lr=params$min_lr, patience_stop=params$patience_stop, tol=params$tol)
  
  
  qrn_pred_train <- QRN_seq_predict(fit_qrn, X_train, y_train, q_level=interm_lvl)
  qrn_pred_valid <- QRN_seq_predict(fit_qrn, X_valid, y_valid, q_level=interm_lvl)[params$seq_len+(1:n_valid)]
  qrn_pred_test <- QRN_seq_predict(fit_qrn, X_test, y_test, q_level=interm_lvl)[params$seq_len+(1:ntest)]
  
  
  train_plot <- training_plot_eqrn(fit_qrn, uncond_losses_interm, title="QR RNN training", y_lab="Quantile loss")
  plot(train_plot)
  if(save_plots){
    ggsave(paste0(params_string, "_training_last.pdf"), plot=train_plot, device="pdf",#substr(sims_codes[1],1,str_length(sims_codes[1])-6)
           path=save_path, width=160, height=120, units="mm", dpi=300)
  }
  
  # Compute losses for desired predicted quantiles
  MAE_loss <- mean_absolute_error(true_quantiles_test, qrn_pred_test)
  MSE_loss <- mean_squared_error(true_quantiles_test, qrn_pred_test)
  q_loss <- quantile_loss(y=y_test[params$seq_len+(1:ntest)],
                          y_hat=qrn_pred_test, q=interm_lvl)
  
  output <- c(unlist(params),
              Train_loss=fit_qrn$train_loss[length(fit_qrn$train_loss)], Valid_loss=fit_qrn$valid_loss[length(fit_qrn$valid_loss)],
              min_valid_loss=min(fit_qrn$valid_loss, na.rm=TRUE), min_valid_e=which.min(fit_qrn$valid_loss),
              test_loss=q_loss, MAE_test=MAE_loss, MSE_test=MSE_loss)
  RESULTS <- list(qrn_pred_train=qrn_pred_train, qrn_pred_valid=qrn_pred_valid, qrn_pred_test=qrn_pred_test,
                  train_loss=fit_qrn$train_loss, valid_loss=fit_qrn$valid_loss, params_string=params_string, output=output)
  try(safe_save_rds(RESULTS, paste0(save_path, "Results_",params_string,".rds")))
  try(EQRN_save(fit_qrn, paste0(save_path, "networks/"), params_string))
  cat("==== End: ", params_string, " ====\n")
  output
}
results_tibble <- tibble::as_tibble(results_grid_fit)
results_tibble <- bind_cols(Y_distr=Y_distr, X_distr=X_distr, n=n, p=ncol(X_train), df=df, interm_lvl=interm_lvl, results_tibble)

filename <- paste0(save_path,"results_QR_RNN_grid_fit_",format(Sys.time(),'%Y%m%d_%H%M%S'),".csv")

write_csv(results_tibble, file=filename)

end_doFuture_strategy()

end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)



