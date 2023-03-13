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
source("R/EQRN_loader.R")


## =============================== PARAMETERS ===============================
save_path <- "Results/Switzerland/intermediate_QR_RNN_gridsearch/"
check_directory(save_path, recursive=TRUE)
parallel_strat <- "multicore"#"sequential", "multisession", "multicore"
n_workers <- 8#(availableCores() - 1)
save_plots <- TRUE

seedR <- 0
seedT <- seedR
set.seed(seedR)

# PARAM: Data
Y_name <- c("station_62")
X_disch_names <- c("station_43")#c()
precip_names <- c("LTB","BRZ","MER","THU","GHS","BEP")#"INT","LTB","BRZ","MER","THU","SCE","GHS","BEP")
use_stl_remainder <- FALSE
use_mean_precip <- FALSE
prop_test <- 2/3
prop_valid <- 1/4

# PARAM: General
interm_lvl = 0.8

# Params: QRN
par_qrn <- list(
  nb_fits=3,
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

disch_names <- c(Y_name, X_disch_names)
sgroup <- c(disch_names, precip_names)

## =============================== Data Processing ===============================

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
X_train <- X[1:n_train, , drop=F]
y_train <- Y[1:n_train]
X_valid <- X[(n_train+1-seq_len):(n_train+n_valid), , drop=F]
y_valid <- Y[(n_train+1-seq_len):(n_train+n_valid)]
X_test <- X[(n_train+n_valid+1-seq_len):(n_train+n_valid+n_test), , drop=F]
y_test <- Y[(n_train+n_valid+1-seq_len):(n_train+n_valid+n_test)]

X_train_all <- X[1:(n_train+n_valid), , drop=F]
y_train_all <- Y[1:(n_train+n_valid)]

data_save <- list(X=X, Y=Y, Dates=Dates, n_train=n_train, n_valid=n_valid, n_test=n_test, seq_len=seq_len)
safe_save_rds(data_save, paste0(save_path, "Data_","backup",".rds"))

suppressWarnings(rm(Disch, Precipitations, Disch_stl_remainder, Disch_Precip, dat_tb, dat))

## =============================== QRN ===============================

# Unconditional losses (for comparison with QRN fit)
uncond_losses_interm <- list(train_loss=quantile_loss(y_train, quantile(y_train, interm_lvl), interm_lvl),
                             valid_loss=quantile_loss(y_valid, quantile(y_train, interm_lvl), interm_lvl))


grid_l <- purrr::cross(par_qrn, .filter = NULL)

`%fun%` <- set_doFuture_strategy(parallel_strat, n_workers=n_workers)

results_grid_fit <- foreach(params=grid_l, .errorhandling="remove", .combine=rbind) %fun% {
  
  params_string <- paste0("qrrnn_", params$rnn_type, "_", params$num_layers, "x", params$hidden_size, "_s", params$seq_len, "_do", params$p_drop*100,
                          "_L2", str_replace(toString(params$L2_pen),"([.])","p"), "_lr", str_replace(toString(params$learning_rate),"([.])","p"))
  cat("======== Start: ", params_string, " ========\n")
  
  # QRN fit
  torch_manual_seed(seedT)
  fit_qrn <- QRN_fit_multiple(X=X_train, y=y_train, q_level=interm_lvl, number_fits=3, hidden_size=params$hidden_size, num_layers=params$num_layers,
                               rnn_type=params$rnn_type, p_drop=params$p_drop, learning_rate=params$learning_rate, L2_pen=params$L2_pen,
                               seq_len=params$seq_len, scale_features=params$scale_features, n_epochs=params$n_epochs,
                               batch_size=params$batch_size, X_valid=X_valid, Y_valid=y_valid, lr_decay=params$lr_decay,
                               patience_decay=params$patience_decay, min_lr=params$min_lr, patience_stop=params$patience_stop, tol=params$tol)
  
  
  qrn_pred_train <- QRN_seq_predict(fit_qrn, X_train, y_train, q_level=interm_lvl)
  qrn_pred_valid <- QRN_seq_predict(fit_qrn, X_valid, y_valid, q_level=interm_lvl)[params$seq_len+(1:n_valid)]
  qrn_pred_test <- QRN_seq_predict(fit_qrn, X_test, y_test, q_level=interm_lvl)[params$seq_len+(1:n_test)]
  
  
  train_plot <- training_plot_eqrn(fit_qrn, uncond_losses_interm, title="QR RNN training", y_lab="Quantile loss")
  plot(train_plot)
  if(save_plots){
    ggsave(paste0(params_string, "_training_last.pdf"), plot=train_plot, device="pdf",
           path=save_path, width=160, height=120, units="mm", dpi=300)
  }
  
  # Compute losses for desired predicted quantiles
  q_loss <- quantile_loss(y=y_test[params$seq_len+(1:n_test)],
                          y_hat=qrn_pred_test, q=interm_lvl)
  
  output <- c(unlist(params),
              Train_loss=fit_qrn$train_loss[length(fit_qrn$train_loss)], Valid_loss=fit_qrn$valid_loss[length(fit_qrn$valid_loss)],
              min_valid_loss=min(fit_qrn$valid_loss, na.rm=TRUE), min_valid_e=which.min(fit_qrn$valid_loss),
              test_loss=q_loss)
  RESULTS <- list(qrn_pred_train=qrn_pred_train, qrn_pred_valid=qrn_pred_valid, qrn_pred_test=qrn_pred_test,
                  train_loss=fit_qrn$train_loss, valid_loss=fit_qrn$valid_loss, params_string=params_string, output=output)
  try(safe_save_rds(RESULTS, paste0(save_path, "Results_",params_string,".rds")))
  try(EQRN_save(fit_qrn, paste0(save_path, "networks/"), params_string))
  cat("==== End: ", params_string, " ====\n")
  output
}
results_tibble <- tibble::as_tibble(results_grid_fit)
results_tibble <- bind_cols(use_stl_remainder=use_stl_remainder, use_mean_precip=use_mean_precip,
                            n_train=n_train, p=ncol(X_train), interm_lvl=interm_lvl, results_tibble)

filename <- paste0(save_path,"results_QR_RNN_grid_fit_",format(Sys.time(),'%Y%m%d_%H%M%S'),".csv")

write_csv(results_tibble, file=filename)

end_doFuture_strategy()

end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)


