library(tidyverse)
library(here)
library(ggpubr)
library(grf)
source("R/gbex_wrappers.R")
source("R/EGAM_wrappers.R")
source("R/EQRN_loader.R")

## =============================== PARAMETERS ===============================

param_folder <- "main/EQRN_ts/parameters/"
param_file_name <- "threshold_analysis_EQRN_bestval.R"
# param_file_name <- "threshold_analysis_EQRN_bestval_nopen.R"

param_file <- paste0(param_folder, param_file_name)
source(param_file)

parallel_strat <- "sequential"#"sequential", "multisession", "multicore"
n_workers <- min(length(interm_lvls), 12)


## =============================== END PARAMETERS ===============================

check_directory(save_path, recursive=TRUE)
check_directory(path_eqrn, recursive=TRUE)

start_time <- Sys.time()

`%fun%` <- set_doFuture_strategy(parallel_strat, n_workers=n_workers)

results_thresh_fit <- foreach(interm_lvl=interm_lvls, .errorhandling="stop", .combine=rbind) %fun% {
  set.seed(seedR)
  itt_time <- Sys.time()
  
  qrn_params_string <- paste0("qrrnn_", par_qrn$rnn_type, "_", par_qrn$num_layers, "x", par_qrn$hidden_size, "_s", par_qrn$seq_len, "_do", par_qrn$p_drop*100,
                              "_L2", str_replace(toString(par_qrn$L2_pen),"([.])","p"), "_lr", str_replace(toString(par_qrn$learning_rate),"([.])","p"),
                              "_ilv",interm_lvl*100)
  
  params_string <- paste0("reqrn_", par_eqrn$rnn_type, "_", par_eqrn$num_layers, "x", par_eqrn$hidden_size, "_",
                          "n"[!par_eqrn$intermediate_q_feature], "u_s", par_eqrn$seq_len, "_do", par_eqrn$p_drop*100,
                          "_L2", str_replace(toString(par_eqrn$L2_pen),"([.])","p"), "_lr", str_replace(toString(par_eqrn$learning_rate),"([.])","p"),
                          "_ilv",interm_lvl*100)
  
  
  cat("==== Start: ", params_string, " ====\n")
  
  # Data
  if(!file.exists(paste0(interm_path,"Data_backup.rds"))){
    stop("Data file not found.")
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
                width = 150, height = 30, cairo=FALSE)
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
  interm_quantiles_all <- thresh_quant_all[1:(n_train_all), , drop=F]#rbind(intermediate_quantiles,valid_quantiles[par_qrn$seq_len:n_valid, , drop=F])
  #Predict intermediate quantiles on X_test
  pred_interm <- thresh_quant_all[(n_train_all+1-par_qrn$seq_len):(n_train_all+ntest), , drop=F]
  
  
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
                                  patience_stop=par_eqrn$patience_stop, tol=par_eqrn$tol, orthogonal_gpd=par_eqrn$orthogonal_gpd, data_type="seq")
    EQRN_save(fit_eqrn, paste0(path_eqrn, "networks/"), params_string)
  }else{
    fit_eqrn <- EQRN_load(paste0(path_eqrn, "networks/"), params_string)
  }
  
  
  cat("\nElapsed time (EQRN fit):\n")
  print(Sys.time() - itt_time)
  
  
  ## ======== TESTING ========
  
  #Final EQRN predictions on X_test
  pred_eqrn_val <- EQRN_predict_seq(fit_eqrn, X_valid, y_valid, quantiles_predict, valid_quantiles, interm_lvl, crop_predictions=TRUE)
  pred_eqrn <- EQRN_predict_seq(fit_eqrn, X_test, y_test, quantiles_predict, pred_interm, interm_lvl, crop_predictions=TRUE)#[params$seq_len+(1:ntest)]
  
  
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
  
  ## ======== GRF GBEX and GPDGAM FITS =========
  
  lagged_train_all <- lagged_features(X=cbind(y_train_all,X_train_all), max_lag=grf_pars$timedep, drop_present=TRUE)
  lagged_test <- lagged_features(X=cbind(y_test,X_test), max_lag=grf_pars$timedep, drop_present=TRUE)
  
  fit_grf <- quantile_forest(lagged_train_all, y_train_all[(grf_pars$timedep+1):length(y_train_all)],
                             num.trees=grf_pars$num.trees, quantiles=grf_pars$quantiles_fit, sample.fraction=grf_pars$sample.fraction,# mtry=grf_pars$mtry,
                             min.node.size=grf_pars$min.node.size, honesty=grf_pars$honesty,
                             honesty.fraction=grf_pars$honesty.fraction, honesty.prune.leaves=grf_pars$honesty.prune.leaves,
                             alpha=grf_pars$alpha, imbalance.penalty=grf_pars$imbalance.penalty, seed=seedGRF)#,compute.oob.predictions = FALSE
  
  if(interm_method_competitors=="qrn"){
    # QRN fit
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
  # fit_egam <- fit_gpd_gam(X=lagged_train_all, y=y_train_all[(grf_pars$timedep+1):length(y_train_all)],
  #                         intermediate_quantiles=lagged_interm_q_trall, interm_lvl=interm_lvl, model_shape=egam_params$model_shape,
  #                         intermediate_q_feature=egam_params$intermediate_q_feature, scale_features=egam_params$scale_features)
  # pred_egam <- predict_gpd_gam(fit_egam, lagged_test, to_predict=quantiles_predict,
  #                              intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
  ## ===========================
  
  
  train_plot <- training_plot_eqrn(fit_eqrn, uncond_losses_interm, show_legend=TRUE)
  plot(train_plot)
  if(save_plots){
    save_myplot(plt=train_plot, plt_nm=paste0(save_path, params_string, "_training.pdf"),
                width = 75, height = 75, cairo=FALSE)
  }
  
  # Compute losses for desired predicted quantiles
  nb_quantiles_predict <- length(quantiles_predict)
  MSE_obj <- multilevel_MSE(pred_true, pred_eqrn, quantiles_predict, prefix="test_", give_names=TRUE, sd=TRUE)
  MAE_obj <- multilevel_MAE(pred_true, pred_eqrn, quantiles_predict, prefix="test_", give_names=TRUE, sd=TRUE)
  MSE_losses <- MSE_obj$MSEs
  MAE_losses <- MAE_obj$MAEs
  MSE_SDs <- MSE_obj$SDs
  MAE_SDs <- MAE_obj$SDs
  for(i in 1:nb_quantiles_predict){
    predplot <- plot_predictions_ts(pred_eqrn[,i], pred_gbex[,i], pred_true[,i], NaN, quantiles_predict[i], name_other="gbex")
    plot(predplot)
    if(save_plots){
      ggsave(paste0(params_string, "_q",quantiles_predict[i]*10000,"_preds.pdf"), plot=predplot, device="pdf",
             path=save_path, width=300, height=200, units="mm", dpi=300)
    }
    
    predplotdiff <- plot_predictions_diff_ts(pred_eqrn[,i], pred_gbex[,i], pred_true[,i], quantiles_predict[i], name_other="gbex")
    plot(predplotdiff)
    if(save_plots){
      ggsave(paste0(params_string, "_q",quantiles_predict[i]*10000,"_preds_diff.pdf"), plot=predplotdiff, device="pdf",
             path=save_path, width=300, height=200, units="mm", dpi=300)
    }
    
    resid_box2 <- residuals_boxplot_eqrn2(pred_eqrn[,i], pred_grf_test[,i], pred_gbex[,i], pred_true[,i], quantiles_predict[i])
    plot(resid_box2)
    if(save_plots){
      ggsave(paste0(params_string, "_q",quantiles_predict[i]*10000,"_residbox2.pdf"), plot=resid_box2, device="pdf",
             path=save_path, width=200, height=150, units="mm", dpi=300)
    }
  }
  pred_eqrn_params <- EQRN_predict_params_seq(fit_eqrn, X_test, y_test, pred_interm, return_parametrization="classical", interm_lvl)
  # test_loss <- loss_GPD(pred_eqrn_params[,1], pred_eqrn_params[,2], (y_test-pred_interm), rescaled=TRUE, interm_lvl=interm_lvl)
  
  ## ===================== predictions vs truth plots ============================
  
  pred_vs_true_plot_grf <- plot_predictions_vs_truth(pred_eqrn, pred_grf_test, pred_semicond$predictions,
                                                     pred_true, quantiles_predict, name_other="GRF")
  pred_vs_true_plot_gbex <- plot_predictions_vs_truth(pred_eqrn, pred_gbex, pred_semicond$predictions,
                                                      pred_true, quantiles_predict, name_other="gbex")
  pred_vs_true_plot <- plot_predictions_vs_truth_solo(pred_eqrn, pred_semicond$predictions,
                                                      pred_true, quantiles_predict)
  # plot(pred_vs_true_plot_grf)
  # plot(pred_vs_true_plot_gbex)
  # plot(pred_vs_true_plot)
  if(save_plots){
    ggsave(paste0(params_string, "_predVStruth_grf.pdf"), plot=pred_vs_true_plot_grf, device="pdf",
           path=save_path, width=300, height=100, units="mm", dpi=300)
    ggsave(paste0(params_string, "_predVStruth_gbex.pdf"), plot=pred_vs_true_plot_gbex, device="pdf",
           path=save_path, width=300, height=100, units="mm", dpi=300)
    ggsave(paste0(params_string, "_predVStruth.pdf"), plot=pred_vs_true_plot, device="pdf",
           path=save_path, width=300, height=100, units="mm", dpi=300)
  }
  # 
  # ppvgbex <- plot_predictions_vs_competitor(pred_eqrn, pred_gbex, quantiles_predict, name_comp="gbex", xaxis="EQRN")
  # plot(ppvgbex)
  # ppvgbex2 <- plot_predictions_vs_competitor(pred_eqrn, pred_gbex, quantiles_predict, name_comp="gbex", xaxis="competitor")
  # plot(ppvgbex2)
  # 
  # ppvevgam <- plot_predictions_vs_competitor(pred_eqrn, pred_egam, quantiles_predict, name_comp="evGAM", xaxis="EQRN")
  # plot(ppvevgam)
  # ppvevgam2 <- plot_predictions_vs_competitor(pred_eqrn, pred_egam, quantiles_predict, name_comp="evGAM", xaxis="competitor")
  # plot(ppvevgam2)
  # 
  # ppvscond <- plot_predictions_vs_competitor(pred_eqrn, pred_semicond$predictions, quantiles_predict, name_comp="Semi-conditional", xaxis="EQRN")
  # plot(ppvscond)
  # ppvscond2 <- plot_predictions_vs_competitor(pred_eqrn, pred_semicond$predictions, quantiles_predict, name_comp="Semi-conditional", xaxis="competitor")
  # plot(ppvscond2)
  # if(save_plots){
  #   ggsave(paste0(params_string, "_eqrnn_vs_gbex.pdf"), plot=ppvgbex, device="pdf",
  #          path=save_path, width=300, height=100, units="mm", dpi=300)
  #   ggsave(paste0(params_string, "_eqrnn_vs_gbex_cx.pdf"), plot=ppvgbex2, device="pdf",
  #          path=save_path, width=300, height=100, units="mm", dpi=300)
  #   ggsave(paste0(params_string, "_eqrnn_vs_evgam.pdf"), plot=ppvevgam, device="pdf",
  #          path=save_path, width=300, height=100, units="mm", dpi=300)
  #   ggsave(paste0(params_string, "_eqrnn_vs_evgam_cx.pdf"), plot=ppvevgam2, device="pdf",
  #          path=save_path, width=300, height=100, units="mm", dpi=300)
  #   ggsave(paste0(params_string, "_eqrnn_vs_semicond.pdf"), plot=ppvscond, device="pdf",
  #          path=save_path, width=300, height=100, units="mm", dpi=300)
  #   ggsave(paste0(params_string, "_eqrnn_vs_semicond_cx.pdf"), plot=ppvscond2, device="pdf",
  #          path=save_path, width=300, height=100, units="mm", dpi=300)
  # }
  # 
  # quant_pred_plot <- c(0.8, 0.85, 0.9, 0.95, 0.99, 0.995, 0.999, 0.9995)
  # quant_pred_plot_smooth <- c(0.8,1-10^c(-1,-1.5,-2,-2.5,-3,-3.5,-4))
  # for (Q_pred_levels in list(list(qs=quant_pred_plot, nm="dq"),
  #                            list(qs=quant_pred_plot_smooth, nm="sq"))){
  #   
  #   true_q_test_plot <- series_theoretical_quantiles(Q_pred_levels$qs, dat, Y_distr=Y_distr)[(n_train_all+1):(n_train_all+ntest), , drop=F]
  #   err_quant_plt <- plot_error_quantile_ts_old(fit_eqrn, fit_grf, fit_gbex, fit_egam, y_train, X_test, y_test, Q_pred_levels$qs, pred_interm,
  #                                               interm_lvl, intermediate_quantiles, true_q_test_plot, factorize=TRUE, test_data="other", pred_interm_c=pred_interm_c)
  #   rmse_quant_plt_np <- plot_rmse_quantile_ts(fit_eqrn, fit_grf, fit_gbex, fit_egam, y_train, X_test, y_test, Q_pred_levels$qs, pred_interm,
  #                                              interm_lvl, intermediate_quantiles, true_q_test_plot, factorize=TRUE, test_data="other", pred_interm_c=pred_interm_c)
  #   plot(err_quant_plt)
  #   plot(rmse_quant_plt_np)
  #   if(save_plots){
  #     ggsave(paste0("deprecated_", params_string, "_errquant_",Q_pred_levels$nm,".pdf"), plot=err_quant_plt, device="pdf",
  #            path=save_path, width=200, height=100, units="mm", dpi=300)
  #     save_myplot(plt=rmse_quant_plt_np, plt_nm=paste0(save_path, params_string, "_rmsequant_",Q_pred_levels$nm,".pdf"),
  #                 width = 75, height = 75, cairo=FALSE)
  #   }
  # }
  # as.character(Sys.time()-itt_time)
  # pred_eqrn_params <- EQRN_predict_params(fit_eqrn, X_test, pred_interm, return_parametrization="classical", interm_lvl)
  # test_loss <- compute_EQRN_GPDLoss(fit_eqrn, X_test, y_test, intermediate_quantiles=pred_interm, interm_lvl=interm_lvl)
  output <- c(interm_lvl=interm_lvl,
              Train_loss=fit_eqrn$train_loss[length(fit_eqrn$train_loss)], Valid_loss=fit_eqrn$valid_loss[length(fit_eqrn$valid_loss)],
              MAE_losses, MAE_SDs, MSE_losses, MSE_SDs,# test_loss=test_loss,
              unlist(par_eqrn), params_string=params_string)
  RESULTS <- list(pred_true=pred_true, pred_eqrn=pred_eqrn, pred_grf_test=pred_grf_test, quantiles_predict=quantiles_predict,
                  train_loss=fit_eqrn$train_loss, valid_loss=fit_eqrn$valid_loss,
                  uncond_losses_interm=uncond_losses_interm, params_string=params_string, qrn_params_string=qrn_params_string,
                  output=output)
  safe_save_rds(RESULTS, paste0(save_path, "Results_",params_string,".rds"))
  resids <- (pred_eqrn-pred_true)
  colnames(resids) <- paste0("test_resid_q",quantiles_predict)
  write_csv(as_tibble(resids), file=paste0(save_path,"Residuals_",params_string,".csv"))
  cat("==== End: ", params_string, " ====\n")
  output
}
results_tibble <- tibble::as_tibble(results_thresh_fit)
results_tibble <- bind_cols(Y_distr=Y_distr, n=n, df=df, results_tibble)

filename <- paste0(save_path,"results_EQRN_threshold_analysis_",format(Sys.time(),'%Y%m%d_%H%M%S'),".csv")

write_csv(results_tibble, file=filename)

end_doFuture_strategy()

results_tibble <- read_csv(paste0(save_path,"results_EQRN_threshold_analysis.csv"))
IntQ_Q_errs <- results_tibble %>%
  tidyr::pivot_longer(cols = starts_with("test_M"), names_to = c(".value", "Prob_lvl"), names_sep = "_q") %>% 
  dplyr::mutate(Prob_lvl=as.double(Prob_lvl)) %>% 
  dplyr::relocate(Prob_lvl, starts_with("test_M"), .after="Valid_loss")

resid_tbl <- foreach(ilvl=interm_lvls, .errorhandling="stop", .combine=bind_rows) %do% {
  params_string <- filter(results_tibble, interm_lvl==ilvl)[["params_string"]]
  residsi <- read_csv(paste0(save_path,"Residuals_",params_string,".csv"))
  residsi$interm_lvl <- ilvl
  residsi
}

resid_tbl_l <- resid_tbl %>% tidyr::pivot_longer(cols=starts_with("test_resid_q"), names_to=c(".value", "Prob_lvl"), names_sep="_q") %>% 
  dplyr::mutate(Prob_lvl=as.double(Prob_lvl)) %>% 
  dplyr::relocate(interm_lvl, Prob_lvl, test_resid)


thresh_mse_plt <- plot_rmse_threshold(IntQ_Q_errs, interm_lvls="interm_lvl", predict_lvls="Prob_lvl", errs="test_MSE",
                                      factorize=TRUE, test_data="other", strip.position=NULL, hide_ylab=FALSE)
thresh_se_bplt <- plot_rse_box_threshold(resid_tbl_l, interm_lvls="interm_lvl", predict_lvls="Prob_lvl", residuals="test_resid",
                                         factorize=TRUE, test_data="other", strip.position=NULL, hide_ylab=FALSE)
plot(thresh_mse_plt)
plot(thresh_se_bplt)
if(!is.null(save_path)){
  save_myplot(plt=thresh_mse_plt, plt_nm=paste0(save_path, "threshold_rmse.pdf"),
              width = 50, height = 50, cairo=FALSE)
  save_myplot(plt=thresh_se_bplt, plt_nm=paste0(save_path, "threshold_squarederrs_box.pdf"),
              width = 50, height = 50, cairo=FALSE)
}



end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)



#Sys.sleep(300)
#system('shutdown -t 30 -s')
