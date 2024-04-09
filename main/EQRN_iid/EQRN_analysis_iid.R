library(tidyverse)
library(grf)
library(here)
library(ggpubr)
source("R/gbex_wrappers.R")
source("R/EGAM_wrappers.R")
source("R/EQRN_loader.R")

## =============================== PARAMETERS ===============================

param_folder <- "main/EQRN_iid/parameters/analysis/"
param_file_names <- c("binorm_sigmoid_bestval.R",
                      "binorm_sigmoid_bestMSE.R",
                      "cosnorm2d_sigmoid_bestval.R",
                      "cosnorm2d_sigmoid_bestMSE.R",
                      "cosnorm_sigmoid_bestval.R",
                      "cosnorm_sigmoid_bestMSE.R")
# binorm_sigmoid, cosnorm2d_sigmoid, cosnorm2d_sigmoid: three choices of data models
# bestval: best validation loss EQRN model from 'EQRN_grid_search'
# bestMSE: best MSE EQRN model (not realistic choice on real data where truth is unknown)


parallel_strat <- "sequential"#"sequential", "multisession", "multicore"
n_workers <- 6


## =============================== END PARAMETERS ===============================

param_files <- paste0(param_folder, param_file_names)


`%fun%` <- set_doFuture_strategy(parallel_strat, n_workers=n_workers)

fe_out <- foreach(param_file=param_files, .errorhandling="remove", .combine=rbind) %fun% {
  source(param_file)
  check_directory(save_path, recursive=TRUE)
  check_directory(path_data, recursive=TRUE)
  check_directory(path_eqrn, recursive=TRUE)
  set.seed(seedR)
  
  if(hid_f_str=="SSNN"){
    stuct_str <- paste0("sc", paste0(net_structure$scale, collapse="_"), "_sh", paste0(net_structure$shape, collapse="_"))
  } else {
    stuct_str <- paste0(net_structure,"h", collapse="_")
  }
  params_string <- paste0(stuct_str, "_", hid_f_str, "_", "n"[!intermediate_q_feature], "u_do", p_drop*100,
                          "_L2", str_replace(toString(L2_pen),"([.])","p"), "_shp",
                          str_replace(toString(shape_penalty),"([.])","p"), "_", lr_str)
  
  start_time <- Sys.time()
  
  ## Generate (or load) Data
  if(!file.exists(paste0(path_data,"Data_backup.rds"))){
    warning("File not found, new data is generated.")
    # Training data
    dat <- generate_joint_distribution(n = n, p = p, model = model, distr = distr, df = df)
    dat_valid <- generate_joint_distribution(n = n_valid, p = p, model = model, distr = distr, df = df)
    
    # Generate test data
    if (test_data == "grid"){
      ntest <- floor(sqrt(ntest)) ** 2
      warning(paste0("Modified ntest to: ", ntest))
    }
    X_test <- generate_test_data(ntest, p, test_data)
    y_test <- generate_conditional_distribution(model, distr, df, X_test)$Y
    
    data_save <- list(dat=dat, dat_valid=dat_valid, X_test=X_test, y_test=y_test,
                      n=n, n_valid=n_valid, ntest=ntest, p=p, model=model, distr=distr, df=df,
                      test_data=test_data)
    safe_save_rds(data_save, paste0(path_data,"Data_backup.rds"))
  }else{
    data_save <- readRDS(paste0(path_data,"Data_backup.rds"))
    dat <- data_save$dat
    dat_valid <- data_save$dat_valid
    X_test <- data_save$X_test
    y_test <- data_save$y_test
    
    verif_dat <- c(n==data_save$n, n_valid==data_save$n_valid, ntest==data_save$ntest, p==data_save$p,
                   model==data_save$model, distr==data_save$distr, df==data_save$df, test_data==data_save$test_data)
    if(any(!verif_dat)){stop("Data mixup: Issue with data generating parameters")}
  }
  rm(data_save)
  X_train_all <- rbind(dat$X,dat_valid$X)
  y_train_all <- c(dat$Y,dat_valid$Y)
  
  
  # GRF fit
  if(!file.exists(paste0(path_data,"GRF_fit.rds")) | force_refit){
    fit_grf <- quantile_forest(dat$X, dat$Y, num.trees=num.trees, quantiles=quantiles_fit, sample.fraction=sample.fraction, mtry=mtry,
                               min.node.size=min.node.size, honesty=honesty, honesty.fraction=honesty.fraction, honesty.prune.leaves=honesty.prune.leaves,
                               alpha=alpha, imbalance.penalty=imbalance.penalty, seed=seedGRF)
    warning("Save not found, new intermediate GRF quantiles are fitted.")
    safe_save_rds(fit_grf, paste0(path_data,"GRF_fit.rds"))
  }else{
    fit_grf <- readRDS(paste0(path_data,"GRF_fit.rds"))
  }
  cat("\nElapsed time (GRF fit):\n")
  print(Sys.time() - start_time)
  
  ## ======= INTERMEDIATE QUANTILES =======
  
  if(intermediate_method=="grf"){
    #Out of bag quantiles prediction on dat$X (for other QR ML methods, would need foldwise construction)
    intermediate_quantiles <- predict(fit_grf, newdata = NULL, quantiles = c(interm_lvl))$predictions
    valid_quantiles <- predict(fit_grf, newdata = dat_valid$X, quantiles = c(interm_lvl))$predictions
    interm_quantiles_all <- rbind(intermediate_quantiles,valid_quantiles)
    #Predict intermediate quantiles on X_test with GRF
    pred_interm <- predict(fit_grf, newdata = X_test, quantiles = c(interm_lvl))$predictions
    pred_interm_fromall <- pred_interm
  }else if(intermediate_method=="oracle"){
    #Intermediate quantiles prediction on dat$X and validation set
    intermediate_quantiles <- generate_theoretical_quantiles(quantiles = c(interm_lvl), X = dat$X, model = model, distr = distr, df = df)
    valid_quantiles <- generate_theoretical_quantiles(quantiles = c(interm_lvl), X = dat_valid$X, model = model, distr = distr, df = df)
    interm_quantiles_all <- rbind(intermediate_quantiles,valid_quantiles)
    #Predict intermediate quantiles on X_test
    pred_interm <- generate_theoretical_quantiles(quantiles = c(interm_lvl), X = X_test, model = model, distr = distr, df = df)
    pred_interm_fromall <- pred_interm
  }
  
  ## ======== EQRN FIT ========
  
  #Fit EQRN with intermediate OOB (or oracle) quantiles
  if((!dir.exists(paste0(path_eqrn, "networks/", params_string,"/"))) | force_refit){
    warning("Save not found, new EQRN model is fitted.")
    torch_manual_seed(seedT)
    fit_eqrn <- EQRN_fit_restart(dat$X, dat$Y, intermediate_quantiles, interm_lvl=interm_lvl, number_fits=nb_fits_eqrn, shape_fixed=shape_fixed,
                                  net_structure=net_structure, hidden_fct=hidden_fct, p_drop=p_drop, intermediate_q_feature=intermediate_q_feature,
                                  learning_rate=learning_rate, L2_pen=L2_pen, shape_penalty=shape_penalty, scale_features=scale_features, n_epochs=n_epochs,
                                  batch_size=batch_size, X_valid=dat_valid$X, y_valid=dat_valid$Y, quant_valid=valid_quantiles,
                                  lr_decay=lr_decay, patience_decay=patience_decay, min_lr=min_lr, patience_stop=patience_stop, orthogonal_gpd=orthogonal_gpd)
    EQRN_save(fit_eqrn, paste0(path_eqrn, "networks/"), params_string)
  }else{
    fit_eqrn <- EQRN_load(paste0(path_eqrn, "networks/"), params_string)
  }
  cat("\nElapsed time (EQRN fit):\n")
  print(Sys.time() - start_time)
  
  
  ## ======== TESTING ========
  
  #High quantile prediction with GRF
  pred_grf_test <- predict(fit_grf, newdata = X_test, quantiles = prob_lvls_predict)$predictions
  
  # UNCONDITIONAL predicted quantile(s) (Y quantile on X_train)
  pred_unc <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = prob_lvls_predict, Y = dat$Y, ntest = ntest)
  
  # SEMI-CONDITIONAL predicted quantiles
  pred_semicond <- predict_GPD_semiconditional(Y=dat$Y, interm_lvl=interm_lvl, thresh_quantiles=intermediate_quantiles,
                                               interm_quantiles_test=pred_interm, prob_lvls_predict=prob_lvls_predict)
  
  # GROUND-TRUTH (y_test)
  pred_true <- generate_theoretical_quantiles(quantiles = prob_lvls_predict, X = X_test, model = model, distr = distr, df = df)
  
  # Unconditional parameters and losses (for comparison with EQRN fit)
  uncond_losses_fixed <- unconditional_train_valid_GPD_loss(Y_train=dat$Y, interm_lvl=interm_lvl, Y_valid=dat_valid$Y)
  uncond_losses_interm <- semiconditional_train_valid_GPD_loss(Y_train=dat$Y, Y_valid=dat_valid$Y,
                                                                               interm_quant_train=intermediate_quantiles,
                                                                               interm_quant_valid=valid_quantiles)
  
  ## ======== GBEX and GPD GAM FITS =========
  fit_gbex <- gbex_fit(X=X_train_all, y=y_train_all, intermediate_quantiles=interm_quantiles_all,
                       interm_lvl=interm_lvl, intermediate_q_feature=gbex_params$intermediate_q_feature, scale_features=gbex_params$scale_features,
                       B=gbex_params$B, lambda=gbex_params$lambda, lambda_ratio=gbex_params$lambda_ratio,
                       lambda_scale=gbex_params$lambda_scale, depth=gbex_params$depth, sf=gbex_params$sf)
  pred_gbex <- gbex_predict(fit_gbex, X_test, to_predict=prob_lvls_predict, intermediate_quantiles=pred_interm, interm_lvl=interm_lvl)
  
  fit_egam <- fit_gpd_gam(X=X_train_all, y=y_train_all, intermediate_quantiles=interm_quantiles_all, interm_lvl=interm_lvl,
                          model_shape=egam_params$model_shape,
                          intermediate_q_feature=egam_params$intermediate_q_feature, scale_features=egam_params$scale_features)
  pred_egam <- predict_gpd_gam(fit_egam, X_test, to_predict=prob_lvls_predict,
                               intermediate_quantiles=pred_interm, interm_lvl=interm_lvl)
  ## ===========================
  
  train_plot <- training_plot_eqrn(fit_eqrn, uncond_losses_interm, show_legend=TRUE)
  valid_plot <- validation_plot_eqrn(fit_eqrn, uncond_losses_interm, show_legend=FALSE)
  plot(train_plot)
  plot(valid_plot)
  if(!is.null(save_path)){
    save_myplot(plt=train_plot, plt_nm=paste0(save_path, params_string, "_training.pdf"),
                width = 75, height = 75, cairo=FALSE)
    save_myplot(plt=valid_plot, plt_nm=paste0(save_path, params_string, "_validation.pdf"),
                width = 75, height = 75, cairo=FALSE)
  }
  
  #Final EQRN predictions on X_test
  pred_eqrn <- EQRN_predict(fit_eqrn, X_test, prob_lvls_predict, pred_interm, interm_lvl)
  
  # Compute losses for desired predicted quantiles
  MSE_losses <- multilevel_MSE(pred_true, pred_eqrn, prob_lvls_predict, prefix="test_", give_names=TRUE)
  MAE_losses <- multilevel_MAE(pred_true, pred_eqrn, prob_lvls_predict, prefix="test_", give_names=TRUE)
  
  nb_prob_lvls_predict <- length(prob_lvls_predict)
  for(i in 1:nb_prob_lvls_predict){
    if(effective_dim_model(model)<2){
      pred_plot <- plot_predictions_1D(pred_eqrn[,i], pred_grf_test[,i], pred_true[,i], X_test, prob_lvls_predict[i])
      plot(pred_plot)
      if(!is.null(save_path)){
        ggsave(paste0(params_string, "_q",prob_lvls_predict[i]*10000,"_predX.pdf"), plot=pred_plot, device="pdf",
               path=save_path, width=300, height=200, units="mm", dpi=300)
      }
    }
    
    # resid_box <- residuals_boxplot_eqrn(pred_eqrn[,i], pred_grf_test[,i], pred_true[,i], prob_lvls_predict[i])
    resid_box2 <- residuals_boxplot_eqrn2(pred_eqrn[,i], pred_grf_test[,i], pred_gbex[,i], pred_true[,i], prob_lvls_predict[i])
    plot(resid_box2)
    if(!is.null(save_path)){
      ggsave(paste0(params_string, "_q",prob_lvls_predict[i]*10000,"_residbox2.pdf"), plot=resid_box2, device="pdf",
             path=save_path, width=200, height=150, units="mm", dpi=300)
    }
  }
  pred_eqrn_params <- EQRN_predict_params(fit_eqrn, X_test, pred_interm, return_parametrization="classical", interm_lvl)
  test_loss <- compute_EQRN_GPDLoss(fit_eqrn, X_test, y_test, intermediate_quantiles=pred_interm, interm_lvl=interm_lvl)
  RESULTS <- list(pred_true=pred_true, pred_eqrn=pred_eqrn, pred_grf_test=pred_grf_test, prob_lvls_predict=prob_lvls_predict,
                  train_loss=fit_eqrn$train_loss, valid_loss=fit_eqrn$valid_loss, test_loss=test_loss,
                  uncond_losses_interm=uncond_losses_interm, X_test=X_test[,1:4], params_string=params_string)
  # if(!is.null(save_path)){
  #   safe_save_rds(RESULTS, paste0(save_path, "Results_",params_string,".rds"))
  # }
  
  ## ===================== predictions vs truth plots ============================
  pred_vs_true_plot_grf <- plot_predictions_vs_truth(pred_eqrn, pred_grf_test, pred_semicond$predictions,
                                                     pred_true, prob_lvls_predict, name_other="GRF")
  plot(pred_vs_true_plot_grf)
  if(!is.null(save_path)){
    ggsave(paste0(params_string, "_predVStruth_grf.pdf"), plot=pred_vs_true_plot_grf, device="pdf",
           path=save_path, width=300, height=100, units="mm", dpi=300)
  }
  pred_vs_true_plot_gbex <- plot_predictions_vs_truth(pred_eqrn, pred_gbex, pred_semicond$predictions,
                                                      pred_true, prob_lvls_predict, name_other="gbex")
  plot(pred_vs_true_plot_gbex)
  if(!is.null(save_path)){
    ggsave(paste0(params_string, "_predVStruth_gbex.pdf"), plot=pred_vs_true_plot_gbex, device="pdf",
           path=save_path, width=300, height=100, units="mm", dpi=300)
  }
  
  #quant_pred_plot <- seq(interm_lvl,0.9995,length.out=100)
  quant_pred_plot <- c(0.8, 0.85, 0.9, 0.95, 0.99, 0.995, 0.999, 0.9995)
  quant_pred_plot_smooth <- c(0.8,1-10^c(-1,-1.5,-2,-2.5,-3,-3.5,-4))
  err_quant_plt <- plot_error_quantile(fit_eqrn, fit_grf, fit_gbex, fit_egam, dat$Y, X_test, quant_pred_plot, intermediate_method,
                                       interm_lvl, intermediate_quantiles, model, distr, df, factorize=TRUE, test_data=test_data)
  plot(err_quant_plt)
  if(!is.null(save_path)){
    ggsave(paste0(params_string, "_errquant.pdf"), plot=err_quant_plt, device="pdf",
           path=save_path, width=200, height=100, units="mm", dpi=300)
  }
  for (Q_pred_levels in list(list(qs=quant_pred_plot, nm="dq"),
                             list(qs=quant_pred_plot_smooth, nm="sq"))){
    
    rmse_quant_plt <- plot_rmse_quantile(fit_eqrn, fit_grf, fit_gbex, fit_egam, dat$Y, X_test, Q_pred_levels$qs, intermediate_method,
                                         interm_lvl, intermediate_quantiles, model, distr, df, factorize=TRUE, test_data=test_data,
                                         legend.position="bottom", crop_obs=0, bias_variance=TRUE, crop_Rbv=c(0,0,0))
    rmse_quant_plt_c <- plot_rmse_quantile(fit_eqrn, fit_grf, fit_gbex, fit_egam, dat$Y, X_test, Q_pred_levels$qs, intermediate_method,
                                           interm_lvl, intermediate_quantiles, model, distr, df, factorize=TRUE, test_data=test_data,
                                           legend.position="bottom", crop_obs=1, bias_variance=TRUE, crop_Rbv=c(2,0,1))
    plot(rmse_quant_plt_c$rmse_plot)
    if(!is.null(save_path)){
      save_myplot(plt=rmse_quant_plt$rmse_plot, plt_nm=paste0(save_path, params_string, "_rmsequant_p_",Q_pred_levels$nm,".pdf"),
                  width = 60, height = 60, cairo=FALSE)
      save_myplot(plt=rmse_quant_plt_c$rmse_plot, plt_nm=paste0(save_path, params_string, "_rmsequant_p_",Q_pred_levels$nm,"_c.pdf"),
                  width = 60, height = 60, cairo=FALSE)
      save_myplot(plt=rmse_quant_plt$Rbv_plot, plt_nm=paste0(save_path, params_string, "_Rbvquant_p_",Q_pred_levels$nm,".pdf"),
                  width = 210, height = 82, cairo=FALSE)
      save_myplot(plt=rmse_quant_plt_c$Rbv_plot, plt_nm=paste0(save_path, params_string, "_Rbvquant_p_",Q_pred_levels$nm,"_c.pdf"),
                  width = 210, height = 82, cairo=FALSE)
    }
  }
  
  
  # Gridplot/heatmap of predictions for "2D" data models
  if(effective_dim_model(model)==2){
    # If test data was not a grid, we need it for 2D predictions gridplot
    if (test_data != "grid"){
      ntest_2 <- floor(sqrt(ntest)) ** 2
      cat(paste0("Modified ntest to: ", ntest_2,"\n"))
      X_test_2 <- generate_test_data(ntest_2, p, "grid")
      
      #Predict intermediate quantiles on X_test_2
      if(intermediate_method=="grf"){
        pred_interm_2 <- predict(fit_grf, newdata = X_test_2, quantiles = c(interm_lvl))$predictions
      }else if(intermediate_method=="oracle"){
        pred_interm_2 <- generate_theoretical_quantiles(quantiles = c(interm_lvl), X = X_test_2, model = model, distr = distr, df = df)
      }
      #High quantile prediction with GRF
      pred_grf_test_2 <- predict(fit_grf, newdata = X_test_2, quantiles = prob_lvls_predict)$predictions
      # UNCONDITIONAL predicted quantile(s) (Y quantile on X_train)
      pred_unc_2 <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = prob_lvls_predict, Y = dat$Y, ntest = ntest_2)
      # SEMI-CONDITIONAL predicted quantiles
      pred_semicond_2 <- predict_GPD_semiconditional(Y=dat$Y, interm_lvl=interm_lvl, thresh_quantiles=intermediate_quantiles,
                                                     interm_quantiles_test=pred_interm_2, prob_lvls_predict=prob_lvls_predict)
      # GROUND-TRUTH (y_test_2)
      pred_true_2 <- generate_theoretical_quantiles(quantiles = prob_lvls_predict, X = X_test_2, model = model, distr = distr, df = df)
      # EQRN predictions on X_test_2
      pred_eqrn_2 <- EQRN_predict(fit_eqrn, X_test_2, prob_lvls_predict, pred_interm_2, interm_lvl)
      # gbex and gpdGAM predictions on X_test_2
      pred_gbex_2 <- gbex_predict(fit_gbex, X_test_2, to_predict=prob_lvls_predict, intermediate_quantiles=pred_interm_2, interm_lvl=interm_lvl)
      pred_egam_2 <- predict_gpd_gam(fit_egam, X_test_2, to_predict=prob_lvls_predict, intermediate_quantiles=pred_interm_2, interm_lvl=interm_lvl)
    } else {
      X_test_2 <- X_test
      pred_interm_2 <- pred_interm
      pred_grf_test_2 <- pred_grf_test
      pred_unc_2 <- pred_unc
      pred_semicond_2 <- pred_semicond
      pred_true_2 <- pred_true
      pred_eqrn_2 <- pred_eqrn
      pred_gbex_2 <- pred_gbex
      pred_egam_2 <- pred_egam
    }
    
    for(i in 1:nb_prob_lvls_predict){
      pred2D <- plot_predictions_2D(pred_eqrn_2[,i], pred_gbex_2[,i], pred_grf_test_2[,i], pred_egam_2[,i],
                                    pred_semicond_2$predictions[,i], pred_true_2[,i], X_test_2, prob_lvls_predict[i])
      plot(pred2D)
      if(!is.null(save_path)){
        # ggsave(paste0(params_string, "_q",prob_lvls_predict[i]*10000,"_pred2D.pdf"), plot=pred2D, device="pdf",
        #        path=save_path, width=200, height=140, units="mm", dpi=300)
        save_myplot(plt=pred2D, plt_nm=paste0(save_path, params_string,"_q",prob_lvls_predict[i]*10000,"_pred2D.pdf"),
                    width = 40, height = 40, cairo=FALSE)
      }
      pred2D_ex <- plot_predictions_2D(pred_eqrn_2[,i]-pred_interm_2, pred_gbex_2[,i]-pred_interm_2,
                                       pred_grf_test_2[,i]-pred_interm_2, pred_egam_2[,i]-pred_interm_2,
                                       pred_semicond_2$predictions[,i]-pred_interm_2, pred_true_2[,i]-pred_interm_2,
                                       X_test_2, prob_lvls_predict[i])
      plot(pred2D_ex)
      if(!is.null(save_path)){
        # ggsave(paste0(params_string, "_q",prob_lvls_predict[i]*10000,"_pred2D_ex.pdf"), plot=pred2D_ex, device="pdf",
        #        path=save_path, width=200, height=140, units="mm", dpi=300)
        save_myplot(plt=pred2D_ex, plt_nm=paste0(save_path, params_string,"_q",prob_lvls_predict[i]*10000,"_pred2D_ex.pdf"),
                    width = 40, height = 40, cairo=FALSE)
      }
    }
  }
  as.character(Sys.time()-start_time)
}

end_doFuture_strategy()

cat("\nRun times:\n")
print(cbind(param_file_names, fe_out))


