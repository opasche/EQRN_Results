
plot_predictions_ts <- function(pred_q, pred_q_other, true_q, y=NaN, prob_lvl_predicted, logscale=FALSE, name_other="GRF"){
  preds_df <- data.frame(t=seq_along(true_q), Other=pred_q_other, EQRN=pred_q, Truth=true_q)
  preds_df <- preds_df %>% rename_with(~str_replace(., "Other", name_other))
  pred_plot <- preds_df %>% tidyr::gather(key="Method", value="Prediction", .data[[name_other]], Truth, EQRN, factor_key=TRUE) %>%
    ggplot(aes(x=t, y=Prediction, group=Method, color=Method)) +
    geom_line() + geom_point(data=data.frame(t=seq_along(true_q), y=y, Method=factor("observations")), aes(x=t, y=y, group=Method, color=Method)) +
    scale_color_manual(values=alpha(c(rgb(57/255,106/255,177/255), "#69b3a2", "black", "black"),c(1,0.5,0.2,0.4))) +
    labs(title=paste0("Predictions (q=",prob_lvl_predicted,")"), x="time", y="Predictions", color=NULL)
  if(logscale){pred_plot <- pred_plot + scale_y_log10()} 
  return(pred_plot)
}

plot_predictions_diff_ts <- function(pred_q, pred_q_other, true_q, prob_lvl_predicted, logscale=FALSE, name_other="GRF"){
  preds_df <- data.frame(t=seq_along(true_q), Other=pred_q_other-true_q, EQRN=pred_q-true_q)
  preds_df <- preds_df %>% rename_with(~str_replace(., "Other", name_other))
  pred_plot <- preds_df %>% tidyr::gather(key="Method", value="Prediction", .data[[name_other]], EQRN, factor_key=TRUE) %>%
    ggplot(aes(x=t, y=Prediction, group=Method, color=Method)) +
    geom_line() + scale_color_manual(values=alpha(c("#69b3a2",rgb(57/255,106/255,177/255)),c(0.5,1))) +
    labs(title=paste0("Predictions (q=",prob_lvl_predicted,")"), x="time", y="Predictions", color=NULL)
  if(logscale){pred_plot <- pred_plot + scale_y_log10()} 
  return(pred_plot)
}


plot_rmse_quantile_ts <- function(fit_eqrn, fit_grf, fit_gbex, fit_egam, Y_train, X_test, y_test, prob_lvls_predict, pred_interm,
                                  interm_lvl, interm_quant_train, true_quantiles_test, factorize=TRUE, test_data="other", pred_interm_c=pred_interm,
                                  legend.position="bottom", crop_obs=0, bias_variance=FALSE, crop_Rbv=c(0,0,0)){
  ntest <- nrow(X_test)
  nb_prob_lvls_predict <- length(prob_lvls_predict)
  #lagged datasets for competitor methods
  lagged_test <- lagged_features(X=cbind(y_test,X_test), max_lag=fit_eqrn$seq_len, drop_present=TRUE)
  laged_interm_q_test <- pred_interm_c[(fit_eqrn$seq_len+1):length(pred_interm_c), , drop=F]
  #High quantile prediction with GRF
  pred_grf_test <- predict(fit_grf, newdata=lagged_test, quantiles = prob_lvls_predict)$predictions
  # UNCONDITIONAL predicted quantile(s) (Y quantile on X_train)
  pred_unc <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = prob_lvls_predict, Y = Y_train, ntest = nrow(true_quantiles_test))
  # SEMI-CONDITIONAL predicted quantiles
  pred_semicond <- predict_GPD_semiconditional(Y=Y_train[(seq_len+1):length(Y_train)], interm_lvl=interm_lvl,
                                               thresh_quantiles=interm_quant_train[(seq_len+1):length(interm_quant_train)],
                                               interm_quantiles_test=pred_interm[(seq_len+1):length(pred_interm)], prob_lvls_predict=prob_lvls_predict)
  # GROUND-TRUTH (y_test)
  pred_true <- true_quantiles_test
  # Prediction GBEX
  pred_gbex <- gbex_predict(fit_gbex, lagged_test, to_predict=prob_lvls_predict, intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
  pred_exgam <- predict_gpd_gam(fit_egam, lagged_test, to_predict=prob_lvls_predict,
                                intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
  
  #Final EQRN predictions on X_test
  pred_eqrnn <- EQRN_predict_seq(fit_eqrn, X_test, y_test, prob_lvls_predict, pred_interm, interm_lvl, crop_predictions=TRUE)
  
  # Compute losses for desired predicted quantiles
  RMSEs_eqrnn <- sqrt(multilevel_MSE(pred_true,pred_eqrnn,prob_lvls_predict,give_names=FALSE))
  RMSEs_grf <- sqrt(multilevel_MSE(pred_true,pred_grf_test,prob_lvls_predict,give_names=FALSE))
  RMSEs_unc <- sqrt(multilevel_MSE(pred_true,pred_unc$predictions,prob_lvls_predict,give_names=FALSE))
  RMSEs_semicond <- sqrt(multilevel_MSE(pred_true,pred_semicond$predictions,prob_lvls_predict,give_names=FALSE))
  RMSEs_gbex <- sqrt(multilevel_MSE(pred_true,pred_gbex,prob_lvls_predict,give_names=FALSE))
  RMSEs_exgam <- sqrt(multilevel_MSE(pred_true,pred_exgam,prob_lvls_predict,give_names=FALSE))
  
  Q_lvls <- if(factorize){as_factor(roundm(prob_lvls_predict,4))}else{prob_lvls_predict}
  
  met_names <- c("Uncond", "Semi-cond", "GRF", "EGAM", "GBEX", "EQRN")
  df_rmse <- data.frame(Quantile=Q_lvls, Uncond=RMSEs_unc, Semi_cond=RMSEs_semicond,
                        GRF=RMSEs_grf, EGAM=RMSEs_exgam, GBEX=RMSEs_gbex, EQRN=RMSEs_eqrnn) %>% 
    rename_with(~str_replace(., "Semi_cond", "Semi-cond")) %>% 
    tidyr::gather(key="Model", value="Error", all_of(met_names), factor_key=TRUE)
  
  y_lab_rmse <- if(test_data=="halton"){"RISE"}else{"RMSE"}
  
  color_scale <- my_palette_methods[met_names]
  linetypes <- c("dotted","dashed","dotdash","longdash","twodash","solid")
  
  mse_plot <- df_rmse %>% ggplot( aes(x=Quantile, y=Error, group=Model, color=Model, linetype=Model)) +
    geom_line(size=1) + scale_linetype_manual(values=linetypes) + scale_color_manual(values=color_scale) +
    scale_y_continuous(expand=c(0.02,0)) + labs(title=NULL, x="Probability level", y=y_lab_rmse, color=NULL, linetype=NULL) +
    theme(legend.position=legend.position)
  
  mse_plot <- mse_plot + if(factorize){scale_x_discrete(expand=c(0.02,0))}else{scale_x_continuous(expand=c(0.02,0))}
  if(crop_obs>0){mse_plot <- mse_plot + coord_cartesian(ylim=c(min(df_rmse$Error), sort(df_rmse$Error, decreasing=TRUE)[[crop_obs+1]]))}
  
  if(!bias_variance){
    return(mse_plot)
    
  }else{
    # Compute metrics for desired predicted quantiles
    bias_eqrnn <- multilevel_pred_bias(pred_true, pred_eqrnn, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_grf <- multilevel_pred_bias(pred_true, pred_grf_test, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_unc <- multilevel_pred_bias(pred_true, pred_unc$predictions, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_semicond <- multilevel_pred_bias(pred_true, pred_semicond$predictions, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_gbex <- multilevel_pred_bias(pred_true, pred_gbex, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_exgam <- multilevel_pred_bias(pred_true, pred_exgam, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    
    rstd_eqrnn <- sqrt(multilevel_resid_var(pred_true, pred_eqrnn, prob_lvls_predict, give_names=FALSE))
    rstd_grf <- sqrt(multilevel_resid_var(pred_true, pred_grf_test, prob_lvls_predict, give_names=FALSE))
    rstd_unc <- sqrt(multilevel_resid_var(pred_true, pred_unc$predictions, prob_lvls_predict, give_names=FALSE))
    rstd_semicond <- sqrt(multilevel_resid_var(pred_true, pred_semicond$predictions, prob_lvls_predict, give_names=FALSE))
    rstd_gbex <- sqrt(multilevel_resid_var(pred_true, pred_gbex, prob_lvls_predict, give_names=FALSE))
    rstd_exgam <- sqrt(multilevel_resid_var(pred_true, pred_exgam, prob_lvls_predict, give_names=FALSE))
    
    Rsqr_eqrnn <- multilevel_R_squared(pred_true, pred_eqrnn, prob_lvls_predict, give_names=FALSE)
    Rsqr_grf <- multilevel_R_squared(pred_true, pred_grf_test, prob_lvls_predict, give_names=FALSE)
    Rsqr_unc <- multilevel_R_squared(pred_true, pred_unc$predictions, prob_lvls_predict, give_names=FALSE)
    Rsqr_semicond <- multilevel_R_squared(pred_true, pred_semicond$predictions, prob_lvls_predict, give_names=FALSE)
    Rsqr_gbex <- multilevel_R_squared(pred_true, pred_gbex, prob_lvls_predict, give_names=FALSE)
    Rsqr_exgam <- multilevel_R_squared(pred_true, pred_exgam, prob_lvls_predict, give_names=FALSE)
    
    Rbv_names <- c("Quantile R squared", "Bias", "Residual standard deviation")
    
    df_Rbv <- data.frame(Quantile=rep(Q_lvls, 3),
                         Uncond=c(Rsqr_unc,bias_unc,rstd_unc),
                         Semi_cond=c(Rsqr_semicond,bias_semicond,rstd_semicond),
                         GRF=c(Rsqr_grf,bias_grf,rstd_grf),
                         EGAM=c(Rsqr_exgam,bias_exgam,rstd_exgam),
                         GBEX=c(Rsqr_gbex,bias_gbex,rstd_gbex),
                         EQRN=c(Rsqr_eqrnn,bias_eqrnn,rstd_eqrnn),
                         metric=factor(rep(Rbv_names, each=nb_prob_lvls_predict), levels=Rbv_names)) %>%
      rename_with(~str_replace(., "Semi_cond", "Semi-cond")) %>%
      tidyr::gather(key="Model", value="Error", all_of(met_names), factor_key=TRUE)
    
    Rbv_lims <- list(c(max(min(sort(df_Rbv$Error[df_Rbv$metric==Rbv_names[1]], decreasing=FALSE)[[crop_Rbv[1]+1]],0),-5.), 1.),
                     (c(-1,1) * sort(abs(df_Rbv$Error[df_Rbv$metric==Rbv_names[2]]), decreasing=TRUE)[[crop_Rbv[2]+1]]),
                     c(min(df_Rbv$Error[df_Rbv$metric==Rbv_names[3]]),
                       sort(df_Rbv$Error[df_Rbv$metric==Rbv_names[3]], decreasing=TRUE)[[crop_Rbv[3]+1]]))
    Rbv_plots <- list()
    for(i in seq_along(Rbv_names)){
      Rbv_plots[[i]] <- df_Rbv %>% filter(metric==Rbv_names[i]) %>%
        ggplot( aes(x=Quantile, y=Error, group=Model, color=Model, linetype=Model)) +
        geom_line(size=1) + scale_linetype_manual(values=linetypes) + scale_color_manual(values=color_scale) +
        scale_y_continuous(expand=c(0.02,0)) +
        labs(title=NULL, x=if(i==2){"Probability level"}else{NULL}, y=Rbv_names[i], color=NULL, linetype=NULL) +
        theme(legend.position=legend.position) + coord_cartesian(ylim=Rbv_lims[[i]]) +
        if(factorize){scale_x_discrete(expand=c(0.02,0))}else{scale_x_continuous(expand=c(0.02,0))}
    }
    #facet_wrap(~metric, scales="free_y", ncol=3, strip.position="left") + 
    Rbv_plot <- ggpubr::ggarrange(plotlist=Rbv_plots, nrow=1, ncol=3, labels=NULL,
                                  common.legend=TRUE, legend="bottom", align="hv")
    
    return(list(rmse_plot=mse_plot, Rbv_plot=Rbv_plot))
  }
}

plot_rmse_quantile_comp_ts <- function(fit_eqrn1, fit_eqrn2, fit_grf, fit_exqar, fit_gbex, fit_egam, Y_train, X_test, y_test, prob_lvls_predict, pred_interm,
                                       interm_lvl, interm_quant_train, true_quantiles_test, factorize=TRUE, test_data="other", pred_interm_c=pred_interm,
                                       names_EQRN=c("EQRN_1","EQRN_2"), legend.position="bottom", crop_obs=0, bias_variance=FALSE, crop_Rbv=c(0,0,0)){
  
  ntest <- nrow(X_test)
  nb_prob_lvls_predict <- length(prob_lvls_predict)
  if(fit_eqrn1$seq_len!=fit_eqrn2$seq_len){
    stop("Unfair comparison due to different eqrnn 'seq_len' in 'plot_error_quantile_comp_ts'.")
  }
  #lagged datasets for competitor methods
  lagged_test <- lagged_features(X=cbind(y_test,X_test), max_lag=fit_eqrn1$seq_len, drop_present=TRUE)
  laged_interm_q_test <- pred_interm_c[(fit_eqrn1$seq_len+1):length(pred_interm_c), , drop=F]
  #High quantile prediction with GRF
  pred_grf_test <- predict(fit_grf, newdata=lagged_test, quantiles = prob_lvls_predict)$predictions
  # UNCONDITIONAL predicted quantile(s) (Y quantile on X_train)
  pred_unc <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = prob_lvls_predict, Y = Y_train, ntest = nrow(true_quantiles_test))
  # SEMI-CONDITIONAL predicted quantiles
  pred_semicond <- predict_GPD_semiconditional(Y=Y_train[(seq_len+1):length(Y_train)], interm_lvl=interm_lvl,
                                               thresh_quantiles=interm_quant_train[(seq_len+1):length(interm_quant_train)],
                                               interm_quantiles_test=pred_interm[(seq_len+1):length(pred_interm)], prob_lvls_predict=prob_lvls_predict)
  # GROUND-TRUTH (y_test)
  pred_true <- true_quantiles_test
  # Prediction competitors
  pred_exqar <- EXQAR_predict(fit_exqar, y_test, X=X_test, prob_lvls_predict, tol=1e-4, min_prop=0.3, return_infos=FALSE)
  pred_gbex <- gbex_predict(fit_gbex, lagged_test, to_predict=prob_lvls_predict, intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
  pred_exgam <- predict_gpd_gam(fit_egam, lagged_test, to_predict=prob_lvls_predict,
                                intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
  
  #Final EQRN predictions on X_test
  pred_eqrnn1 <- EQRN_predict_seq(fit_eqrn1, X_test, y_test, prob_lvls_predict, pred_interm, interm_lvl, crop_predictions=TRUE)
  pred_eqrnn2 <- EQRN_predict_seq(fit_eqrn2, X_test, y_test, prob_lvls_predict, pred_interm, interm_lvl, crop_predictions=TRUE)
  
  # Compute losses for desired predicted quantiles
  RMSEs_eqrnn1 <- sqrt(multilevel_MSE(pred_true,pred_eqrnn1,prob_lvls_predict,give_names=FALSE))
  RMSEs_eqrnn2 <- sqrt(multilevel_MSE(pred_true,pred_eqrnn2,prob_lvls_predict,give_names=FALSE))
  RMSEs_grf <- sqrt(multilevel_MSE(pred_true,pred_grf_test,prob_lvls_predict,give_names=FALSE))
  RMSEs_unc <- sqrt(multilevel_MSE(pred_true,pred_unc$predictions,prob_lvls_predict,give_names=FALSE))
  RMSEs_semicond <- sqrt(multilevel_MSE(pred_true,pred_semicond$predictions,prob_lvls_predict,give_names=FALSE))
  RMSEs_exqar <- sqrt(multilevel_MSE(pred_true,pred_exqar,prob_lvls_predict,give_names=FALSE))
  RMSEs_gbex <- sqrt(multilevel_MSE(pred_true,pred_gbex,prob_lvls_predict,give_names=FALSE))
  RMSEs_exgam <- sqrt(multilevel_MSE(pred_true,pred_exgam,prob_lvls_predict,give_names=FALSE))
  
  Q_lvls <- if(factorize){as_factor(roundm(prob_lvls_predict,4))}else{prob_lvls_predict}
  
  met_names <- c("Uncond", "Semi-cond", "GRF", "EXQAR", "EGAM", "GBEX", names_EQRN[2], names_EQRN[1])
  
  df_rmse <- data.frame(Quantile=Q_lvls, Uncond=RMSEs_unc, Semi_cond=RMSEs_semicond,
                        GRF=RMSEs_grf, EXQAR=RMSEs_exqar, EGAM=RMSEs_exgam, GBEX=RMSEs_gbex,
                        EQRN2=RMSEs_eqrnn2, EQRN1=RMSEs_eqrnn1) %>% 
    rename_with(~str_replace(., "Semi_cond", "Semi-cond")) %>%
    rename_with(~str_replace(., "EQRN2", names_EQRN[2])) %>%
    rename_with(~str_replace(., "EQRN1", names_EQRN[1])) %>% 
    tidyr::gather(key="Model", value="Error", all_of(met_names), factor_key=TRUE)
  
  y_lab_rmse <- if(test_data=="halton"){"RISE"}else{"RMSE"}
  
  linetypes <- c("dotted","dashed","dotdash","dotdash","longdash","twodash","dashed","solid")
  color_scale <- my_palette_methods[c("Uncond", "Semi-cond", "GRF", "EXQAR", "EGAM", "GBEX", "EQRN2", "EQRN")]
  names(linetypes) <- met_names
  names(color_scale) <- met_names
  
  mse_plot <- df_rmse %>% ggplot( aes(x=Quantile, y=Error, group=Model, color=Model, linetype=Model)) +
    geom_line(size=1) + scale_linetype_manual(values=linetypes) + scale_color_manual(values=color_scale) + 
    scale_y_continuous(expand=c(0.02,0)) + labs(title=NULL, x="Probability level", y=y_lab_rmse, color=NULL, linetype=NULL) +
    theme(legend.position=legend.position)
  
  mse_plot <- mse_plot + if(factorize){scale_x_discrete(expand=c(0.02,0))}else{scale_x_continuous(expand=c(0.02,0))}
  if(crop_obs>0){mse_plot <- mse_plot + coord_cartesian(ylim=c(min(df_rmse$Error), sort(df_rmse$Error, decreasing=TRUE)[[crop_obs+1]]))}
  
  if(!bias_variance){
    return(mse_plot)
    
  }else{
    # Compute metrics for desired predicted quantiles
    bias_eqrnn1 <- multilevel_pred_bias(pred_true, pred_eqrnn1, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_eqrnn2 <- multilevel_pred_bias(pred_true, pred_eqrnn2, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_grf <- multilevel_pred_bias(pred_true, pred_grf_test, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_unc <- multilevel_pred_bias(pred_true, pred_unc$predictions, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_semicond <- multilevel_pred_bias(pred_true, pred_semicond$predictions, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_exqar <- multilevel_pred_bias(pred_true, pred_exqar, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_gbex <- multilevel_pred_bias(pred_true, pred_gbex, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_exgam <- multilevel_pred_bias(pred_true, pred_exgam, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    
    rstd_eqrnn1 <- sqrt(multilevel_resid_var(pred_true, pred_eqrnn1, prob_lvls_predict, give_names=FALSE))
    rstd_eqrnn2 <- sqrt(multilevel_resid_var(pred_true, pred_eqrnn2, prob_lvls_predict, give_names=FALSE))
    rstd_grf <- sqrt(multilevel_resid_var(pred_true, pred_grf_test, prob_lvls_predict, give_names=FALSE))
    rstd_unc <- sqrt(multilevel_resid_var(pred_true, pred_unc$predictions, prob_lvls_predict, give_names=FALSE))
    rstd_semicond <- sqrt(multilevel_resid_var(pred_true, pred_semicond$predictions, prob_lvls_predict, give_names=FALSE))
    rstd_exqar <- sqrt(multilevel_resid_var(pred_true, pred_exqar, prob_lvls_predict, give_names=FALSE))
    rstd_gbex <- sqrt(multilevel_resid_var(pred_true, pred_gbex, prob_lvls_predict, give_names=FALSE))
    rstd_exgam <- sqrt(multilevel_resid_var(pred_true, pred_exgam, prob_lvls_predict, give_names=FALSE))
    
    Rsqr_eqrnn1 <- multilevel_R_squared(pred_true, pred_eqrnn1, prob_lvls_predict, give_names=FALSE)
    Rsqr_eqrnn2 <- multilevel_R_squared(pred_true, pred_eqrnn2, prob_lvls_predict, give_names=FALSE)
    Rsqr_grf <- multilevel_R_squared(pred_true, pred_grf_test, prob_lvls_predict, give_names=FALSE)
    Rsqr_unc <- multilevel_R_squared(pred_true, pred_unc$predictions, prob_lvls_predict, give_names=FALSE)
    Rsqr_semicond <- multilevel_R_squared(pred_true, pred_semicond$predictions, prob_lvls_predict, give_names=FALSE)
    Rsqr_exqar <- multilevel_R_squared(pred_true, pred_exqar, prob_lvls_predict, give_names=FALSE)
    Rsqr_gbex <- multilevel_R_squared(pred_true, pred_gbex, prob_lvls_predict, give_names=FALSE)
    Rsqr_exgam <- multilevel_R_squared(pred_true, pred_exgam, prob_lvls_predict, give_names=FALSE)
    
    Rbv_names <- c("Quantile R squared", "Bias", "Residual standard deviation")
    
    df_Rbv <- data.frame(Quantile=rep(Q_lvls, 3),
                         Uncond=c(Rsqr_unc,bias_unc,rstd_unc),
                         Semi_cond=c(Rsqr_semicond,bias_semicond,rstd_semicond),
                         GRF=c(Rsqr_grf,bias_grf,rstd_grf),
                         EXQAR=c(Rsqr_exqar,bias_exqar,rstd_exqar),
                         EGAM=c(Rsqr_exgam,bias_exgam,rstd_exgam),
                         GBEX=c(Rsqr_gbex,bias_gbex,rstd_gbex),
                         EQRN2=c(Rsqr_eqrnn2,bias_eqrnn2,rstd_eqrnn2),
                         EQRN1=c(Rsqr_eqrnn1,bias_eqrnn1,rstd_eqrnn1),
                         metric=factor(rep(Rbv_names, each=nb_prob_lvls_predict), levels=Rbv_names)) %>%
      rename_with(~str_replace(., "Semi_cond", "Semi-cond")) %>%
      rename_with(~str_replace(., "EQRN2", names_EQRN[2])) %>%
      rename_with(~str_replace(., "EQRN1", names_EQRN[1])) %>%
      tidyr::gather(key="Model", value="Error", all_of(met_names), factor_key=TRUE)
    
    Rbv_lims <- list(c(max(min(sort(df_Rbv$Error[df_Rbv$metric==Rbv_names[1]], decreasing=FALSE)[[crop_Rbv[1]+1]],0),-5.), 1.),
                     (c(-1,1) * sort(abs(df_Rbv$Error[df_Rbv$metric==Rbv_names[2]]), decreasing=TRUE)[[crop_Rbv[2]+1]]),
                     c(min(df_Rbv$Error[df_Rbv$metric==Rbv_names[3]]),
                       sort(df_Rbv$Error[df_Rbv$metric==Rbv_names[3]], decreasing=TRUE)[[crop_Rbv[3]+1]]))
    Rbv_plots <- list()
    for(i in seq_along(Rbv_names)){
      Rbv_plots[[i]] <- df_Rbv %>% filter(metric==Rbv_names[i]) %>%
        ggplot( aes(x=Quantile, y=Error, group=Model, color=Model, linetype=Model)) +
        geom_line(size=1) + scale_linetype_manual(values=linetypes) + scale_color_manual(values=color_scale) +
        scale_y_continuous(expand=c(0.02,0)) +
        labs(title=NULL, x=if(i==2){"Probability level"}else{NULL}, y=Rbv_names[i], color=NULL, linetype=NULL) +
        theme(legend.position=legend.position) + coord_cartesian(ylim=Rbv_lims[[i]]) +
        if(factorize){scale_x_discrete(expand=c(0.02,0))}else{scale_x_continuous(expand=c(0.02,0))}
    }
    Rbv_plot <- ggpubr::ggarrange(plotlist=Rbv_plots, nrow=1, ncol=3, labels=NULL,
                                  common.legend=TRUE, legend="bottom", align="hv")
    
    return(list(rmse_plot=mse_plot, Rbv_plot=Rbv_plot))
  }
}



plot_exceedence_proba_ts <- function(exceedence_proba, exceedence_level, Date_index=NULL,
                                     int_index=seq_along(exceedence_proba), type_probas=c("probability","years"),
                                     proba_ratio=FALSE, log_scaling=FALSE, date_breaks="3 month", date_minor_breaks="1 day",
                                     event_dates=NULL, p_lab=if(proba_ratio){"Probability ratio"}else{"Exceedence probability"},
                                     color_met=my_palette_methods[["EQRN"]]){
  type_probas <- match.arg(type_probas)
  if(proba_ratio & log_scaling){stop("Only one of 'proba_ratio' or 'log_scaling' is allowed in 'plot_exceedence_proba_ts'.")}
  ep <- as.numeric(exceedence_proba)
  exceedence_thresh <- if(type_probas=="probability"){1-exceedence_level}else{exceedence_level}
  logp1_ticks <- c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)
  trans_logp1 <- scales::trans_new(name="logp1", transform=function(x){log(x*100+0.01)}, inverse=function(x){(exp(x)-0.01)/100}, domain=c(-0,Inf),
                                   breaks=function(a)logp1_ticks, format=scales::number_format(accuracy=0.0001))
  if(proba_ratio){
    ep <- ep/exceedence_thresh
  }
  lims <- c(0, max(ep,exceedence_thresh))
  if(!is.null(Date_index)){
    df <- data.frame(id=int_index, Date=Date_index, x_ind=Date_index, Proba=ep)
    ep_plot <- df %>% ggplot( aes(x=x_ind, y=Proba)) + geom_line(color=color_met) + labs(y=p_lab, x="Date") +
      scale_x_date(date_breaks=date_breaks, date_minor_breaks = date_minor_breaks, date_labels = "%Y %b", expand=c(0.01,0))
  }else{
    df <- data.frame(id=int_index, x_ind=int_index, Proba=ep)
    ep_plot <- df %>% ggplot( aes(x=x_ind, y=Proba)) + geom_line(color=color_met) + labs(y=p_lab, x="id")
  }
  ep_plot <- ep_plot + geom_hline(yintercept = as.numeric(exceedence_thresh), linetype=4, colour=my_palette$red) +
    geom_point(data=filter(df, x_ind %in% event_dates), shape=17, size=2, color=color_met) +
    expand_limits(y=lims) + theme(panel.grid.major.x=element_line(size=0.5), axis.ticks.x=element_line(size=0.5))
  ep_plot <- ep_plot + if(log_scaling){scale_y_continuous(trans=trans_logp1, labels=logp1_ticks, expand=c(0.02,0))}else{scale_y_continuous(expand=c(0.02,0))}
  if(!is.null(event_dates)){ep_plot <- ep_plot + geom_vline(xintercept=event_dates, linetype="dashed")}
  return(ep_plot)
}

plot_exceedence_proba_box <- function(exceedence_proba, exceedence_level, type_probas=c("probability","years"),
                                      log_scaling=FALSE, outlier_p=1){
  type_probas <- match.arg(type_probas)
  ep <- as.numeric(exceedence_proba)
  exceedence_thresh <- if(type_probas=="probability"){1-exceedence_level}else{exceedence_level}
  logp1_ticks <- c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)
  trans_logp1 <- scales::trans_new(name="logp1", transform=function(x){log(x*100+0.01)}, inverse=function(x){(exp(x)-0.01)/100}, domain=c(-0,Inf),
                                   breaks=function(a)logp1_ticks, format=scales::number_format(accuracy=0.0001))
  lims <- c(0, max(ep,exceedence_thresh))
  ep_plot <- data.frame(x="EQRN", Proba=ep) %>% ggplot( aes(x=x, y=Proba)) +
    geom_hline(yintercept = as.numeric(exceedence_thresh), linetype=4, colour=my_palette$red) + geom_violin(width=1, color="grey", alpha=0.5) +
    stat_boxplot(geom = "errorbar", width=0.5) + geom_boxplot(width=0.75, color="black", alpha=0.8)+
    stat_summary(fun=mean, geom="point", shape=4, size=3, color="black", fill="black") +
    labs(x=NULL, y="Exceedence probability (EQRN)") +
    theme(legend.position="none", panel.grid.major.y = element_blank()) + scale_x_discrete(breaks=NULL, labels=NULL) + expand_limits(y=lims)
  if(log_scaling){
    ep_plot <- ep_plot +
      scale_y_continuous(limits=c(lims[1], max(quantile(ep,outlier_p,na.rm=TRUE),exceedence_thresh)), trans=trans_logp1, labels=logp1_ticks) +
      coord_flip()
  }else{
    ep_plot <- ep_plot +
      scale_y_continuous(limits=c(lims[1], max(quantile(ep,outlier_p,na.rm=TRUE),exceedence_thresh))) + coord_flip()
  }
  return(ep_plot)
}

plot_pred_vs_return_lvl_ts <- function(preds, return_level, PoT_quantile, empirical_quantile,
                                       y_lab=latex2exp::TeX(r'(Discharge $\[m^3s^{-1}\]$)'), y_obs=NULL,
                                       Date_index=NULL, int_index=seq_along(preds), date_breaks="3 month", date_minor_breaks="1 day",
                                       event_dates=NULL, legend.position="bottom", color_met=my_palette_methods[["EQRN"]]){
  ablines <- data.frame(intercepts=c(return_level, PoT_quantile, empirical_quantile), slopes=c(0,0,0),
                        cols=c("Block-Maxima","PoT","Empirical"))
  lims <- c(min(preds, return_level, PoT_quantile, empirical_quantile),
            max(preds, return_level, PoT_quantile, empirical_quantile, y_obs))
  df <- data.frame(id=int_index, preds=preds, Method="Prediction", Met_y="Observations")
  if(!is.null(y_obs)){
    y_masked <- y_obs
    y_masked[y_masked<min(return_level, PoT_quantile, empirical_quantile)] <- NA
    df$y_ex <- y_obs
  }
  if(!is.null(Date_index)){
    df$Date <- Date_index
    df$x_ind <- Date_index
    ep_plot <- df %>% ggplot( aes(x=x_ind, y=preds, group=Method, color=Method, linetype=Method)) + geom_line() +
      labs(y=y_lab, x="Date", color=NULL, linetype=NULL) +
      scale_x_date(date_breaks=date_breaks, date_minor_breaks = date_minor_breaks, date_labels = "%Y %b", expand=c(0.01,0))
  }else{
    df$x_ind <- int_index
    ep_plot <- df %>% ggplot( aes(x=x_ind, y=preds, group=Method, color=Method, linetype=Method)) + geom_line() +
      labs(y=y_lab, x="id", color=NULL, linetype=NULL)
  }
  ep_plot <- ep_plot + scale_color_manual(values = c("Prediction"=color_met,"Block-Maxima"=my_palette$red,
                                                     "PoT"=my_palette_methods[["Uncond"]],"Empirical"=my_palette$light_blue)) + #, guide=guide_legend(reverse=TRUE)) +
    scale_linetype_manual(values=c("Prediction"="solid","Block-Maxima"="twodash","PoT"="dotted","Empirical"="dotted")) +
    geom_abline(data=ablines, aes(intercept=intercepts, slope=slopes, col=cols, linetype=cols), size=1, show.legend=FALSE) +
    geom_point(data=filter(df, x_ind %in% event_dates), shape=17, size=2) +
    expand_limits(y=lims) + scale_y_continuous(expand=c(0.02,0)) +
    theme(panel.grid.major.x=element_line(size=0.5), axis.ticks.x=element_line(size=0.5), legend.position=legend.position)
  if(!is.null(event_dates)){ep_plot <- ep_plot + geom_vline(xintercept=event_dates, linetype="dashed")}
  if(!is.null(y_obs)){ep_plot <- ep_plot + geom_point(aes(y=y_ex), show.legend=FALSE, color="darkslategrey", alpha=0.8)}
  return(ep_plot)
}


plot_pred_quant_risk_ts <- function(preds, exceedence_proba, exceedence_level, return_level, PoT_quantile, empirical_quantile,
                                    y_lab=latex2exp::TeX(r'(Discharge $\[m^3s^{-1}\]$)'), y_obs=NULL,
                                    Date_index=NULL, int_index=seq_along(preds), type_probas=c("probability","years"),
                                    proba_ratio=FALSE, date_breaks="3 month", date_minor_breaks="1 day",
                                    event_dates=NULL, legend.position="bottom", p_lab=if(proba_ratio){"'Probability ratio'"}else{"'Probability'"},
                                    log_scaling=FALSE, color_met=my_palette_methods[["EQRN"]]){
  type_probas <- match.arg(type_probas)
  ep <- as.numeric(exceedence_proba)
  exceedence_thresh <- if(type_probas=="probability"){1-exceedence_level}else{exceedence_level}
  logp1_ticks <- c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)
  trans_logp1 <- scales::trans_new(name="logp1", transform=function(x){log(x*100+0.01)}, inverse=function(x){(exp(x)-0.01)/100}, domain=c(-0,Inf),
                                   breaks=function(a)logp1_ticks, format=scales::number_format(accuracy=0.0001))
  ablines <- data.frame(intercepts=c(return_level, PoT_quantile, empirical_quantile, as.numeric(exceedence_thresh)),
                        slopes=c(0,0,0,0), cols=c("Block-Maxima","PoT","Empirical","B-M threshold"),
                        Metric=factor(c("Quantile","Quantile","Quantile","EP"),
                                      levels=c("Quantile","EP")))
  if(proba_ratio){
    ep <- ep/exceedence_thresh
  }
  df <- data.frame(id=int_index, Quantile=preds, EP=ep, Method="Prediction", Met_y="Observations")
  if(!is.null(y_obs)){
    df$y_ex <- y_obs
  }
  df <- df %>% tidyr::gather(key="Metric", value="Prediction", Quantile, EP, factor_key=TRUE)
  df[df$Metric=="EP", ]$y_ex <- NA
  levels(df$Metric) <- c("Quantile"=y_lab, "EP"=p_lab)
  levels(ablines$Metric) <- c("Quantile"=y_lab, "EP"=p_lab)
  
  if(!is.null(Date_index)){
    df$Date <- Date_index
    df$x_ind <- Date_index
    ep_plot <- df %>% ggplot( aes(x=x_ind, y=Prediction, group=Method, color=Method, linetype=Method)) + geom_line() +
      labs(y=NULL, x="Date", color=NULL, linetype=NULL) +
      scale_x_date(date_breaks=date_breaks, date_minor_breaks = date_minor_breaks, date_labels = "%Y %b", expand=c(0.01,0))
  }else{
    df$x_ind <- int_index
    ep_plot <- df %>% ggplot( aes(x=x_ind, y=Prediction, group=Method, color=Method, linetype=Method)) + geom_line() +
      labs(y=NULL, x="id", color=NULL, linetype=NULL)
  }
  ep_plot <- ep_plot +
    scale_color_manual(values = c("Prediction"=color_met,"Block-Maxima"=my_palette$red,"B-M threshold"=my_palette$red,
                                  "PoT"=my_palette_methods[["Uncond"]],"Empirical"=my_palette$light_blue)) + 
    scale_linetype_manual(values=c("Prediction"="solid","Block-Maxima"="twodash","B-M threshold"="blank","PoT"="dotted","Empirical"="dotted")) +
    geom_abline(data=ablines, aes(intercept=intercepts, slope=slopes, col=cols, linetype=cols), size=1, show.legend=FALSE) +
    geom_point(data=filter(df, x_ind %in% event_dates), shape=17, size=2) +
    theme(panel.grid.major.x=element_line(size=0.5), axis.ticks.x=element_line(size=0.5), legend.position=legend.position) +
    facet_grid(Metric ~ ., scales="free_y", switch="y", labeller=label_parsed)
  ep_plot <- ep_plot +
    if(log_scaling){stop("'log_scaling' not yet possible in 'plot_pred_quant_risk_ts'.")}else{scale_y_continuous(expand=c(0.02,0))}
  if(!is.null(event_dates)){ep_plot <- ep_plot + geom_vline(xintercept=event_dates, linetype="dashed")}
  if(!is.null(y_obs)){ep_plot <- ep_plot + geom_point(aes(y=y_ex), show.legend=FALSE, color="darkslategrey", alpha=0.8)}
  return(ep_plot)
}

plot_features_for_pred_ts <- function(Y, X, Z, predictions, return_level, event_ind, Date_index=NULL, int_index=seq_along(Y), seq_len=10,
                                      var_names=c("Y","X","Z"), date_breaks="1 month", date_minor_breaks="1 day", legend.position="bottom"){
  break_f = function(x) c(seq.Date(from = as.Date(paste0(year(min(x)),"-01-01")), to = max(x), by = date_breaks), event_ind)
  mbreak_f = function(x) seq.Date(from = as.Date(paste0(year(min(x)),"-01-01")), to = max(x), by = date_minor_breaks)
  ablines <- data.frame(intercepts=c(return_level), slopes=c(0), cols=c("Block-Maxima"),
                        Station=factor(c("Y"), levels=c("Y","X","Z")))
  levels(ablines$Station) <- c(Y = var_names[1], X = var_names[2], Z = var_names[3])
  if(!is.null(Date_index)){
    pred <- predictions[which(Date_index==event_ind)]
    df <- data.frame(id=c(int_index,NA), Date=c(Date_index, event_ind),
                     Y=c(Y,pred), X=c(X,NA), Z=c(Z,NA), Features=c(rep(NA,length(Y)), "Pred. Quantile"))
    df[df$Date<event_ind & df$Date>=event_ind-seq_len,]$Features <- "Covariate"
    df <- df %>% tidyr::gather(key="Station", value="Val", Y, X, Z, factor_key=TRUE)
    levels(df$Station) <- c(Y = var_names[1], X = var_names[2], Z = var_names[3])
    
    preds_df <- data.frame(id=int_index, Date=Date_index, Prediction=predictions, Station=factor(c("Y"), levels=c("Y","X","Z")))
    levels(preds_df$Station) <- c(Y = var_names[1], X = var_names[2], Z = var_names[3])
    
    sc_plot <- df %>% ggplot(aes(x=Date, y=Val, color=Features)) + 
      geom_abline(data=ablines, aes(intercept=intercepts, slope=slopes), color=my_palette$red, linetype="twodash", size=1, show.legend=FALSE) +
      geom_line(data=preds_df, aes(x=Date, y=Prediction), color=alpha(my_palette_methods[["EQRN"]], 0.4), size=1, inherit.aes=FALSE) +
      geom_point(aes(shape=Features, size=Features)) + labs(y=NULL, x="Date", color=NULL, linetype=NULL, shape=NULL) +
      scale_x_date(breaks=break_f, minor_breaks=mbreak_f, date_labels="%Y %b %d",
                   limits=c(min(Date_index), max(Date_index)), expand=c(0.02,0)) +
      scale_color_manual(values=c("Covariate"=my_palette$green, "Pred. Quantile"=my_palette$blue), na.value="black", guide="legend") +
      scale_shape_manual(values=c("Covariate"=18, "Pred. Quantile"=17), na.value=19, guide="legend") +
      scale_size_manual(values=c("Covariate"=2.5, "Pred. Quantile"=2), na.value=1.5, guide="legend") +
      scale_y_continuous(expand=c(0.05,0)) + facet_grid(Station ~ ., scales="free_y", switch="y", labeller=label_parsed) +
      theme(panel.grid.major.x=element_line(size=0.5), axis.ticks.x=element_line(size=0.5), legend.position=legend.position)
  }else{
    pred <- predictions[which(int_index==event_ind)]
    df <- data.frame(id=c(int_index, event_ind),
                     Y=c(Y,pred), X=c(X,NA), Z=c(Z,NA), Features=c(rep(NA,length(Y)), "Pred. Quantile"))
    df[df$id<event_ind & df$id>=event_ind-seq_len,]$Features <- "Covariate"
    df <- df %>% tidyr::gather(key="Station", value="Val", Y, X, Z, factor_key=TRUE)
    levels(df$Station) <- c(Y = var_names[1], X = var_names[2], Z = var_names[3])
    
    preds_df <- data.frame(id=int_index, Date=Date_index, Prediction=predictions, Station=factor(c("Y"), levels=c("Y","X","Z")))
    levels(preds_df$Station) <- c(Y = var_names[1], X = var_names[2], Z = var_names[3])
    
    sc_plot <- df %>% ggplot( aes(x=id, y=Val, color=Features)) + 
      geom_abline(data=ablines, aes(intercept=intercepts, slope=slopes), color=my_palette$red, linetype="twodash", size=1, show.legend=FALSE) +
      geom_line(data=preds_df, aes(x=id, y=Prediction), color=alpha(my_palette_methods[["EQRN"]], 0.4), size=1, inherit.aes=FALSE) +
      geom_point(aes(shape=Features, size=Features)) + labs(y=NULL, x="Date", color=NULL, linetype=NULL, shape=NULL) +
      scale_x_continuous(expand=c(0.02,0)) +
      scale_color_manual(values=c("Covariate"=my_palette$green, "Pred. Quantile"=my_palette$blue), na.value="black", guide="legend") +
      scale_shape_manual(values=c("Covariate"=18, "Pred. Quantile"=17), na.value=19, guide="legend") +
      scale_size_manual(values=c("Covariate"=2.5, "Pred. Quantile"=2), na.value=1.5, guide="legend") +
      scale_y_continuous(expand=c(0.05,0)) + facet_grid(Station ~ ., scales="free_y", switch="y", labeller=label_parsed) +
      theme(panel.grid.major.x=element_line(size=0.5), axis.ticks.x=element_line(size=0.5), legend.position=legend.position)
  }
  return(sc_plot)
}


plot_series_comparison2 <- function(Y, X, Date_index=NULL, int_index=seq_along(Y),
                                    var_names=c("Y","X"), date_breaks="3 month", date_minor_breaks="1 day",
                                    event_dates=NULL, show_top=0){
  
  break_f = function(x) seq.Date(from = as.Date(paste0(year(min(x)),"-01-01")), to = max(x), by = date_breaks)
  mbreak_f = function(x) seq.Date(from = as.Date(paste0(year(min(x)),"-01-01")), to = max(x), by = date_minor_breaks)
  if(!is.null(Date_index)){
    df <- data.frame(id=int_index, Date=Date_index, Y=Y, X=X) %>% tidyr::gather(key="Station", value="Val", Y, X, factor_key=TRUE)
    levels(df$Station) <- c(Y = var_names[1], X = var_names[2])
    sc_plot <- df %>% ggplot( aes(x=Date, y=Val)) + geom_line() + labs(y=NULL, x="Date") +
      scale_x_date(breaks=break_f, minor_breaks=mbreak_f, date_labels="%Y %b",
                   limits=c(min(Date_index), max(Date_index)), expand=c(0.01,0)) +
      scale_y_continuous(expand=c(0.02,0)) + facet_grid(Station ~ ., scales="free_y", switch="y", labeller=label_parsed) +
      theme(panel.grid.major.x=element_line(size=0.5), axis.ticks.x=element_line(size=0.5))
  }else{
    df <- data.frame(id=int_index, Y=Y, X=X) %>% tidyr::gather(key="Variable", value="Val", Y, X, factor_key=TRUE)
    levels(df$Variable) <- c(Y = var_names[1], X = var_names[2])
    sc_plot <- df %>% ggplot( aes(x=id, y=Val)) + geom_line() + labs(y=NULL, x=NULL) +
      scale_x_continuous(expand=c(0.01,0)) +
      scale_y_continuous(expand=c(0.02,0)) + facet_wrap(~Variable, scales="free_y", ncol=1, strip.position="left", labeller=label_parsed)
  }
  if(!is.null(event_dates)){
    sc_plot <- sc_plot + geom_vline(xintercept=event_dates, linetype="dashed")
  }
  if(show_top>0){
    sc_plot <- sc_plot + geom_vline(xintercept=as.numeric((df %>% arrange(desc(Val)) %>% head(show_top))[["Date"]]), linetype=4)
  }
  return(sc_plot)
}

plot_series_comparison3 <- function(Y, X, Z, Date_index=NULL, int_index=seq_along(Y),
                                    var_names=c("'Y'","'X'","'Z'"), date_breaks="3 month", date_minor_breaks="1 day",
                                    event_dates=NULL, show_top=0){
  break_f = function(x) seq.Date(from = as.Date(paste0(year(min(x)),"-01-01")), to = max(x), by = date_breaks)
  mbreak_f = function(x) seq.Date(from = as.Date(paste0(year(min(x)),"-01-01")), to = max(x), by = date_minor_breaks)
  if(!is.null(Date_index)){
    df <- data.frame(id=int_index, Date=Date_index, Y=Y, X=X, Z=Z) %>% tidyr::gather(key="Station", value="Val", Y, X, Z, factor_key=TRUE)
    levels(df$Station) <- c(Y = var_names[1], X = var_names[2], Z = var_names[3])
    sc_plot <- df %>% ggplot( aes(x=Date, y=Val)) + geom_line() + labs(y=NULL, x="Date") +
      scale_x_date(breaks=break_f, minor_breaks=mbreak_f, date_labels="%Y %b",
                   limits=c(min(Date_index), max(Date_index)), expand=c(0.01,0)) +
      scale_y_continuous(expand=c(0.02,0)) + facet_grid(Station ~ ., scales="free_y", switch="y", labeller=label_parsed) +
      theme(panel.grid.major.x=element_line(size=0.5), axis.ticks.x=element_line(size=0.5))
  }else{
    df <- data.frame(id=int_index, Y=Y, X=X, Z=Z) %>% tidyr::gather(key="Station", value="Val", X, Y, Z, factor_key=TRUE)
    levels(df$Station) <- c(Y = var_names[1], X = var_names[2], Z = var_names[3])
    sc_plot <- df %>% ggplot( aes(x=id, y=Val)) + geom_line() + labs(y=NULL, x=NULL) +
      scale_x_continuous(expand=c(0.01,0)) +
      scale_y_continuous(expand=c(0.02,0)) + facet_wrap(~Station, scales="free_y", ncol=1, strip.position="left", labeller=label_parsed)
  }
  if(!is.null(event_dates)){
    sc_plot <- sc_plot + geom_vline(xintercept=event_dates, linetype="dashed")
  }
  if(show_top>0){
    sc_plot <- sc_plot + geom_vline(xintercept=as.numeric((df %>% arrange(desc(Val)) %>% head(show_top))[["Date"]]), linetype=4)
  }
  return(sc_plot)
}

plot_series_comparison4 <- function(Y, X, Z, W, Date_index=NULL, int_index=seq_along(Y),
                                    var_names=c("'Y'","'X'","'Z'","'W'"), date_breaks="3 month", date_minor_breaks="1 day",
                                    event_dates=NULL, show_top=0){
  break_f = function(x) seq.Date(from = as.Date(paste0(year(min(x)),"-01-01")), to = max(x), by = date_breaks)
  mbreak_f = function(x) seq.Date(from = as.Date(paste0(year(min(x)),"-01-01")), to = max(x), by = date_minor_breaks)
  if(!is.null(Date_index)){
    df <- data.frame(id=int_index, Date=Date_index, Y=Y, X=X, Z=Z, W=W) %>%
      tidyr::gather(key="Station", value="Val", -id, -Date, factor_key=TRUE)
    levels(df$Station) <- c(Y = var_names[1], X = var_names[2], Z = var_names[3], W = var_names[4])
    sc_plot <- df %>% ggplot( aes(x=Date, y=Val)) + geom_line() + labs(y=NULL, x="Date") +
      scale_x_date(breaks=break_f, minor_breaks=mbreak_f, date_labels="%Y %b",
                   limits=c(min(Date_index), max(Date_index)), expand=c(0.01,0)) +
      scale_y_continuous(expand=c(0.02,0)) + facet_grid(Station ~ ., scales="free_y", switch="y", labeller=label_parsed) +
      theme(panel.grid.major.x=element_line(size=0.5), axis.ticks.x=element_line(size=0.5))
  }else{
    df <- data.frame(id=int_index, Y=Y, X=X, Z=Z, W=W) %>% tidyr::gather(key="Station", value="Val", -id, factor_key=TRUE)
    levels(df$Station) <- c(Y = var_names[1], X = var_names[2], Z = var_names[3], W = var_names[4])
    sc_plot <- df %>% ggplot( aes(x=id, y=Val)) + geom_line() + labs(y=NULL, x=NULL) +
      scale_x_continuous(expand=c(0.01,0)) +
      scale_y_continuous(expand=c(0.02,0)) + facet_wrap(~Station, scales="free_y", ncol=1, strip.position="left", labeller=label_parsed)
  }
  if(!is.null(event_dates)){
    sc_plot <- sc_plot + geom_vline(xintercept=event_dates, linetype="dashed")
  }
  if(show_top>0){
    sc_plot <- sc_plot + geom_vline(xintercept=as.numeric((df %>% arrange(desc(Val)) %>% head(show_top))[["Date"]]), linetype=4)
  }
  return(sc_plot)
}

# ================================= DEPRECATED =================================



plot_series_comparison <- function(df, col1, col2, datecol="Date"){
  ps1 <- data.frame(df) %>% ggplot( aes_string(x=datecol, y=col1)) + geom_line() + labs(y=col1, x=NULL) +
    scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 day", date_labels = "%Y %b", expand=c(0.01,0)) +
    geom_vline(xintercept = as.numeric((df %>% arrange(desc(.data[[col1]])) %>% head(5))[["Date"]]), linetype=4) +
    scale_y_continuous(expand=c(0.02,0)) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  ps2 <- data.frame(df) %>% ggplot( aes_string(x=datecol, y=col2)) + geom_line() + labs(y=col2, x="Date") +
    scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 day", date_labels = "%Y %b", expand=c(0.01,0)) +
    scale_y_continuous(expand=c(0.02,0)) + 
    geom_vline(xintercept = as.numeric((df %>% arrange(desc(.data[[col2]])) %>% head(5))[["Date"]]), linetype=4)
  sc_plot <- ggpubr::ggarrange(ps1, ps2, labels=NULL, ncol=1, nrow=2, align="v")
  return(sc_plot)
}


plot_series_comparison_sim <- function(Y, X){
  ps1 <- data.frame(id=seq_along(Y), Y=Y, X=X) %>% ggplot( aes(x=id, y=Y)) + geom_line() + labs(y="Y", x=NULL) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  ps2 <- data.frame(id=seq_along(Y), Y=Y, X=X) %>% ggplot( aes(x=id, y=X)) + geom_line() + labs(y="X", x=NULL)
  sc_plot <- ggpubr::ggarrange(ps1, ps2, labels=NULL, ncol=1, nrow=2, align="v")
  return(sc_plot)
}



plot_error_quantile_ts_old <- function(fit_eqrn, fit_grf, fit_gbex, fit_egam, Y_train, X_test, y_test, prob_lvls_predict, pred_interm,
                                       interm_lvl, interm_quant_train, true_quantiles_test, factorize=TRUE, test_data="other", pred_interm_c=pred_interm){
  
  ntest <- nrow(X_test)
  nb_prob_lvls_predict <- length(prob_lvls_predict)
  #lagged datasets for competitor methods
  lagged_test <- lagged_features(X=cbind(y_test,X_test), max_lag=fit_eqrn$seq_len, drop_present=TRUE)
  laged_interm_q_test <- pred_interm_c[(fit_eqrn$seq_len+1):length(pred_interm_c), , drop=F]
  #High quantile prediction with GRF
  pred_grf_test <- predict(fit_grf, newdata=lagged_test, quantiles = prob_lvls_predict)$predictions
  # UNCONDITIONAL predicted quantile(s) (Y quantile on X_train)
  pred_unc <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = prob_lvls_predict, Y = Y_train, ntest = nrow(true_quantiles_test))
  # SEMI-CONDITIONAL predicted quantiles
  pred_semicond <- predict_GPD_semiconditional(Y=Y_train[(seq_len+1):length(Y_train)], interm_lvl=interm_lvl,
                                               thresh_quantiles=interm_quant_train[(seq_len+1):length(interm_quant_train)],
                                               interm_quantiles_test=pred_interm[(seq_len+1):length(pred_interm)], prob_lvls_predict=prob_lvls_predict)
  # GROUND-TRUTH (y_test)
  pred_true <- true_quantiles_test
  # Prediction GBEX
  pred_gbex <- gbex_predict(fit_gbex, lagged_test, to_predict=prob_lvls_predict, intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
  pred_exgam <- predict_gpd_gam(fit_egam, lagged_test, to_predict=prob_lvls_predict,
                                intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
  
  #Final EQRN predictions on X_test
  pred_eqrnn <- EQRN_predict_seq(fit_eqrn, X_test, y_test, prob_lvls_predict, pred_interm, interm_lvl, crop_predictions=TRUE)
  
  # Compute losses for desired predicted quantiles
  MAEs_eqrnn <- multilevel_MAE(pred_true,pred_eqrnn,prob_lvls_predict,give_names=FALSE)
  RMSEs_eqrnn <- sqrt(multilevel_MSE(pred_true,pred_eqrnn,prob_lvls_predict,give_names=FALSE))
  MAEs_grf <- multilevel_MAE(pred_true,pred_grf_test,prob_lvls_predict,give_names=FALSE)
  RMSEs_grf <- sqrt(multilevel_MSE(pred_true,pred_grf_test,prob_lvls_predict,give_names=FALSE))
  MAEs_unc <- multilevel_MAE(pred_true,pred_unc$predictions,prob_lvls_predict,give_names=FALSE)
  RMSEs_unc <- sqrt(multilevel_MSE(pred_true,pred_unc$predictions,prob_lvls_predict,give_names=FALSE))
  MAEs_semicond <- multilevel_MAE(pred_true,pred_semicond$predictions,prob_lvls_predict,give_names=FALSE)
  RMSEs_semicond <- sqrt(multilevel_MSE(pred_true,pred_semicond$predictions,prob_lvls_predict,give_names=FALSE))
  MAEs_gbex <- multilevel_MAE(pred_true,pred_gbex,prob_lvls_predict,give_names=FALSE)
  RMSEs_gbex <- sqrt(multilevel_MSE(pred_true,pred_gbex,prob_lvls_predict,give_names=FALSE))
  MAEs_exgam <- multilevel_MAE(pred_true,pred_exgam,prob_lvls_predict,give_names=FALSE)
  RMSEs_exgam <- sqrt(multilevel_MSE(pred_true,pred_exgam,prob_lvls_predict,give_names=FALSE))
  
  if(factorize){
    df_rmse <- data.frame(Quantile=as_factor(roundm(prob_lvls_predict,4)), Uncond=RMSEs_unc, Semi_cond=RMSEs_semicond,
                          GRF=RMSEs_grf, gpd_GAM=RMSEs_exgam, gbex=RMSEs_gbex, EQRN=RMSEs_eqrnn)
    df_mae <- data.frame(Quantile=as_factor(roundm(prob_lvls_predict,4)), Uncond=MAEs_unc, Semi_cond=MAEs_semicond,
                         GRF=MAEs_grf, gpd_GAM=MAEs_exgam, gbex=MAEs_gbex, EQRN=MAEs_eqrnn)
  } else {
    df_rmse <- data.frame(Quantile=prob_lvls_predict, Uncond=RMSEs_unc, Semi_cond=RMSEs_semicond,
                          GRF=RMSEs_grf, gpd_GAM=RMSEs_exgam, gbex=RMSEs_gbex, EQRN=RMSEs_eqrnn)
    df_mae <- data.frame(Quantile=prob_lvls_predict, Uncond=MAEs_unc, Semi_cond=MAEs_semicond,
                         GRF=MAEs_grf, gpd_GAM=MAEs_exgam, gbex=MAEs_gbex, EQRN=MAEs_eqrnn)
  }
  
  y_lab_rmse <- if(test_data=="halton"){"RISE"}else{"RMSE"}
  
  linetypes <- c("dotted","dotted","dashed","solid","solid","solid")
  color_scale <- c("#ff7f00","#984ea3","#4daf4a","#69b3a2","#e41a1c","#377eb8")
  
  mse_plot <- df_rmse %>% tidyr::gather(key="Model", value="Error", Uncond, Semi_cond, GRF, gpd_GAM, gbex, EQRN, factor_key=TRUE) %>%
    ggplot( aes(x=Quantile, y=Error, group=Model, color=Model, linetype=Model)) +
    geom_line(size=1) + scale_linetype_manual(values=linetypes) + scale_color_manual(values=color_scale) + 
    scale_y_continuous(expand=c(0.02,0)) + labs(title=NULL, x="Probability level", y=y_lab_rmse, color=NULL, linetype=NULL)
  mae_plot <- df_mae %>% tidyr::gather(key="Model", value="Error", Uncond, Semi_cond, GRF, gpd_GAM, gbex, EQRN, factor_key=TRUE) %>%
    ggplot( aes(x=Quantile, y=Error, group=Model, color=Model, linetype=Model)) +
    geom_line(size=1) + scale_linetype_manual(values=linetypes) + scale_color_manual(values=color_scale) +
    scale_y_continuous(expand=c(0.02,0)) + labs(title=NULL, x="Probability level", y="MAE", color=NULL, linetype=NULL)
  if(factorize){
    mse_plot <- mse_plot + scale_x_discrete(expand=c(0.02,0))
    mae_plot <- mae_plot + scale_x_discrete(expand=c(0.02,0))
  }else{
    mse_plot <- mse_plot + scale_x_continuous(expand=c(0.02,0))
    mae_plot <- mae_plot + scale_x_continuous(expand=c(0.02,0))
  }
  err_plot <- ggpubr::ggarrange(mse_plot, mae_plot, labels=NULL, ncol=2, nrow=1,
                                common.legend=TRUE,legend="bottom", align="hv")
  return(err_plot)
}



plot_rmse_quantile_comp_ts_old <- function(fit_eqrn1, fit_eqrn2, fit_grf, fit_gbex, fit_egam, Y_train, X_test, y_test, prob_lvls_predict, pred_interm,
                                           interm_lvl, interm_quant_train, true_quantiles_test, factorize=TRUE, test_data="other", pred_interm_c=pred_interm,
                                           names_EQRN=c("EQRN_1","EQRN_2"), legend.position="bottom", crop_obs=0, bias_variance=FALSE, crop_Rbv=c(0,0,0)){
  
  ntest <- nrow(X_test)
  nb_prob_lvls_predict <- length(prob_lvls_predict)
  if(fit_eqrn1$seq_len!=fit_eqrn2$seq_len){
    stop("Unfair comparison due to different eqrnn 'seq_len' in 'plot_error_quantile_comp_ts'.")
  }
  #lagged datasets for competitor methods
  lagged_test <- lagged_features(X=cbind(y_test,X_test), max_lag=fit_eqrn1$seq_len, drop_present=TRUE)
  laged_interm_q_test <- pred_interm_c[(fit_eqrn1$seq_len+1):length(pred_interm_c), , drop=F]
  #High quantile prediction with GRF
  pred_grf_test <- predict(fit_grf, newdata=lagged_test, quantiles = prob_lvls_predict)$predictions
  # UNCONDITIONAL predicted quantile(s) (Y quantile on X_train)
  pred_unc <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = prob_lvls_predict, Y = Y_train, ntest = nrow(true_quantiles_test))
  # SEMI-CONDITIONAL predicted quantiles
  pred_semicond <- predict_GPD_semiconditional(Y=Y_train[(seq_len+1):length(Y_train)], interm_lvl=interm_lvl,
                                               thresh_quantiles=interm_quant_train[(seq_len+1):length(interm_quant_train)],
                                               interm_quantiles_test=pred_interm[(seq_len+1):length(pred_interm)], prob_lvls_predict=prob_lvls_predict)
  # GROUND-TRUTH (y_test)
  pred_true <- true_quantiles_test
  # Prediction GBEX
  pred_gbex <- gbex_predict(fit_gbex, lagged_test, to_predict=prob_lvls_predict, intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
  pred_exgam <- predict_gpd_gam(fit_egam, lagged_test, to_predict=prob_lvls_predict,
                                intermediate_quantiles=laged_interm_q_test, interm_lvl=interm_lvl)
  
  #Final EQRN predictions on X_test
  pred_eqrnn1 <- EQRN_predict_seq(fit_eqrn1, X_test, y_test, prob_lvls_predict, pred_interm, interm_lvl, crop_predictions=TRUE)
  pred_eqrnn2 <- EQRN_predict_seq(fit_eqrn2, X_test, y_test, prob_lvls_predict, pred_interm, interm_lvl, crop_predictions=TRUE)
  
  # Compute losses for desired predicted quantiles
  RMSEs_eqrnn1 <- sqrt(multilevel_MSE(pred_true,pred_eqrnn1,prob_lvls_predict,give_names=FALSE))
  RMSEs_eqrnn2 <- sqrt(multilevel_MSE(pred_true,pred_eqrnn2,prob_lvls_predict,give_names=FALSE))
  RMSEs_grf <- sqrt(multilevel_MSE(pred_true,pred_grf_test,prob_lvls_predict,give_names=FALSE))
  RMSEs_unc <- sqrt(multilevel_MSE(pred_true,pred_unc$predictions,prob_lvls_predict,give_names=FALSE))
  RMSEs_semicond <- sqrt(multilevel_MSE(pred_true,pred_semicond$predictions,prob_lvls_predict,give_names=FALSE))
  RMSEs_gbex <- sqrt(multilevel_MSE(pred_true,pred_gbex,prob_lvls_predict,give_names=FALSE))
  RMSEs_exgam <- sqrt(multilevel_MSE(pred_true,pred_exgam,prob_lvls_predict,give_names=FALSE))
  
  Q_lvls <- if(factorize){as_factor(roundm(prob_lvls_predict,4))}else{prob_lvls_predict}
  
  met_names <- c("Uncond", "Semi-cond", "GRF", "EGAM", "GBEX", names_EQRN[2], names_EQRN[1])
  df_rmse <- data.frame(Quantile=Q_lvls, Uncond=RMSEs_unc, Semi_cond=RMSEs_semicond,
                        GRF=RMSEs_grf, EGAM=RMSEs_exgam, GBEX=RMSEs_gbex,
                        EQRN2=RMSEs_eqrnn2, EQRN1=RMSEs_eqrnn1) %>% 
    rename_with(~str_replace(., "Semi_cond", "Semi-cond")) %>%
    rename_with(~str_replace(., "EQRN2", names_EQRN[2])) %>%
    rename_with(~str_replace(., "EQRN1", names_EQRN[1])) %>% 
    tidyr::gather(key="Model", value="Error", all_of(met_names), factor_key=TRUE)
  
  y_lab_rmse <- if(test_data=="halton"){"RISE"}else{"RMSE"}
  
  linetypes <- c("dotted","dashed","dotdash","longdash","twodash","dashed","solid")
  color_scale <- my_palette_methods[c("Uncond", "Semi-cond", "GRF", "EGAM", "GBEX", "EQRN2", "EQRN")]
  names(linetypes) <- met_names
  names(color_scale) <- met_names
  
  mse_plot <- df_rmse %>% ggplot( aes(x=Quantile, y=Error, group=Model, color=Model, linetype=Model)) +
    geom_line(size=1) + scale_linetype_manual(values=linetypes) + scale_color_manual(values=color_scale) + 
    scale_y_continuous(expand=c(0.02,0)) + labs(title=NULL, x="Probability level", y=y_lab_rmse, color=NULL, linetype=NULL) +
    theme(legend.position=legend.position)
  
  mse_plot <- mse_plot + if(factorize){scale_x_discrete(expand=c(0.02,0))}else{scale_x_continuous(expand=c(0.02,0))}
  if(crop_obs>0){mse_plot <- mse_plot + coord_cartesian(ylim=c(min(df_rmse$Error), sort(df_rmse$Error, decreasing=TRUE)[[crop_obs+1]]))}
  
  if(!bias_variance){
    return(mse_plot)
    
  }else{
    # Compute metrics for desired predicted quantiles
    bias_eqrnn1 <- multilevel_pred_bias(pred_true, pred_eqrnn1, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_eqrnn2 <- multilevel_pred_bias(pred_true, pred_eqrnn2, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_grf <- multilevel_pred_bias(pred_true, pred_grf_test, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_unc <- multilevel_pred_bias(pred_true, pred_unc$predictions, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_semicond <- multilevel_pred_bias(pred_true, pred_semicond$predictions, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_gbex <- multilevel_pred_bias(pred_true, pred_gbex, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    bias_exgam <- multilevel_pred_bias(pred_true, pred_exgam, prob_lvls_predict, square_bias=FALSE, give_names=FALSE)
    
    rstd_eqrnn1 <- sqrt(multilevel_resid_var(pred_true, pred_eqrnn1, prob_lvls_predict, give_names=FALSE))
    rstd_eqrnn2 <- sqrt(multilevel_resid_var(pred_true, pred_eqrnn2, prob_lvls_predict, give_names=FALSE))
    rstd_grf <- sqrt(multilevel_resid_var(pred_true, pred_grf_test, prob_lvls_predict, give_names=FALSE))
    rstd_unc <- sqrt(multilevel_resid_var(pred_true, pred_unc$predictions, prob_lvls_predict, give_names=FALSE))
    rstd_semicond <- sqrt(multilevel_resid_var(pred_true, pred_semicond$predictions, prob_lvls_predict, give_names=FALSE))
    rstd_gbex <- sqrt(multilevel_resid_var(pred_true, pred_gbex, prob_lvls_predict, give_names=FALSE))
    rstd_exgam <- sqrt(multilevel_resid_var(pred_true, pred_exgam, prob_lvls_predict, give_names=FALSE))
    
    Rsqr_eqrnn1 <- multilevel_R_squared(pred_true, pred_eqrnn1, prob_lvls_predict, give_names=FALSE)
    Rsqr_eqrnn2 <- multilevel_R_squared(pred_true, pred_eqrnn2, prob_lvls_predict, give_names=FALSE)
    Rsqr_grf <- multilevel_R_squared(pred_true, pred_grf_test, prob_lvls_predict, give_names=FALSE)
    Rsqr_unc <- multilevel_R_squared(pred_true, pred_unc$predictions, prob_lvls_predict, give_names=FALSE)
    Rsqr_semicond <- multilevel_R_squared(pred_true, pred_semicond$predictions, prob_lvls_predict, give_names=FALSE)
    Rsqr_gbex <- multilevel_R_squared(pred_true, pred_gbex, prob_lvls_predict, give_names=FALSE)
    Rsqr_exgam <- multilevel_R_squared(pred_true, pred_exgam, prob_lvls_predict, give_names=FALSE)
    
    Rbv_names <- c("Quantile R squared", "Bias", "Residual standard deviation")
    
    df_Rbv <- data.frame(Quantile=rep(Q_lvls, 3),
                         Uncond=c(Rsqr_unc,bias_unc,rstd_unc),
                         Semi_cond=c(Rsqr_semicond,bias_semicond,rstd_semicond),
                         GRF=c(Rsqr_grf,bias_grf,rstd_grf),
                         EGAM=c(Rsqr_exgam,bias_exgam,rstd_exgam),
                         GBEX=c(Rsqr_gbex,bias_gbex,rstd_gbex),
                         EQRN2=c(Rsqr_eqrnn2,bias_eqrnn2,rstd_eqrnn2),
                         EQRN1=c(Rsqr_eqrnn1,bias_eqrnn1,rstd_eqrnn1),
                         metric=factor(rep(Rbv_names, each=nb_prob_lvls_predict), levels=Rbv_names)) %>%
      rename_with(~str_replace(., "Semi_cond", "Semi-cond")) %>%
      rename_with(~str_replace(., "EQRN2", names_EQRN[2])) %>%
      rename_with(~str_replace(., "EQRN1", names_EQRN[1])) %>%
      tidyr::gather(key="Model", value="Error", all_of(met_names), factor_key=TRUE)
    
    Rbv_lims <- list(c(max(min(sort(df_Rbv$Error[df_Rbv$metric==Rbv_names[1]], decreasing=FALSE)[[crop_Rbv[1]+1]],0),-5.), 1.),
                     (c(-1,1) * sort(abs(df_Rbv$Error[df_Rbv$metric==Rbv_names[2]]), decreasing=TRUE)[[crop_Rbv[2]+1]]),
                     c(min(df_Rbv$Error[df_Rbv$metric==Rbv_names[3]]),
                       sort(df_Rbv$Error[df_Rbv$metric==Rbv_names[3]], decreasing=TRUE)[[crop_Rbv[3]+1]]))
    Rbv_plots <- list()
    for(i in seq_along(Rbv_names)){
      Rbv_plots[[i]] <- df_Rbv %>% filter(metric==Rbv_names[i]) %>%
        ggplot( aes(x=Quantile, y=Error, group=Model, color=Model, linetype=Model)) +
        geom_line(size=1) + scale_linetype_manual(values=linetypes) + scale_color_manual(values=color_scale) +
        scale_y_continuous(expand=c(0.02,0)) +
        labs(title=NULL, x=if(i==2){"Probability level"}else{NULL}, y=Rbv_names[i], color=NULL, linetype=NULL) +
        theme(legend.position=legend.position) + coord_cartesian(ylim=Rbv_lims[[i]]) +
        if(factorize){scale_x_discrete(expand=c(0.02,0))}else{scale_x_continuous(expand=c(0.02,0))}
    }
    #facet_wrap(~metric, scales="free_y", ncol=3, strip.position="left") + 
    Rbv_plot <- ggpubr::ggarrange(plotlist=Rbv_plots, nrow=1, ncol=3, labels=NULL,
                                  common.legend=TRUE, legend="bottom", align="hv")
    
    return(list(rmse_plot=mse_plot, Rbv_plot=Rbv_plot))
  }
}

