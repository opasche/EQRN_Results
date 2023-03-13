
training_plot_eqrn <- function(fit_eqrn, uncond_losses_interm, title=NULL, y_lab="GPD likelihood loss", show_legend=TRUE){
  
  ablines <- data.frame(intercepts=c(uncond_losses_interm$train_loss, uncond_losses_interm$valid_loss), slopes=c(0,0),
                        cols=c("Semi-cond train","Semi-cond valid"))
  
  train_plot <- data.frame(Epoch=seq_along(fit_eqrn$train_loss), Train=fit_eqrn$train_loss, Validation=fit_eqrn$valid_loss) %>%
    tidyr::gather(key="Type", value="Value", Train, Validation, factor_key=TRUE) %>%
    ggplot( aes(x=Epoch, y=Value, group=Type, color=Type, linetype=Type)) +
    geom_line(size=1) + scale_color_manual(values = c(my_palette$red,my_palette$blue,my_palette$red,my_palette$blue), guide=guide_legend(reverse=TRUE)) +
    scale_linetype_manual(values=c("dashed","dashed","solid","solid"), guide=guide_legend(reverse=TRUE)) +
    geom_abline(data=ablines, aes(intercept=intercepts, slope=slopes, col=cols, linetype=cols), size=1, show.legend=FALSE) +
    scale_y_continuous(expand=c(0.02,0)) + scale_x_continuous(expand=c(0.02,0)) +
    labs(title=title, x="Epoch", y=y_lab, color=NULL, linetype=NULL)
  if(show_legend){
    train_plot <- train_plot + theme(
      legend.position = c(.99, .99),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(5, 5, 5, 5)
    )
  }else{
    train_plot <- train_plot + theme(legend.position="none")
  }
  return(train_plot)
}

validation_plot_eqrn <- function(fit_eqrn, uncond_losses_interm, title=NULL, y_lab="GPD likelihood loss", show_legend=FALSE){
  
  ablines <- data.frame(intercepts=c(uncond_losses_interm$valid_loss), slopes=c(0),
                        cols=c("Semi-cond valid"))
  lims <- c(min(fit_eqrn$valid_loss,uncond_losses_interm$valid_loss),max(fit_eqrn$valid_loss,uncond_losses_interm$valid_loss))
  train_plot <- data.frame(Epoch=seq_along(fit_eqrn$valid_loss), Validation=fit_eqrn$valid_loss) %>%
    tidyr::gather(key="Type", value="Value", Validation, factor_key=TRUE) %>%
    ggplot( aes(x=Epoch, y=Value, group=Type, color=Type, linetype=Type)) +
    geom_line(size=1) + scale_color_manual(values = c(my_palette_methods[["Semi-cond"]],my_palette_methods[["EQRN"]]), guide=guide_legend(reverse=TRUE)) +
    scale_linetype_manual(values=c("dashed","solid"), guide=guide_legend(reverse=TRUE)) +
    geom_abline(data=ablines, aes(intercept=intercepts, slope=slopes, col=cols, linetype=cols), size=1, show.legend=FALSE) +
    scale_y_continuous(expand=c(0.02,0)) + scale_x_continuous(expand=c(0.02,0)) + expand_limits(y=lims)
    labs(title=title, x="Epoch", y=y_lab, color=NULL, linetype=NULL)
  if(show_legend){
    train_plot <- train_plot + theme(
      legend.position = c(.99, .99),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(5, 5, 5, 5)
    )
  }else{
    train_plot <- train_plot + theme(legend.position="none")
  }
  return(train_plot)
}

plot_predictions_1D <- function(pred_eqrnn_test, pred_grf_test, pred_true, X_test, quantile_predicted){
  
  pred_plots <- list()
  for (j in 1:4) {
    pred_plots[[j]] <- data.frame(X=X_test[,j], GRF=pred_grf_test, EQRNN=pred_eqrnn_test, Truth=pred_true) %>%
      tidyr::gather(key="Method", value="Prediction", GRF, EQRNN, Truth, factor_key=TRUE) %>%
      ggplot(aes(x=X, y=Prediction, group=Method, color=Method)) +
      geom_point() + 
      scale_color_manual(values=alpha(c("#69b3a2",rgb(57/255,106/255,177/255), "black"),c(0.2,0.2,1))) +
      labs(x=paste0("X",j), y=NULL, color=NULL)
  }
  pred_plot <- ggpubr::ggarrange(plotlist=pred_plots, labels="AUTO", ncol=min(2,length(pred_plots)), nrow=2,
                                 common.legend=TRUE,legend="bottom", vjust=1.5, align="hv") %>%
    ggpubr::annotate_figure(top=text_grob(paste0("Predictions (q=",quantile_predicted,")"), face="bold", size=14),
                            left=text_grob("Predictions", rot = 90))
  return(pred_plot)
}

residuals_boxplot_eqrn <- function(pred_eqrnn_test, pred_grf_test, pred_true, quantile_predicted){
  
  res_box <- data.frame(GRF=c(pred_grf_test-pred_true), EQRNN=c(pred_eqrnn_test-pred_true)) %>%
    tidyr::gather(key="Method", value="Value", factor_key=TRUE) %>% ggplot( aes(x=Method, y=Value, fill=Method)) +
    geom_abline(intercept=0, slope=0, col="white", size=1) + geom_violin(width=1, color="grey", alpha=0.5) +
    stat_boxplot(geom = "errorbar", width=0.5) + geom_boxplot(width=0.75, color="black", alpha=0.8)+ coord_flip() +
    stat_summary(fun=mean, geom="point", shape=4, size=3, color="black", fill="black") +
    labs(x=NULL, y="Residuals") +
    theme(legend.position="none", panel.grid.major.y = element_blank())
  abs_box <- data.frame(GRF=abs(c(pred_grf_test-pred_true)), EQRNN=abs(c(pred_eqrnn_test-pred_true))) %>%
    tidyr::gather(key="Method", value="Value", factor_key=TRUE) %>% ggplot( aes(x=Method, y=Value, fill=Method)) +
    geom_abline(intercept=0, slope=0, col="white", size=1) + geom_violin(width=1, color="grey", alpha=0.5) +
    stat_boxplot(geom = "errorbar", width=0.5) + geom_boxplot(width=0.75, color="black", alpha=0.8)+ coord_flip() +
    stat_summary(fun=mean, geom="point", shape=4, size=3, color="black", fill="black") +
    labs(x=NULL, y="Absolute Residuals") +
    theme(legend.position="none", panel.grid.major.y = element_blank())
  sq_box <- data.frame(GRF=(c(pred_grf_test-pred_true))^2, EQRNN=(c(pred_eqrnn_test-pred_true))^2) %>%
    tidyr::gather(key="Method", value="Value", factor_key=TRUE) %>% ggplot( aes(x=Method, y=Value, fill=Method)) +
    geom_abline(intercept=0, slope=0, col="white", size=1) + geom_violin(width=1, color="grey", alpha=0.5) +
    stat_boxplot(geom = "errorbar", width=0.5) + geom_boxplot(width=0.75, color="black", alpha=0.8)+ coord_flip() +
    stat_summary(fun=mean, geom="point", shape=4, size=3, color="black", fill="black") +
    labs(x=NULL, y="Squared Residuals") +
    theme(legend.position="none", panel.grid.major.y = element_blank())
  resid_boxes <- ggpubr::ggarrange(res_box, abs_box, sq_box, labels=NULL, ncol=1, nrow=3) %>%
    ggpubr::annotate_figure(top=text_grob(paste0("Quantile Residuals (q=",quantile_predicted,")"), face="bold", size=14))
  return(resid_boxes)
}

residuals_boxplot_eqrn2 <- function(pred_eqrnn_test, pred_grf_test, pred_gbex_test, pred_true, quantile_predicted){
  
  res_box <- data.frame(GRF=c(pred_grf_test-pred_true), gbex=c(pred_gbex_test-pred_true),
                        EQRNN=c(pred_eqrnn_test-pred_true)) %>%
    tidyr::gather(key="Method", value="Value", factor_key=TRUE) %>% ggplot( aes(x=Method, y=Value, fill=Method)) +
    geom_abline(intercept=0, slope=0, col="white", size=1) + geom_violin(width=1, color="grey", alpha=0.5) +
    stat_boxplot(geom = "errorbar", width=0.5) + geom_boxplot(width=0.75, color="black", alpha=0.8)+ coord_flip() +
    stat_summary(fun=mean, geom="point", shape=4, size=3, color="black", fill="black") +
    labs(x=NULL, y="Residuals") +
    theme(legend.position="none", panel.grid.major.y = element_blank())
  abs_box <- data.frame(GRF=abs(c(pred_grf_test-pred_true)), gbex=abs(c(pred_gbex_test-pred_true)),
                        EQRNN=abs(c(pred_eqrnn_test-pred_true))) %>%
    tidyr::gather(key="Method", value="Value", factor_key=TRUE) %>% ggplot( aes(x=Method, y=Value, fill=Method)) +
    geom_abline(intercept=0, slope=0, col="white", size=1) + geom_violin(width=1, color="grey", alpha=0.5) +
    stat_boxplot(geom = "errorbar", width=0.5) + geom_boxplot(width=0.75, color="black", alpha=0.8)+ coord_flip() +
    stat_summary(fun=mean, geom="point", shape=4, size=3, color="black", fill="black") +
    labs(x=NULL, y="Absolute Residuals") +
    theme(legend.position="none", panel.grid.major.y = element_blank())
  sq_box <- data.frame(GRF=(c(pred_grf_test-pred_true))^2, gbex=(c(pred_gbex_test-pred_true))^2,
                       EQRNN=(c(pred_eqrnn_test-pred_true))^2) %>%
    tidyr::gather(key="Method", value="Value", factor_key=TRUE) %>% ggplot( aes(x=Method, y=Value, fill=Method)) +
    geom_abline(intercept=0, slope=0, col="white", size=1) + geom_violin(width=1, color="grey", alpha=0.5) +
    stat_boxplot(geom = "errorbar", width=0.5) + geom_boxplot(width=0.75, color="black", alpha=0.8)+ coord_flip() +
    stat_summary(fun=mean, geom="point", shape=4, size=3, color="black", fill="black") +
    labs(x=NULL, y="Squared Residuals") +
    theme(legend.position="none", panel.grid.major.y = element_blank())
  resid_boxes <- ggpubr::ggarrange(res_box, abs_box, sq_box, labels=NULL, ncol=1, nrow=3) %>%
    ggpubr::annotate_figure(top=text_grob(paste0("Quantile Residuals (q=",quantile_predicted,")"), face="bold", size=14))
  return(resid_boxes)
}


plot_predictions_vs_truth <- function(pred_eqrnn, pred_other, pred_semicond, pred_true, quantiles_predicted, name_other="GRF",
                                      legend.position="bottom"){
  
  nb_quantiles_predict <- length(quantiles_predicted)
  pred_vs_true_plots <- list()
  col_other <- if(is.na(my_palette_methods[name_other])){my_palette$red}else{my_palette_methods[[name_other]]}
  for(i in 1:nb_quantiles_predict){
    preds_df <- data.frame(Truth=pred_true[,i], Semi_cond=pred_semicond[,i], Other=pred_other[,i], EQRNN=pred_eqrnn[,i])
    preds_df <- preds_df %>% rename_with(~str_replace(., "Other", name_other))
    lims <- c(min(preds_df), max(preds_df))
    pred_vs_true_plots[[i]] <- preds_df %>%
      tidyr::gather(key="Method", value="Prediction", Semi_cond, .data[[name_other]], EQRNN, factor_key=TRUE) %>%
      ggplot(aes(x=Truth, y=Prediction, group=Method, color=Method, shape=Method)) +
      geom_abline(intercept=0, slope=1, col="black", size=1) + geom_point() + coord_fixed() +
      theme(aspect.ratio=1, legend.position=legend.position) + 
      scale_color_manual(values=alpha(c(my_palette_methods[["Semi_cond"]], col_other, my_palette_methods[["EQRNN"]]),c(0.4,0.4,0.4))) +
      scale_shape_manual(values=c(4, 18, 16)) +
      expand_limits(x=lims, y=lims) + scale_y_continuous(expand=c(0.02,0)) + scale_x_continuous(expand=c(0.02,0))
    if(legend.position=="none"){
      pred_vs_true_plots[[i]] <- pred_vs_true_plots[[i]] + labs(title=NULL, x=NULL, y=NULL, color=NULL)
    }else{
      pred_vs_true_plots[[i]] <- pred_vs_true_plots[[i]] +
        labs(title=paste0(quantiles_predicted[i]*100,"% quantile"), x=NULL, y=NULL, color=NULL)
    }
  }
  pred_vs_true_plot <- ggpubr::ggarrange(plotlist=pred_vs_true_plots, labels=NULL, ncol=min(3,length(pred_vs_true_plots)),
                                         common.legend=TRUE,legend=legend.position, align="hv") %>%
    ggpubr::annotate_figure(left=text_grob("Predicted quantiles", rot = 90), bottom=text_grob("True quantiles"))
  return(pred_vs_true_plot)
}


plot_predictions_vs_truth_solo <- function(pred_eqrnn, pred_semicond, pred_true, quantiles_predicted,
                                           legend.position="bottom"){
  
  nb_quantiles_predict <- length(quantiles_predicted)
  pred_vs_true_plots <- list()
  for(i in 1:nb_quantiles_predict){
    preds_df <- data.frame(Truth=pred_true[,i], Semi_cond=pred_semicond[,i], EQRNN=pred_eqrnn[,i])
    lims <- c(min(preds_df), max(preds_df))
    pred_vs_true_plots[[i]] <- preds_df %>%
      tidyr::gather(key="Method", value="Prediction", Semi_cond, EQRNN, factor_key=TRUE) %>%
      ggplot(aes(x=Truth, y=Prediction, group=Method, color=Method, shape=Method)) +
      geom_abline(intercept=0, slope=1, col="black", size=1) + geom_point() + coord_fixed() +
      theme(aspect.ratio=1, legend.position=legend.position) + 
      scale_color_manual(values=alpha(c(my_palette_methods[["Semi_cond"]], my_palette_methods[["EQRNN"]]),c(0.5,0.5))) +
      scale_shape_manual(values=c(4, 16)) +
      expand_limits(x=lims, y=lims) + scale_y_continuous(expand=c(0.02,0)) + scale_x_continuous(expand=c(0.02,0))
    if(legend.position=="none"){
      pred_vs_true_plots[[i]] <- pred_vs_true_plots[[i]] + labs(title=NULL, x=NULL, y=NULL, color=NULL)
    }else{
      pred_vs_true_plots[[i]] <- pred_vs_true_plots[[i]] +
        labs(title=paste0(quantiles_predicted[i]*100,"% quantile"), x=NULL, y=NULL, color=NULL)
    }
  }
  pred_vs_true_plot <- ggpubr::ggarrange(plotlist=pred_vs_true_plots, labels=NULL, ncol=min(3,length(pred_vs_true_plots)),
                                         common.legend=TRUE,legend=legend.position, align="hv") %>%
    ggpubr::annotate_figure(left=text_grob("Predicted quantiles", rot = 90), bottom=text_grob("True quantiles"))
  return(pred_vs_true_plot)
}

plot_predictions_vs_truth_comp <- function(pred_eqrnn1, pred_eqrnn2, pred_other, pred_true,
                                           names_EQRN=c("EQRNN_1","EQRNN_2"), name_other="Semi_cond", legend.position="bottom"){
  
  n <- length(pred_eqrnn1)
  preds_df <- data.frame(Truth=c(pred_true,pred_true), Other=c(pred_other,pred_other), EQRNN=c(pred_eqrnn1, pred_eqrnn2),
                         EQRNN_mod=factor(c(rep(names_EQRN[1], n),rep(names_EQRN[2], n)), levels=names_EQRN))
  lims <- c(min(pred_eqrnn1, pred_eqrnn2, pred_other, pred_true, na.rm=TRUE),
            max(pred_eqrnn1, pred_eqrnn2, pred_other, pred_true, na.rm=TRUE))
  preds_df <- preds_df %>% rename_with(~str_replace(., "Other", name_other))
  col_other <- if(is.na(my_palette_methods[name_other])){my_palette$red}else{my_palette_methods[[name_other]]}
  pred_vs_true_plot <- preds_df %>%
    tidyr::gather(key="Method", value="Prediction", .data[[name_other]], EQRNN, factor_key=TRUE) %>%
    ggplot(aes(x=Truth, y=Prediction, group=Method, color=Method, shape=Method)) +
    geom_abline(intercept=0, slope=1, col="black", size=1) + geom_point() + coord_fixed() +
    theme(aspect.ratio=1, legend.position=legend.position) + 
    scale_color_manual(values=alpha(c(col_other, my_palette_methods[["EQRNN"]]),c(0.5,0.5))) +
    scale_shape_manual(values=c(4, 16)) +
    expand_limits(x=lims, y=lims) + scale_y_continuous(expand=c(0.02,0)) + scale_x_continuous(expand=c(0.02,0)) +
    labs(title=NULL, x="True quantiles", y="Predicted quantiles", color=NULL) +
    facet_wrap(~EQRNN_mod, scales="fixed", ncol=2, strip.position="top")
  if(legend.position=="none"){
    pred_vs_true_plot <- pred_vs_true_plot + theme(strip.text.x = element_blank())
  }
  return(pred_vs_true_plot)
}


plot_predictions_vs_competitor <- function(pred_eqrnn, pred_comp, quantiles_predicted, name_comp="Competitor", xaxis=c("EQRN","competitor")){
  # 
  xaxis <- match.arg(xaxis)
  nb_quantiles_predict <- length(quantiles_predicted)
  pred_vs_comp_plots <- list()
  for(i in 1:nb_quantiles_predict){
    preds_df <- data.frame(competitor=pred_comp[,i], EQRNN=pred_eqrnn[,i])
    lims <- c(min(preds_df), max(preds_df))
    if(xaxis=="competitor"){
      tmp_plt <- preds_df %>% ggplot(aes(x=competitor, y=EQRNN))
    }else{
      tmp_plt <- preds_df %>% ggplot(aes(x=EQRNN, y=competitor))
    }
    pred_vs_comp_plots[[i]] <- tmp_plt +
      geom_abline(intercept=0, slope=1, col="black", size=1) + geom_point() + coord_fixed() + theme(aspect.ratio=1) + 
      expand_limits(x=lims, y=lims) + scale_y_continuous(expand=c(0.02,0)) + scale_x_continuous(expand=c(0.02,0)) +
      labs(title=paste0("Predicted ",quantiles_predicted[i]*100,"% quantile"), x=NULL, y=NULL, color=NULL)
  }
  tmp_plt <- ggpubr::ggarrange(plotlist=pred_vs_comp_plots, labels="AUTO", ncol=min(3,length(pred_vs_comp_plots)),
                               common.legend=TRUE,legend="bottom", vjust=1.5, align="hv")
  if(xaxis=="competitor"){
    pred_vs_true_plot <- tmp_plt %>% ggpubr::annotate_figure(left=text_grob("EQRNN", rot = 90), bottom=text_grob(name_comp))
  }else{
    pred_vs_true_plot <- tmp_plt %>% ggpubr::annotate_figure(left=text_grob(name_comp, rot = 90), bottom=text_grob("EQRNN"))
  }
  return(pred_vs_true_plot)
}

plot_rmse_quantile <- function(fit_eqrn, fit_grf, fit_gbex, fit_egam, Y_train, X_test, quantiles_predict, intermediate_method=c("grf","oracle"),
                               interm_lvl, interm_quant_train, model, distr, df, factorize=TRUE, test_data="other",
                               legend.position="bottom", crop_obs=0, hide_ylab=FALSE, bias_variance=FALSE, crop_Rbv=c(0,0,0)){
  #X_test should be halton to get ISE
  
  intermediate_method <- match.arg(intermediate_method)
  ntest <- nrow(X_test)
  nb_quantiles_predict <- length(quantiles_predict)
  #Predict intermediate quantiles on X_test
  if(intermediate_method=="grf"){
    pred_interm <- predict(fit_grf, newdata = X_test, quantiles = c(interm_lvl))$predictions
  }else if(intermediate_method=="oracle"){
    pred_interm <- generate_theoretical_quantiles(quantiles = c(interm_lvl), X = X_test, model = model, distr = distr, df = df)
  }
  #High quantile prediction with GRF
  pred_grf_test <- predict(fit_grf, newdata = X_test, quantiles = quantiles_predict)$predictions
  # UNCONDITIONAL predicted quantile(s) (Y quantile on X_train)
  pred_unc <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = quantiles_predict, Y = Y_train, ntest = ntest)
  # SEMI-CONDITIONAL predicted quantiles
  pred_semicond <- predict_GPD_semiconditional(Y=Y_train, interm_lvl=interm_lvl, thresh_quantiles=interm_quant_train,
                                               interm_quantiles_test=pred_interm, quantiles_predict=quantiles_predict)
  # GROUND-TRUTH (y_test)
  pred_true <- generate_theoretical_quantiles(quantiles = quantiles_predict, X = X_test, model = model, distr = distr, df = df)
  # Prediction GBEX
  pred_gbex <- gbex_predict(fit_gbex, X_test, to_predict=quantiles_predict, intermediate_quantiles=pred_interm, interm_lvl=interm_lvl)
  pred_exgam <- predict_gpd_gam(fit_egam, X_test, to_predict=quantiles_predict,
                                intermediate_quantiles=pred_interm, interm_lvl=interm_lvl)
  
  #Final EQRN predictions on X_test
  pred_eqrnn <- EQRN_predict(fit_eqrn, X_test, quantiles_predict, pred_interm, interm_lvl)
  
  # Compute losses for desired predicted quantiles
  RMSEs_eqrnn <- sqrt(multilevel_MSE(pred_true,pred_eqrnn,quantiles_predict,give_names=FALSE))
  RMSEs_grf <- sqrt(multilevel_MSE(pred_true,pred_grf_test,quantiles_predict,give_names=FALSE))
  RMSEs_unc <- sqrt(multilevel_MSE(pred_true,pred_unc$predictions,quantiles_predict,give_names=FALSE))
  RMSEs_semicond <- sqrt(multilevel_MSE(pred_true,pred_semicond$predictions,quantiles_predict,give_names=FALSE))
  RMSEs_gbex <- sqrt(multilevel_MSE(pred_true,pred_gbex,quantiles_predict,give_names=FALSE))
  RMSEs_exgam <- sqrt(multilevel_MSE(pred_true,pred_exgam,quantiles_predict,give_names=FALSE))
  
  Q_lvls <- if(factorize){as_factor(roundm(quantiles_predict,4))}else{quantiles_predict}
  
  met_names <- c("Uncond", "Semi-cond", "GRF", "EGAM", "GBEX", "EQRN")
  df_rmse <- data.frame(Quantile=Q_lvls, Uncond=RMSEs_unc, Semi_cond=RMSEs_semicond,
                        GRF=RMSEs_grf, EGAM=RMSEs_exgam, GBEX=RMSEs_gbex, EQRN=RMSEs_eqrnn) %>%
    rename_with(~str_replace(., "Semi_cond", "Semi-cond")) %>% 
    tidyr::gather(key="Model", value="Error", all_of(met_names), factor_key=TRUE)
  
  y_lab_rmse <- if(test_data=="halton"){"RISE"}else{"RMSE"}
  if(hide_ylab){y_lab_rmse <- NULL}
  
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
    bias_eqrnn <- multilevel_pred_bias(pred_true, pred_eqrnn, quantiles_predict, square_bias=FALSE, give_names=FALSE)
    bias_grf <- multilevel_pred_bias(pred_true, pred_grf_test, quantiles_predict, square_bias=FALSE, give_names=FALSE)
    bias_unc <- multilevel_pred_bias(pred_true, pred_unc$predictions, quantiles_predict, square_bias=FALSE, give_names=FALSE)
    bias_semicond <- multilevel_pred_bias(pred_true, pred_semicond$predictions, quantiles_predict, square_bias=FALSE, give_names=FALSE)
    bias_gbex <- multilevel_pred_bias(pred_true, pred_gbex, quantiles_predict, square_bias=FALSE, give_names=FALSE)
    bias_exgam <- multilevel_pred_bias(pred_true, pred_exgam, quantiles_predict, square_bias=FALSE, give_names=FALSE)
    
    rstd_eqrnn <- sqrt(multilevel_resid_var(pred_true, pred_eqrnn, quantiles_predict, give_names=FALSE))
    rstd_grf <- sqrt(multilevel_resid_var(pred_true, pred_grf_test, quantiles_predict, give_names=FALSE))
    rstd_unc <- sqrt(multilevel_resid_var(pred_true, pred_unc$predictions, quantiles_predict, give_names=FALSE))
    rstd_semicond <- sqrt(multilevel_resid_var(pred_true, pred_semicond$predictions, quantiles_predict, give_names=FALSE))
    rstd_gbex <- sqrt(multilevel_resid_var(pred_true, pred_gbex, quantiles_predict, give_names=FALSE))
    rstd_exgam <- sqrt(multilevel_resid_var(pred_true, pred_exgam, quantiles_predict, give_names=FALSE))
    
    Rsqr_eqrnn <- multilevel_R_squared(pred_true, pred_eqrnn, quantiles_predict, give_names=FALSE)
    Rsqr_grf <- multilevel_R_squared(pred_true, pred_grf_test, quantiles_predict, give_names=FALSE)
    Rsqr_unc <- multilevel_R_squared(pred_true, pred_unc$predictions, quantiles_predict, give_names=FALSE)
    Rsqr_semicond <- multilevel_R_squared(pred_true, pred_semicond$predictions, quantiles_predict, give_names=FALSE)
    Rsqr_gbex <- multilevel_R_squared(pred_true, pred_gbex, quantiles_predict, give_names=FALSE)
    Rsqr_exgam <- multilevel_R_squared(pred_true, pred_exgam, quantiles_predict, give_names=FALSE)
    
    Rbv_names <- c("Quantile R squared", "Bias", "Residual standard deviation")
    
    df_Rbv <- data.frame(Quantile=rep(Q_lvls, 3),
                          Uncond=c(Rsqr_unc,bias_unc,rstd_unc),
                          Semi_cond=c(Rsqr_semicond,bias_semicond,rstd_semicond),
                          GRF=c(Rsqr_grf,bias_grf,rstd_grf),
                          EGAM=c(Rsqr_exgam,bias_exgam,rstd_exgam),
                          GBEX=c(Rsqr_gbex,bias_gbex,rstd_gbex),
                          EQRN=c(Rsqr_eqrnn,bias_eqrnn,rstd_eqrnn),
                          metric=factor(rep(Rbv_names, each=nb_quantiles_predict), levels=Rbv_names)) %>%
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
    # Rbv_plot <- df_Rbv %>% ggplot( aes(x=Quantile, y=Error, group=Model, color=Model, linetype=Model)) +
    #   geom_line(size=1) + scale_linetype_manual(values=linetypes) + scale_color_manual(values=color_scale) +
    #   scale_y_continuous(expand=c(0.02,0)) + labs(title=NULL, x="Probability level", y=NULL, color=NULL, linetype=NULL) +
    #   facet_wrap(~metric, scales="free_y", ncol=3, strip.position="left") + theme(legend.position=legend.position)
    # Rbv_plot <- Rbv_plot + if(factorize){scale_x_discrete(expand=c(0.02,0))}else{scale_x_continuous(expand=c(0.02,0))}
    # Rbv_plot <- Rbv_plots[[1]] + Rbv_plots[[2]] + Rbv_plots[[3]] +
    #   patchwork::plot_layout(ncol=3, nrow=1, guides='collect') & theme(legend.position=legend.position)
    Rbv_plot <- ggpubr::ggarrange(plotlist=Rbv_plots, nrow=1, ncol=3, labels=NULL,
                                  common.legend=TRUE, legend="bottom", align="hv")
    
    return(list(rmse_plot=mse_plot, Rbv_plot=Rbv_plot))
  }
}



plot_predictions_2D <- function(pred_eqrnn_test, pred_gbex_test, pred_grf_test, pred_exgam_test,
                                pred_semicond, ground_truth, X_test, quantile_predicted, 
                                legend.position="bottom", show_title=FALSE){
  
  title_q <- if(show_title){paste0("Predictions (q=",quantile_predicted,")")}else{NULL}
  
  df <- data.frame(X1=X_test[,1], X2=X_test[,2], EQRNN=pred_eqrnn_test,
                   gbex=pred_gbex_test, GRF=pred_grf_test, ExGAM=pred_exgam_test,
                   Semi_cond=pred_semicond, Truth=ground_truth) %>% 
    tidyr::gather(key="Model", value="Quantile", EQRNN, gbex, GRF, ExGAM, Semi_cond, Truth, factor_key=TRUE)
  
  pred2D <- df %>% ggplot(aes(X1, X2, z = Quantile)) + theme_minimal() +
    geom_raster(aes(fill = Quantile)) +
    scale_fill_continuous(name = latex2exp::TeX(paste0(r"(Quantile ($\tau=)",quantile_predicted,"$):\t"))) +
    coord_fixed() +
    theme(aspect.ratio=1, legend.position=legend.position, axis.ticks=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), legend.background=element_blank(),
          strip.background=element_blank(),strip.placement="outside",
          strip.text=element_text(size=12), strip.switch.pad.grid=unit(0, "pt"),
          strip.switch.pad.wrap=unit(0, "pt"),
          plot.caption=element_text(size=7.5,hjust=0,margin=margin(t=15)),
          text=element_text(size=11), axis.text=element_text(size=10),#
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.box.spacing=margin(0.5),
          plot.margin=margin(0,0,0,0, unit="pt")) +
    scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0)) +
    labs(title=title_q, x=latex2exp::TeX("$X_1$"), y=latex2exp::TeX("$X_2$")) +
    facet_wrap(~Model, scales="fixed", ncol=3, nrow=2, strip.position="top")
  return(pred2D)
}


plot_validation_vs_resid <- function(results_tibble, groupby="hidden_fct", add_title="", loss="Valid_loss", MSAE_prefix="test_"){
  y_label <- if(loss=="Valid_loss"){"Validation loss"}else{str_replace_all(loss,"_"," ")}
  val_res_plt <- list()
  results_tibble <- results_tibble %>% mutate_at(vars("test_MSE_q0.995","test_MSE_q0.999","test_MSE_q0.9995"), list(~ sqrt(.)))
  val_res_plt[[1]] <- results_tibble %>% ggplot(aes_string(x=paste0(MSAE_prefix,"MSE_q0.995"), y=loss, group=groupby, color=groupby)) +
    geom_point() + labs(y=NULL, x=paste0(str_replace_all(MSAE_prefix,"_"," "),"RMSE (q=0.995)"))
  val_res_plt[[2]] <- results_tibble %>% ggplot(aes_string(x=paste0(MSAE_prefix,"MSE_q0.999"), y=loss, group=groupby, color=groupby)) +
    geom_point() + labs(y=NULL, x=paste0(str_replace_all(MSAE_prefix,"_"," "),"RMSE (q=0.999)"))
  val_res_plt[[3]] <- results_tibble %>% ggplot(aes_string(x=paste0(MSAE_prefix,"MSE_q0.9995"), y=loss, group=groupby, color=groupby)) +
    geom_point() + labs(y=NULL, x=paste0(str_replace_all(MSAE_prefix,"_"," "),"RMSE (q=0.9995)"))
  
  val_res_plt[[4]] <- results_tibble %>% ggplot(aes_string(x=paste0(MSAE_prefix,"MAE_q0.995"), y=loss, group=groupby, color=groupby)) +
    geom_point() + labs(y=NULL, x=paste0(str_replace_all(MSAE_prefix,"_"," "),"MAE (q=0.995)"))
  val_res_plt[[5]] <- results_tibble %>% ggplot(aes_string(x=paste0(MSAE_prefix,"MAE_q0.999"), y=loss, group=groupby, color=groupby)) +
    geom_point() + labs(y=NULL, x=paste0(str_replace_all(MSAE_prefix,"_"," "),"MAE (q=0.999)"))
  val_res_plt[[6]] <- results_tibble %>% ggplot(aes_string(x=paste0(MSAE_prefix,"MAE_q0.9995"), y=loss, group=groupby, color=groupby)) +
    geom_point() + labs(y=NULL, x=paste0(str_replace_all(MSAE_prefix,"_"," "),"MAE (q=0.9995)"))
  
  validation_vs_resid <- ggpubr::ggarrange(plotlist=val_res_plt, labels=NULL, ncol=3, nrow=2,
                                           common.legend=TRUE,legend="bottom", vjust=1.5, align="hv") %>%
    ggpubr::annotate_figure(top=text_grob(paste0("Validation loss-error relationship ", add_title)), left=text_grob(y_label, rot = 90))
  return(validation_vs_resid)
}

YvsX_scatter <- function(results_tibble, x_name, y_name, groupby=NULL, title=NULL){
  YvsX_plt  <- results_tibble %>% ggplot(aes_string(x=x_name, y=y_name, group=groupby, color=groupby)) +
    geom_point() + labs(y=y_name, x=x_name, title=title)
  return(YvsX_plt)
}

YvsX_scatter_vect <- function(x, y, x_name, y_name, groupby=NULL, title=NULL){
  YvsX_plt  <- data.frame(x=x, y=y) %>% ggplot(aes_string(x="x", y="y", group=groupby, color=groupby)) +
    geom_point() + labs(y=y_name, x=x_name, title=title)
  return(YvsX_plt)
}

plot_exceedences_quantile <- function(pred_eqrnn, y, quantile_levels, log_scaling=FALSE, factorize=TRUE, show_legend=FALSE){
  
  n <- nrow(pred_eqrnn)
  nb_quantiles_predict <- length(quantile_levels)
  if(n!=length(y)){stop("Observation number mismatch between predictions and y in 'plot_exceedences_quantile'.")}
  if(nb_quantiles_predict!=ncol(pred_eqrnn)){stop("Quantile levels mismatch between predictions and quantile_levels in 'plot_exceedences_quantile'.")}
  
  if(factorize){
    df <- data.frame(Quantile=as_factor(quantile_levels), Expected=(1-quantile_levels)*n, EQRNN=colSums(c(y)>pred_eqrnn))
  } else {
    df <- data.frame(Quantile=quantile_levels, Expected=(1-quantile_levels)*n, EQRNN=colSums(c(y)>pred_eqrnn))
  }
  logp1_ticks <- c(1,2,5,10,20,50,100,200,500,1000,2000,5000)
  trans_logp1 <- scales::trans_new(name="logp1", transform=function(x){log(x*100+1)}, inverse=function(x){(exp(x)-1)/100}, domain=c(-0,Inf),
                                   breaks=function(a)logp1_ticks, format=scales::number_format(accuracy=0.1))
  
  linetypes <- c(Expected="solid",EQRNN="dotted")
  color_scale <- c(Expected="black",EQRNN=my_palette_methods[["EQRNN"]])
  
  nex_plt <- df %>% tidyr::gather(key="Exceedences", value="Number", Expected, EQRNN, factor_key=TRUE) %>%
    ggplot( aes(x=Quantile, y=Number, group=Exceedences, color=Exceedences, linetype=Exceedences)) +
    geom_line(size=1) + geom_point() + scale_linetype_manual(values=linetypes) + scale_color_manual(values=color_scale) +
    labs(title=NULL, x="Quantile level", y="Number of exceedences", color=NULL, linetype=NULL) + 
    if(factorize){scale_x_discrete(expand=c(0.02,0))}else{scale_x_continuous(expand=c(0.02,0))}
  if(log_scaling){nex_plt <- nex_plt + scale_y_continuous(trans=trans_logp1, labels=logp1_ticks, expand=c(0.02,0))}
  else{nex_plt <- nex_plt + scale_y_continuous(expand=c(0.02,0))}
  if(!show_legend){nex_plt <- nex_plt + theme(legend.position="none")}
  return(nex_plt)
}

plot_rmse_threshold <- function(result_tbl, interm_lvls="interm_lvl", predict_lvls="Prob_lvl", errs="test_MSE", err_sd=paste0(errs,"_sd"),
                                factorize=TRUE, test_data="other", strip.position="top", hide_ylab=FALSE){
  
  if(factorize){result_tbl[interm_lvls] <- as_factor(roundm(result_tbl[[interm_lvls]],4))}
  
  result_tbl$sd_max <- sqrt(result_tbl[[errs]] + result_tbl[[err_sd]])
  result_tbl$sd_min <- sqrt(pmax(result_tbl[[errs]] - result_tbl[[err_sd]], 0))
  result_tbl[errs] <- sqrt(result_tbl[[errs]])
  
  y_lab_rmse <- if(test_data=="halton"){"RISE"}else{"RMSE"}
  if(hide_ylab){y_lab_rmse <- NULL}
  
  rmse_plot <- result_tbl %>% 
    ggplot( aes_string(x=interm_lvls, y=errs)) +
    geom_errorbar(aes(ymin=sd_min, ymax=sd_max), color=my_palette_methods["EQRN"], width=.2) + #, position=position_dodge(0.05)
    geom_point(size=2, color=my_palette_methods["EQRN"]) + geom_line(size=1, color=my_palette_methods["EQRN"]) +
    scale_y_continuous(expand=c(0.02,0)) + labs(title=NULL, x="Intermediate probability level", y=y_lab_rmse) +
    facet_wrap(~.data[[predict_lvls]], scales="free_y", ncol=3, strip.position=strip.position)# + theme(legend.position=legend.position)
  
  rmse_plot <- rmse_plot + if(factorize){scale_x_discrete(expand=c(0.02,0))}else{scale_x_continuous(expand=c(0.02,0))}
  # if(crop_obs>0){rmse_plot <- rmse_plot + coord_cartesian(ylim=c(min(select(result_tbl, errs)), sort(unlist(select(result_tbl,errs)), decreasing=TRUE)[[2]]))}
  if(is.null(strip.position)){rmse_plot <- rmse_plot + theme(strip.text.x = element_blank())}
  return(rmse_plot)
}

plot_rse_box_threshold <- function(residuals_tbl, interm_lvls="interm_lvl", predict_lvls="Prob_lvl", residuals="test_resid",
                                   factorize=TRUE, test_data="other", strip.position="top", hide_ylab=FALSE){
  
  if(factorize){residuals_tbl[interm_lvls] <- as_factor(roundm(residuals_tbl[[interm_lvls]],4))}
  
  # residuals_tbl$SEs <- residuals_tbl[[residuals]]^2
  residuals_tbl$AbsRes <- abs(residuals_tbl[[residuals]])
  
  rootmeansquared <- function(x, ...){sqrt(mean(x^2, ...))}
  
  y_lab_rmse <- if(test_data=="halton"){"RISE"}else{"RMSE"}
  if(hide_ylab){y_lab_rmse <- NULL}
  
  boxplot_col <- my_palette_methods["EQRN"]
  rmse_col <- my_palette$red
  
  rse_plot <- residuals_tbl %>% 
    ggplot(aes_string(x=interm_lvls, y="AbsRes", group=interm_lvls)) +
    stat_boxplot(geom = "errorbar", width=0.5, color=boxplot_col) +
    geom_boxplot(width=0.75, color=boxplot_col, alpha=0.6) +
    stat_summary(fun=rootmeansquared, geom="point", shape=18, size=4, color=rmse_col, fill=rmse_col) +
    scale_y_continuous(expand=c(0.02,0)) + labs(title=NULL, x="Intermediate probability level", y=y_lab_rmse) +
    facet_wrap(~.data[[predict_lvls]], scales="free_y", ncol=3, strip.position=strip.position)# + theme(legend.position=legend.position)
  
  rse_plot <- rse_plot + if(factorize){scale_x_discrete(expand=c(0.02,0))}else{scale_x_continuous(expand=c(0.02,0))}
  # if(crop_obs>0){rse_plot <- rse_plot + coord_cartesian(ylim=c(min(select(residuals_tbl, residuals)), sort(unlist(select(residuals_tbl,residuals)), decreasing=TRUE)[[2]]))}
  if(is.null(strip.position)){rse_plot <- rse_plot + theme(strip.text.x = element_blank())}
  return(rse_plot)
}


# ================================= DEPRECATED =================================


plot_error_quantile <- function(fit_eqrn, fit_grf, fit_gbex, fit_egam, Y_train, X_test, quantiles_predict, intermediate_method=c("grf","oracle"),
                                interm_lvl, interm_quant_train, model, distr, df, factorize=TRUE, test_data="other"){
  #X_test should ba halton to get ISE
  
  intermediate_method <- match.arg(intermediate_method)
  ntest <- nrow(X_test)
  nb_quantiles_predict <- length(quantiles_predict)
  #Predict intermediate quantiles on X_test
  if(intermediate_method=="grf"){
    pred_interm <- predict(fit_grf, newdata = X_test, quantiles = c(interm_lvl))$predictions
  }else if(intermediate_method=="oracle"){
    pred_interm <- generate_theoretical_quantiles(quantiles = c(interm_lvl), X = X_test, model = model, distr = distr, df = df)
  }
  #High quantile prediction with GRF
  pred_grf_test <- predict(fit_grf, newdata = X_test, quantiles = quantiles_predict)$predictions
  # UNCONDITIONAL predicted quantile(s) (Y quantile on X_train)
  pred_unc <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = quantiles_predict, Y = Y_train, ntest = ntest)
  # SEMI-CONDITIONAL predicted quantiles
  pred_semicond <- predict_GPD_semiconditional(Y=Y_train, interm_lvl=interm_lvl, thresh_quantiles=interm_quant_train,
                                               interm_quantiles_test=pred_interm, quantiles_predict=quantiles_predict)
  # GROUND-TRUTH (y_test)
  pred_true <- generate_theoretical_quantiles(quantiles = quantiles_predict, X = X_test, model = model, distr = distr, df = df)
  # Prediction GBEX
  pred_gbex <- gbex_predict(fit_gbex, X_test, to_predict=quantiles_predict, intermediate_quantiles=pred_interm, interm_lvl=interm_lvl)
  pred_exgam <- predict_gpd_gam(fit_egam, X_test, to_predict=quantiles_predict,
                                intermediate_quantiles=pred_interm, interm_lvl=interm_lvl)
  
  #Final EQRNN predictions on X_test
  pred_eqrnn <- EQRN_predict(fit_eqrn, X_test, quantiles_predict, pred_interm, interm_lvl)
  
  # Compute losses for desired predicted quantiles
  MAEs_eqrnn <- multilevel_MAE(pred_true,pred_eqrnn,quantiles_predict,give_names=FALSE)
  RMSEs_eqrnn <- sqrt(multilevel_MSE(pred_true,pred_eqrnn,quantiles_predict,give_names=FALSE))
  MAEs_grf <- multilevel_MAE(pred_true,pred_grf_test,quantiles_predict,give_names=FALSE)
  RMSEs_grf <- sqrt(multilevel_MSE(pred_true,pred_grf_test,quantiles_predict,give_names=FALSE))
  MAEs_unc <- multilevel_MAE(pred_true,pred_unc$predictions,quantiles_predict,give_names=FALSE)
  RMSEs_unc <- sqrt(multilevel_MSE(pred_true,pred_unc$predictions,quantiles_predict,give_names=FALSE))
  MAEs_semicond <- multilevel_MAE(pred_true,pred_semicond$predictions,quantiles_predict,give_names=FALSE)
  RMSEs_semicond <- sqrt(multilevel_MSE(pred_true,pred_semicond$predictions,quantiles_predict,give_names=FALSE))
  MAEs_gbex <- multilevel_MAE(pred_true,pred_gbex,quantiles_predict,give_names=FALSE)
  RMSEs_gbex <- sqrt(multilevel_MSE(pred_true,pred_gbex,quantiles_predict,give_names=FALSE))
  MAEs_exgam <- multilevel_MAE(pred_true,pred_exgam,quantiles_predict,give_names=FALSE)
  RMSEs_exgam <- sqrt(multilevel_MSE(pred_true,pred_exgam,quantiles_predict,give_names=FALSE))
  
  Q_lvls <- if(factorize){as_factor(roundm(quantiles_predict,4))}else{quantiles_predict}
  df_rmse <- data.frame(Quantile=Q_lvls, Uncond=RMSEs_unc, Semi_cond=RMSEs_semicond,
                        GRF=RMSEs_grf, gpd_GAM=RMSEs_exgam, gbex=RMSEs_gbex, EQRNN=RMSEs_eqrnn)
  df_mae <- data.frame(Quantile=Q_lvls, Uncond=MAEs_unc, Semi_cond=MAEs_semicond,
                       GRF=MAEs_grf, gpd_GAM=MAEs_exgam, gbex=MAEs_gbex, EQRNN=MAEs_eqrnn)
  
  y_lab_rmse <- if(test_data=="halton"){"RISE"}else{"RMSE"}
  
  linetypes <- c("dotted","dotted","dashed","solid","solid","solid")
  color_scale <- c("#ff7f00","#984ea3","#4daf4a","#69b3a2","#e41a1c","#377eb8")
  
  mse_plot <- df_rmse %>% tidyr::gather(key="Model", value="Error", Uncond, Semi_cond, GRF, gpd_GAM, gbex, EQRNN, factor_key=TRUE) %>%
    ggplot( aes(x=Quantile, y=Error, group=Model, color=Model, linetype=Model)) +
    geom_line(size=1) + scale_linetype_manual(values=linetypes) + scale_color_manual(values=color_scale) +
    scale_y_continuous(expand=c(0.02,0)) + labs(title=NULL, x="Probability level", y=y_lab_rmse, color=NULL, linetype=NULL)
  mae_plot <- df_mae %>% tidyr::gather(key="Model", value="Error", Uncond, Semi_cond, GRF, gpd_GAM, gbex, EQRNN, factor_key=TRUE) %>%
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

plot_error_distrib_results <- function(results_tibble, quantile_predicted){
  resid_box <- data.frame(Uncond_MSEs=results_tibble[[paste0("Uncond_test_MSE_q", quantile_predicted)]],
                          SemiCond_MSEs=results_tibble[[paste0("Semicond_test_MSE_q", quantile_predicted)]],
                          GRF_MSEs=results_tibble[[paste0("GRF_test_MSE_q", quantile_predicted)]],
                          EQRNN_MSEs=results_tibble[[paste0("EQRNN_test_MSE_q", quantile_predicted)]],
                          Uncond_MAEs=results_tibble[[paste0("Uncond_test_MAE_q", quantile_predicted)]],
                          SemiCond_MAEs=results_tibble[[paste0("Semicond_test_MAE_q", quantile_predicted)]],
                          GRF_MAEs=results_tibble[[paste0("GRF_test_MAE_q", quantile_predicted)]],
                          EQRNN_MAEs=results_tibble[[paste0("EQRNN_test_MAE_q", quantile_predicted)]]) %>%
    tidyr::gather(key="Method", value="Value", factor_key=TRUE) %>% ggplot( aes(x=Method, y=Value, group=Method, fill=Method)) +
    geom_abline(intercept=0, slope=0, col="white", size=1) + geom_violin(width=1, color="grey", alpha=0.5) +
    stat_boxplot(geom = "errorbar", width=0.5) + geom_boxplot(width=0.75, color="black", alpha=0.8)+ coord_flip() +
    stat_summary(fun=mean, geom="point", shape=4, size=3, color="black", fill="black") +
    labs(title = paste0("Quantile Residuals (q=",quantile_predicted,")"), x=NULL, y=NULL) +
    theme(legend.position="none", panel.grid.major.y = element_blank())
  return(resid_box)
}

