library(tidyverse)
library(grf)
library(here)
library(ggpubr)
library(evd)
library(ismev)
library(torch)
library(future)
library(doFuture)
source("R/EQRN_loader.R")

## =============================== PARAMETERS ===============================

# param_file <- "main/EQRN_iid/parameters/grid_search/binorm_sigmoid_oracle.R"
param_file <- "main/EQRN_iid/parameters/grid_search/cosnorm2d_sigmoid_oracle.R"
# param_file <- "main/EQRN_iid/parameters/grid_search/cosnorm2d_sigmoid_oracle_SSNN.R"
# param_file <- "main/EQRN_iid/parameters/grid_search/cosnorm_sigmoid_oracle.R"


parallel_strat <- "multisession"#"sequential", "multisession", "multicore"
n_workers <- 12#(availableCores() - 1)#PARAM

## =============================== END PARAMETERS ===============================

source(param_file)
check_directory(save_path, recursive=TRUE)
set.seed(seedR)

start_time <- Sys.time()

# Training data
dat <- generate_joint_distribution(n = n, p = p, model = model, distr = distr, df = df)
dat_valid <- generate_joint_distribution(n = n_valid, p = p, model = model, distr = distr, df = df)
X_train_all <- rbind(dat$X,dat_valid$X)
y_train_all <- c(dat$Y,dat_valid$Y)

# Generate test data
if (test_data == "grid"){
  ntest <- floor(sqrt(ntest)) ** 2
  warning(paste0("Modified ntest to: ", ntest))
}
X_test <- generate_test_data(ntest, p, test_data)
y_test <- generate_conditional_distribution(model, distr, df, X_test)$Y

# GRF fit
fit_grf <- quantile_forest(dat$X, dat$Y, num.trees=num.trees, quantiles=quantiles_fit, sample.fraction=sample.fraction, mtry=mtry,
                           min.node.size=min.node.size, honesty=honesty, honesty.fraction=honesty.fraction, honesty.prune.leaves=honesty.prune.leaves,
                           alpha=alpha, imbalance.penalty=imbalance.penalty, seed=seedGRF)

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
  #Intermediate quantiles prediction on dat$X
  intermediate_quantiles <- generate_theoretical_quantiles(quantiles = c(interm_lvl), X = dat$X, model = model, distr = distr, df = df)
  valid_quantiles <- generate_theoretical_quantiles(quantiles = c(interm_lvl), X = dat_valid$X, model = model, distr = distr, df = df)
  interm_quantiles_all <- rbind(intermediate_quantiles,valid_quantiles)
  #Predict intermediate quantiles on X_test
  pred_interm <- generate_theoretical_quantiles(quantiles = c(interm_lvl), X = X_test, model = model, distr = distr, df = df)
  pred_interm_fromall <- pred_interm
}

## ======== TESTING (GRF and Ground truth) ========

#High quantile prediction with GRF
pred_grf_test <- predict(fit_grf, newdata = X_test, quantiles = prob_lvls_predict)$predictions

# UNCONDITIONAL predicted quantile(s) (Y quantile on X_train)
pred_unc <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = prob_lvls_predict, Y = dat$Y, ntest = ntest)

# SEMI-CONDITIONAL predicted quantiles
pred_semicond <- predict_GPD_semiconditional(Y=dat$Y, interm_lvl=interm_lvl, thresh_quantiles=intermediate_quantiles,
                                             interm_quantiles_test=pred_interm, prob_lvls_predict=prob_lvls_predict)

# GROUND-TRUTH (y_test)
pred_true <- generate_theoretical_quantiles(quantiles = prob_lvls_predict, X = X_test, model = model, distr = distr, df = df)

# Unconditional parameters and losses
uncond_losses_fixed <- unconditional_train_valid_GPD_loss(Y_train=dat$Y, interm_lvl=interm_lvl, Y_valid=dat_valid$Y)
uncond_losses_interm <- semiconditional_train_valid_GPD_loss(Y_train=dat$Y, Y_valid=dat_valid$Y,
                                                                             interm_quant_train=intermediate_quantiles,
                                                                             interm_quant_valid=valid_quantiles)

## ======== EQRN FIT ========

grid_l <- purrr::cross(params_list, .filter = NULL)

`%fun%` <- set_doFuture_strategy(parallel_strat, n_workers=n_workers)

results_grid_fit <- foreach(params=grid_l, .errorhandling="remove", .combine=rbind) %fun% {
  source("R/EQRN_loader.R")
  
  stuct_str <- paste0(params$net_structure,"h", collapse="_")
  if(is.character(params$hidden_fct)){
    if(params$hidden_fct=="SSNN"){
      stuct_str <- paste0("sc", paste0(params$net_structure$scale, collapse="_"), "_sh", paste0(params$net_structure$shape, collapse="_"))
    }
  }
  params_string <- paste0(stuct_str, "_", hid_f_str, "_", "n"[!params$intermediate_q_feature], "u_do", params$p_drop*100,
                          "_L2", str_replace(toString(params$L2_pen),"([.])","p"), "_shp",
                          str_replace(toString(params$shape_penalty),"([.])","p"), "_", lr_str)
  cat("======== Start: ", params_string, " ========\n")
  
  #Fit EQRN with intermediate OOB quantiles
  torch_manual_seed(seedT)
  fit_eqrn <- EQRN_fit_restart(dat$X, dat$Y, intermediate_quantiles, interm_lvl=interm_lvl, number_fits=nb_fits_eqrn, shape_fixed=params$shape_fixed,
                                net_structure=params$net_structure, hidden_fct=params$hidden_fct, p_drop=params$p_drop, intermediate_q_feature=params$intermediate_q_feature,
                                learning_rate=params$learning_rate, L2_pen=params$L2_pen, shape_penalty=params$shape_penalty, scale_features=params$scale_features,
                                n_epochs=params$n_epochs, batch_size=params$batch_size, X_valid=dat_valid$X, y_valid=dat_valid$Y, quant_valid=valid_quantiles,
                                lr_decay=params$lr_decay, patience_decay=params$patience_decay, min_lr=params$min_lr, patience_stop=params$patience_stop,
                                orthogonal_gpd=params$orthogonal_gpd)
  
  train_plot <- training_plot_eqrn(fit_eqrn, uncond_losses_interm)
  if(save_plots){
    ggsave(paste0(params_string, "_training.png"), plot=train_plot, device="png",
           path=save_path, width=200, height=150, units="mm", dpi=300)
  }
  
  ## ======== TESTING ========
  
  #Final EQRN predictions on X_test
  pred_eqrn <- EQRN_predict(fit_eqrn, X_test, prob_lvls_predict, pred_interm, interm_lvl)
  
  # Compute losses for desired predicted quantiles
  MSE_losses <- multilevel_MSE(pred_true, pred_eqrn, prob_lvls_predict, prefix="test_", give_names=TRUE)
  MAE_losses <- multilevel_MAE(pred_true, pred_eqrn, prob_lvls_predict, prefix="test_", give_names=TRUE)
  
  nb_prob_lvls_predict <- length(prob_lvls_predict)
  for(i in 1:nb_prob_lvls_predict){
    resid_box <- residuals_boxplot_eqrn(pred_eqrn[,i], pred_grf_test[,i], pred_true[,i], prob_lvls_predict[i])
    if(save_plots){
      ggsave(paste0(params_string, "_q",prob_lvls_predict[i]*10000,"_residbox.png"), plot=resid_box, device="png",
             path=save_path, width=250, height=100, units="mm", dpi=300)
    }
  }
  
  test_loss <- compute_EQRN_GPDLoss(fit_eqrn, X_test, y_test, intermediate_quantiles=pred_interm, interm_lvl=interm_lvl)
  output <- c(shape_fixed=params$shape_fixed,
              net_structure=stuct_str,
              hidden_fct=hid_f_str,
              p_drop=params$p_drop,
              intermediate_q_feature=params$intermediate_q_feature,
              L2_pen=params$L2_pen,
              shape_penalty=params$shape_penalty,
              Train_loss=fit_eqrn$train_loss[length(fit_eqrn$train_loss)], Valid_loss=fit_eqrn$valid_loss[length(fit_eqrn$valid_loss)],
              test_loss=test_loss, MAE_losses, MSE_losses,
              learning_rate=params$learning_rate,
              n_epochs=params$n_epochs,
              lr_decay=params$lr_decay,
              patience_decay=params$patience_decay,
              min_lr=params$min_lr,
              patience_stop=params$patience_stop,
              orthogonal_gpd=params$orthogonal_gpd)
  RESULTS <- list(pred_true=pred_true, pred_eqrn=pred_eqrn, pred_grf_test=pred_grf_test, prob_lvls_predict=prob_lvls_predict,
                  train_loss=fit_eqrn$train_loss, valid_loss=fit_eqrn$valid_loss,
                  uncond_losses_interm=uncond_losses_interm, params_string=params_string, output=output)
  safe_save_rds(RESULTS, paste0(save_path, "Results_",params_string,".rds"))
  cat("==== End: ", params_string, " ====\n")
  output
}
results_tibble <- tibble::as_tibble(results_grid_fit)
results_tibble <- bind_cols(model=model, distr=distr, n=n, p=p, df=df, interm_lvl=interm_lvl, results_tibble)

filename <- paste0(save_path,"results_EQRN_grid_fit_",format(Sys.time(),'%Y%m%d_%H%M%S'),".csv")

write_csv(results_tibble, file=filename)

end_doFuture_strategy()

end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)


#Sys.sleep(300)
#system('shutdown -t 30 -s')

