library(tidyverse)
library(grf)
library(here)
library(ggpubr)
library(evd)
library(ismev)
library(future)
library(doFuture)
source("R/gbex_wrappers.R")
source("R/EQRN_loader.R")

save_path <- "Results/gbex/cosnorm_sigmoid_n7e3_p10_q80_oracle_Qfeat/"
check_directory(save_path, recursive=TRUE)
parallel_strat <- "multisession"#"sequential", "multisession", "multicore"
n_workers <- 7#(availableCores() - 1)#PARAM
err_handling <- "remove"#"stop", "remove", "pass"

seedR <- 0
seedgbex <- 1
set.seed(seedR)

# PARAM: Data
n <- 5e3
p <- 10
ntest <- 5e3
n_valid <- 2e3
model = "cosnorm_sigmoid"#"binorm_sigmoid", "cosnorm2d_sigmoid", "cosnorm_sigmoid"
distr = "student_t"
df <- 4
test_data = "halton"

# PARAM: General
interm_lvl = 0.8
prob_lvls_predict = c(0.995,0.999,0.9995)

# PARAM: gbex
intermediate_q_feature=TRUE
scale_features=FALSE
num_folds=5
Bmax=1000
grid_lambda_ratio=c(5,6,7,8,9,10)
grid_depth=list(c(1,0),c(1,1),c(2,1),c(2,2),c(3,1),c(3,2),c(3,3))
stratified=TRUE
lambda_scale=0.01
min_leaf_size=NULL
sf=0.75

start_time <- Sys.time()

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


#Intermediate quantiles prediction on dat$X
intermediate_quantiles <- generate_theoretical_quantiles(quantiles = c(interm_lvl), X = dat$X, model = model, distr = distr, df = df)

valid_quantiles <- generate_theoretical_quantiles(quantiles = c(interm_lvl), X = dat_valid$X, model = model, distr = distr, df = df)

## ======== TESTING (GRF and Ground truth) ========

#Predict intermediate quantiles on X_test
pred_interm <- generate_theoretical_quantiles(quantiles = c(interm_lvl), X = X_test, model = model, distr = distr, df = df)

CV_results <- gbex_CV(X=rbind(dat$X,dat_valid$X), y=c(dat$Y,dat_valid$Y), intermediate_quantiles=rbind(intermediate_quantiles,valid_quantiles),
                      interm_lvl=interm_lvl, intermediate_q_feature=intermediate_q_feature, scale_features=scale_features,
                      num_folds=num_folds, Bmax=Bmax, grid_lambda_ratio=grid_lambda_ratio, grid_depth=grid_depth,
                      stratified=stratified, lambda_scale=lambda_scale, min_leaf_size=min_leaf_size, sf=sf,
                      parallel_strat=parallel_strat, n_workers=n_workers, seed=seedgbex, err_handling=err_handling)


filename <- paste0(save_path,"results_gbex_CV_",format(Sys.time(),'%Y%m%d_%H%M%S'),".csv")

write_csv(CV_results$results_tibble, file=filename)

cat(paste0("\n Best params:\n",CV_results$best_params,"\n\n"))

end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)

