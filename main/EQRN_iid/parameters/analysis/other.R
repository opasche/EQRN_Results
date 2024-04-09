
## =============================== PARAMETERS ===============================
save_path <- "Results/iid/EQRN_analysis/test_new/"

seedR <- 0
seedGRF <- 1
seedT <- seedR
options(torch.old_seed_behavior=TRUE)
force_refit <- FALSE

# PARAM: Data
n <- 5e3
p <- 10
ntest <- 5e3
n_valid <- 2e3
model = "cosnorm2d_sigmoid"# "binorm_sigmoid", "cosnorm2d_sigmoid", "cosnorm_sigmoid"
distr = "student_t"
df <- 4
test_data = "halton"
path_data <- paste0("data/simulations/iid_data_interm/",model,"/")

# PARAM: General
intermediate_method <- "oracle"
interm_lvl = 0.8
prob_lvls_predict = c(0.995,0.999,0.9995)#c(interm_lvl,0.995,0.999,0.9995)

# Params: GRF
num.trees = 5e3
quantiles_fit = c(0.1, 0.5, 0.9)
sample.fraction = 0.5
mtry = min(ceiling(sqrt(p) + 20), p)
min.node.size = 5
honesty = TRUE
honesty.fraction = 0.5
honesty.prune.leaves = TRUE
alpha = 0.05
imbalance.penalty = 0

# Params: EQRN
path_eqrn <- paste0("Results/iid/EQRN_best_models_saves/",model,"/",intermediate_method,"/")#save_path
nb_fits_eqrn <- 3
shape_fixed=FALSE
net_structure=c(5,3,3)
hidden_fct=torch_tanh
p_drop=0
intermediate_q_feature=TRUE
learning_rate=1e-4
L2_pen=0
shape_penalty=0
scale_features=TRUE
orthogonal_gpd=TRUE
n_epochs=1000
batch_size=256
lr_decay=0.4
patience_decay=5
min_lr=1e-5
patience_stop=20

hid_f_str <- "tanh"
lr_str <- "lrd4"

# Params: gbex
gbex_params <- list(
  B=125,
  lambda=NULL,
  lambda_ratio=7,
  lambda_scale=0.01,
  depth=c(3, 1),
  sf=0.75,
  intermediate_q_feature=TRUE,
  scale_features=FALSE)

# Params: egam
egam_params <- list(
  model_shape=TRUE,
  intermediate_q_feature=gbex_params$intermediate_q_feature,
  scale_features=gbex_params$scale_features)

## =============================== END PARAMETERS ===============================
