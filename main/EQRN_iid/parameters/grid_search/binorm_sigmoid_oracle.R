
## =============================== PARAMETERS ===============================
save_path <- "Results/iid/EQRN_grid_search/binorm_sigmoid_n5e3_p10_q80_oracle_ortho_tanh/"
save_plots <- TRUE

seedR <- 0
seedGRF <- 1
seedT <- seedR
options(torch.old_seed_behavior=TRUE)

# PARAM: Data
n <- 5e3
p <- 10
ntest <- 5e3
n_valid <- 2e3
model = "binorm_sigmoid"# "binorm_sigmoid", "cosnorm2d_sigmoid", "cosnorm_sigmoid"
distr = "student_t"
df <- 4
test_data = "halton"

# PARAM: General
intermediate_method <- "oracle"
interm_lvl = 0.8
quantiles_predict = c(0.995,0.999,0.9995)

# PARAM: GRF
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

# PARAM: EQRN
nb_fits_eqrn <- 3
params_list <- list(
  shape_fixed=FALSE,
  net_structure=list(c(5,3,3),c(20,10),c(10,5,5),c(10,10,10),c(20,10,10),c(64,64,64),c(128,128,128),c(32,32,32,32),c(256,256,256,256)),
  hidden_fct=list(torch_tanh),
  p_drop=list(0),
  intermediate_q_feature=list(TRUE),
  learning_rate=list(1e-3,1e-4),
  L2_pen=list(0,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3),
  shape_penalty=list(0),
  scale_features=TRUE,
  orthogonal_gpd=TRUE,
  n_epochs=1000,
  batch_size=256,
  lr_decay=0.4,
  patience_decay=5,
  min_lr=1e-5,
  patience_stop=20
)
hid_f_str <- "tanh"
lr_str <- "lrd3or4"

## =============================== END PARAMETERS ===============================

