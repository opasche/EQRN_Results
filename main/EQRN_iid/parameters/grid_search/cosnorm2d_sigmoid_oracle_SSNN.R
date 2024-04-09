
## =============================== PARAMETERS ===============================
save_path <- "Results/iid/EQRN_grid_search/cosnorm2d_sigmoid_n5e3_p10_q80_oracle_ortho_SSNN/"
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
model = "cosnorm2d_sigmoid"# "binorm_sigmoid", "cosnorm2d_sigmoid", "cosnorm_sigmoid"
distr = "student_t"
df <- 4
test_data = "halton"

# PARAM: General
intermediate_method <- "oracle"
interm_lvl = 0.8
prob_lvls_predict = c(0.995,0.999,0.9995)

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
  net_structure=purrr::cross(list(scale=list(c(64,64,64),c(64,64,64,64),c(64,64,64,64,64,64)),
                                  shape=list(c(1),c(5,3),c(5,3,3))), .filter = NULL),
  hidden_fct=list("SSNN"),
  p_drop=list(0,0.01,0.02,0.05),
  intermediate_q_feature=list(TRUE),
  learning_rate=list(1e-4,1e-5),
  L2_pen=list(0),
  shape_penalty=list(0),
  scale_features=TRUE,
  orthogonal_gpd=TRUE,
  n_epochs=1000,
  batch_size=256,
  lr_decay=0.4,
  patience_decay=5,
  min_lr=1e-6,
  patience_stop=20
)
hid_f_str <- "SSNN"
lr_str <- "lrd4or5"

## =============================== END PARAMETERS ===============================

