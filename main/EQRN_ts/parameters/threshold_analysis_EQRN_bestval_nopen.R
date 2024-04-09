## =============================== PARAMETERS ===============================
save_path <- "Results/ts/EQRN_threshold_analysis/bestval_nopen/"
save_plots <- TRUE
force_refit <- FALSE

seedR <- 0
seedGRF <- 1
seedT <- seedR

# PARAM: Data
n <- 5e3
ntest <- 5e3
n_valid <- 2e3
X_distr="foldnormal"
Y_distr="foldnormal"
df <- 4
alphas=c(.2,.1,.1,.1,.1)
betas=c(0)
sX=c(.3,.2,.1,.1,.1)
S0=1
AR=c(0)
MA=c(0)
muX=c(0)
mu0=0
ARX=c(0.4)
seasonal_hetero=0

# PARAM: General
intermediate_method <- "qrn"#"qrn""oracle"
interm_method_competitors <- intermediate_method#"qrn""grf""oracle"
interm_lvls = c(0.7,0.75,0.8,0.85,0.9,0.95)
prob_lvls_predict = c(0.995,0.999,0.9995)#c(interm_lvl,0.995,0.999,0.9995)

# Params: QRNN
interm_path <- "data/simulations/ts_intermediate_quantile_threshold_analysis/"
par_qrn <- list(
  nb_fits=3,
  rnn_type="lstm",
  num_layers=1,
  hidden_size=128,
  p_drop=0,
  L2_pen=0,
  seq_len=10,
  learning_rate=1e-3,
  n_epochs=1e3,
  batch_size=256,
  lr_decay=0.4,
  patience_decay=5,
  min_lr=1e-5,
  patience_stop=20,
  scale_features=TRUE,
  tol=1e-4
)

# ============== Params: EQRN ==============
path_eqrn <- "Results/ts/EQRN_threshold_analysis/network_saves/"#save_path
par_eqrn <- list(
  nb_fits = 3,
  shape_fixed=TRUE,
  rnn_type="lstm",
  num_layers=2,
  hidden_size=128,
  p_drop=0,
  intermediate_q_feature=TRUE,
  L2_pen=0,
  shape_penalty=0,
  seq_len=par_qrn$seq_len,
  learning_rate=1e-3,
  scale_features=TRUE,
  orthogonal_gpd=TRUE,
  n_epochs=1000,
  batch_size=256,
  lr_decay=0.4,
  patience_decay=5,
  min_lr=1e-5,
  patience_stop=20,
  tol=1e-5
)

hid_f_str <- "tanh"
lr_str <- "lrd4"

# Params: GRF
grf_pars <- list(
  timedep = par_qrn$seq_len,
  num.trees = 5e3,
  quantiles_fit = c(0.1, 0.5, 0.9),
  sample.fraction = 0.5,
  #mtry = min(ceiling(sqrt(p) + 20), p),
  min.node.size = 5,
  honesty = TRUE,
  honesty.fraction = 0.5,
  honesty.prune.leaves = TRUE,
  alpha = 0.05,
  imbalance.penalty = 0
)

# Params: gbex
gbex_params <- list(
  B=556,
  lambda=NULL,
  lambda_ratio=6,
  lambda_scale=0.01,
  depth=c(1, 0),
  sf=0.75,
  intermediate_q_feature=TRUE,
  scale_features=FALSE)

# Params: egam
egam_params <- list(
  model_shape=FALSE,
  intermediate_q_feature=gbex_params$intermediate_q_feature,
  scale_features=gbex_params$scale_features)

## =============================== END PARAMETERS ===============================
