generate_joint_distribution <- function(n, p,
                                        model = c("iid",
                                                  "step", 
                                                  "step2d",
                                                  "linear",
                                                  "quadratic",
                                                  "mixture",
                                                  "step_step", 
                                                  "step2d_step",
                                                  "linear_step",
                                                  "quadratic_step",
                                                  "binorm_sigmoid",
                                                  "cosnorm",
                                                  "cosnorm_step",
                                                  "cosnorm_widestep",
                                                  "cosnorm_sigmoid",
                                                  "cosnorm2d_sigmoid"),
                                        distr = c("gaussian", "student_t"), df){
  ## integer (x2) character (x2) numeric -> list
  ## generate n iid observations of (X, Y), where X is p-dimensional predictor
  ## and Y is the response following the given model with the
  ## given distribution.
  ## Returns a list with:
  ## - X, nxp matrix, p-dimensional predictor
  ## - Y, vector with n elements, response variable
  
  generate_joint_distribution_internal(n, p, model, distr, df)
}

generate_conditional_distribution <- function(model = c("iid",
                                                        "step", 
                                                        "step2d",
                                                        "linear",
                                                        "quadratic",
                                                        "mixture",
                                                        "step_step", 
                                                        "step2d_step",
                                                        "linear_step",
                                                        "quadratic_step",
                                                        "binorm_sigmoid",
                                                        "cosnorm",
                                                        "cosnorm_step",
                                                        "cosnorm_widestep",
                                                        "cosnorm_sigmoid",
                                                        "cosnorm2d_sigmoid"),
                                              distr = c("gaussian", "student_t"),
                                              df, X){
  ## integer (x2) character (x2) numeric numeric_matrix -> list
  ## generate n iid observations of (Y| X = x), where X is p-dimensional predictor
  ## and Y is the response following the given model with the
  ## given distribution.
  ## Returns a list with:
  ## - X, nxp matrix, p-dimensional predictor
  ## - Y, vector with n elements, response variable
  n <- nrow(X)
  p <- ncol(X)
  generate_joint_distribution_internal(n, p, model, distr, df, X)
  
}

generate_theoretical_quantiles <- function(quantiles, X,
                                           model = c("iid",
                                                     "step", 
                                                     "step2d",
                                                     "linear",
                                                     "quadratic",
                                                     "mixture",
                                                     "step_step", 
                                                     "step2d_step",
                                                     "linear_step",
                                                     "quadratic_step",
                                                     "binorm_sigmoid",
                                                     "cosnorm",
                                                     "cosnorm_step",
                                                     "cosnorm_widestep",
                                                     "cosnorm_sigmoid",
                                                     "cosnorm2d_sigmoid"),
                                           distr = c("gaussian", "student_t"),
                                           df){
  ## numeric_vector numeric_matrix character (x2) numeric -> numeric_matrix
  ## produce theoretical quantiles for the given model and distribution
  ## for the different observations (rows of X)
  
  model <- match.arg(model)
  distr <- match.arg(distr)
  n <- nrow(X)
  p <- ncol(X)
  
  if (!missing(df) & distr == "gaussian"){
    warning("df is not considered when distr = 'gaussian'")
  }
  
  switch(model,
         "iid" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_iid, df_constant)
         },
         "step" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_step, df_constant)
         },
         "step2d" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_step2d, df_constant)
         },
         "linear" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_linear, df_constant)
         },
         "quadratic" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_quadratic, df_constant)
         },
         "mixture" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_mixture, df_constant)
         },
         "step_step" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_step, df_step)
         },
         "step2d_step" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_step2d, df_step)
         },
         "linear_step" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_linear, df_step)
         },
         "quadratic_step" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_quadratic, df_step)
         },
         "binorm_sigmoid" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_binorm, df_sigmoid)
         },
         "cosnorm" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_cosnorm, df_constant)
         },
         "cosnorm_step" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_cosnorm, df_step)
         },
         "cosnorm_widestep" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_cosnorm, df_widestep)
         },
         "cosnorm_sigmoid" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_cosnorm, df_sigmoid)
         },
         "cosnorm2d_sigmoid" = {
           q_vec <- generate_quantile_model(quantiles, X, distr, df,
                                            sigma_cosnorm2d, df_sigmoid)
         }
  )
  
  return(q_vec)
}

generate_Y_model <- function(X, distr = c("gaussian", "student_t"),
                             df, sigma_fun, df_fun){
  ## numeric_matrix character numeric function function -> numeric_vector
  ## generate random response Y for the given model
  
  distr <- match.arg(distr)
  n <- nrow(X)
  p <- ncol(X)
  
  switch(distr,
         "gaussian" = {
           Y_tilde <- rnorm(n)
         },
         "student_t" = {
           Y_tilde <- rt(n, df = df_fun(X, df))
         })
  
  sigma_x <- sigma_fun(X)
  
  sigma_x * Y_tilde
}

generate_quantile_model <- function(quantiles, X, 
                                    distr = c("gaussian", "student_t"),
                                    df, sigma_fun, df_fun){
  ## numeric_vector numeric_matrix character numeric -> numeric_matrix
  ## produce theoretical quantiles for the given model
  
  distr <- match.arg(distr)
  n <- nrow(X)
  p <- ncol(X)
  l <- length(quantiles)
  
  switch(distr,
         "gaussian" = {
           q_tilde <- quantiles_gaussian(quantiles, n)
         },
         "student_t" = {
           q_tilde <- quantiles_student_t(quantiles, df_fun(X, df))
           
         })
  
  sigma_x <- sigma_fun(X) %>% rep_vector2matrix(nrep = l, dim = "col")
  
  sigma_x * q_tilde
}

quantiles_student_t <- function(quantiles, df){
  ## numeric_vector numeric_vector -> numeric_matrix
  ## produce matrix with student_t quantiles where nrow = length(df) and
  ## ncol = length(quantiles)
  
  purrr::map(quantiles, function(q){qt(q, df = df)}) %>% 
    list2matrix(dim = "col")
}

cdf_student_t <- function(q, df){
  ## numeric_vector numeric_vector -> numeric_matrix
  ## produce matrix with student_t quantiles where nrow = length(df) and
  ## ncol = length(quantiles)
  
  pt(q, df = df)
}

quantiles_pareto <- function(quantiles, df){
  ## numeric_vector numeric_vector -> numeric_matrix
  ## produce matrix with pareto quantiles where nrow = length(df) and
  ## ncol = length(quantiles)
  
  purrr::map(quantiles, function(q){EnvStats::qpareto(q, 1, df)}) %>% 
    list2matrix(dim = "col")
}

cdf_pareto <- function(q, df){
  ## numeric_vector numeric_vector -> numeric_matrix
  ## produce matrix with pareto quantiles where nrow = length(df) and
  ## ncol = length(quantiles)
  
  EnvStats::ppareto(q, 1, df)
}

quantiles_gaussian <- function(quantiles, n){
  ## numeric_vector integer -> numeric_matrix
  ## produce matrix with gaussian quantiles where nrow = n and 
  ## ncol = length(quantiles)
  
  rep_vector2matrix(qnorm(quantiles), nrep = n,
                    dim = "row")
}

cdf_gaussian <- function(q){
  ## numeric_vector integer -> numeric_matrix
  ## produce matrix with gaussian quantiles where nrow = n and 
  ## ncol = length(quantiles)
  
  pnorm(q)
}

quantiles_foldnorm <- function(quantiles, n, mean = 0, sd = 1, a1 = 1, a2 = 1){
  ## numeric_vector integer -> numeric_matrix
  ## produce matrix with folded-normal quantiles where nrow = n and 
  ## ncol = length(quantiles)
  
  rep_vector2matrix(VGAM::qfoldnorm(quantiles, mean=mean, sd=sd, a1=a1, a2=a2), nrep = n,
                    dim = "row")
}

cdf_foldnorm <- function(q, mean = 0, sd = 1, a1 = 1, a2 = 1){
  ## numeric_vector integer -> numeric_matrix
  ## produce matrix with folded-normal quantiles where nrow = n and 
  ## ncol = length(quantiles)
  
  VGAM::pfoldnorm(q, mean=mean, sd=sd, a1=a1, a2=a2)
}

# sigma functions

sigma_iid <- function(X){
  ## numeric_matrix -> numeric_vector
  ## produce scale function: iid data not depending on X's
  
  rep(1, times = nrow(X))
}

sigma_step <- function(X){
  ## numeric_matrix -> numeric_vector
  ## produce scale function: scale(X1) = step function
  
  sigma_x <- 1 + 1 * (X[, 1] > 0)
  
  return(sigma_x)
  
}

sigma_step2d <- function(X){
  ## numeric_matrix -> numeric_vector
  ## produce scale function: scale(X1, X2) = step function
  
  sigma_x <- numeric(nrow(X))
  sigma_x[X[, 1] < 0] <- 1
  sigma_x[X[, 1] >= 0 & X[, 2] < 0.5] <- 3
  sigma_x[X[, 1] >= 0 & X[, 2] >= 0.5] <- 2
  
  return(sigma_x)
}

sigma_linear <- function(X){
  ## numeric_matrix -> numeric_vector
  ## produce scale function: linear model scale(X1) = 2 + 1 * X1
  
  X1 <- X[, 1]
  2 + 1 * X1
}

sigma_quadratic <- function(X){
  ## numeric_matrix -> numeric_vector
  ## produce scale function: quadratic model scale(X1) = 3 - 2 * X1 ^2
  
  X1 <- X[, 1]
  3 - 2 * X1 ^2
}

sigma_mixture <- function(X){
  ## numeric_matrix -> numeric_vector
  ## produce scale function: scale(X1, X2) = mixture of 2 Gaussians
  
  mu_1 <- c(-.5, .5)
  sigma_1 <- rbind(c(1/8, 0), c(0, 1/8))
  w_1 <- 6
  mu_2 <- c(.5, -.5)
  sigma_2 <- rbind(c(1/8, 0), c(0, 1/8))
  w_2 <- 6
  
  sigma_x <- 1 +
    w_1 * mvtnorm::dmvnorm(X[, c(1, 2)], mean = mu_1, sigma = sigma_1) +
    w_2 * mvtnorm::dmvnorm(X[, c(1, 2)], mean = mu_2, sigma = sigma_2)
  
  return(sigma_x)
}

sigma_binorm <- function(X){
  ## numeric_matrix -> numeric_vector
  ## produce scale function: scale(X1, X2) = 1 + 6*binorm(X1,X2;0.9)
  
  mu <- c(0, 0)
  sigma <- rbind(c(1, 0.9), c(0.9, 1))
  sigma_x <- 1 + 6 * mvtnorm::dmvnorm(X[, c(1, 2)], mean = mu, sigma = sigma)
  
  return(sigma_x)
}

sigma_cosnorm <- function(X){
  ## numeric_matrix -> numeric_vector
  ## produce scale function: scale(X) = 4 + 3*cos(6*||X||_2^2 + 3.5)
  
  sigma_x <- 4 + 3*cos(6*sqrt(rowSums(X^2))+3.5)
  
  return(sigma_x)
}

sigma_cosnorm2d <- function(X){
  ## numeric_matrix -> numeric_vector
  ## produce scale function: scale(X) = 4 + 3*cos(7*||(X1,X2)||_2^2 + 3)
  
  sigma_x <- 4 + 3*cos(7*sqrt(rowSums(X[, c(1, 2)]^2))+3)
  
  return(sigma_x)
}

# shape functions (expressed as degrees of freedom (df), where df = 1 / shape)

df_constant <- function(X, df){
  ## numeric_matrix numeric -> numeric_vector
  ## produce degree of freedom function: constant equal to df
  
  rep(df, times = nrow(X))
}

df_step <- function(X, df){
  ## numeric_matrix numeric -> numeric_vector
  ## produce degree of freedom function: shape(X2) = step function
  
  df_x <- df + 2 * df * (X[, 2] > 0)
  
  return(df_x)
}

df_sigmoid <- function(X, df){
  ## numeric_matrix numeric -> numeric_vector
  ## produce degree of freedom function: shape(X1) = 7*(1 + exp(4*X1 + 1.2))^(-1) + 3
  
  df_x <- 7/(1 + exp(4*X[,1] + 1.2)) + 3
  
  return(df_x)
}

df_widestep <- function(X, df){
  ## numeric_matrix numeric -> numeric_vector
  ## produce degree of freedom function: shape(X2) = df + (2*X2+0.5) * (X2>-0.25 & X2<0.25) + 1 * (X2>0.25)
  
  df_x <- df + df*((2*X[,1]+0.5) * (X[, 1]>-0.25 & X[, 1]<0.25) + 1 * (X[, 1]>=0.25))
  
  return(df_x)
}


# other functions

generate_test_data <- function(ntest, p, method=c("halton", "grid", "uniform"),
                               warn=FALSE){
  ## integer (x2) character -> tibble
  ## generate the predictor test data using "halton" or "grid" method
  
  method <- match.arg(method)
  if (method == "halton"){
    X_test <- randtoolbox::halton(ntest, p) * 2 - 1
    
  } else if (method == "grid"){
    ntest_mod <- floor(sqrt(ntest))
    if(warn){warning(paste0("Modified ntest to: ", ntest_mod ** 2))}
    X_test <- expand_grid(X1 = seq(-1, 1, length.out = ntest_mod),
                          X2 = seq(-1, 1, length.out = ntest_mod)) %>%
      as.matrix()
    
    if (p > 2){
      X_test <- cbind(X_test,
                      matrix(runif(ntest_mod ** 2, min = -1, max = 1),
                             ncol = p - 2,
                             nrow = ntest_mod ** 2))
    }
  } else if (method == "uniform"){
    X_test <- matrix(runif(ntest * p, min = -1, max = 1), ntest, p)
  }
  
  colnames(X_test) <- paste0("X", 1:p)
  return(X_test)
}

generate_iid_data <- function(n, distr = c("gaussian", "log_normal", "student_t"),
                              df = 2){
  ## integer character numeric -> numeric_vector
  ## generate n random observations from distribution distr (with dof = df)
  
  distr <- match.arg(distr)
  
  if (distr == "gaussian"){
    
    Y <- rnorm(n)
    
  } else if (distr == "log_normal"){
    
    Y <- rlnorm(n)
    
  } else if (distr == "student_t"){
    
    Y <- rt(n, df)
    
  }
  
  return(Y)
}

generate_iid_quantile <- function(quantiles_predict,
                                  distr = c("gaussian", "log_normal", "student_t"),
                                  df = 2){
  ## numeric_vector character numeric -> numeric_vector
  ## generate theoretical quantiles of distribution distr (with dof = df)
  
  distr <- match.arg(distr)
  
  if (distr == "gaussian"){
    
    q <- qnorm(p = quantiles_predict)
    
  } else if (distr == "log_normal"){
    
    q <- qlnorm(p = quantiles_predict)
    
  } else if (distr == "student_t"){
    
    q <- qt(p = quantiles_predict, df = df)
    
  }
  
  return(q)
}  

generate_joint_distribution_internal <- function(n, p,
                                                 model = c("iid",
                                                           "step", 
                                                           "step2d",
                                                           "linear",
                                                           "quadratic",
                                                           "mixture",
                                                           "step_step", 
                                                           "step2d_step",
                                                           "linear_step",
                                                           "quadratic_step",
                                                           "binorm_sigmoid",
                                                           "cosnorm",
                                                           "cosnorm_step",
                                                           "cosnorm_widestep",
                                                           "cosnorm_sigmoid",
                                                           "cosnorm2d_sigmoid"),
                                                 distr = c("gaussian", "student_t"),
                                                 df, X = NULL){
  ## integer (x2) character (x2) numeric numeric_matrix|NULL-> list
  ## generate n iid observations of (X, Y), where X is p-dimensional predictor
  ## and Y is the response following the given model with the
  ## given distribution.
  ## Returns a list with:
  ## - X, nxp matrix, p-dimensional predictor
  ## - Y, vector with n elements, response variable
  
  model <- match.arg(model)
  distr <- match.arg(distr)
  
  if (is.null(X)){
    X <- matrix(runif(n * p, min = -1, max = 1), n, p)
  } else {
    check_X_matrix(X, n, p)
  }
  colnames(X) <- paste0("X", 1:p)
  
  if (!missing(df) & distr == "gaussian"){
    warning("df is not considered when distr = 'gaussian'")
  }
  
  switch(model,
         "iid" = {
           Y <- generate_Y_model(X, distr, df, sigma_iid, df_constant)
         },
         "step" = {
           Y <- generate_Y_model(X, distr, df, sigma_step, df_constant)
         },
         "step2d" = {
           Y <- generate_Y_model(X, distr, df, sigma_step2d, df_constant)
         },
         "linear" = {
           Y <- generate_Y_model(X, distr, df, sigma_linear, df_constant)
         },
         "quadratic" = {
           Y <- generate_Y_model(X, distr, df, sigma_quadratic, df_constant)
         },
         "mixture" = {
           Y <- generate_Y_model(X, distr, df, sigma_mixture, df_constant)
         },
         "step_step" = {
           Y <- generate_Y_model(X, distr, df, sigma_step, df_step)
         },
         "step2d_step" = {
           Y <- generate_Y_model(X, distr, df, sigma_step2d, df_step)
         },
         "linear_step" = {
           Y <- generate_Y_model(X, distr, df, sigma_linear, df_step)
         },
         "quadratic_step" = {
           Y <- generate_Y_model(X, distr, df, sigma_quadratic, df_step)
         },
         "binorm_sigmoid" = {
           Y <- generate_Y_model(X, distr, df, sigma_binorm, df_sigmoid)
         },
         "cosnorm" = {
           Y <- generate_Y_model(X, distr, df, sigma_cosnorm, df_constant)
         },
         "cosnorm_step" = {
           Y <- generate_Y_model(X, distr, df, sigma_cosnorm, df_step)
         },
         "cosnorm_widestep" = {
           Y <- generate_Y_model(X, distr, df, sigma_cosnorm, df_widestep)
         },
         "cosnorm_sigmoid" = {
           Y <- generate_Y_model(X, distr, df, sigma_cosnorm, df_sigmoid)
         },
         "cosnorm2d_sigmoid" = {
           Y <- generate_Y_model(X, distr, df, sigma_cosnorm2d, df_sigmoid)
         }
         
  )
  
  return(list(X = X, Y = Y))
}

effective_dim_model <- function(model = c("iid",
                                          "step", 
                                          "step2d",
                                          "linear",
                                          "quadratic",
                                          "mixture",
                                          "step_step", 
                                          "step2d_step",
                                          "linear_step",
                                          "quadratic_step",
                                          "binorm_sigmoid",
                                          "cosnorm",
                                          "cosnorm_step",
                                          "cosnorm_widestep",
                                          "cosnorm_sigmoid",
                                          "cosnorm2d_sigmoid"),
                                p=10){
  ## character -> int
  ## Returns the effective dimension of the model, i.e. the number of X features having an effect on Y.
  
  model <- match.arg(model)
  
  switch(model,
         "iid" = {
           return(0)
         },
         "step" = {
           return(1)
         },
         "step2d" = {
           return(2)
         },
         "linear" = {
           return(1)
         },
         "quadratic" = {
           return(1)
         },
         "mixture" = {
           return(2)
         },
         "step_step" = {
           return(2)
         },
         "step2d_step" = {
           return(2)
         },
         "linear_step" = {
           return(2)
         },
         "quadratic_step" = {
           return(2)
         },
         "binorm_sigmoid" = {
           return(2)
         },
         "cosnorm" = {
           return(p)
         },
         "cosnorm_step" = {
           return(p)
         },
         "cosnorm_widestep" = {
           return(p)
         },
         "cosnorm_sigmoid" = {
           return(p)
         },
         "cosnorm2d_sigmoid" = {
           return(2)
         })
}
