
generate_series_model <- function(n, df, AR, MA, muX, mu0=0, alphas, betas, sX, S0=1,
                                  ARX=c(0), X_distr=c("foldnormal","normal","uniform"),
                                  Y_distr=c("t","pareto","foldnormal","normal"), presample_n=2*time_dep,
                                  seasonal_hetero=0){
  
  X_distr <- match.arg(X_distr)
  Y_distr <- match.arg(Y_distr)
  sX <- vec2mat(sX)
  muX <- vec2mat(muX)
  ARX <- matrix(ARX, nrow=1)
  
  p <- max(nrow(sX), nrow(muX))
  
  if(sum(abs(c(alphas, betas)))>=1){stop("Please provide alpha and beta coefficients absolutily summing below unity in 'generate_series_model'.")}
  if(any(abs(polyroot(c(1, -AR)))<=1)){stop("root(s) of AR process polynomial lie in the unit circle in 'generate_series_model'.")}
  if(any(abs(polyroot(c(1, MA)))<=1)){stop("root(s) of MA process polynomial lie in the unit circle in 'generate_series_model'.")}
  
  time_dep <- max(sapply(list(t(alphas), t(betas), sX, t(AR), t(MA), muX, ARX), ncol))
  alphas <- pad_params(alphas, time_dep)
  betas <- pad_params(betas, time_dep)
  sX <- pad_params(sX, time_dep)
  sX <- pad_params(sX, p, "row")
  AR <- pad_params(AR, time_dep)
  MA <- pad_params(MA, time_dep)
  muX <- pad_params(muX, time_dep)
  muX <- pad_params(muX, p, "row")
  ARX <- pad_params(ARX, time_dep)
  X <- matrix(as.double(NA), nrow=(n+presample_n), ncol=p)
  X_eps <- matrix(if(X_distr=="uniform") runif(n=p*(n+presample_n),min=-1,max=1)
                  else if(X_distr=="foldnormal") VGAM::rfoldnorm(n=p*(n+presample_n), mean=0, sd=1, a1=1, a2=1)
                  else rnorm(n=p*(n+presample_n), mean=0, sd=1),
                  nrow=(n+presample_n), ncol=p)
  Z <-  if(Y_distr=="t") rt(n=(n+presample_n), df=df)
        else if(Y_distr=="pareto") (EnvStats::rpareto(n=(n+presample_n), 1, df) - df/(df-1))
        else if(Y_distr=="foldnormal") VGAM::rfoldnorm(n=(n+presample_n), mean=0, sd=1, a1=1, a2=1)
        else rnorm(n=(n+presample_n), mean=0, sd=1)
  S <- rep(as.double(NA), n)
  eps <- rep(as.double(NA), n)
  trend <- rep(as.double(NA), n)
  Y <- rep(as.double(NA), n)
  X[1,] <- X_eps[1, , drop=F]
  trend[1] <- mu0
  S[1] <- sqrt(S0)
  eps[1] <- S[1]*Z[1]
  Y[1] <- trend[1] + eps[1]
  
  for(i in 2:(n+presample_n)){
    tds <- 1:min((i-1),time_dep)
    X[i,] <- ARX[,tds, drop=F] %*% X[i-tds, , drop=F] + X_eps[i, , drop=F]
    trend[i] <- mu0 + sum(AR[tds]*Y[i-tds]) + sum(MA[tds]*eps[i-tds]) + sum(muX[,tds, drop=F]%*%X[i-tds, , drop=F])
    S[i] <- sqrt(S0 + sum(alphas[tds]*eps[i-tds]^2) + sum(betas[tds]*S[i-tds]^2) + sum(sX[,tds, drop=F] %*% X[i-tds, , drop=F]^2))
    if(seasonal_hetero>0){S[i] <- S[i] + sin(2*pi*(i-1)/seasonal_hetero)+1}
    eps[i] <- S[i]*Z[i]
    Y[i] <- trend[i] + eps[i]
  }
  if(seasonal_hetero>0){X <- cbind(X,sin(2*pi*(1:(n+presample_n))/seasonal_hetero)+1)}
  out_range <- (presample_n+1):(n+presample_n)
  out <- list(Y=Y[out_range], X=X[out_range, , drop=F],
              trend=trend[out_range], S=S[out_range], Z=Z[out_range],
              df=df)
}

series_theoretical_quantiles <- function(quantiles, gen, Y_distr=c("t","pareto","foldnormal","normal")) {
  Y_distr <- match.arg(Y_distr)
  l <- length(quantiles)
  n <- length(gen$S)
  
  sigma <- gen$S %>% rep_vector2matrix(nrep = l, dim = "col")
  
  trend <- gen$trend %>% rep_vector2matrix(nrep = l, dim = "col")
  
  if(length(gen$df)==1 & n>1){
    df <- rep(gen$df, n)
  }else{
    df <- gen$df
  }
  trend + sigma * if(Y_distr=="t") quantiles_student_t(quantiles, df)
                  else if(Y_distr=="pareto") (quantiles_pareto(quantiles, df) - df/(df-1))
                  else if(Y_distr=="foldnormal") quantiles_foldnorm(quantiles, n=n, mean=0, sd=1, a1=1, a2=1)
                  else quantiles_gaussian(quantiles, n=n)
}

series_theoretical_cdf <- function(val, gen, Y_distr=c("t","pareto","foldnormal","normal"), lower.tail=TRUE){
  Y_distr <- match.arg(Y_distr)
  l <- length(val)
  n <- length(gen$S)
  
  if(l!=1 & l!=n){stop("Input length issue in 'series_theoretical_cdf'.")}
  
  sigma <- gen$S
  trend <- gen$trend
  
  if(length(gen$df)==1 & n>1){
    df <- rep(gen$df, n)
  }else{
    df <- gen$df
  }
  cdf <- if(Y_distr=="t") cdf_student_t((val-trend)/sigma, df)
  else if(Y_distr=="pareto") cdf_pareto((val-trend)/sigma + df/(df-1), df)
  else if(Y_distr=="foldnormal") cdf_foldnorm((val-trend)/sigma, mean=0, sd=1, a1=1, a2=1)
  else cdf_gaussian((val-trend)/sigma)
  if(lower.tail){
    return(cdf)
  }else{
    return(1-cdf)
  }
}

pad_params <- function(params, size, axis=c("col","row")){
  axis <- match.arg(axis)
  if(is.null(dim(params))){
    return( c(params, rep(0, size-length(params))) )
  }else{
    if(axis=="col"){
      return( cbind(params, matrix(0, nrow=nrow(params), ncol=size-ncol(params))) )
    }else{
      return( rbind(params, matrix(0, nrow=size-nrow(params), ncol=ncol(params))) )
    }
  }
}

