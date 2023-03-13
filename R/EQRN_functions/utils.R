
#' Check directory existence
#' 
#' @description Checks if the desired directory exists. If not, the desired directory is created.
#'
#' @param dir_name Path to the desired directory, as a string.
#' @param recursive Should elements of the path other than the last be created? 
#' If \code{TRUE}, behaves like the Unix command \code{mkdir -p}.
#' @param no_warning Whether to cancel the warning issued if a directory is created (bool).
#' 
#' @export
#'
#' @examples
check_directory <- function(dir_name, recursive=TRUE, no_warning=FALSE){
  if (!dir.exists(dir_name)){
    dir.create(dir_name, recursive=recursive)
    if(!no_warning){
      warning(paste0("The following given directory did not exist and was created by 'check_directory': ", dir_name))
    }
  }
}

#' Safe RDS save
#' 
#' @description Safe version of \code{\link{saveRDS}}. 
#' If the given save path (i.e. \code{dirname(file_path)}) does not exist, it is created instead of raising an error.
#'
#' @param object R variable or object to save on disk.
#' @param file_path Path and name of the save file, as a string.
#' @param recursive Should elements of the path other than the last be created? 
#' If \code{TRUE}, behaves like the Unix command \code{mkdir -p}.
#' @param no_warning Whether to cancel the warning issued if a directory is created (bool).
#'
#' @export
#'
#' @examples
safe_save_rds <- function(object, file_path, recursive=TRUE, no_warning=FALSE){
  dir_name <- dirname(file_path)
  check_directory(dir_name, recursive=recursive, no_warning=no_warning)
  
  saveRDS(object, file = file_path)
  
}

#' Last element of a vector
#'
#' @param x Vector.
#' 
#' @description Returns the last element of the given vector in the most efficient way.
#'
#' @return The last element in the vector \code{x}.
#' 
#' @details The last element is obtained using code{x[length(x)]}, which is done in \code{O(1)} and faster than, for example, any of
#' \code{Rcpp::mylast(x)}, \code{tail(x, n=1)}, \code{dplyr::last(x)}, \code{x[end(x)[1]]]}, and \code{rev(x)[1]}.
#' @export
#'
#' @examples
last_elem <- function(x){
  x[length(x)]
}

#' Mathematical number rounding
#' 
#' @description This function rounds numbers in the mathematical sense, 
#' as opposed to the base \code{R} function \code{\link{round}} that rounds 'to the even digit'.
#'
#' @param x Vector of numerical values to round.
#' @param decimals Integer indicating the number of decimal places to be used.
#'
#' @return A vector containing the entries of \code{x}, rounded to \code{decimals} decimals.
#' @export
#'
#' @examples
roundm = function(x, decimals=0){
  posneg <- sign(x)
  z <- abs(x)*10^decimals
  z <- z + 0.5 + sqrt(.Machine$double.eps)
  z <- trunc(z)
  z <- z/10^decimals
  z*posneg
}

#' Convert a vector to a matrix
#'
#' @param v Vector.
#' @param axis One of \code{"col"} (default) or \code{"row"}.
#'
#' @return The vector \code{v} as a matrix. 
#' If \code{axis=="col"} (default) the column vector \code{v} is returned as a \code{length(v)} times \code{1} matrix. 
#' If \code{axis=="row"}, the vector \code{v} is returned as a transposed \code{1} times \code{length(v)} matrix.
#' @export
#'
#' @examples
vec2mat <- function(v, axis=c("col","row")){
  axis <- match.arg(axis)
  if (is.null(dim(v))) {
    v <- if(axis=="col"){matrix(v, nrow=1)}else{matrix(v, ncol=1)}
  }
  return(v)
}

#' Tibble replicatior
#'
#' @param tbl A tibble::tibble.
#' @param m An integer.
#'
#' @return The tibble is replicated \code{m} times and colums names appended with \code{rep_id = 1:m}.
#'
#' @examples
rep_tibble <- function(tbl, m){
  tbl <-  tbl %>% tibble::rownames_to_column()
  
  tidyr::expand_grid(rep_id = 1:m, rowname = tbl$rowname) %>% 
    left_join(tbl, by = "rowname") %>% 
    select(-rowname)
}

#' Replicated vector to matrix
#'
#' @param vec Vector.
#' @param nrep Number of repetitions.
#' @param dim One of \code{"row"} (default) or \code{"col"}.
#'
#' @return Matrix of replicated vector.
#'
#' @examples
rep_vector2matrix <- function(vec, nrep, dim = c("row", "col")){
  ## vector integer character -> matrix
  ## stack nrep of vec in (row|col) of a matrix
  
  dim <- match.arg(dim)
  l <- length(vec)
  
  if (l == 0){
    stop("vec must contain at least one element.")
  }
  
  if (dim == "col"){
    matrix(vec, nrow = l, ncol = nrep)
  } else {
    matrix(vec, nrow = nrep, ncol = l, byrow = TRUE)
  }
}

#' Convert a list to a matrix
#'
#' @param lst A list.
#' @param dim One of \code{"row"} (default) or \code{"col"}.
#'
#' @return The list converted to a matrix, by stacking the elements of \code{lst} in the rows or columns of a matrix.
#'
#' @examples
list2matrix <- function(lst, dim = c("row", "col")){
  dim <- match.arg(dim)
  l <- length(lst)
  
  if (l == 0){
    stop("lst must contain at least one element.")
  }
  
  if (dim == "col"){
    matrix(unlist(lst), ncol = l)
  } else {
    matrix(unlist(lst), nrow = l, byrow = TRUE)
  }
}

#' Convert a matrix to a list
#'
#' @param mat A matrix.
#'
#' @return A list with elements corresponding to rows of \code{mat}.
#'
#' @examples
matrix2list <- function(mat){
  split(mat, rep(1:nrow(mat), times = ncol(mat)))
}

#' Check the simulation X matrix
#'
#' @param X Covariate matrix.
#' @param n Number of observations.
#' @param p Number of covariates.
#' 
#' @return Returns TRUE if X is a matrix with dimension n * p. Otherwise an error is raised.
#'
#' @examples
check_X_matrix <- function(X, n, p){
  cond_1 <- is.matrix(X)
  
  if (cond_1){
    cond_2 <- all.equal(dim(X), c(n, p))
  } else {
    cond_2 <- FALSE
  }
  
  if (cond_1 & cond_2){
    return(TRUE)
  } else {
    stop(paste0("X must be a matrix with ", deparse(substitute(n)),
                " rows and ", deparse(substitute(p)), " columns."))
  }
}

#' Create cross-validation folds
#' 
#' @description Utility function to create folds of data, used in cross-validation proceidures. 
#' The implementation is from the \code{gbex} \code{R} package
#'
#' @param y Numerical vector of observations
#' @param num_folds Number of folds to create.
#' @param stratified Logical value. If \code{TRUE}, the folds are stratified along \code{rank(y)}.
#'
#' @return Vector of indices of the assigned folds for each observation.
#' @export
#'
#' @examples
make_folds <- function(y, num_folds, stratified=FALSE){
  n = length(y)
  if(stratified) {
    folds_matrix <- sapply(1:ceiling(n/num_folds), function(i) {
      sample(1:num_folds)
    })
    folds_vector <- folds_matrix[1:n]
    folds <- folds_vector[rank(-y)]
  } else {
    index_shuffled = sample(1:n)
    folds = cut(seq(1, length(index_shuffled)), breaks = num_folds, 
                labels = F)[order(index_shuffled)]
  }
  return(folds)
}

#' Covariate lagged replication for temporal dependence
#'
#' @param X Covariate matrix.
#' @param max_lag Integer giving the maximum lag (i.e. the number of temporal dependence steps).
#' @param drop_present Whether to drop the "present" features (bool).
#'
#' @return Matrix with the original columns replicated, and shifted by \code{1:max_lag} if \code{drop_present==TRUE} (default) 
#' or by \code{0:max_lag} if \code{drop_present==FALSE}.
#' @export
#'
#' @examples
lagged_features <- function(X, max_lag, drop_present=TRUE){
  n <- nrow(X)
  p <- ncol(X)
  Xl <- matrix(as.double(NA), nrow=n-max_lag, ncol=p*(max_lag+1))
  for(i in 0:max_lag){
    Xl[, (p*i+(1:p))] <- X[(max_lag+1-i):(n-i), , drop=F]
  }
  if(drop_present){
    Xl <- Xl[, (p+1):(p*(max_lag+1)), drop=F]
  }
  return(Xl)
}


# ==== Parallel helpers ====

#' Get doFuture operator
#'
#' @param strategy One of \code{"sequential"} (default), \code{"multisession"}, \code{"multicore"}, or \code{"mixed"}.
#'
#' @return Returns the appropriate operator to use in a \code{\link{foreach::foreach}} loop. 
#' The \code{`%do%`} operator is returned if \code{strategy=="sequential"}. 
#' Otherwise, the \code{`%dopar%`} operator is returned.
#' @export
#'
#' @examples
get_doFuture_operator <- function(strategy=c("sequential", "multisession", "multicore", "mixed")){
  ## character integer -> ___
  ## get doFuture operator
  
  strategy <- match.arg(strategy)
  
  if(strategy == "sequential"){
    return(`%do%`)
  } else {
    return(`%dopar%`)
  }
}

#' Set a doFuture execution strategy
#'
#' @param strategy One of \code{"sequential"} (default), \code{"multisession"}, \code{"multicore"}, or \code{"mixed"}.
#' @param n_workers A positive numeric scalar or a function specifying the maximum number of parallel futures 
#' that can be active at the same time before blocking. 
#' If a function, it is called without arguments when the future is created and its value is used to configure the workers. 
#' The function should return a numeric scalar. 
#' Defaults to \code{\link{future::availableCores}()-1} if \code{NULL} (default), with \code{"multicore"} constraint in the relevant case. 
#' Ignored if \code{strategy=="sequential"}.
#'
#' @return The corresponding \code{\link{get_doFuture_operator}} operator to use in a \code{\link{foreach::foreach}} loop.
#' @export
#'
#' @examples
set_doFuture_strategy <- function(strategy=c("sequential", "multisession", "multicore", "mixed"),
                                  n_workers=NULL){
  strategy <- match.arg(strategy)
  
  doFuture::registerDoFuture()
  if(strategy == "sequential"){
    future::plan(future::sequential)
    
  } else if (strategy == "multisession"){
    if(is.null(n_workers)){
      n_workers <- max(future::availableCores() - 1, 1)
    }
    future::plan(future::multisession, workers = n_workers)
    
  } else if (strategy == "multicore"){
    if(is.null(n_workers)){
      n_workers <- max(future::availableCores(constraints = "multicore") - 1, 1)
    }
    future::plan(future::multicore, workers = n_workers)
    
  } else if (strategy == "mixed"){
    if(is.null(n_workers)){
      n_workers <- max(future::availableCores() - 1, 1)
    }
    strategy_1 <- future::tweak(future::sequential)
    strategy_2 <- future::tweak(future::multisession, workers = n_workers)
    future::plan(list(strategy_1, strategy_2))
  }
  return(get_doFuture_operator(strategy))
}

#' End the currently set doFuture strategy
#'
#' @description Resets the default strategy using \code{future::plan("default")}.
#' 
#' @export
#'
#' @examples
end_doFuture_strategy <- function(){
  ## ends the doFuture execution strategy
  
  future::plan("default")
}

#' Start a doParallel execution strategy
#'
#' @param strategy One of \code{"sequential"} (default) or \code{"parallel"}.
#' @param n_workers Number of parallel workers as an integer.
#' Defaults to \code{\link{parallell::detectCores}()-1} if \code{NULL} (default). 
#' Ignored if \code{strategy=="sequential"}.
#'
#' @return A named list containing: 
#' \item{par_operator}{the relevant \code{\link{foreach::foreach}} loop operator,}
#' \item{cl}{the cluster object.}
#'
#' @examples
start_doParallel_strategy <- function(strategy=c("sequential", "parallel"),
                                      n_workers=NULL){
  
  strategy <- match.arg(strategy)
  
  if(is.null(n_workers)){
    n_workers <- max(parallell::detectCores() - 1, 1)
  }
  if(strategy=="parallel"){
    cl <- parallel::makeCluster(n_workers)
    doParallel::registerDoParallel(cl)
    `%fun%` <- `%dopar%`
  } else {
    cl <- NULL
    `%fun%` <- `%do%`
  }
  return(list(par_operator=`%fun%`, cl=cl))
}

#' Stop the current doParallel strategy
#' 
#' @description Stops the given cluster, using \code{\link{parallel::stopCluster}(cl)}, if \code{strategy=="parallel"}.
#'
#' @param strategy One of \code{"sequential"} (default) or \code{"parallel"}.
#' @param cl Cluster object, returned by \code{\link{start_doParallel_strategy}(strategy, ...)}.
#' 
#' @examples
stop_doParallel_strategy <- function(strategy=c("sequential", "parallel"), cl){
  ## character cluster -> ___
  ## closes the doParallel execution strategy
  
  strategy <- match.arg(strategy)
  if(strategy=="parallel"){
    parallel::stopCluster(cl)
  }
}
