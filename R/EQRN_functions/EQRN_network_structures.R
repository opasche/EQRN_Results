
#' MLP module for GPD parameter prediction
#'
#' @description
#' A fully-connected network (or multi-layer perception) as a \code{\link{torch::nn_module}}, 
#' designed for generalized Pareto distribution parameter prediction.
#' 
#' @export
#'
#' @details
#' The constructor allows specifying:
#' \item{D_in}{the input size (i.e. the number of features),}
#' \item{Hidden_vect}{a vector of integers whose length determines the number of layers in the neural network 
#' and entries the number of neurons in each corresponding successive layer,}
#' \item{activation}{the activation function for the hidden layers 
#' (should be either a callable function, preferably from the \code{torch} library),}
#' \item{p_drop}{probability parameter for dropout before each hidden layer for regularization during training,}
#' \item{shape_fixed}{whether the shape estimate depends on the covariates or not (bool).}
FC_GPD_net <- torch::nn_module(
  "FC_GPD_net",
  initialize = function(D_in, Hidden_vect=c(5,5,5), activation=nnf_sigmoid, p_drop=0,
                        shape_fixed=FALSE) {
    self$activ <- activation
    
    self$dropout <- nn_dropout(p=p_drop)
    self$shape_fixed <- shape_fixed
    
    dims_in <- c(D_in, Hidden_vect)
    
    layers <- list()
    for(i in seq_along(Hidden_vect)){
      layers[[i]] <- nn_linear(dims_in[i], Hidden_vect[i])
    }
    self$linears <- nn_module_list(layers)
    
    if(self$shape_fixed){
      self$linSIGMA <- nn_linear(last_elem(Hidden_vect), 1)
      self$linXI <- nn_linear(1, 1, bias = FALSE)
    }else{
      self$linOUT <- nn_linear(last_elem(Hidden_vect), 2)
    }
  },
  forward = function(x) {
    for (i in seq_along(self$linears)){
      x <- self$linears[[i]](x)
      x <- self$activ(x)
      x <- self$dropout(x)
    }
    
    if(self$shape_fixed){
      xi <- torch_ones(c(x$shape[1],1), device=device)
      x <- x %>% self$linSIGMA()
      xi <- xi %>% self$linXI()
      x <- torch_cat(list(x$exp(), xi), dim = 2)
    }else{
      x <- x %>% self$linOUT()
      s <- torch_split(x, 1, dim = 2)
      x <- torch_cat(list(s[[1]]$exp(), torch_tanh(s[[2]]) * 0.6 + 0.1), dim = 2)
    }
    x
  }
)

#' Dropout module
#'
#' @description
#' A dropout layer as a \code{\link{torch::nn_module}}.
#'
#' @details
#' The constructor allows specifying:
#' \item{p}{probability of an element to be zeroed (default is 0.5),}
#' \item{inplace}{if set to TRUE, will do the operation in-place (default is FALSE).}
nn_dropout_nd <- torch::nn_module(
  "nn_dropout_nd",
  initialize = function(p = 0.5, inplace = FALSE) {
    if (p < 0 || p > 1)
      value_error("dropout probability has to be between 0 and 1 but got {p}")
    
    self$p <- p
    self$inplace <- inplace
    
  }
)

#' Alpha-dropout module
#'
#' @description
#' An alpha-dropout layer as a \code{\link{torch::nn_module}}, used in self-normalizing networks.
#' 
#' @export
#'
#' @details
#' The constructor allows specifying:
#' \item{p}{probability of an element to be zeroed (default is 0.5),}
#' \item{inplace}{if set to TRUE, will do the operation in-place (default is FALSE).}
#' 
#' @reference
#' Gunter Klambauer, Thomas Unterthiner, Andreas Mayr, Sepp Hochreiter. Self-Normalizing Neural Networks. 
#' Advances in Neural Information Processing Systems 30 (NIPS 2017), 2017.
nn_alpha_dropout <- torch::nn_module(
  "nn_alpha_dropout",
  inherit = nn_dropout_nd,
  forward = function(input) {
    torch::nnf_alpha_dropout(input, self$p, self$training, self$inplace)
  }
)

#' Self-normalized fully-connected network module for GPD parameter prediction
#'
#' @description
#' A fully-connected self-normalizing network (or multi-layer perception) as a \code{\link{torch::nn_module}}, 
#' designed for generalized Pareto distribution parameter prediction.
#' 
#' @export
#'
#' @details
#' The constructor allows specifying:
#' \item{D_in}{the input size (i.e. the number of features),}
#' \item{Hidden_vect}{a vector of integers whose length determines the number of layers in the neural network 
#' and entries the number of neurons in each corresponding successive layer,}
#' \item{p_drop}{probability parameter for the \code{alpha-dropout} before each hidden layer for regularization during training.}
#' 
#' @reference
#' Gunter Klambauer, Thomas Unterthiner, Andreas Mayr, Sepp Hochreiter. Self-Normalizing Neural Networks. 
#' Advances in Neural Information Processing Systems 30 (NIPS 2017), 2017.
FC_GPD_SNN <- torch::nn_module(
  "FC_GPD_SNN",
  initialize = function(D_in, Hidden_vect=c(64,64,64), p_drop=0.01) {
    self$activ <- nn_selu()
    self$negalphaprime <- 1.0507009873554804934193349852947 * 1.6732632423543772848170429916718
    
    #self$dropout <- nn_alpha_dropout(p=p_drop)
    self$p_drop <- p_drop
    self$dropout <- nnf_alpha_dropout
    
    dims_in <- c(D_in, Hidden_vect)
    
    layers <- list()
    for(i in seq_along(Hidden_vect)){
      layers[[i]] <- nn_linear(dims_in[i], Hidden_vect[i])
      nn_init_kaiming_normal_(layers[[i]]$weight, mode="fan_in", nonlinearity="linear")
      nn_init_zeros_(layers[[i]]$bias)
    }
    self$linears <- nn_module_list(layers)
    
    self$linOUT <- nn_linear(last_elem(Hidden_vect), 2)
    nn_init_kaiming_normal_(self$linOUT$weight, mode="fan_in", nonlinearity="linear")
    nn_init_zeros_(self$linOUT$bias)
  },
  forward = function(x) {
    for (i in seq_along(self$linears)){
      x <- self$linears[[i]](x)
      x <- self$activ(x)
      x <- self$dropout(x, self$p_drop, self$training)
    }
    
    x <- x %>%  self$linOUT()
    s <- torch_split(x, 1, dim = 2)
    x <- torch_cat(list(self$activ(s[[1]]) + self$negalphaprime, torch_tanh(s[[2]]*0.1) * 0.6 + 0.1), dim = 2)
    x
  }
)

#' Self-normalized separated network module for GPD parameter prediction
#'
#' @description
#' A parameter-separated self-normalizing network (or multi-layer perception) as a \code{\link{torch::nn_module}}, 
#' designed for generalized Pareto distribution parameter prediction.
#' 
#' @export
#'
#' @details
#' The constructor allows specifying:
#' \item{D_in}{the input size (i.e. the number of features),}
#' \item{Hidden_vect_scale}{a vector of integers whose length determines the number of layers in the sub-network 
#' for the scale parameter and entries the number of neurons in each corresponding successive layer,}
#' \item{Hidden_vect_shape}{a vector of integers whose length determines the number of layers in the sub-network 
#' for the shape parameter and entries the number of neurons in each corresponding successive layer,}
#' \item{p_drop}{probability parameter for the \code{alpha-dropout} before each hidden layer for regularization during training.}
#' 
#' @reference
#' Gunter Klambauer, Thomas Unterthiner, Andreas Mayr, Sepp Hochreiter. Self-Normalizing Neural Networks. 
#' Advances in Neural Information Processing Systems 30 (NIPS 2017), 2017.
Separated_GPD_SNN <- torch::nn_module(
  "Separated_GPD_SNN",
  initialize = function(D_in, Hidden_vect_scale=c(64,64,64), Hidden_vect_shape=c(5,3), p_drop=0.01) {
    # Parameters and activations
    self$activ_scale <- nn_selu()
    self$activ_shape <- torch_tanh
    self$negalphaprime <- 1.0507009873554804934193349852947 * 1.6732632423543772848170429916718
    
    # Dropout
    #self$dropout <- nn_alpha_dropout(p=p_drop)
    self$p_drop <- p_drop
    self$dropout <- nnf_alpha_dropout
    
    # Scale attributes
    
    dims_in_scale <- c(D_in, Hidden_vect_scale)
    
    layers_scale <- list()
    for(i in seq_along(Hidden_vect_scale)){
      layers_scale[[i]] <- nn_linear(dims_in_scale[i], Hidden_vect_scale[i])
      nn_init_kaiming_normal_(layers_scale[[i]]$weight, mode="fan_in", nonlinearity="linear")
      nn_init_zeros_(layers_scale[[i]]$bias)
    }
    self$linears_scale <- nn_module_list(layers_scale)
    
    self$linOUT_scale <- nn_linear(last_elem(Hidden_vect_scale), 1)
    nn_init_kaiming_normal_(self$linOUT_scale$weight, mode="fan_in", nonlinearity="linear")
    nn_init_zeros_(self$linOUT_scale$bias)
    
    # Shape attributes
    
    dims_in_shape <- c(D_in, Hidden_vect_shape)
    
    layers_shape <- list()
    for(i in seq_along(Hidden_vect_shape)){
      layers_shape[[i]] <- nn_linear(dims_in_shape[i], Hidden_vect_shape[i])
    }
    self$linears_shape <- nn_module_list(layers_shape)
    
    self$linOUT_shape <- nn_linear(last_elem(Hidden_vect_shape), 1)
  },
  forward = function(x) {
    # Scale
    scale_activ <- x
    
    for (i in seq_along(self$linears_scale)){
      scale_activ <- self$linears_scale[[i]](scale_activ)
      scale_activ <- self$activ_scale(scale_activ)
      scale_activ <- self$dropout(scale_activ, self$p_drop, self$training)
    }
    
    scale_activ <- scale_activ %>%  self$linOUT_scale()
    
    # Shape
    shape_activ <- x
    
    for (i in seq_along(self$linears_shape)){
      shape_activ <- self$linears_shape[[i]](shape_activ)
      shape_activ <- self$activ_shape(shape_activ)
    }
    
    shape_activ <- shape_activ %>%  self$linOUT_shape()
    
    # OUT
    x <- torch_cat(list(self$activ_scale(scale_activ) + self$negalphaprime, torch_tanh(shape_activ*0.1) * 0.6 + 0.1), dim = 2)
    #x <- torch_cat(list(self$activ_scale(scale_activ) + self$negalphaprime, torch_tanh(shape_activ) * 0.6 + 0.1), dim = 2)
    x
  }
)

