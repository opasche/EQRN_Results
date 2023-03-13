
#' Recurrent network module for GPD parameter prediction
#'
#' @description
#' A recurrent neural network as a \code{\link{torch::nn_module}}, 
#' designed for generalized Pareto distribution parameter prediction, with sequential dependence.
#' 
#' @export
#'
#' @details
#' The constructor allows specifying:
#' \item{type}{the type of recurrent architecture, can be one of \code{"lstm"} (default) or \code{"gru"},}
#' \item{nb_input_features}{the input size (i.e. the number of features),}
#' \item{hidden_size}{the dimension of the hidden latent state variables in the recurrent network,}
#' \item{num_layers}{the number of recurrent layers,}
#' \item{dropout}{probability parameter for dropout before each hidden layer for regularization during training,}
#' \item{shape_fixed}{whether the shape estimate depends on the covariates or not (bool).}
Recurrent_GPD_net <- nn_module(
  "Recurrent_GPD_net",
  initialize = function(type=c("lstm","gru"), nb_input_features, hidden_size,
                        num_layers=1, dropout=0, shape_fixed=FALSE) {
    
    self$type <- match.arg(type)
    self$num_layers <- num_layers
    self$shape_fixed <- shape_fixed
    
    self$rnn <- if (self$type == "gru") {
      nn_gru(input_size = nb_input_features, hidden_size = hidden_size,
             num_layers = num_layers, dropout = dropout, batch_first = TRUE)
    } else {
      nn_lstm(input_size = nb_input_features, hidden_size = hidden_size,
              num_layers = num_layers, dropout = dropout, batch_first = TRUE)
    }
    
    if(self$shape_fixed){
      self$FC_scale <- nn_linear(hidden_size, 1)
      self$FC_shape <- nn_linear(1, 1, bias = FALSE)
    }else{
      self$FC_out <- nn_linear(hidden_size, 2)
    }
  },
  
  forward = function(x) {
    
    # list of [output, hidden(, cell)]
    # we keep the output, which is of size (batch_size, n_timesteps, hidden_size)
    x <- self$rnn(x)[[1]]
    
    # from the output, we only want the final timestep
    x <- x[ , dim(x)[2], ]# shape now is (batch_size, hidden_size)
    
    
    if(self$shape_fixed){
      x <- x %>% self$FC_scale()
      xi <- torch_ones(c(x$shape[1],1), device=device)
      xi <- xi %>% self$FC_shape()
      
      x <- torch_cat(list(x$exp(), xi), dim = 2)
      #x <- torch_cat(list(x$exp(), torch_tanh(xi) * 0.6 + 0.1), dim = 2)
    }else{
      # final shape then is (batch_size, 2)
      x <- x %>% self$FC_out()
      
      s <- torch_split(x, 1, dim = 2)
      x <- torch_cat(list(s[[1]]$exp(), torch_tanh(s[[2]]) * 0.6 + 0.1), dim = 2)
    }
    x
  }
)



#' Recurrent quantile regression neural network module
#'
#' @description
#' A recurrent neural network as a \code{\link{torch::nn_module}}, 
#' designed for quantile regression.
#' 
#' @export
#'
#' @details
#' The constructor allows specifying:
#' \item{type}{the type of recurrent architecture, can be one of \code{"lstm"} (default) or \code{"gru"},}
#' \item{nb_input_features}{the input size (i.e. the number of features),}
#' \item{hidden_size}{the dimension of the hidden latent state variables in the recurrent network,}
#' \item{num_layers}{the number of recurrent layers,}
#' \item{dropout}{probability parameter for dropout before each hidden layer for regularization during training.}
QRNN_RNN_net <- nn_module(
  "QRNN_RNN_net",
  initialize = function(type=c("lstm","gru"), nb_input_features, hidden_size, num_layers=1, dropout=0) {
    
    self$type <- match.arg(type)
    self$num_layers <- num_layers
    
    self$rnn <- if (self$type == "gru") {
      nn_gru(input_size = nb_input_features, hidden_size = hidden_size,
             num_layers = num_layers, dropout = dropout, batch_first = TRUE)
    } else {
      nn_lstm(input_size = nb_input_features, hidden_size = hidden_size,
              num_layers = num_layers, dropout = dropout, batch_first = TRUE)
    }
    
    self$FC_out <- nn_linear(hidden_size, 1)
    
  },
  
  forward = function(x) {
    
    # list of [output, hidden(, cell)]
    # we keep the output, which is of size (batch_size, n_timesteps, hidden_size)
    x <- self$rnn(x)[[1]]
    
    # from the output, we only want the final timestep
    x <- x[ , dim(x)[2], ]# shape now is (batch_size, hidden_size)
    
    # feed this to a single output neuron
    # final shape then is (batch_size, 1)
    x <- x %>% self$FC_out()
    x
  }
)

