# This script loads recursively the EQRN dependencies and routines for the './main/' scipts
# including the snapshot of the EQRN package, simulation functions and plotting helpers.
# Olivier PASCHE, 2022

library(tidyverse)
library(evd)
library(ismev)
library(torch)

library(future)
library(doFuture)

library(lubridate)
library(tsibble) # Tidy Temporal Data Frames and Tools
library(feasts) # Feature Extraction and Statistics for Time Series
library(EnvStats)
library(VGAM)

library(latex2exp)
library(egg)
#library(coro)
library(ggpubr)
library(scales)

#coro, mvtnorm, randtoolbox, # stats, tools

source("R/EQRN_functions/utils.R")
source("R/EQRN_functions/accuracy_metrics.R")

source("R/EQRN_functions/EVT_utils.R")

source("R/EQRN_functions/EQRN_network_structures.R")
source("R/EQRN_functions/EQRN.R")

source("R/EQRN_functions/EQRN_seq_network_structures.R")
source("R/EQRN_functions/QRN.R")
source("R/EQRN_functions/EQRN_seq.R")

source("R/simulation_functions/model_functions.R")
source("R/simulation_functions/series_model_functions.R")

source("R/plotting_functions/plot_utils.R")
source("R/plotting_functions/plot_helpers.R")
source("R/plotting_functions/plot_helpers_ts.R")


# Override for default torch device
#default_device <- function(){
#  if(torch::cuda_is_available()) {
#    device <- torch::torch_device("cuda")
#  } else {
#    device <- torch::torch_device("cpu")
#  }
#  return(device)
#}

# Legacy for default torch device
device <- default_device()

