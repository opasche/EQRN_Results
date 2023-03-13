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

if (cuda_is_available()) {
  device <- torch_device("cuda")
} else {
  device <- torch_device("cpu")
}

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
