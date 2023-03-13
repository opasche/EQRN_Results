

scripts <- c("main/other_models/grid_serach_gbex_oracle.R",
             "main/grid_run_EQRNN_oracle.R",
             "main/grid_run_EQRNN_3.R")

for (i in 1:length(scripts)) {
  sink(paste0("perform_sims_log_script",i,".txt"))
  cat(paste0("\n================================ SCRIPT ", i," ================================\n\n"))
  try(source(scripts[i]), outFile=paste0("perform_sims_errors_last_script_",i,".txt"))
  sink()
}

