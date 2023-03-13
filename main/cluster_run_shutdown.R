
# Olivier C. Pasche, University of Geneva, Jan 2021 (Olivier.Pasche@unige.ch).

library(optparse)
library(here)

# Specify the Rscript command call options in a list
option_list <- list(
  optparse::make_option(c("-s", "--script_path"), action="store", type="character", default=NULL,
                        help="Path of the R script file to run (with respect to the project root directory)",
                        metavar="character")
)
# Parse the desired option
parsed_options <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
script_path <- parsed_options$script_path

# Path and name of the log files containing the Rscript output and eventual errors
logfile_path <- here("last_cluster_run_shutdown_log.txt")
errfile_path <- here("last_cluster_run_shutdown_err.txt")
# Clean logs from last run
suppressWarnings(file.remove(logfile_path,errfile_path))

custom_error_raise <- function(condition, err_msg){
  if(condition){
    sink(errfile_path)
    cat("\n", err_msg, "\n")
    sink()
    stop(err_msg)
  }
}

custom_error_raise(is.null(script_path),
                   "Must specify script_path to 'cluster_run_shutdown' via the 'Rscript' command option flag -s or --script_path.")

# Path is relative to the project home directory
script_path <- here::here(script_path)

custom_error_raise(!file.exists(script_path),
                   "Specified script_path passed to 'cluster_run_shutdown' via the 'Rscript' command does not exist.")

# Try to run the desired script, logging its output and eventual errors
sink(logfile_path)
start_time_crs <- Sys.time()
try(source(script_path), outFile=errfile_path)
end_time_crs <- Sys.time()
cat("\nRun time (cluster_run_shutdown.R):\n")
print(end_time_crs - start_time_crs)
sink()

# Wait for 2 mins and then shutdown the machine
Sys.sleep(120)
#system('shutdown -t 30 -s')# For Windows
system('sudo poweroff')# For Linux (with non-root sudo user)
