script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("This script must be run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_argument[1L]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

source(file.path(project_root, "config", "simulation_settings.R"))
source(file.path(project_root, "R", "method.R"))
source(file.path(project_root, "R", "data_generators.R"))
source(file.path(project_root, "R", "simulation_helpers.R"))

arguments <- parse_release_arguments()
settings <- centrality_settings(arguments$profile)
run_all_simulations(
  settings,
  project_root,
  workers = arguments$workers,
  force = arguments$force
)
message("Simulation results written under results/", settings$profile, ".")
