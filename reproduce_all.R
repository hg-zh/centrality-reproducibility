script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("This script must be run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_argument[1L]), mustWork = TRUE)
project_root <- dirname(script_path)

source(file.path(project_root, "config", "simulation_settings.R"))
source(file.path(project_root, "R", "method.R"))
source(file.path(project_root, "R", "data_generators.R"))
source(file.path(project_root, "R", "simulation_helpers.R"))
source(file.path(project_root, "R", "plotting_helpers.R"))

arguments <- parse_release_arguments()
settings <- centrality_settings(arguments$profile)

reproduce_figures_1_2(settings, project_root)
run_all_simulations(
  settings,
  project_root,
  workers = arguments$workers,
  force = arguments$force
)
plot_all_simulations(settings, project_root)
message("End-to-end reproduction completed for profile '", settings$profile, "'.")
