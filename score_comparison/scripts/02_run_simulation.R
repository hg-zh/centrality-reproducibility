script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("This script must be run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
module_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
clean_root <- normalizePath(file.path(module_root, ".."), mustWork = TRUE)

source(file.path(clean_root, "R", "method.R"))
source(file.path(module_root, "R", "utilities.R"))
source(file.path(module_root, "R", "distributions.R"))
source(file.path(module_root, "R", "profiles.R"))
source(file.path(module_root, "R", "score_proposed.R"))
source(file.path(module_root, "R", "score_transport_rank.R"))
source(file.path(module_root, "R", "score_transport_cost.R"))
source(file.path(module_root, "R", "score_metric_spatial_depth.R"))
source(file.path(module_root, "R", "score_metric_halfspace_depth.R"))
source(file.path(module_root, "R", "set_geometry.R"))
source(file.path(module_root, "R", "set_length.R"))
source(file.path(module_root, "R", "prediction_sets.R"))
source(file.path(module_root, "R", "metrics.R"))
source(file.path(module_root, "R", "simulation.R"))

arguments <- parse_comparison_arguments()
config <- load_comparison_config(module_root, arguments$config)
output_dir <- resolve_from_root(module_root, arguments$output_dir, must_work = FALSE)
run_comparison_simulation(
  config,
  output_dir,
  workers = arguments$workers,
  force = arguments$force
)
message("Finite-sample simulation completed: ", arguments$output_dir)
