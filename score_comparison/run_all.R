script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("This script must be run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
module_root <- normalizePath(dirname(script_path), mustWork = TRUE)
clean_root <- normalizePath(file.path(module_root, ".."), mustWork = TRUE)

source(file.path(clean_root, "R", "method.R"))
for (file in c(
  "utilities.R", "distributions.R", "profiles.R", "score_proposed.R",
  "score_transport_rank.R", "score_transport_cost.R",
  "score_metric_spatial_depth.R", "score_metric_halfspace_depth.R",
  "set_geometry.R", "set_length.R", "prediction_sets.R", "metrics.R",
  "simulation.R", "population.R", "aggregation.R", "markdown_tables.R"
)) {
  source(file.path(module_root, "R", file))
}

arguments <- parse_comparison_arguments()
config <- load_comparison_config(module_root, arguments$config)
output_dir <- resolve_from_root(module_root, arguments$output_dir, must_work = FALSE)
run_population_checks(config, output_dir, arguments$workers)
run_comparison_simulation(
  config,
  output_dir,
  workers = arguments$workers,
  force = arguments$force
)
aggregate_comparison_results(config, output_dir, output_dir)
make_results_markdown(
  config,
  output_dir,
  file.path(output_dir, "results_table.md")
)
message("The complete score comparison finished successfully.")
