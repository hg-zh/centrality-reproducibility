script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("This script must be run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
module_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
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

training <- c(-2, -1, 0, 1, 2)
tau_grid <- c(0, 1, 4)
fit <- fit_profile_scores(training, tau_grid, chunk_size = 2L)
stopifnot(
  identical(dim(fit$training_profiles), c(5L, 5L)),
  all(fit$training_profiles[, 1L] == 0),
  max(abs(colSums(fit$weights) - c(
    1,
    1 - exp(-1),
    (1 - exp(-4)) / 4
  ))) < 1e-14
)

query <- c(-1.5, 0, 0.75, 4)
profiles <- profiles_for_fit(query, fit)
c0 <- proposed_scores_from_profiles(profiles, fit, tau = 0)[, 1L]
rank <- transport_rank_scores_from_profiles(profiles, fit)
stopifnot(max(abs(c0 - rank)) < 1e-14)

c0_raw <- proposed_raw_scores_from_profiles(profiles, fit, tau = 0)[, 1L]
published_raw <- published_transport_rank_raw(query, training)
stopifnot(diff(range(published_raw - c0_raw)) < 1e-14)

fast_cost <- transport_cost_from_profiles(profiles, fit)
naive_cost <- vapply(seq_len(nrow(profiles)), function(index) {
  mean(abs(sweep(
    fit$training_profiles,
    2L,
    profiles[index, ],
    "-"
  )))
}, numeric(1))
stopifnot(max(abs(fast_cost - naive_cost)) < 1e-14)

scaled_fit <- fit_profile_scores(3 * training, tau_grid, chunk_size = 2L)
scaled_profiles <- profiles_for_fit(3 * query, scaled_fit)
stopifnot(
  max(abs(
    proposed_raw_scores_from_profiles(scaled_profiles, scaled_fit) -
      3 * proposed_raw_scores_from_profiles(profiles, fit)
  )) < 1e-13,
  max(abs(
    transport_cost_from_profiles(scaled_profiles, scaled_fit) -
      3 * fast_cost
  )) < 1e-13
)

calibration <- seq(0.1, 0.9, length.out = 20L)
threshold <- manuscript_threshold(calibration, alpha = 0.1, n_training = 20L)
probability <- floor((length(calibration) + 1) * 0.1) / length(calibration)
stopifnot(identical(
  unname(threshold),
  unname(stats::quantile(calibration, probability, type = 1L, names = FALSE))
))
stopifnot(is.infinite(manuscript_threshold(1:10, alpha = 0.05, n_training = 10L)))

symmetric_grid <- seq(-4, 4, by = 0.01)
symmetric_density <- stats::dnorm(symmetric_grid)
symmetric_hdr <- oracle_hdr_membership(
  symmetric_grid, symmetric_density, probability = 0.95
)
stopifnot(identical(symmetric_hdr, rev(symmetric_hdr)))

pilot_environment <- new.env(parent = baseenv())
sys.source(file.path(module_root, "config", "pilot_config.R"), pilot_environment)
tiny <- pilot_environment$comparison_config()
tiny$profile <- "smoke"
tiny$replicates <- 2L
tiny$alpha <- 0.1
tiny$n_train <- 30L
tiny$n_calibration <- 30L
tiny$n_test <- 80L
tiny$candidate_grid <- seq(-6, 6, by = 0.1)
tiny$tau_grid <- c(0, 2, 6)
tiny$chunk_size <- 100L
tiny$outside_margin_points <- 3L
tiny$refinement_replicates <- c(1L, 2L)
tiny$refinement_factor <- 2L
tiny$length_tolerance <- 0.11
tiny$scenarios <- tiny$scenarios[
  tiny$scenarios$id %in% c("normal", "mixture_250"), , drop = FALSE
]
tiny$population <- list(
  primary = list(
    reference_n = 81L,
    tuning_reference_n = 31L,
    grid_limits = c(-6, 6),
    grid_step = 0.1
  ),
  refined = list(
    reference_n = 101L,
    tuning_reference_n = 41L,
    grid_limits = c(-6, 6),
    grid_step = 0.05
  )
)
validate_comparison_config(tiny)

temporary_output <- tempfile("score-comparison-smoke-")
dir.create(temporary_output)
on.exit(unlink(temporary_output, recursive = TRUE, force = TRUE), add = TRUE)
run_population_checks(tiny, temporary_output, workers = 1L)
run_comparison_simulation(tiny, temporary_output, workers = 1L, force = TRUE)
aggregate_comparison_results(tiny, temporary_output, temporary_output)
make_results_markdown(
  tiny,
  temporary_output,
  file.path(temporary_output, "results_table.md")
)

expected <- c(
  "population_summary.csv", "population_refinement.csv",
  "population_manifest.rds",
  "raw_replicates.csv", "tau_selection.csv", "equivalence_checks.csv",
  "depth_geometry_diagnostics.csv", "depth_geometry_summary.csv",
  "aggregate_summary.csv", "tau_summary.csv", "validation_checks.csv",
  "results_table.md"
)
stopifnot(all(file.exists(file.path(temporary_output, expected))))
validation <- utils::read.csv(file.path(temporary_output, "validation_checks.csv"))
stopifnot(all(validation$passed))
raw <- utils::read.csv(file.path(temporary_output, "raw_replicates.csv"))
stopifnot(
  nrow(raw) == nrow(tiny$scenarios) * tiny$replicates * 5L,
  all(table(raw$scenario, raw$replicate) == 5L)
)
depth_diagnostics <- utils::read.csv(
  file.path(temporary_output, "depth_geometry_diagnostics.csv")
)
stopifnot(
  nrow(depth_diagnostics) ==
    nrow(tiny$scenarios) * tiny$replicates * 2L,
  all(depth_diagnostics$grid_exact_component_equal)
)

message("All score-comparison unit and reduced end-to-end checks passed.")
