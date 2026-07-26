script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("This script must be run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_argument[1L]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

source(file.path(project_root, "config", "simulation_settings.R"))
source(file.path(project_root, "R", "method.R"))
source(file.path(project_root, "R", "data_generators.R"))
source(file.path(project_root, "R", "simulation_helpers.R"))
source(file.path(project_root, "R", "plotting_helpers.R"))

set_reproducible_rng(11L)
training <- matrix(c(-2, -1, 0, 1, 2), ncol = 1L)
distance_grid <- seq(0, 12, length.out = 81L)
probability_grid <- make_probability_grid(41L, include_endpoints = TRUE)
fit <- fit_centrality_score(
  training,
  "euclidean",
  distance_grid,
  probability_grid
)
stopifnot(
  inherits(fit, "centrality_fit"),
  identical(dim(fit$training_quantiles), c(5L, 41L)),
  all(diff(apply(fit$training_quantiles, 1L, identity)[, 1L]) >= 0)
)

candidate_scores <- centrality_scores(matrix(c(0, 10), ncol = 1L), fit, c(0, 4))
stopifnot(
  identical(dim(candidate_scores), c(2L, 2L)),
  all(is.finite(candidate_scores)),
  all(candidate_scores > 0 & candidate_scores < 1),
  all(candidate_scores[1L, ] > candidate_scores[2L, ])
)

shifted_fit <- fit_centrality_score(
  training + 7,
  "euclidean",
  distance_grid,
  probability_grid
)
shifted_scores <- centrality_scores(
  matrix(c(7, 17), ncol = 1L),
  shifted_fit,
  c(0, 4)
)
stopifnot(isTRUE(all.equal(candidate_scores, shifted_scores, tolerance = 1e-12)))

sphere_points <- rbind(c(1, 0, 0), c(0, 1, 0), c(-1, 0, 0))
sphere_distances <- pairwise_distance(sphere_points, sphere_points, "sphere")
stopifnot(
  max(abs(diag(sphere_distances))) < 2e-8,
  isTRUE(all.equal(sphere_distances, t(sphere_distances), tolerance = 1e-12)),
  abs(sphere_distances[1L, 3L] - pi) < 1e-12
)

u <- seq(0.01, 0.99, length.out = 101L)
normal_a <- stats::qnorm(u, mean = 0, sd = 1)
normal_b <- stats::qnorm(u, mean = 2, sd = 1)
w_distance <- pairwise_distance(normal_a, normal_b, "wasserstein")
stopifnot(abs(w_distance[1L, 1L] - 2) < 1e-12)

calibration <- matrix(seq(0.1, 0.9, length.out = 20L), ncol = 1L)
threshold <- calibration_thresholds(
  calibration,
  alpha = 0.05,
  n_training = 20L,
  method = "legacy_figure_code",
  quantile_type = 7L
)
expected <- stats::quantile(
  calibration[, 1L],
  0.05 * (1 - 1 / 20),
  names = FALSE,
  type = 7L
)
stopifnot(identical(unname(threshold), unname(expected)))

smoke_settings <- centrality_settings("smoke")
changed_alpha_settings <- smoke_settings
changed_alpha_settings$alpha <- 0.1
stopifnot(
  !identical(
    simulation_cache_configuration(smoke_settings, "euclidean"),
    simulation_cache_configuration(changed_alpha_settings, "euclidean")
  )
)
reproduce_figures_1_2(smoke_settings, project_root)
run_all_simulations(
  smoke_settings,
  project_root,
  workers = 1L,
  force = TRUE
)
plot_all_simulations(smoke_settings, project_root)

expected_files <- c(
  file.path("figures", "smoke", "Fig2_new_DP.pdf"),
  file.path("figures", "smoke", "illu_2d_2plots.pdf"),
  file.path("figures", "smoke", "mix_2d_1.pdf"),
  file.path("figures", "smoke", "mix_2d_heatmap_2x3.pdf"),
  file.path("figures", "smoke", "S2_W2_lines_2x2.pdf"),
  file.path("figures", "smoke", "S2W2_heat_2_by_1.pdf"),
  file.path("results", "smoke", "figure3_euclidean.rds"),
  file.path("results", "smoke", "figure4_sphere.rds"),
  file.path("results", "smoke", "figure4_wasserstein.rds"),
  file.path("results", "smoke", "manifest.rds")
)
stopifnot(all(file.exists(file.path(project_root, expected_files))))

message("All method checks and reduced end-to-end checks passed.")
