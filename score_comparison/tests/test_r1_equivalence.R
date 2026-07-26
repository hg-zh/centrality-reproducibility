script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("This script must be run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
module_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
clean_root <- normalizePath(file.path(module_root, ".."), mustWork = TRUE)

source(file.path(clean_root, "R", "method.R"))
source(file.path(module_root, "R", "score_metric_spatial_depth.R"))
source(file.path(module_root, "R", "score_metric_halfspace_depth.R"))
source(file.path(module_root, "R", "prediction_sets.R"))

probability <- seq(0.001, 0.999, length.out = 401L)
spatial_population <- metric_spatial_depth_continuous_r1(probability)
halfspace_population <- metric_halfspace_depth_continuous_r1(probability)
stopifnot(max(abs(
  spatial_population -
    4 * halfspace_population * (1 - halfspace_population)
)) < 1e-14)

reference_probability <- (seq_len(99L) - 0.5) / 99
reference <- stats::qnorm(reference_probability)
calibration <- seq(-2.87, 2.93, length.out = 100L) + 1e-7
testing <- seq(-3.11, 3.07, length.out = 151L) + 3e-7
grid <- seq(-4.03, 4.07, length.out = 501L) + 5e-7

spatial_calibration <- metric_spatial_depth(calibration, reference)
halfspace_calibration <- metric_halfspace_depth(calibration, reference)
spatial_testing <- metric_spatial_depth(testing, reference)
halfspace_testing <- metric_halfspace_depth(testing, reference)
spatial_grid <- metric_spatial_depth(grid, reference)
halfspace_grid <- metric_halfspace_depth(grid, reference)
stopifnot(
  max(abs(
    spatial_grid - 4 * halfspace_grid * (1 - halfspace_grid)
  )) < 1e-14,
  identical(
    rank(-spatial_grid, ties.method = "min"),
    rank(-halfspace_grid, ties.method = "min")
  )
)

spatial_evaluation <- evaluate_oriented_scores(
  spatial_calibration,
  spatial_testing,
  spatial_grid,
  alpha = 0.05,
  n_training = length(reference)
)
halfspace_evaluation <- evaluate_oriented_scores(
  halfspace_calibration,
  halfspace_testing,
  halfspace_grid,
  alpha = 0.05,
  n_training = length(reference)
)
stopifnot(
  identical(
    spatial_evaluation$test_membership,
    halfspace_evaluation$test_membership
  ),
  identical(
    spatial_evaluation$grid_membership,
    halfspace_evaluation$grid_membership
  )
)

spatial_nonconformity <- evaluate_nonconformity_scores(
  depth_to_nonconformity(spatial_calibration),
  depth_to_nonconformity(spatial_testing),
  depth_to_nonconformity(spatial_grid),
  alpha = 0.05
)
stopifnot(
  identical(
    spatial_evaluation$test_membership,
    spatial_nonconformity$test_membership
  ),
  identical(
    spatial_evaluation$grid_membership,
    spatial_nonconformity$grid_membership
  )
)

tied_reference <- c(0, 0, 0, 1)
stopifnot(
  metric_spatial_depth(0, tied_reference) == 15 / 16,
  metric_halfspace_depth(0, tied_reference) == 3 / 4
)

message("Continuous R1 ranking, region, tie, and orientation checks passed.")

