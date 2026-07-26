script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("This script must be run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
module_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
clean_root <- normalizePath(file.path(module_root, ".."), mustWork = TRUE)

source(file.path(clean_root, "R", "method.R"))
source(file.path(module_root, "R", "score_metric_spatial_depth.R"))
source(file.path(module_root, "R", "score_metric_halfspace_depth.R"))
source(file.path(module_root, "R", "prediction_sets.R"))
source(file.path(module_root, "R", "set_geometry.R"))
source(file.path(module_root, "R", "set_length.R"))

normal_reference <- stats::qnorm((seq_len(101L) - 0.5) / 101)
normal_points <- c(0, 8)
normal_spatial <- metric_spatial_depth(normal_points, normal_reference)
normal_halfspace <- metric_halfspace_depth(normal_points, normal_reference)
stopifnot(
  normal_spatial[[1L]] > normal_spatial[[2L]],
  normal_halfspace[[1L]] > normal_halfspace[[2L]],
  depth_to_nonconformity(normal_spatial)[[1L]] <
    depth_to_nonconformity(normal_spatial)[[2L]]
)

calibration_scores <- seq(0, 1, length.out = 200L)
conformity_index <- floor((length(calibration_scores) + 1) * 0.05)
nonconformity_index <- ceiling(
  (length(calibration_scores) + 1) * (1 - 0.05)
)
conformity_threshold <- manuscript_threshold(
  calibration_scores,
  alpha = 0.05,
  n_training = 200L
)
nonconformity_cutoff <- nonconformity_threshold(
  depth_to_nonconformity(calibration_scores),
  alpha = 0.05
)
stopifnot(
  conformity_index == 10L,
  nonconformity_index == 191L,
  conformity_threshold == sort(calibration_scores)[[conformity_index]],
  -nonconformity_cutoff == conformity_threshold,
  identical(
    manuscript_threshold(1:10, alpha = 0.05, n_training = 10L),
    -Inf
  ),
  identical(nonconformity_threshold(1:10, alpha = 0.05), Inf)
)

tied_reference <- c(-2, -1, 0, 0, 3)
gap_probes <- c(-3, -1.5, -0.5, 1.5, 4)
support <- sort(unique(tied_reference))
for (method in c("spatial", "halfspace")) {
  cells <- empirical_depth_cells_r1(tied_reference, method)
  expected_gap <- if (method == "spatial") {
    metric_spatial_depth(gap_probes, tied_reference)
  } else {
    metric_halfspace_depth(gap_probes, tied_reference)
  }
  expected_atom <- if (method == "spatial") {
    metric_spatial_depth(support, tied_reference)
  } else {
    metric_halfspace_depth(support, tied_reference)
  }
  stopifnot(
    identical(cells$kind, c(rep("gap", 5L), rep("atom", 4L))),
    max(abs(cells$score[cells$kind == "gap"] - expected_gap)) < 1e-14,
    max(abs(cells$score[cells$kind == "atom"] - expected_atom)) < 1e-14
  )
}

for (method in c("spatial", "halfspace")) {
  below_tail <- empirical_depth_intervals_r1(
    tied_reference, -1e-12, method
  )
  at_tail <- empirical_depth_intervals_r1(tied_reference, 0, method)
  above_tail <- empirical_depth_intervals_r1(
    tied_reference, 1e-12, method
  )
  above_maximum <- empirical_depth_intervals_r1(
    tied_reference, 1 + 1e-12, method
  )
  stopifnot(
    identical(below_tail, data.frame(left = -Inf, right = Inf)),
    identical(at_tail, data.frame(left = -Inf, right = Inf)),
    all(is.finite(unlist(above_tail))),
    !nrow(above_maximum),
    is.infinite(total_interval_length(at_tail))
  )
}

halfspace_singleton <- empirical_depth_intervals_r1(
  tied_reference, 0.6, "halfspace"
)
spatial_open_endpoint_closure <- empirical_depth_intervals_r1(
  tied_reference, 0.96, "spatial"
)
stopifnot(
  identical(halfspace_singleton, data.frame(left = 0, right = 0)),
  identical(
    spatial_open_endpoint_closure,
    data.frame(left = -1, right = 0)
  ),
  total_interval_length(spatial_open_endpoint_closure) == 1
)

two_point_spatial <- empirical_depth_intervals_r1(c(0, 1), 1, "spatial")
stopifnot(
  identical(two_point_spatial, data.frame(left = 0, right = 1)),
  total_interval_length(two_point_spatial) == 1,
  metric_spatial_depth(c(0, 0.5, 1), c(0, 1)) ==
    c(0.75, 1, 0.75)
)

all_equal <- rep(2, 5L)
for (method in c("spatial", "halfspace")) {
  stopifnot(
    identical(
      empirical_depth_intervals_r1(all_equal, 0.5, method),
      data.frame(left = 2, right = 2)
    ),
    !nrow(empirical_depth_intervals_r1(all_equal, 1 + 1e-12, method))
  )
}

lattice_reference <- seq_len(25L)
lattice_threshold <- 7 / 25
lattice_interval <- empirical_depth_intervals_r1(
  lattice_reference,
  lattice_threshold,
  "halfspace"
)
stopifnot(
  identical(
    lattice_interval,
    data.frame(left = 7, right = 19)
  )
)

continuous_reference <- stats::qnorm((seq_len(200L) - 0.5) / 200)
halfspace_threshold <- 10 / 200
spatial_threshold <- 1 - ((200 - 2 * 10) / 200)^2
halfspace_interval <- empirical_depth_intervals_r1(
  continuous_reference,
  halfspace_threshold,
  "halfspace"
)
spatial_interval <- empirical_depth_intervals_r1(
  continuous_reference,
  spatial_threshold,
  "spatial"
)
continuous_grid <- seq(-4, 4, length.out = 801L) + 1e-7
stopifnot(
  identical(halfspace_interval, spatial_interval),
  total_interval_length(halfspace_interval) ==
    total_interval_length(spatial_interval),
  identical(
    metric_halfspace_depth(continuous_grid, continuous_reference) >=
      halfspace_threshold,
    metric_spatial_depth(continuous_grid, continuous_reference) >=
      spatial_threshold
  ),
  identical(
    points_in_intervals(continuous_grid, halfspace_interval),
    points_in_intervals(continuous_grid, spatial_interval)
  )
)

two_rays <- data.frame(
  left = c(-Inf, 1),
  right = c(-1, Inf)
)
empty <- data.frame(left = numeric(), right = numeric())
stopifnot(
  nrow(two_rays) == 2L,
  is.infinite(total_interval_length(two_rays)),
  total_interval_length(empty) == 0
)

message(
  "Exact depth cells, orientation, quantile, tails, ties, and length checks passed."
)
