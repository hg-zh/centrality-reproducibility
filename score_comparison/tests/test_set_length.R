script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("This script must be run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
module_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

source(file.path(module_root, "R", "set_geometry.R"))
source(file.path(module_root, "R", "set_length.R"))

grid <- seq(-3, 3, by = 0.25)
score <- pmax(1 - abs(grid + 1), 1 - abs(grid - 1))
intervals <- score_level_intervals(
  grid,
  score,
  threshold = 0.5,
  outside_margin_points = 2L
)
stopifnot(
  nrow(intervals) == 2L,
  max(abs(intervals$left - c(-1.5, 0.5))) < 1e-14,
  max(abs(intervals$right - c(-0.5, 1.5))) < 1e-14,
  abs(total_interval_length(intervals) - 2) < 1e-14
)

single_interval <- data.frame(left = -1.5, right = 1.5)
stopifnot(
  abs(interval_overlap_length(intervals, single_interval) - 2) < 1e-14,
  abs(
    interval_symmetric_difference_length(intervals, single_interval) - 1
  ) < 1e-14
)

normal_grid <- seq(-8, 8, by = 0.005)
normal_geometry <- oracle_hdr_geometry(
  normal_grid,
  stats::dnorm(normal_grid),
  probability = 0.95,
  cdf_function = stats::pnorm,
  outside_margin_points = 10L
)
expected_normal_length <- 2 * stats::qnorm(0.975)
stopifnot(
  normal_geometry$components == 1L,
  abs(normal_geometry$probability_content - 0.95) < 1e-7,
  abs(normal_geometry$total_length - expected_normal_length) < 0.002
)

expanded <- ensure_score_domain(
  initial_grid = seq(-1, 1, by = 0.1),
  score_function = function(point) {
    answer <- matrix(-point^2, ncol = 1L)
    colnames(answer) <- "quadratic"
    answer
  },
  thresholds = c(quadratic = -4),
  outside_margin_points = 3L,
  expansion_fraction = 0.25,
  max_expansions = 8L
)
stopifnot(
  expanded$expansions > 0L,
  min(expanded$grid) < -2,
  max(expanded$grid) > 2,
  expanded$left_outside_points >= 3L,
  expanded$right_outside_points >= 3L
)

fine_grid <- seq(-3, 3, by = 0.125)
fine_scores <- matrix(
  pmax(1 - abs(fine_grid + 1), 1 - abs(fine_grid - 1)),
  ncol = 1L,
  dimnames = list(NULL, "two_peaks")
)
refinement <- refined_geometry_comparison(
  coarse_geometry = list(
    grid = grid,
    lengths = c(two_peaks = total_interval_length(intervals)),
    components = c(two_peaks = nrow(intervals))
  ),
  fine_grid = fine_grid,
  fine_scores = fine_scores,
  thresholds = c(two_peaks = 0.5),
  outside_margin_points = 4L
)
stopifnot(
  refinement$absolute_length_change < 1e-14,
  refinement$component_count_stable
)

unbounded_intervals <- data.frame(left = -Inf, right = Inf)
stopifnot(
  is.infinite(total_interval_length(unbounded_intervals)),
  is.infinite(interval_symmetric_difference_length(
    unbounded_intervals,
    single_interval
  ))
)

message("Interpolated interval, total-length, HDR, and refinement checks passed.")
