script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("This script must be run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
module_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

source(file.path(module_root, "R", "score_metric_spatial_depth.R"))
source(file.path(module_root, "R", "score_metric_halfspace_depth.R"))

reference <- c(-2, -1, 0, 0, 3)
query <- c(-3, -1, 0, 0.5, 4)
spatial_fast <- metric_spatial_depth(query, reference, chunk_size = 2L)
spatial_direct <- metric_spatial_depth_direct(query, reference)
halfspace_fast <- metric_halfspace_depth(query, reference, chunk_size = 2L)
halfspace_direct <- metric_halfspace_depth_direct(query, reference)
stopifnot(
  max(abs(spatial_fast - spatial_direct)) < 1e-14,
  max(abs(halfspace_fast - halfspace_direct)) < 1e-14
)

counts <- empirical_side_counts(0, reference)
stopifnot(
  counts$less == 2L,
  counts$equal == 2L,
  counts$greater == 1L,
  abs(metric_spatial_depth(0, reference) - (1 - (1 / 5)^2)) < 1e-14,
  abs(metric_halfspace_depth(0, reference) - 3 / 5) < 1e-14
)

atom_reference <- c(0, 0, 0, 1)
atom_spatial <- metric_spatial_depth(0, atom_reference)
atom_halfspace <- metric_halfspace_depth(0, atom_reference)
stopifnot(
  abs(atom_spatial - 15 / 16) < 1e-14,
  abs(atom_halfspace - 3 / 4) < 1e-14,
  abs(atom_spatial - 4 * atom_halfspace * (1 - atom_halfspace)) > 0.1
)

all_equal <- rep(2, 5)
stopifnot(
  metric_spatial_depth(2, all_equal) == 1,
  metric_halfspace_depth(2, all_equal) == 1
)

full_spatial <- metric_spatial_training_depth(reference)
loo_spatial <- metric_spatial_training_depth(reference, leave_one_out = TRUE)
full_halfspace <- metric_halfspace_training_depth(reference)
loo_halfspace <- metric_halfspace_training_depth(reference, leave_one_out = TRUE)
stopifnot(
  length(full_spatial) == length(reference),
  length(loo_spatial) == length(reference),
  length(full_halfspace) == length(reference),
  length(loo_halfspace) == length(reference),
  all(is.finite(c(
    full_spatial, loo_spatial, full_halfspace, loo_halfspace
  )))
)

message("Metric-depth direct, shortcut, atom, tie, and LOO checks passed.")

