nearest_grid_value <- function(grid, target) {
  which.min(abs(grid - target))
}

normalized_score <- function(score) {
  score <- as.numeric(score)
  extent <- diff(range(score))
  if (!is.finite(extent) || extent <= 100 * .Machine$double.eps) {
    return(rep(0.5, length(score)))
  }
  (score - min(score)) / extent
}

score_curve_metrics <- function(grid, oriented_score, density, modes) {
  score <- as.numeric(oriented_score)
  standardized <- normalized_score(score)
  center_score <- stats::approx(
    grid, standardized, xout = 0, rule = 2, ties = "ordered"
  )$y
  mode_scores <- stats::approx(
    grid, standardized, xout = modes, rule = 2, ties = "ordered"
  )$y
  optimum_indices <- which(score >= max(score) - 1e-12 * max(1, abs(max(score))))
  optimum_distance <- min(vapply(
    grid[optimum_indices],
    function(point) min(abs(point - modes)),
    numeric(1)
  ))
  list(
    density_spearman = suppressWarnings(stats::cor(
      score, density, method = "spearman"
    )),
    center_loss = max(standardized) - center_score,
    mode_minus_center = mean(mode_scores) - center_score,
    optimum_to_mode = optimum_distance,
    optimum_location = grid[optimum_indices[[which.min(abs(grid[optimum_indices]))]]]
  )
}

prediction_set_metrics <- function(
    intervals,
    test_membership,
    oracle_intervals,
    mode_membership,
    center_membership,
    antimode_membership = NA) {
  mode_membership <- as.logical(mode_membership)
  center_membership <- as.logical(center_membership)
  antimode_membership <- as.logical(antimode_membership)
  set_length <- total_interval_length(intervals)
  oracle_length <- total_interval_length(oracle_intervals)
  is_unbounded <- is.infinite(set_length)
  hdr_difference <- interval_symmetric_difference_length(
    intervals, oracle_intervals
  )
  list(
    coverage = mean(test_membership),
    total_set_length = set_length,
    finite_replication_total_length =
      if (is_unbounded) NA_real_ else set_length,
    set_unbounded = as.numeric(is_unbounded),
    components = nrow(intervals),
    hdr_symmetric_difference_length = hdr_difference,
    finite_replication_hdr_symmetric_difference_length =
      if (is_unbounded) NA_real_ else hdr_difference,
    oracle_hdr_length = oracle_length,
    excess_length = set_length - oracle_length,
    center_included = as.numeric(center_membership),
    valley_excluded = if (is.na(antimode_membership)) {
      NA_real_
    } else {
      as.numeric(!antimode_membership)
    },
    all_modes_included = as.numeric(all(mode_membership)),
    peak_inclusion_fraction = mean(mode_membership)
  )
}

summarize_with_mcse <- function(values) {
  values <- as.numeric(values)
  values <- values[!is.na(values)]
  if (!length(values)) return(c(estimate = NA_real_, mcse = NA_real_))
  if (any(is.infinite(values))) {
    return(c(estimate = Inf, mcse = NA_real_))
  }
  if (any(!is.finite(values))) stop("Monte Carlo summaries require finite values.")
  c(
    estimate = mean(values),
    mcse = if (length(values) > 1L) stats::sd(values) / sqrt(length(values)) else NA_real_
  )
}
