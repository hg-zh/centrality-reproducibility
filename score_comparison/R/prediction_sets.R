manuscript_threshold <- function(scores, alpha, n_training) {
  scores <- as.numeric(scores)
  if (any(!is.finite(scores))) stop("Calibration scores must be finite.")
  order_index <- floor((length(scores) + 1) * alpha)
  if (order_index == 0L) return(-Inf)
  as.numeric(calibration_thresholds(
    matrix(scores, ncol = 1L),
    alpha = alpha,
    n_training = n_training,
    method = "manuscript_empirical"
  ))
}

depth_to_nonconformity <- function(depth) {
  depth <- as.numeric(depth)
  if (any(!is.finite(depth))) stop("Depth values must be finite.")
  -depth
}

nonconformity_threshold <- function(scores, alpha) {
  scores <- as.numeric(scores)
  if (any(!is.finite(scores))) stop("Calibration scores must be finite.")
  index <- ceiling((length(scores) + 1) * (1 - alpha))
  if (index > length(scores)) return(Inf)
  sort(scores, partial = index)[[index]]
}

evaluate_nonconformity_scores <- function(
    calibration_scores,
    test_scores,
    grid_scores,
    alpha) {
  threshold <- nonconformity_threshold(calibration_scores, alpha)
  list(
    threshold = threshold,
    test_membership = as.numeric(test_scores) <= threshold,
    grid_membership = as.numeric(grid_scores) <= threshold
  )
}

calibration_threshold_vector <- function(score_matrix, alpha, n_training) {
  score_matrix <- as.matrix(score_matrix)
  answer <- vapply(seq_len(ncol(score_matrix)), function(column) {
    manuscript_threshold(score_matrix[, column], alpha, n_training)
  }, numeric(1))
  names(answer) <- colnames(score_matrix)
  answer
}

select_tau_by_midpoints <- function(
    fit,
    calibration_profiles,
    alpha,
    penalty = 0.0001) {
  calibration_scores <- proposed_scores_from_profiles(calibration_profiles, fit)
  thresholds <- vapply(
    seq_len(ncol(calibration_scores)),
    function(column) manuscript_threshold(
      calibration_scores[, column],
      alpha,
      length(fit$training)
    ),
    numeric(1)
  )

  pairs <- utils::combn(fit$training, 2L)
  midpoints <- colMeans(pairs)
  midpoint_profiles <- profiles_for_fit(midpoints, fit)
  midpoint_scores <- proposed_scores_from_profiles(midpoint_profiles, fit)
  proportions <- colMeans(sweep(midpoint_scores, 2L, thresholds, ">="))
  criterion <- proportions + penalty * fit$tau_grid
  selected_index <- which.min(criterion)

  list(
    tau = fit$tau_grid[[selected_index]],
    index = selected_index,
    proportion = proportions,
    criterion = criterion,
    thresholds = thresholds,
    midpoint_count = length(midpoints)
  )
}

evaluate_oriented_scores <- function(
    calibration_scores,
    test_scores,
    grid_scores,
    alpha,
    n_training) {
  threshold <- manuscript_threshold(calibration_scores, alpha, n_training)
  list(
    threshold = threshold,
    test_membership = as.numeric(test_scores) >= threshold,
    grid_membership = as.numeric(grid_scores) >= threshold
  )
}

count_grid_components <- function(membership) {
  membership <- as.logical(membership)
  if (!length(membership) || !any(membership)) return(0L)
  as.integer(sum(membership & !c(FALSE, membership[-length(membership)])))
}

oracle_hdr_membership <- function(grid, density, probability = 0.95) {
  grid <- as.numeric(grid)
  density <- as.numeric(density)
  if (length(grid) != length(density) || length(grid) < 3L ||
      any(!is.finite(density)) || any(density < 0)) {
    stop("The oracle density grid is invalid.")
  }
  spacing <- diff(grid)
  if (max(abs(spacing - spacing[[1L]])) > 1e-10) {
    stop("The oracle HDR requires an equally spaced grid.")
  }
  cell_mass <- density * spacing[[1L]]
  if (sum(cell_mass) < probability) {
    stop("The candidate grid contains less mass than the requested HDR.")
  }
  ordering <- order(density, decreasing = TRUE)
  crossing <- which(cumsum(cell_mass[ordering]) >= probability)[[1L]]
  cutoff <- density[ordering[[crossing]]]
  tolerance <- 100 * .Machine$double.eps * max(1, max(density))
  density >= cutoff - tolerance
}
