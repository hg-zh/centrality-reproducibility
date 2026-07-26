empirical_side_counts <- function(query, reference, chunk_size = 1000L) {
  query <- as.numeric(query)
  reference <- as.numeric(reference)
  if (!length(query) || !length(reference) ||
      any(!is.finite(query)) || any(!is.finite(reference))) {
    stop("Depth inputs must be nonempty and finite.")
  }
  chunk_size <- max(1L, as.integer(chunk_size))
  less <- equal <- greater <- integer(length(query))
  starts <- seq.int(1L, length(query), by = chunk_size)
  for (start in starts) {
    end <- min(length(query), start + chunk_size - 1L)
    difference <- outer(query[start:end], reference, "-")
    if (end == start) difference <- matrix(difference, nrow = 1L)
    less[start:end] <- rowSums(difference > 0)
    equal[start:end] <- rowSums(difference == 0)
    greater[start:end] <- rowSums(difference < 0)
  }
  list(less = less, equal = equal, greater = greater)
}

metric_spatial_depth <- function(query, reference, chunk_size = 1000L) {
  counts <- empirical_side_counts(query, reference, chunk_size)
  size <- length(reference)
  1 - ((counts$greater - counts$less) / size)^2
}

metric_spatial_depth_direct <- function(query, reference) {
  query <- as.numeric(query)
  reference <- as.numeric(reference)
  if (!length(query) || !length(reference) ||
      any(!is.finite(query)) || any(!is.finite(reference))) {
    stop("Depth inputs must be nonempty and finite.")
  }
  vapply(query, function(candidate) {
    distance_one <- abs(reference - candidate)
    numerator <- outer(distance_one^2, distance_one^2, "+") -
      outer(reference, reference, function(x, y) (x - y)^2)
    denominator <- outer(distance_one, distance_one, "*")
    active <- outer(reference != candidate, reference != candidate, "&")
    h <- matrix(0, nrow = length(reference), ncol = length(reference))
    h[active] <- numerator[active] / denominator[active]
    1 - mean(h) / 2
  }, numeric(1))
}

metric_spatial_depth_population_r1 <- function(
    probability_below,
    probability_above) {
  probability_below <- as.numeric(probability_below)
  probability_above <- as.numeric(probability_above)
  if (length(probability_below) != length(probability_above) ||
      any(!is.finite(probability_below)) ||
      any(!is.finite(probability_above)) ||
      any(probability_below < 0 | probability_above < 0) ||
      any(probability_below + probability_above > 1 + 1e-12)) {
    stop("The supplied one-dimensional probabilities are invalid.")
  }
  1 - (probability_above - probability_below)^2
}

metric_spatial_depth_continuous_r1 <- function(cdf) {
  cdf <- as.numeric(cdf)
  if (any(!is.finite(cdf)) || any(cdf < 0 | cdf > 1)) {
    stop("CDF values must lie in [0,1].")
  }
  4 * cdf * (1 - cdf)
}

metric_spatial_training_depth <- function(reference, leave_one_out = FALSE) {
  reference <- as.numeric(reference)
  if (!leave_one_out) return(metric_spatial_depth(reference, reference))
  if (length(reference) < 2L) {
    stop("Leave-one-out depth requires at least two observations.")
  }
  vapply(seq_along(reference), function(index) {
    metric_spatial_depth(reference[[index]], reference[-index])
  }, numeric(1))
}

