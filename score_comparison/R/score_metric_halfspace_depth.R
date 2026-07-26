metric_halfspace_depth <- function(query, reference, chunk_size = 1000L) {
  counts <- empirical_side_counts(query, reference, chunk_size)
  size <- length(reference)
  pmin(counts$less + counts$equal, counts$greater + counts$equal) / size
}

metric_halfspace_depth_direct <- function(query, reference) {
  query <- as.numeric(query)
  reference <- as.numeric(reference)
  if (!length(query) || !length(reference) ||
      any(!is.finite(query)) || any(!is.finite(reference))) {
    stop("Depth inputs must be nonempty and finite.")
  }
  span <- max(1, diff(range(c(query, reference))))
  vapply(query, function(candidate) {
    anchors <- unique(c(reference, candidate - span, candidate + span))
    pairs <- expand.grid(
      first = anchors,
      second = anchors,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    contains <- abs(candidate - pairs$first) <=
      abs(candidate - pairs$second)
    probabilities <- vapply(which(contains), function(index) {
      mean(
        abs(reference - pairs$first[[index]]) <=
          abs(reference - pairs$second[[index]])
      )
    }, numeric(1))
    min(probabilities)
  }, numeric(1))
}

metric_halfspace_depth_population_r1 <- function(cdf, cdf_left = cdf) {
  cdf <- as.numeric(cdf)
  cdf_left <- as.numeric(cdf_left)
  if (length(cdf) != length(cdf_left) ||
      any(!is.finite(cdf)) || any(!is.finite(cdf_left)) ||
      any(cdf < 0 | cdf > 1) || any(cdf_left < 0 | cdf_left > 1) ||
      any(cdf_left > cdf + 1e-12)) {
    stop("CDF and left-limit values are invalid.")
  }
  pmin(cdf, 1 - cdf_left)
}

metric_halfspace_depth_continuous_r1 <- function(cdf) {
  metric_halfspace_depth_population_r1(cdf, cdf)
}

metric_halfspace_training_depth <- function(reference, leave_one_out = FALSE) {
  reference <- as.numeric(reference)
  if (!leave_one_out) return(metric_halfspace_depth(reference, reference))
  if (length(reference) < 2L) {
    stop("Leave-one-out depth requires at least two observations.")
  }
  vapply(seq_along(reference), function(index) {
    metric_halfspace_depth(reference[[index]], reference[-index])
  }, numeric(1))
}

