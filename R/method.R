set_reproducible_rng <- function(seed) {
  RNGkind(
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  set.seed(as.integer(seed))
  invisible(seed)
}

make_probability_grid <- function(size, include_endpoints = TRUE) {
  size <- as.integer(size)
  if (include_endpoints) {
    if (size < 3L) stop("An endpoint-inclusive probability grid needs at least 3 points.")
    return(c(0, sort(stats::rbeta(size - 2L, 0.5, 0.5)), 1))
  }
  if (size < 1L) stop("The probability-grid size must be positive.")
  sort(stats::rbeta(size, 0.5, 0.5))
}

as_observation_matrix <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.vector(x)) x <- matrix(x, nrow = 1L)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("Observations must be supplied as a numeric matrix.")
  }
  storage.mode(x) <- "double"
  x
}

pairwise_distance <- function(query, reference, metric) {
  query <- as_observation_matrix(query)
  reference <- as_observation_matrix(reference)
  if (ncol(query) != ncol(reference)) {
    stop("Query and reference observations have different dimensions.")
  }

  if (identical(metric, "euclidean")) {
    squared <- outer(rowSums(query^2), rowSums(reference^2), "+") -
      2 * tcrossprod(query, reference)
    return(sqrt(pmax(squared, 0)))
  }

  if (identical(metric, "sphere")) {
    inner <- tcrossprod(query, reference)
    return(acos(pmin(pmax(inner, -1), 1)))
  }

  if (identical(metric, "wasserstein")) {
    squared <- outer(rowSums(query^2), rowSums(reference^2), "+") -
      2 * tcrossprod(query, reference)
    return(sqrt(pmax(squared, 0) / ncol(query)))
  }

  stop("Unknown metric: ", metric)
}

inverse_distance_profile <- function(distances, distance_grid, probability_grid) {
  if (anyNA(distances) || any(!is.finite(distances))) {
    stop("Distances must be finite.")
  }
  empirical_cdf <- findInterval(distance_grid, sort(distances)) / length(distances)
  stats::approx(
    x = empirical_cdf,
    y = distance_grid,
    xout = probability_grid,
    rule = 2,
    ties = min
  )$y
}

distance_quantile_profiles <- function(
    query,
    training,
    metric,
    distance_grid,
    probability_grid,
    chunk_size = 250L) {
  query <- as_observation_matrix(query)
  training <- as_observation_matrix(training)
  chunk_size <- max(1L, as.integer(chunk_size))
  answer <- matrix(
    NA_real_,
    nrow = nrow(query),
    ncol = length(probability_grid)
  )

  starts <- seq.int(1L, nrow(query), by = chunk_size)
  for (start in starts) {
    end <- min(nrow(query), start + chunk_size - 1L)
    distances <- pairwise_distance(query[start:end, , drop = FALSE], training, metric)
    block <- t(vapply(
      seq_len(nrow(distances)),
      function(i) inverse_distance_profile(
        distances[i, ],
        distance_grid,
        probability_grid
      ),
      numeric(length(probability_grid))
    ))
    answer[start:end, ] <- block
  }
  answer
}

fit_centrality_score <- function(
    training,
    metric,
    distance_grid,
    probability_grid,
    chunk_size = 250L) {
  training <- as_observation_matrix(training)
  quantiles <- distance_quantile_profiles(
    training,
    training,
    metric,
    distance_grid,
    probability_grid,
    chunk_size
  )
  structure(
    list(
      training = training,
      metric = metric,
      distance_grid = distance_grid,
      probability_grid = probability_grid,
      training_quantiles = quantiles,
      mean_quantile = colMeans(quantiles)
    ),
    class = "centrality_fit"
  )
}

centrality_scores_from_quantiles <- function(quantile_profiles, fit, tau) {
  if (!inherits(fit, "centrality_fit")) stop("fit must be a centrality_fit object.")
  quantile_profiles <- as_observation_matrix(quantile_profiles)
  if (ncol(quantile_profiles) != length(fit$probability_grid)) {
    stop("The quantile-profile grid does not match the fitted object.")
  }
  tau <- as.numeric(tau)
  weights <- vapply(
    tau,
    function(value) exp(-value * fit$probability_grid),
    numeric(length(fit$probability_grid))
  )
  if (is.null(dim(weights))) weights <- matrix(weights, ncol = 1L)

  # Equation (5) in the manuscript: average the weighted difference between
  # training and candidate inverse distance profiles, then apply expit.
  candidate_component <- quantile_profiles %*% weights
  training_component <- as.numeric(crossprod(fit$mean_quantile, weights))
  raw <- sweep(-candidate_component, 2L, training_component, "+") /
    length(fit$probability_grid)
  scores <- stats::plogis(raw)
  colnames(scores) <- format(tau, trim = TRUE, scientific = FALSE)
  scores
}

centrality_scores <- function(candidates, fit, tau, chunk_size = 250L) {
  profiles <- distance_quantile_profiles(
    candidates,
    fit$training,
    fit$metric,
    fit$distance_grid,
    fit$probability_grid,
    chunk_size
  )
  centrality_scores_from_quantiles(profiles, fit, tau)
}

calibration_thresholds <- function(
    calibration_scores,
    alpha,
    n_training,
    method = "legacy_figure_code",
    quantile_type = 7L) {
  calibration_scores <- as_observation_matrix(calibration_scores)
  if (identical(method, "legacy_figure_code")) {
    probability <- alpha * (1 - 1 / n_training)
    return(vapply(
      seq_len(ncol(calibration_scores)),
      function(j) stats::quantile(
        calibration_scores[, j],
        probs = probability,
        names = FALSE,
        type = quantile_type
      ),
      numeric(1)
    ))
  }
  if (identical(method, "manuscript_empirical")) {
    probability <- floor((nrow(calibration_scores) + 1) * alpha) /
      nrow(calibration_scores)
    return(vapply(
      seq_len(ncol(calibration_scores)),
      function(j) stats::quantile(
        calibration_scores[, j],
        probs = probability,
        names = FALSE,
        type = 1L
      ),
      numeric(1)
    ))
  }
  stop("Unknown calibration-threshold method: ", method)
}

evaluate_centrality_region <- function(
    fit,
    calibration,
    testing,
    search_grid,
    tau,
    alpha,
    size_scale,
    threshold_method = "legacy_figure_code",
    quantile_type = 7L,
    chunk_size = 250L) {
  calibration_q <- distance_quantile_profiles(
    calibration, fit$training, fit$metric,
    fit$distance_grid, fit$probability_grid, chunk_size
  )
  testing_q <- distance_quantile_profiles(
    testing, fit$training, fit$metric,
    fit$distance_grid, fit$probability_grid, chunk_size
  )
  search_q <- distance_quantile_profiles(
    search_grid, fit$training, fit$metric,
    fit$distance_grid, fit$probability_grid, chunk_size
  )

  calibration_scores <- centrality_scores_from_quantiles(calibration_q, fit, tau)
  testing_scores <- centrality_scores_from_quantiles(testing_q, fit, tau)
  search_scores <- centrality_scores_from_quantiles(search_q, fit, tau)
  thresholds <- calibration_thresholds(
    calibration_scores,
    alpha,
    nrow(fit$training),
    threshold_method,
    quantile_type
  )
  inside_test <- sweep(testing_scores, 2L, thresholds, ">=")
  inside_search <- sweep(search_scores, 2L, thresholds, ">=")

  list(
    coverage = colMeans(inside_test),
    size = colMeans(inside_search) * size_scale,
    thresholds = thresholds,
    search_membership = inside_search
  )
}
