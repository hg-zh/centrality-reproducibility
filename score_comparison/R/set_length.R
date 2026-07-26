total_interval_length <- function(intervals) {
  if (!nrow(intervals)) return(0)
  widths <- intervals$right - intervals$left
  if (anyNA(widths) || any(widths < -1e-12)) {
    stop("Prediction-set intervals are invalid.")
  }
  if (any(is.infinite(widths))) return(Inf)
  sum(pmax(widths, 0))
}

ensure_score_domain <- function(
    initial_grid,
    score_function,
    thresholds,
    outside_margin_points = 10L,
    expansion_fraction = 0.25,
    max_expansions = 8L,
    unbounded_methods = character()) {
  spacing <- validate_regular_grid(initial_grid)
  grid <- as.numeric(initial_grid)
  threshold_names <- names(thresholds)
  thresholds <- as.numeric(thresholds)
  names(thresholds) <- threshold_names
  unbounded <- names(thresholds) %in% unbounded_methods
  bounded_columns <- which(!unbounded)
  margin <- max(1L, as.integer(outside_margin_points))
  max_expansions <- max(0L, as.integer(max_expansions))
  expansion_fraction <- as.numeric(expansion_fraction)
  if (!is.finite(expansion_fraction) || expansion_fraction <= 0) {
    stop("The domain expansion fraction must be positive.")
  }

  for (iteration in 0:max_expansions) {
    scores <- as.matrix(score_function(grid))
    if (nrow(scores) != length(grid) ||
        ncol(scores) != length(thresholds) ||
        any(!is.finite(scores))) {
      stop("The score callback returned an invalid matrix.")
    }
    membership <- sweep(scores, 2L, thresholds, ">=")
    margin_here <- min(margin, floor((length(grid) - 1L) / 2L))
    left_ok <- rep(TRUE, ncol(membership))
    right_ok <- rep(TRUE, ncol(membership))
    if (length(bounded_columns)) {
      left_ok[bounded_columns] <- colSums(
        membership[
          seq_len(margin_here),
          bounded_columns,
          drop = FALSE
        ]
      ) == 0L
    }
    right_indices <- seq.int(length(grid) - margin_here + 1L, length(grid))
    if (length(bounded_columns)) {
      right_ok[bounded_columns] <- colSums(
        membership[
          right_indices,
          bounded_columns,
          drop = FALSE
        ]
      ) == 0L
    }
    if (all(left_ok) && all(right_ok)) {
      colnames(scores) <- names(thresholds)
      return(list(
        grid = grid,
        scores = scores,
        expansions = iteration,
        left_outside_points = min(vapply(
          bounded_columns,
          function(column) {
            included <- which(membership[, column])
            if (!length(included)) length(grid) else min(included) - 1L
          },
          integer(1)
        )),
        right_outside_points = min(vapply(
          bounded_columns,
          function(column) {
            included <- which(membership[, column])
            if (!length(included)) length(grid) else length(grid) - max(included)
          },
          integer(1)
        )),
        unbounded_methods = names(thresholds)[unbounded]
      ))
    }
    if (iteration == max_expansions) break
    block <- max(2L * margin, ceiling(length(grid) * expansion_fraction))
    if (any(!left_ok)) {
      grid <- c(grid[[1L]] - spacing * rev(seq_len(block)), grid)
    }
    if (any(!right_ok)) {
      grid <- c(grid, tail(grid, 1L) + spacing * seq_len(block))
    }
  }
  stop("The score domain could not be expanded to an outside margin.")
}

population_level_set_geometry <- function(
    grid,
    oriented_score,
    probability,
    cdf_function,
    outside_margin_points = 10L,
    root_tolerance = 1e-10) {
  validate_regular_grid(grid)
  score <- as.numeric(oriented_score)
  if (length(score) != length(grid) || any(!is.finite(score)) ||
      probability <= 0 || probability >= 1) {
    stop("The population level-set inputs are invalid.")
  }
  content_at <- function(threshold) {
    intervals <- score_level_intervals(
      grid, score, threshold,
      outside_margin_points = 0L,
      enforce_margin = FALSE
    )
    interval_probability_content(intervals, cdf_function)
  }
  lower <- min(score)
  upper <- max(score)
  lower_content <- content_at(lower)
  upper_content <- content_at(upper)
  if (lower_content < probability || upper_content > probability) {
    stop("The population score grid does not bracket the requested content.")
  }
  threshold <- stats::uniroot(
    function(value) content_at(value) - probability,
    interval = c(lower, upper),
    tol = root_tolerance
  )$root
  intervals <- score_level_intervals(
    grid, score, threshold,
    outside_margin_points = outside_margin_points,
    enforce_margin = TRUE
  )
  list(
    threshold = threshold,
    intervals = intervals,
    probability_content = interval_probability_content(intervals, cdf_function),
    total_length = total_interval_length(intervals),
    components = nrow(intervals)
  )
}

oracle_hdr_geometry <- function(
    grid,
    density,
    probability,
    cdf_function,
    outside_margin_points = 10L) {
  answer <- population_level_set_geometry(
    grid = grid,
    oriented_score = density,
    probability = probability,
    cdf_function = cdf_function,
    outside_margin_points = outside_margin_points
  )
  answer$density_cutoff <- answer$threshold
  answer
}

refined_geometry_comparison <- function(
    coarse_geometry,
    fine_grid,
    fine_scores,
    thresholds,
    outside_margin_points = 10L,
    unbounded_methods = character(),
    exact_intervals = list()) {
  fine_scores <- as.matrix(fine_scores)
  rows <- lapply(seq_len(ncol(fine_scores)), function(column) {
    method <- colnames(fine_scores)[[column]]
    is_unbounded <- method %in% unbounded_methods
    fine_intervals <- if (method %in% names(exact_intervals)) {
      exact_intervals[[method]]
    } else if (is_unbounded) {
      data.frame(left = -Inf, right = Inf)
    } else {
      score_level_intervals(
        fine_grid,
        fine_scores[, column],
        thresholds[[column]],
        outside_margin_points = outside_margin_points,
        enforce_margin = TRUE
      )
    }
    fine_length <- if (is_unbounded) {
      Inf
    } else {
      total_interval_length(fine_intervals)
    }
    absolute_change <- if (
      is.infinite(fine_length) &&
        is.infinite(coarse_geometry$lengths[[column]])
    ) {
      0
    } else {
      abs(fine_length - coarse_geometry$lengths[[column]])
    }
    data.frame(
      method = method,
      coarse_step = validate_regular_grid(coarse_geometry$grid),
      fine_step = validate_regular_grid(fine_grid),
      coarse_length = coarse_geometry$lengths[[column]],
      fine_length = fine_length,
      absolute_length_change = absolute_change,
      coarse_components = coarse_geometry$components[[column]],
      fine_components = nrow(fine_intervals),
      component_count_stable =
        coarse_geometry$components[[column]] == nrow(fine_intervals),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
