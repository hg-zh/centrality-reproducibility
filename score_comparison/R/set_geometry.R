validate_regular_grid <- function(grid) {
  grid <- as.numeric(grid)
  if (length(grid) < 3L || any(!is.finite(grid)) ||
      is.unsorted(grid, strictly = TRUE)) {
    stop("A geometry grid must contain at least three increasing finite points.")
  }
  spacing <- diff(grid)
  tolerance <- 1e-10 * max(1, abs(spacing[[1L]]))
  if (max(abs(spacing - spacing[[1L]])) > tolerance) {
    stop("Set geometry requires an equally spaced grid.")
  }
  spacing[[1L]]
}

threshold_crossing <- function(x_left, x_right, score_left, score_right, threshold) {
  residual_left <- score_left - threshold
  residual_right <- score_right - threshold
  denominator <- residual_right - residual_left
  if (!is.finite(denominator) || abs(denominator) <= .Machine$double.eps) {
    return(mean(c(x_left, x_right)))
  }
  weight <- -residual_left / denominator
  weight <- min(1, max(0, weight))
  x_left + weight * (x_right - x_left)
}

score_level_intervals <- function(
    grid,
    oriented_score,
    threshold,
    outside_margin_points = 0L,
    enforce_margin = TRUE) {
  validate_regular_grid(grid)
  score <- as.numeric(oriented_score)
  if (length(score) != length(grid) || any(!is.finite(score)) ||
      length(threshold) != 1L || !is.finite(threshold)) {
    stop("The score curve or threshold is invalid.")
  }
  membership <- score >= threshold
  margin <- max(0L, as.integer(outside_margin_points))
  if (enforce_margin && margin > 0L) {
    if (2L * margin >= length(grid)) {
      stop("The requested outside margin is too large for the grid.")
    }
    if (any(head(membership, margin)) || any(tail(membership, margin))) {
      stop("The prediction set does not have the required outside margin.")
    }
  }
  if (!any(membership)) {
    return(data.frame(left = numeric(), right = numeric()))
  }

  runs <- rle(membership)
  run_ends <- cumsum(runs$lengths)
  run_starts <- c(1L, head(run_ends, -1L) + 1L)
  included <- which(runs$values)
  rows <- lapply(included, function(run_index) {
    start <- run_starts[[run_index]]
    end <- run_ends[[run_index]]
    left <- if (start == 1L) {
      grid[[1L]]
    } else {
      threshold_crossing(
        grid[[start - 1L]], grid[[start]],
        score[[start - 1L]], score[[start]], threshold
      )
    }
    right <- if (end == length(grid)) {
      grid[[length(grid)]]
    } else {
      threshold_crossing(
        grid[[end]], grid[[end + 1L]],
        score[[end]], score[[end + 1L]], threshold
      )
    }
    data.frame(left = left, right = right)
  })
  intervals <- do.call(rbind, rows)
  rownames(intervals) <- NULL
  intervals
}

empirical_depth_score_from_counts <- function(
    less,
    equal,
    greater,
    reference_size,
    method = c("spatial", "halfspace")) {
  method <- match.arg(method)
  if (identical(method, "spatial")) {
    return(1 - ((greater - less) / reference_size)^2)
  }
  pmin(less + equal, greater + equal) / reference_size
}

merge_closed_intervals <- function(intervals) {
  if (!nrow(intervals)) {
    return(data.frame(left = numeric(), right = numeric()))
  }
  intervals <- intervals[order(intervals$left, intervals$right), , drop = FALSE]
  merged <- vector("list", nrow(intervals))
  used <- 1L
  current_left <- intervals$left[[1L]]
  current_right <- intervals$right[[1L]]
  if (nrow(intervals) > 1L) {
    for (index in 2:nrow(intervals)) {
      if (intervals$left[[index]] <= current_right) {
        current_right <- max(current_right, intervals$right[[index]])
      } else {
        merged[[used]] <- data.frame(
          left = current_left,
          right = current_right
        )
        used <- used + 1L
        current_left <- intervals$left[[index]]
        current_right <- intervals$right[[index]]
      }
    }
  }
  merged[[used]] <- data.frame(left = current_left, right = current_right)
  answer <- do.call(rbind, merged[seq_len(used)])
  rownames(answer) <- NULL
  answer
}

empirical_depth_cells_r1 <- function(
    reference,
    method = c("spatial", "halfspace")) {
  reference <- as.numeric(reference)
  method <- match.arg(method)
  if (!length(reference) || any(!is.finite(reference))) {
    stop("The exact empirical-depth reference is invalid.")
  }

  support <- sort(unique(reference))
  frequency <- tabulate(match(reference, support), nbins = length(support))
  reference_size <- length(reference)
  cumulative <- cumsum(frequency)

  # The open cells are (-Inf,u_1), (u_1,u_2), ..., (u_K,Inf).
  # Within each cell there are no ties, so the depth is constant.
  cell_less <- c(0L, cumulative)
  cell_greater <- reference_size - cell_less
  cell_score <- empirical_depth_score_from_counts(
    less = cell_less,
    equal = integer(length(cell_less)),
    greater = cell_greater,
    reference_size = reference_size,
    method = method
  )

  # Scores at support points retain the exact method-specific tie convention.
  point_less <- c(0L, head(cumulative, -1L))
  point_greater <- reference_size - point_less - frequency
  point_score <- empirical_depth_score_from_counts(
    less = point_less,
    equal = frequency,
    greater = point_greater,
    reference_size = reference_size,
    method = method
  )

  gap_rows <- lapply(seq_along(cell_score), function(cell_index) {
    passed_support_points <- cell_index - 1L
    data.frame(
      kind = "gap",
      left = if (passed_support_points == 0L) {
        -Inf
      } else {
        support[[passed_support_points]]
      },
      right = if (passed_support_points == length(support)) {
        Inf
      } else {
        support[[passed_support_points + 1L]]
      },
      score = cell_score[[cell_index]],
      stringsAsFactors = FALSE
    )
  })
  atom_rows <- lapply(seq_along(support), function(point_index) {
    data.frame(
      kind = "atom",
      left = support[[point_index]],
      right = support[[point_index]],
      score = point_score[[point_index]],
      stringsAsFactors = FALSE
    )
  })
  answer <- do.call(rbind, c(gap_rows, atom_rows))
  rownames(answer) <- NULL
  answer
}

empirical_depth_intervals_r1 <- function(
    reference,
    threshold,
    method = c("spatial", "halfspace")) {
  threshold <- as.numeric(threshold)
  method <- match.arg(method)
  if (length(threshold) != 1L || is.na(threshold)) {
    stop("The exact empirical-depth threshold is invalid.")
  }
  cells <- empirical_depth_cells_r1(reference, method)
  accepted <- cells[cells$score >= threshold, , drop = FALSE]
  if (!nrow(accepted)) {
    return(data.frame(left = numeric(), right = numeric()))
  }
  # Existing metric code represents sets as closed intervals. Taking the
  # closure of each accepted gap preserves exact Lebesgue length and component
  # geometry. Exact membership at training atoms is evaluated from the score
  # itself, so spatial-depth open endpoints are not silently accepted.
  merge_closed_intervals(accepted[c("left", "right")])
}

points_in_intervals <- function(points, intervals, tolerance = 1e-12) {
  points <- as.numeric(points)
  if (!nrow(intervals)) return(rep(FALSE, length(points)))
  vapply(points, function(point) {
    any(
      point >= intervals$left - tolerance &
        point <= intervals$right + tolerance
    )
  }, logical(1))
}

interval_probability_content <- function(intervals, cdf_function) {
  if (!nrow(intervals)) return(0)
  probabilities <- cdf_function(intervals$right) - cdf_function(intervals$left)
  if (any(!is.finite(probabilities)) || any(probabilities < -1e-10)) {
    stop("The interval probability calculation failed.")
  }
  sum(pmax(probabilities, 0))
}

interval_overlap_length <- function(first, second) {
  if (!nrow(first) || !nrow(second)) return(0)
  total <- 0
  for (first_index in seq_len(nrow(first))) {
    for (second_index in seq_len(nrow(second))) {
      total <- total + max(
        0,
        min(first$right[[first_index]], second$right[[second_index]]) -
          max(first$left[[first_index]], second$left[[second_index]])
      )
    }
  }
  total
}

interval_symmetric_difference_length <- function(first, second) {
  total_interval_length(first) + total_interval_length(second) -
    2 * interval_overlap_length(first, second)
}
