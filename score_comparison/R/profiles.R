sorted_distance_profiles <- function(query, training, chunk_size = 1000L) {
  query <- as.numeric(query)
  training <- as.numeric(training)
  if (!length(query) || !length(training) ||
      any(!is.finite(query)) || any(!is.finite(training))) {
    stop("Profile inputs must be nonempty and finite.")
  }
  chunk_size <- max(1L, as.integer(chunk_size))
  answer <- matrix(NA_real_, nrow = length(query), ncol = length(training))
  starts <- seq.int(1L, length(query), by = chunk_size)
  for (start in starts) {
    end <- min(length(query), start + chunk_size - 1L)
    distances <- abs(outer(query[start:end], training, "-"))
    if (end == start) {
      answer[start, ] <- sort(distances)
    } else {
      answer[start:end, ] <- t(apply(distances, 1L, sort))
    }
  }
  answer
}

profile_integral_weights <- function(profile_size, tau) {
  profile_size <- as.integer(profile_size)
  tau <- as.numeric(tau)
  if (profile_size < 1L || any(!is.finite(tau)) || any(tau < 0)) {
    stop("Profile size and tau values are invalid.")
  }
  lower <- (seq_len(profile_size) - 1) / profile_size
  upper <- seq_len(profile_size) / profile_size
  weights <- vapply(
    tau,
    function(value) {
      if (value == 0) {
        rep(1 / profile_size, profile_size)
      } else {
        exp(-value * lower) * (-expm1(-value / profile_size)) / value
      }
    },
    numeric(profile_size)
  )
  if (is.null(dim(weights))) weights <- matrix(weights, ncol = 1L)
  colnames(weights) <- format(tau, trim = TRUE, scientific = FALSE)
  weights
}

fit_profile_scores <- function(training, tau_grid, chunk_size = 1000L) {
  training <- as.numeric(training)
  profiles <- sorted_distance_profiles(training, training, chunk_size)
  weights <- profile_integral_weights(length(training), tau_grid)
  structure(
    list(
      training = training,
      training_profiles = profiles,
      mean_training_profile = colMeans(profiles),
      tau_grid = as.numeric(tau_grid),
      weights = weights,
      chunk_size = as.integer(chunk_size)
    ),
    class = "profile_score_fit"
  )
}

profiles_for_fit <- function(query, fit) {
  if (!inherits(fit, "profile_score_fit")) {
    stop("fit must be a profile_score_fit object.")
  }
  sorted_distance_profiles(query, fit$training, fit$chunk_size)
}
