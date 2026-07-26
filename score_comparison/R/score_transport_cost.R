mean_absolute_profile_cost <- function(profiles, fit) {
  if (!inherits(fit, "profile_score_fit")) {
    stop("fit must be a profile_score_fit object.")
  }
  profiles <- as.matrix(profiles)
  n_anchor <- nrow(fit$training_profiles)
  profile_size <- ncol(fit$training_profiles)
  if (ncol(profiles) != profile_size) {
    stop("Candidate and fitted profile sizes differ.")
  }

  total <- numeric(nrow(profiles))
  for (column in seq_len(profile_size)) {
    reference <- sort(fit$training_profiles[, column])
    cumulative <- c(0, cumsum(reference))
    value <- profiles[, column]
    below <- findInterval(value, reference)
    lower_sum <- cumulative[below + 1L]
    upper_sum <- cumulative[n_anchor + 1L] - lower_sum
    deviation <- below * value - lower_sum +
      upper_sum - (n_anchor - below) * value
    total <- total + deviation / n_anchor
  }
  total / profile_size
}

transport_cost_from_profiles <- function(profiles, fit) {
  mean_absolute_profile_cost(profiles, fit)
}

transport_cost_centrality_from_profiles <- function(profiles, fit) {
  # The primary paper defines a cost, so lower is more central. Negating it
  # gives the common higher-is-more-central orientation used for calibration.
  -transport_cost_from_profiles(profiles, fit)
}

transport_cost <- function(query, fit) {
  transport_cost_from_profiles(profiles_for_fit(query, fit), fit)
}

transport_cost_centrality <- function(query, fit) {
  -transport_cost(query, fit)
}

