proposed_raw_scores_from_profiles <- function(profiles, fit, tau = fit$tau_grid) {
  if (!inherits(fit, "profile_score_fit")) {
    stop("fit must be a profile_score_fit object.")
  }
  profiles <- as.matrix(profiles)
  if (ncol(profiles) != length(fit$training)) {
    stop("Candidate and fitted profile sizes differ.")
  }
  index <- match(as.numeric(tau), fit$tau_grid)
  if (anyNA(index)) stop("Requested tau is absent from the fitted tau grid.")
  weights <- fit$weights[, index, drop = FALSE]
  candidate <- profiles %*% weights
  baseline <- as.numeric(crossprod(fit$mean_training_profile, weights))
  sweep(-candidate, 2L, baseline, "+")
}

proposed_scores_from_profiles <- function(profiles, fit, tau = fit$tau_grid) {
  scores <- stats::plogis(proposed_raw_scores_from_profiles(profiles, fit, tau))
  colnames(scores) <- format(as.numeric(tau), trim = TRUE, scientific = FALSE)
  scores
}

proposed_scores <- function(query, fit, tau = fit$tau_grid) {
  proposed_scores_from_profiles(profiles_for_fit(query, fit), fit, tau)
}

