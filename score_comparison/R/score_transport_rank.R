transport_rank_raw_from_profiles <- function(profiles, fit) {
  if (!inherits(fit, "profile_score_fit")) {
    stop("fit must be a profile_score_fit object.")
  }
  profiles <- as.matrix(profiles)
  if (ncol(profiles) != length(fit$training)) {
    stop("Candidate and fitted profile sizes differ.")
  }
  baseline <- mean(fit$training_profiles)
  baseline - rowMeans(profiles)
}

transport_rank_scores_from_profiles <- function(profiles, fit) {
  stats::plogis(transport_rank_raw_from_profiles(profiles, fit))
}

transport_rank_scores <- function(query, fit) {
  transport_rank_scores_from_profiles(profiles_for_fit(query, fit), fit)
}

published_transport_rank_raw <- function(query, training) {
  query <- as.numeric(query)
  training <- as.numeric(training)
  if (length(training) < 2L) stop("Leave-one-out rank needs at least two anchors.")
  pairwise <- abs(outer(training, training, "-"))
  anchor_baseline <- sum(pairwise) / (length(training) * (length(training) - 1))
  anchor_baseline - rowMeans(abs(outer(query, training, "-")))
}

published_transport_rank_scores <- function(query, training) {
  stats::plogis(published_transport_rank_raw(query, training))
}

