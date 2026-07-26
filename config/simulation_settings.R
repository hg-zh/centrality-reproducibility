centrality_settings <- function(profile = c("manuscript", "smoke")) {
  profile <- match.arg(profile)

  means <- list(c(2, 2), c(-2, -2))
  covariances <- list(
    matrix(c(0.5, 0.3, 0.3, 0.3), nrow = 2),
    matrix(c(0.3, 0, 0, 0.3), nrow = 2)
  )

  settings <- list(
    profile = profile,
    rng = list(
      kind = "Mersenne-Twister",
      normal_kind = "Inversion",
      sample_kind = "Rejection",
      replicate_seed_rule = "seed equals Monte Carlo replicate number"
    ),
    alpha = 0.05,
    threshold = list(
      method = "legacy_figure_code",
      quantile_type = 7
    ),
    figure1 = list(
      seed = 2024L,
      n = 200L,
      mean = c(0, 0),
      # This is the value in the script that produced the archived PDF.
      # The manuscript caption instead reports the second variance as 0.5.
      covariance = matrix(c(2, 0, 0, 1), nrow = 2),
      selected_points = rbind(center = c(0, 0), middle = c(2, 0), far = c(4, 0)),
      distance_grid = seq(0, 10, length.out = 501L),
      probability_grid_size = 501L,
      output = "Fig2_new_DP.pdf"
    ),
    figure2 = list(
      # The archived plotting script did not set a seed. This release adds one.
      seed = 2024L,
      n_train = 200L,
      n_calibration = 200L,
      mixing = 0.5,
      means = means,
      covariances = covariances,
      tau = c(0, 4),
      distance_grid = seq(0, 15, length.out = 501L),
      probability_grid_size = 501L,
      grid_n = 100L,
      grid_limits = c(-5, 5),
      output = "illu_2d_2plots.pdf"
    ),
    euclidean = list(
      replicates = 200L,
      sample_sizes = c(200L, 400L, 800L),
      n_test = 2000L,
      mixing = 0.5,
      heat_mixing = c(0.5, 0.8),
      # The archived generator visited these settings in order. Retaining the
      # order preserves the random stream for the displayed a=0.8 scenario.
      heat_generation_mixing = c(0.5, 0.6, 0.7, 0.8),
      means = means,
      covariances = covariances,
      tau = seq(0, 100, by = 2.5),
      heat_tau = c(0, 2, 4),
      distance_grid = seq(0, 15, length.out = 201L),
      probability_grid_size = 201L,
      grid_n = 100L,
      grid_limits = c(-5, 5),
      # The legacy code uses this normalization although the plotted grid
      # spans a square of area 100.
      size_scale = 25
    ),
    sphere = list(
      replicates = 200L,
      sample_sizes = c(200L, 400L, 800L),
      n_test = 2000L,
      mixing = 0.5,
      noise_sd = 0.2,
      tau = seq(0, 20, by = 1),
      heat_tau = 10,
      distance_grid = seq(0, 3.5, length.out = 101L),
      probability_grid_size = 101L,
      angular_grid_n = 100L,
      size_scale = 4 * pi
    ),
    wasserstein = list(
      replicates = 200L,
      sample_sizes = c(200L, 400L, 800L),
      n_test = 2000L,
      mean_location = 2,
      mean_sd = 0.2,
      sd_location = 0.5,
      sd_sd = 0.1,
      tau = seq(0, 20, by = 2),
      heat_tau = 10,
      distance_grid = seq(0, 3, length.out = 501L),
      probability_grid_size = 499L,
      mu_grid_n = 101L,
      sigma_grid_n = 101L,
      # The final plotting code displays ten times the stored grid proportion.
      size_scale = 10
    )
  )

  if (identical(profile, "smoke")) {
    settings$figure1$n <- 60L
    settings$figure1$distance_grid <- seq(0, 10, length.out = 121L)
    settings$figure1$probability_grid_size <- 61L

    settings$figure2$n_train <- 30L
    settings$figure2$n_calibration <- 30L
    settings$figure2$distance_grid <- seq(0, 15, length.out = 61L)
    settings$figure2$probability_grid_size <- 41L
    settings$figure2$grid_n <- 19L

    settings$euclidean$replicates <- 2L
    settings$euclidean$sample_sizes <- c(16L, 24L, 32L)
    settings$euclidean$n_test <- 60L
    settings$euclidean$tau <- c(0, 2, 4, 8)
    settings$euclidean$heat_tau <- c(0, 2, 4)
    settings$euclidean$distance_grid <- seq(0, 15, length.out = 41L)
    settings$euclidean$probability_grid_size <- 31L
    settings$euclidean$grid_n <- 17L

    settings$sphere$replicates <- 2L
    settings$sphere$sample_sizes <- c(16L, 24L, 32L)
    settings$sphere$n_test <- 60L
    settings$sphere$tau <- c(0, 2, 4, 10)
    settings$sphere$distance_grid <- seq(0, 3.5, length.out = 41L)
    settings$sphere$probability_grid_size <- 31L
    settings$sphere$angular_grid_n <- 12L

    settings$wasserstein$replicates <- 2L
    settings$wasserstein$sample_sizes <- c(16L, 24L, 32L)
    settings$wasserstein$n_test <- 60L
    settings$wasserstein$tau <- c(0, 2, 4, 10)
    settings$wasserstein$distance_grid <- seq(0, 3, length.out = 41L)
    settings$wasserstein$probability_grid_size <- 31L
    settings$wasserstein$mu_grid_n <- 17L
    settings$wasserstein$sigma_grid_n <- 13L
  }

  settings
}
