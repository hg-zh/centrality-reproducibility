generate_gaussian_mixture <- function(n, mixing, means, covariances) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for Gaussian-mixture generation.")
  }
  answer <- matrix(NA_real_, nrow = n, ncol = length(means[[1L]]))
  for (i in seq_len(n)) {
    component <- stats::rbinom(1L, 1L, mixing) + 1L
    answer[i, ] <- MASS::mvrnorm(1L, means[[component]], covariances[[component]])
  }
  answer
}

generate_gaussian_mixture_pair <- function(n, mixing, means, covariances) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for Gaussian-mixture generation.")
  }
  training <- matrix(NA_real_, nrow = n, ncol = length(means[[1L]]))
  calibration <- training
  for (i in seq_len(n)) {
    component <- stats::rbinom(1L, 1L, mixing) + 1L
    training[i, ] <- MASS::mvrnorm(1L, means[[component]], covariances[[component]])
    component <- stats::rbinom(1L, 1L, mixing) + 1L
    calibration[i, ] <- MASS::mvrnorm(1L, means[[component]], covariances[[component]])
  }
  list(training = training, calibration = calibration)
}

make_euclidean_search_grid <- function(grid_n, limits = c(-5, 5)) {
  grid <- expand.grid(
    x = seq(limits[1L], limits[2L], length.out = grid_n),
    y = seq(limits[1L], limits[2L], length.out = grid_n)
  )
  list(data = grid, observations = as.matrix(grid))
}

spherical_to_cartesian <- function(theta, phi) {
  c(
    cos(phi) * sin(theta),
    sin(phi) * sin(theta),
    cos(theta)
  )
}

cartesian_to_spherical <- function(points) {
  points <- as_observation_matrix(points)
  radius <- sqrt(rowSums(points^2))
  phi <- atan2(points[, 2L], points[, 1L])
  phi[phi < 0] <- phi[phi < 0] + 2 * pi
  theta <- acos(pmax(-1, pmin(1, points[, 3L] / radius)))
  data.frame(phi = phi, theta = theta)
}

sphere_component_means <- function() {
  rbind(
    spherical_to_cartesian(pi / 4, pi / 2),
    spherical_to_cartesian(3 * pi / 4, 3 * pi / 2)
  )
}

generate_sphere_mixture <- function(n, mixing = 0.5, noise_sd = 0.2) {
  means <- sphere_component_means()
  answer <- matrix(NA_real_, nrow = n, ncol = 3L)
  for (i in seq_len(n)) {
    shared <- stats::rnorm(1L, 0, noise_sd)
    tangent <- c(stats::rnorm(1L, 0, noise_sd), shared, -shared)
    tangent_norm <- sqrt(sum(tangent^2))
    component <- stats::rbinom(1L, 1L, mixing) + 1L
    answer[i, ] <- cos(tangent_norm) * means[component, ] +
      sin(tangent_norm) * tangent / tangent_norm
  }
  answer
}

generate_sphere_mixture_pair <- function(n, mixing = 0.5, noise_sd = 0.2) {
  means <- sphere_component_means()
  one_draw <- function() {
    shared <- stats::rnorm(1L, 0, noise_sd)
    tangent <- c(stats::rnorm(1L, 0, noise_sd), shared, -shared)
    tangent_norm <- sqrt(sum(tangent^2))
    component <- stats::rbinom(1L, 1L, mixing) + 1L
    cos(tangent_norm) * means[component, ] +
      sin(tangent_norm) * tangent / tangent_norm
  }
  training <- matrix(NA_real_, nrow = n, ncol = 3L)
  calibration <- training
  for (i in seq_len(n)) {
    training[i, ] <- one_draw()
    calibration[i, ] <- one_draw()
  }
  list(training = training, calibration = calibration)
}

make_sphere_search_grid <- function(angular_grid_n) {
  theta_search <- seq(0, pi, length.out = angular_grid_n)
  phi_search <- seq(0, 2 * pi, length.out = angular_grid_n)
  x <- outer(cos(theta_search), cos(phi_search), "*")
  y <- outer(cos(theta_search), sin(phi_search), "*")
  z <- outer(sin(theta_search), rep(1, length(phi_search)), "*")
  observations <- rbind(
    cbind(c(x), c(y), c(z)),
    cbind(c(x), c(y), -c(z))
  )
  list(
    data = cbind(
      as.data.frame(observations),
      cartesian_to_spherical(observations)
    ),
    observations = observations
  )
}

normal_parameters_to_quantiles <- function(parameters, probability_grid) {
  parameters <- as.data.frame(parameters)
  z <- stats::qnorm(probability_grid)
  outer(parameters$sigma, z, "*") + parameters$mu
}

generate_wasserstein_parameters <- function(
    n,
    mean_location,
    mean_sd,
    sd_location,
    sd_sd) {
  data.frame(
    mu = stats::rnorm(n, mean_location, mean_sd),
    sigma = abs(stats::rnorm(n, sd_location, sd_sd))
  )
}

generate_wasserstein_pair <- function(n, settings, probability_grid) {
  training_parameters <- generate_wasserstein_parameters(
    n,
    settings$mean_location,
    settings$mean_sd,
    settings$sd_location,
    settings$sd_sd
  )
  calibration_parameters <- generate_wasserstein_parameters(
    n,
    settings$mean_location,
    settings$mean_sd,
    settings$sd_location,
    settings$sd_sd
  )
  list(
    training = normal_parameters_to_quantiles(training_parameters, probability_grid),
    calibration = normal_parameters_to_quantiles(calibration_parameters, probability_grid)
  )
}

make_wasserstein_search_grid <- function(mu_grid_n, sigma_grid_n, probability_grid) {
  parameters <- expand.grid(
    mu = seq(0, 5, length.out = mu_grid_n),
    sigma = seq(0, 1, length.out = sigma_grid_n + 1L)[-1L]
  )
  list(
    data = parameters,
    observations = normal_parameters_to_quantiles(parameters, probability_grid)
  )
}
