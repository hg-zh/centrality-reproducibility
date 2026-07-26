scenario_row <- function(scenarios, id) {
  match_index <- match(id, scenarios$id)
  if (is.na(match_index)) stop("Unknown scenario: ", id)
  scenarios[match_index, , drop = FALSE]
}

scenario_cdf <- function(x, scenario) {
  family <- scenario$family[[1L]]
  parameter <- as.numeric(scenario$parameter[[1L]])
  if (identical(family, "normal")) {
    return(stats::pnorm(x))
  }
  if (identical(family, "laplace")) {
    scale <- 1 / sqrt(2)
    return(ifelse(
      x < 0,
      0.5 * exp(x / scale),
      1 - 0.5 * exp(-x / scale)
    ))
  }
  if (identical(family, "student_t")) {
    return(stats::pt(x * sqrt(3), df = 3))
  }
  if (identical(family, "normal_mixture")) {
    return(
      0.5 * stats::pnorm(x, mean = -parameter, sd = 1) +
        0.5 * stats::pnorm(x, mean = parameter, sd = 1)
    )
  }
  stop("Unknown scenario family: ", family)
}

scenario_density <- function(x, scenario) {
  family <- scenario$family[[1L]]
  parameter <- as.numeric(scenario$parameter[[1L]])
  if (identical(family, "normal")) {
    return(stats::dnorm(x))
  }
  if (identical(family, "laplace")) {
    scale <- 1 / sqrt(2)
    return(exp(-abs(x) / scale) / (2 * scale))
  }
  if (identical(family, "student_t")) {
    return(stats::dt(x * sqrt(3), df = 3) * sqrt(3))
  }
  if (identical(family, "normal_mixture")) {
    return(
      0.5 * stats::dnorm(x, mean = -parameter, sd = 1) +
        0.5 * stats::dnorm(x, mean = parameter, sd = 1)
    )
  }
  stop("Unknown scenario family: ", family)
}

scenario_quantile <- function(probability, scenario) {
  if (any(probability <= 0 | probability >= 1)) {
    stop("Deterministic quantile probabilities must lie in (0,1).")
  }
  family <- scenario$family[[1L]]
  parameter <- as.numeric(scenario$parameter[[1L]])
  if (identical(family, "normal")) {
    return(stats::qnorm(probability))
  }
  if (identical(family, "laplace")) {
    scale <- 1 / sqrt(2)
    return(ifelse(
      probability < 0.5,
      scale * log(2 * probability),
      -scale * log(2 * (1 - probability))
    ))
  }
  if (identical(family, "student_t")) {
    return(stats::qt(probability, df = 3) / sqrt(3))
  }
  if (identical(family, "normal_mixture")) {
    lower <- -parameter - 12
    upper <- parameter + 12
    return(vapply(
      probability,
      function(value) stats::uniroot(
        function(point) scenario_cdf(point, scenario) - value,
        interval = c(lower, upper),
        tol = 1e-12
      )$root,
      numeric(1)
    ))
  }
  stop("Unknown scenario family: ", family)
}

draw_scenario <- function(size, scenario) {
  family <- scenario$family[[1L]]
  parameter <- as.numeric(scenario$parameter[[1L]])
  if (identical(family, "normal")) {
    return(stats::rnorm(size))
  }
  if (identical(family, "laplace")) {
    return(scenario_quantile(stats::runif(size), scenario))
  }
  if (identical(family, "student_t")) {
    return(stats::rt(size, df = 3) / sqrt(3))
  }
  if (identical(family, "normal_mixture")) {
    component <- ifelse(stats::runif(size) < 0.5, -parameter, parameter)
    return(stats::rnorm(size, mean = component, sd = 1))
  }
  stop("Unknown scenario family: ", family)
}

deterministic_scenario_sample <- function(size, scenario, offset = 0.5) {
  probability <- (seq_len(size) - offset) / size
  if (any(probability <= 0 | probability >= 1)) {
    stop("The deterministic probability offset is invalid.")
  }
  scenario_quantile(probability, scenario)
}

scenario_numerical_modes <- function(scenario, grid = seq(-8, 8, by = 0.001)) {
  density <- scenario_density(grid, scenario)
  local <- which(
    density >= c(-Inf, density[-length(density)]) &
      density >= c(density[-1L], -Inf)
  )
  local <- local[grid[local] > min(grid) & grid[local] < max(grid)]
  if (!length(local)) stop("No interior density mode was detected.")
  peak <- max(density)
  local <- local[density[local] >= peak * (1 - 1e-7)]
  positions <- grid[local]
  if (length(positions) > 1L) {
    groups <- cumsum(c(TRUE, diff(local) > 1L))
    positions <- vapply(
      split(positions, groups),
      mean,
      numeric(1)
    )
  }
  as.numeric(positions)
}

scenario_true_modes <- function(scenario) {
  family <- scenario$family[[1L]]
  parameter <- as.numeric(scenario$parameter[[1L]])
  if (!identical(family, "normal_mixture") || parameter <= 1) {
    return(0)
  }
  positive <- stats::optimize(
    function(point) -scenario_density(point, scenario),
    interval = c(.Machine$double.eps^0.25, parameter + 4),
    tol = 1e-14
  )$minimum
  c(-positive, positive)
}

scenario_true_antimode <- function(scenario) {
  family <- scenario$family[[1L]]
  parameter <- as.numeric(scenario$parameter[[1L]])
  if (identical(family, "normal_mixture") && parameter > 1) 0 else NA_real_
}

scenario_geometry <- function(scenario) {
  modes <- scenario_true_modes(scenario)
  list(
    mode_count = length(modes),
    modes = modes,
    antimode = scenario_true_antimode(scenario),
    median = 0
  )
}

scenario_modes <- function(scenario, grid = NULL) {
  scenario_true_modes(scenario)
}

verify_scenario_modality <- function(scenario) {
  numerical_modes <- scenario_numerical_modes(scenario)
  geometry <- scenario_geometry(scenario)
  expected <- if (
    identical(scenario$family[[1L]], "normal_mixture") &&
      as.numeric(scenario$parameter[[1L]]) > 1
  ) 2L else 1L
  if (length(numerical_modes) != expected ||
      geometry$mode_count != expected ||
      max(abs(numerical_modes - geometry$modes)) > 0.002) {
    stop(
      "Numerical modality check failed for ", scenario$id[[1L]],
      ": expected ", expected, " mode(s), found ",
      length(numerical_modes), "."
    )
  }
  list(
    count = geometry$mode_count,
    locations = geometry$modes,
    antimode = geometry$antimode,
    median = geometry$median,
    numerical_locations = numerical_modes
  )
}
