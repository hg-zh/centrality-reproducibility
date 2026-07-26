parse_release_arguments <- function(
    args = commandArgs(trailingOnly = TRUE),
    default_profile = "manuscript") {
  result <- list(
    profile = default_profile,
    workers = 1L,
    force = FALSE
  )
  for (argument in args) {
    if (grepl("^--profile=", argument)) {
      result$profile <- sub("^--profile=", "", argument)
    } else if (grepl("^--workers=", argument)) {
      result$workers <- as.integer(sub("^--workers=", "", argument))
    } else if (identical(argument, "--force")) {
      result$force <- TRUE
    } else {
      stop("Unknown command-line argument: ", argument)
    }
  }
  if (!result$profile %in% c("manuscript", "smoke")) {
    stop("--profile must be 'manuscript' or 'smoke'.")
  }
  if (is.na(result$workers) || result$workers < 1L) {
    stop("--workers must be a positive integer.")
  }
  result
}

release_output_directory <- function(project_root, kind, profile) {
  path <- file.path(project_root, kind, profile)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

save_rds_atomic <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(basename(path), ".temporary-"),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) {
    stop("Could not move completed result into place: ", path)
  }
  invisible(path)
}

metric_rows <- function(replicate, sample_size, tau, evaluation) {
  data.frame(
    replicate = as.integer(replicate),
    sample_size = as.integer(sample_size),
    tau = as.numeric(tau),
    coverage = as.numeric(evaluation$coverage),
    size = as.numeric(evaluation$size)
  )
}

analysis_settings <- function(settings) {
  list(
    alpha = settings$alpha,
    threshold = settings$threshold
  )
}

simulation_cache_configuration <- function(settings, scenario) {
  list(
    profile = settings$profile,
    scenario = scenario,
    settings = settings[[scenario]],
    analysis = analysis_settings(settings),
    rng = settings$rng
  )
}

simulate_euclidean_replicate <- function(replicate, settings) {
  cfg <- settings$euclidean
  set_reproducible_rng(replicate)
  probability_grid <- make_probability_grid(
    cfg$probability_grid_size,
    include_endpoints = TRUE
  )
  search <- make_euclidean_search_grid(cfg$grid_n, cfg$grid_limits)
  testing <- generate_gaussian_mixture(
    cfg$n_test, cfg$mixing, cfg$means, cfg$covariances
  )

  rows <- vector("list", length(cfg$sample_sizes))
  for (index in seq_along(cfg$sample_sizes)) {
    n <- cfg$sample_sizes[index]
    pair <- generate_gaussian_mixture_pair(
      n, cfg$mixing, cfg$means, cfg$covariances
    )
    fit <- fit_centrality_score(
      pair$training,
      "euclidean",
      cfg$distance_grid,
      probability_grid
    )
    evaluation <- evaluate_centrality_region(
      fit = fit,
      calibration = pair$calibration,
      testing = testing,
      search_grid = search$observations,
      tau = cfg$tau,
      alpha = settings$alpha,
      size_scale = cfg$size_scale,
      threshold_method = settings$threshold$method,
      quantile_type = settings$threshold$quantile_type
    )
    rows[[index]] <- metric_rows(replicate, n, cfg$tau, evaluation)
  }

  heat <- vector("list", length(cfg$heat_mixing))
  set_reproducible_rng(replicate)
  for (mixing in cfg$heat_generation_mixing) {
    heat_probability_grid <- make_probability_grid(
      cfg$probability_grid_size,
      include_endpoints = TRUE
    )
    invisible(generate_gaussian_mixture(
      cfg$n_test, mixing, cfg$means, cfg$covariances
    ))
    pairs <- lapply(cfg$sample_sizes, function(n) {
      generate_gaussian_mixture_pair(n, mixing, cfg$means, cfg$covariances)
    })
    display_index <- match(mixing, cfg$heat_mixing)
    if (is.na(display_index)) next
    pair <- pairs[[1L]]
    fit <- fit_centrality_score(
      pair$training,
      "euclidean",
      cfg$distance_grid,
      heat_probability_grid
    )
    calibration_q <- distance_quantile_profiles(
      pair$calibration,
      fit$training,
      fit$metric,
      fit$distance_grid,
      fit$probability_grid
    )
    search_q <- distance_quantile_profiles(
      search$observations,
      fit$training,
      fit$metric,
      fit$distance_grid,
      fit$probability_grid
    )
    calibration_scores <- centrality_scores_from_quantiles(
      calibration_q, fit, cfg$heat_tau
    )
    search_scores <- centrality_scores_from_quantiles(
      search_q, fit, cfg$heat_tau
    )
    thresholds <- calibration_thresholds(
      calibration_scores,
      settings$alpha,
      nrow(fit$training),
      settings$threshold$method,
      settings$threshold$quantile_type
    )
    heat[[display_index]] <- list(
      mixing = mixing,
      tau = cfg$heat_tau,
      membership = sweep(search_scores, 2L, thresholds, ">=")
    )
  }

  list(metrics = do.call(rbind, rows), heat = heat)
}

simulate_sphere_replicate <- function(replicate, settings) {
  cfg <- settings$sphere
  set_reproducible_rng(replicate)
  probability_grid <- make_probability_grid(
    cfg$probability_grid_size,
    include_endpoints = TRUE
  )
  search <- make_sphere_search_grid(cfg$angular_grid_n)
  testing <- generate_sphere_mixture(
    cfg$n_test, cfg$mixing, cfg$noise_sd
  )

  rows <- vector("list", length(cfg$sample_sizes))
  heat_membership <- NULL
  for (index in seq_along(cfg$sample_sizes)) {
    n <- cfg$sample_sizes[index]
    pair <- generate_sphere_mixture_pair(n, cfg$mixing, cfg$noise_sd)
    fit <- fit_centrality_score(
      pair$training,
      "sphere",
      cfg$distance_grid,
      probability_grid
    )
    evaluation <- evaluate_centrality_region(
      fit = fit,
      calibration = pair$calibration,
      testing = testing,
      search_grid = search$observations,
      tau = cfg$tau,
      alpha = settings$alpha,
      size_scale = cfg$size_scale,
      threshold_method = settings$threshold$method,
      quantile_type = settings$threshold$quantile_type
    )
    rows[[index]] <- metric_rows(replicate, n, cfg$tau, evaluation)
    if (index == 1L) {
      heat_column <- match(cfg$heat_tau, cfg$tau)
      if (is.na(heat_column)) stop("Sphere heat-map tau is absent from the tau grid.")
      heat_membership <- evaluation$search_membership[, heat_column]
    }
  }

  list(
    metrics = do.call(rbind, rows),
    heat = list(tau = cfg$heat_tau, membership = heat_membership)
  )
}

simulate_wasserstein_replicate <- function(replicate, settings) {
  cfg <- settings$wasserstein
  set_reproducible_rng(replicate)
  probability_grid <- make_probability_grid(
    cfg$probability_grid_size,
    include_endpoints = FALSE
  )
  search <- make_wasserstein_search_grid(
    cfg$mu_grid_n,
    cfg$sigma_grid_n,
    probability_grid
  )
  testing_parameters <- generate_wasserstein_parameters(
    cfg$n_test,
    cfg$mean_location,
    cfg$mean_sd,
    cfg$sd_location,
    cfg$sd_sd
  )
  testing <- normal_parameters_to_quantiles(testing_parameters, probability_grid)

  rows <- vector("list", length(cfg$sample_sizes))
  heat_membership <- NULL
  for (index in seq_along(cfg$sample_sizes)) {
    n <- cfg$sample_sizes[index]
    pair <- generate_wasserstein_pair(n, cfg, probability_grid)
    fit <- fit_centrality_score(
      pair$training,
      "wasserstein",
      cfg$distance_grid,
      probability_grid
    )
    evaluation <- evaluate_centrality_region(
      fit = fit,
      calibration = pair$calibration,
      testing = testing,
      search_grid = search$observations,
      tau = cfg$tau,
      alpha = settings$alpha,
      size_scale = cfg$size_scale,
      threshold_method = settings$threshold$method,
      quantile_type = settings$threshold$quantile_type
    )
    rows[[index]] <- metric_rows(replicate, n, cfg$tau, evaluation)
    if (index == 1L) {
      heat_column <- match(cfg$heat_tau, cfg$tau)
      if (is.na(heat_column)) {
        stop("Wasserstein heat-map tau is absent from the tau grid.")
      }
      heat_membership <- evaluation$search_membership[, heat_column]
    }
  }

  list(
    metrics = do.call(rbind, rows),
    heat = list(tau = cfg$heat_tau, membership = heat_membership)
  )
}

run_cached_replicates <- function(
    scenario,
    settings,
    project_root,
    runner,
    workers = 1L,
    force = FALSE) {
  scenario_settings <- settings[[scenario]]
  cache_configuration <- simulation_cache_configuration(settings, scenario)
  cache_directory <- file.path(
    release_output_directory(project_root, "results", settings$profile),
    "cache",
    scenario
  )
  dir.create(cache_directory, recursive = TRUE, showWarnings = FALSE)
  ids <- seq_len(scenario_settings$replicates)

  run_one <- function(replicate) {
    path <- file.path(
      cache_directory,
      sprintf("replicate_%04d.rds", replicate)
    )
    if (!force && file.exists(path)) {
      cached <- readRDS(path)
      if (
        identical(cached$schema_version, 2L) &&
        identical(cached$scenario, scenario) &&
        identical(cached$profile, settings$profile) &&
        identical(cached$replicate, as.integer(replicate)) &&
        identical(cached$configuration, cache_configuration)
      ) {
        return(cached$result)
      }
    }
    result <- runner(replicate, settings)
    save_rds_atomic(
      list(
        schema_version = 2L,
        scenario = scenario,
        profile = settings$profile,
        replicate = as.integer(replicate),
        configuration = cache_configuration,
        result = result
      ),
      path
    )
    result
  }

  workers <- min(as.integer(workers), length(ids))
  if (workers > 1L && .Platform$OS.type == "unix") {
    results <- parallel::mclapply(
      ids,
      run_one,
      mc.cores = workers,
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
  } else {
    if (workers > 1L) {
      message("Parallel workers are unavailable on this platform; using one worker.")
    }
    results <- lapply(ids, run_one)
  }
  results
}

summarize_metric_rows <- function(results) {
  raw <- do.call(rbind, lapply(results, `[[`, "metrics"))
  summary <- stats::aggregate(
    cbind(coverage, size) ~ sample_size + tau,
    data = raw,
    FUN = mean
  )
  list(raw = raw, summary = summary)
}

aggregate_euclidean_results <- function(results, settings) {
  cfg <- settings$euclidean
  metrics <- summarize_metric_rows(results)
  search <- make_euclidean_search_grid(cfg$grid_n, cfg$grid_limits)
  heat_rows <- list()
  row_id <- 1L
  for (mix_index in seq_along(cfg$heat_mixing)) {
    for (tau_index in seq_along(cfg$heat_tau)) {
      count <- numeric(nrow(search$data))
      for (result in results) {
        count <- count + result$heat[[mix_index]]$membership[, tau_index]
      }
      heat_rows[[row_id]] <- data.frame(
        search$data,
        mixing = cfg$heat_mixing[mix_index],
        tau = cfg$heat_tau[tau_index],
        probability = count / length(results)
      )
      row_id <- row_id + 1L
    }
  }
  list(
    schema_version = 2L,
    profile = settings$profile,
    scenario = "euclidean",
    settings = cfg,
    analysis = analysis_settings(settings),
    rng = settings$rng,
    metric_replicates = metrics$raw,
    metric_summary = metrics$summary,
    heatmap = do.call(rbind, heat_rows)
  )
}

aggregate_sphere_results <- function(results, settings) {
  cfg <- settings$sphere
  metrics <- summarize_metric_rows(results)
  search <- make_sphere_search_grid(cfg$angular_grid_n)
  count <- Reduce(
    "+",
    lapply(results, function(result) as.numeric(result$heat$membership))
  )
  heatmap <- cbind(
    search$data[, c("phi", "theta")],
    probability = count / length(results),
    tau = cfg$heat_tau
  )
  list(
    schema_version = 2L,
    profile = settings$profile,
    scenario = "sphere",
    settings = cfg,
    analysis = analysis_settings(settings),
    rng = settings$rng,
    metric_replicates = metrics$raw,
    metric_summary = metrics$summary,
    heatmap = heatmap
  )
}

aggregate_wasserstein_results <- function(results, settings) {
  cfg <- settings$wasserstein
  metrics <- summarize_metric_rows(results)
  probability_grid <- seq(0.01, 0.99, length.out = cfg$probability_grid_size)
  search <- make_wasserstein_search_grid(
    cfg$mu_grid_n,
    cfg$sigma_grid_n,
    probability_grid
  )
  count <- Reduce(
    "+",
    lapply(results, function(result) as.numeric(result$heat$membership))
  )
  heatmap <- cbind(
    search$data,
    probability = count / length(results),
    tau = cfg$heat_tau
  )
  list(
    schema_version = 2L,
    profile = settings$profile,
    scenario = "wasserstein",
    settings = cfg,
    analysis = analysis_settings(settings),
    rng = settings$rng,
    metric_replicates = metrics$raw,
    metric_summary = metrics$summary,
    heatmap = heatmap
  )
}

run_all_simulations <- function(
    settings,
    project_root,
    workers = 1L,
    force = FALSE) {
  output_directory <- release_output_directory(
    project_root, "results", settings$profile
  )

  message("Running Euclidean simulations (", settings$profile, ").")
  euclidean_runs <- run_cached_replicates(
    "euclidean", settings, project_root,
    simulate_euclidean_replicate, workers, force
  )
  euclidean <- aggregate_euclidean_results(euclidean_runs, settings)
  save_rds_atomic(
    euclidean,
    file.path(output_directory, "figure3_euclidean.rds")
  )

  message("Running sphere simulations (", settings$profile, ").")
  sphere_runs <- run_cached_replicates(
    "sphere", settings, project_root,
    simulate_sphere_replicate, workers, force
  )
  sphere <- aggregate_sphere_results(sphere_runs, settings)
  save_rds_atomic(
    sphere,
    file.path(output_directory, "figure4_sphere.rds")
  )

  message("Running Wasserstein simulations (", settings$profile, ").")
  wasserstein_runs <- run_cached_replicates(
    "wasserstein", settings, project_root,
    simulate_wasserstein_replicate, workers, force
  )
  wasserstein <- aggregate_wasserstein_results(wasserstein_runs, settings)
  save_rds_atomic(
    wasserstein,
    file.path(output_directory, "figure4_wasserstein.rds")
  )

  manifest <- list(
    schema_version = 2L,
    profile = settings$profile,
    result_files = c(
      "figure3_euclidean.rds",
      "figure4_sphere.rds",
      "figure4_wasserstein.rds"
    ),
    completed_replicates = c(
      euclidean = length(euclidean_runs),
      sphere = length(sphere_runs),
      wasserstein = length(wasserstein_runs)
    ),
    scenario_settings = list(
      euclidean = settings$euclidean,
      sphere = settings$sphere,
      wasserstein = settings$wasserstein
    ),
    analysis = analysis_settings(settings),
    rng = settings$rng
  )
  save_rds_atomic(manifest, file.path(output_directory, "manifest.rds"))
  invisible(manifest)
}
