population_level_result <- function(scenario, config, level_name, level) {
  modality <- verify_scenario_modality(scenario)
  tuning_training <- deterministic_scenario_sample(
    level$tuning_reference_n, scenario, offset = 0.5
  )
  tuning_calibration <- deterministic_scenario_sample(
    level$tuning_reference_n + 1L, scenario, offset = 0.5
  )
  tuning_fit <- fit_profile_scores(
    tuning_training,
    config$tau_grid,
    config$chunk_size
  )
  tuning_profiles <- profiles_for_fit(tuning_calibration, tuning_fit)
  tau_selection <- select_tau_by_midpoints(
    tuning_fit,
    tuning_profiles,
    config$alpha,
    config$tau_penalty
  )

  reference <- deterministic_scenario_sample(
    level$reference_n, scenario, offset = 0.5
  )
  grid <- seq(
    level$grid_limits[[1L]],
    level$grid_limits[[2L]],
    by = level$grid_step
  )
  fit <- fit_profile_scores(reference, config$tau_grid, config$chunk_size)
  profiles <- profiles_for_fit(grid, fit)
  proposed_score <- stats::plogis(proposed_raw_scores_from_profiles(
    profiles, fit, tau_selection$tau
  )[, 1L])
  rank_score <- stats::plogis(transport_rank_raw_from_profiles(profiles, fit))
  cost_oriented <- transport_cost_centrality_from_profiles(profiles, fit)
  cdf <- scenario_cdf(grid, scenario)
  spatial_depth <- metric_spatial_depth_continuous_r1(cdf)
  halfspace_depth <- metric_halfspace_depth_continuous_r1(cdf)
  density <- scenario_density(grid, scenario)
  cdf_function <- function(point) scenario_cdf(point, scenario)

  curves <- list(
    "Proposed C-score" = list(score = proposed_score, tau = tau_selection$tau),
    "Transport rank" = list(score = rank_score, tau = 0),
    "Transport cost" = list(score = cost_oriented, tau = NA_real_),
    "Metric spatial depth" = list(score = spatial_depth, tau = NA_real_),
    "Metric halfspace (Tukey) depth" =
      list(score = halfspace_depth, tau = NA_real_)
  )
  target_points <- unique(c(
    0,
    modality$locations,
    if (is.finite(modality$antimode)) modality$antimode else numeric()
  ))
  target_profiles <- profiles_for_fit(target_points, fit)
  target_cdf <- scenario_cdf(target_points, scenario)
  target_scores <- list(
    "Proposed C-score" = proposed_scores_from_profiles(
      target_profiles, fit, tau_selection$tau
    )[, 1L],
    "Transport rank" =
      transport_rank_scores_from_profiles(target_profiles, fit),
    "Transport cost" =
      transport_cost_centrality_from_profiles(target_profiles, fit),
    "Metric spatial depth" =
      metric_spatial_depth_continuous_r1(target_cdf),
    "Metric halfspace (Tukey) depth" =
      metric_halfspace_depth_continuous_r1(target_cdf)
  )
  target_index <- function(value) match(value, target_points)
  oracle <- oracle_hdr_geometry(
    grid,
    density,
    probability = 1 - config$alpha,
    cdf_function = cdf_function,
    outside_margin_points = config$outside_margin_points
  )
  geometries <- lapply(curves, function(curve) {
    population_level_set_geometry(
      grid,
      curve$score,
      probability = 1 - config$alpha,
      cdf_function = cdf_function,
      outside_margin_points = config$outside_margin_points
    )
  })

  summary_rows <- lapply(names(curves), function(method) {
    curve <- curves[[method]]
    geometry <- geometries[[method]]
    metrics <- score_curve_metrics(
      grid, curve$score, density, modality$locations
    )
    mode_membership <- vapply(modality$locations, function(mode) {
      target_scores[[method]][target_index(mode)] >= geometry$threshold
    }, logical(1))
    center_membership <-
      target_scores[[method]][target_index(0)] >= geometry$threshold
    valley_excluded <- if (is.finite(modality$antimode)) {
      as.numeric(
        target_scores[[method]][target_index(modality$antimode)] <
          geometry$threshold
      )
    } else {
      NA_real_
    }
    data.frame(
      level = level_name,
      scenario = scenario$id[[1L]],
      scenario_label = scenario$label[[1L]],
      regime = scenario$regime[[1L]],
      method = method,
      tau = curve$tau,
      mode_count = modality$count,
      mode_locations = paste(
        formatC(modality$locations, digits = 8L, format = "fg"),
        collapse = ","
      ),
      density_spearman = metrics$density_spearman,
      center_loss = metrics$center_loss,
      mode_minus_center = metrics$mode_minus_center,
      optimum_to_mode = metrics$optimum_to_mode,
      optimum_location = metrics$optimum_location,
      population_score_threshold = geometry$threshold,
      population_probability_content = geometry$probability_content,
      population_total_length = geometry$total_length,
      population_components = geometry$components,
      population_center_included = as.numeric(center_membership),
      population_peak_inclusion_fraction = mean(mode_membership),
      population_all_modes_included = as.numeric(all(mode_membership)),
      population_valley_excluded = valley_excluded,
      oracle_density_cutoff = oracle$density_cutoff,
      oracle_hdr_length = oracle$total_length,
      oracle_hdr_components = oracle$components,
      reference_n = level$reference_n,
      grid_step = level$grid_step,
      stringsAsFactors = FALSE
    )
  })
  curve_rows <- do.call(rbind, lapply(names(curves), function(method) {
    data.frame(
      level = level_name,
      scenario = scenario$id[[1L]],
      method = method,
      tau = curves[[method]]$tau,
      candidate = grid,
      oriented_score = curves[[method]]$score,
      normalized_score = normalized_score(curves[[method]]$score),
      density = density,
      stringsAsFactors = FALSE
    )
  }))
  tau_rows <- data.frame(
    level = level_name,
    scenario = scenario$id[[1L]],
    tau = config$tau_grid,
    midpoint_proportion = tau_selection$proportion,
    criterion = tau_selection$criterion,
    selected = seq_along(config$tau_grid) == tau_selection$index,
    tuning_reference_n = level$tuning_reference_n,
    midpoint_count = tau_selection$midpoint_count,
    stringsAsFactors = FALSE
  )
  modality_row <- data.frame(
    level = level_name,
    scenario = scenario$id[[1L]],
    numerical_mode_count = modality$count,
    numerical_mode_locations = paste(
      formatC(modality$numerical_locations, digits = 8L, format = "fg"),
      collapse = ","
    ),
    true_mode_locations = paste(
      formatC(modality$locations, digits = 12L, format = "fg"),
      collapse = ","
    ),
    true_antimode = modality$antimode,
    stringsAsFactors = FALSE
  )
  spatial_geometry <- geometries[["Metric spatial depth"]]
  halfspace_geometry <- geometries[["Metric halfspace (Tukey) depth"]]
  depth_equivalence <- data.frame(
    level = level_name,
    scenario = scenario$id[[1L]],
    transform_max_abs = max(abs(
      spatial_depth - 4 * halfspace_depth * (1 - halfspace_depth)
    )),
    rank_correlation = suppressWarnings(stats::cor(
      spatial_depth, halfspace_depth, method = "spearman"
    )),
    region_symmetric_difference_length =
      interval_symmetric_difference_length(
        spatial_geometry$intervals,
        halfspace_geometry$intervals
      ),
    component_count_equal =
      spatial_geometry$components == halfspace_geometry$components,
    all_modes_membership_equal = identical(
      target_scores[["Metric spatial depth"]][
        vapply(modality$locations, target_index, integer(1))
      ] >= spatial_geometry$threshold,
      target_scores[["Metric halfspace (Tukey) depth"]][
        vapply(modality$locations, target_index, integer(1))
      ] >= halfspace_geometry$threshold
    ),
    stringsAsFactors = FALSE
  )
  list(
    summary = do.call(rbind, summary_rows),
    curves = curve_rows,
    tau = tau_rows,
    modality = modality_row,
    depth_equivalence = depth_equivalence
  )
}

scale_equivariance_check <- function(config, scale = 3) {
  scenario <- config$scenarios[
    match("normal", config$scenarios$id), , drop = FALSE
  ]
  if (!nrow(scenario)) stop("The scale check requires the normal scenario.")
  size <- min(401L, config$population$primary$reference_n)
  training <- deterministic_scenario_sample(size, scenario)
  query <- seq(-3, 3, length.out = 301L)
  fit <- fit_profile_scores(training, config$tau_grid, config$chunk_size)
  profiles <- profiles_for_fit(query, fit)
  raw <- proposed_raw_scores_from_profiles(profiles, fit)
  rank <- transport_rank_raw_from_profiles(profiles, fit)
  cost <- transport_cost_from_profiles(profiles, fit)
  spatial <- metric_spatial_depth(query, training, config$chunk_size)
  halfspace <- metric_halfspace_depth(query, training, config$chunk_size)

  scaled_fit <- fit_profile_scores(
    scale * training, config$tau_grid, config$chunk_size
  )
  scaled_profiles <- profiles_for_fit(scale * query, scaled_fit)
  scaled_raw <- proposed_raw_scores_from_profiles(
    scaled_profiles, scaled_fit
  )
  scaled_rank <- transport_rank_raw_from_profiles(
    scaled_profiles, scaled_fit
  )
  scaled_cost <- transport_cost_from_profiles(scaled_profiles, scaled_fit)
  scaled_spatial <- metric_spatial_depth(
    scale * query, scale * training, config$chunk_size
  )
  scaled_halfspace <- metric_halfspace_depth(
    scale * query, scale * training, config$chunk_size
  )

  minimum_correlation <- function(original, transformed) {
    original <- as.matrix(original)
    transformed <- as.matrix(transformed)
    min(vapply(seq_len(ncol(original)), function(column) {
      suppressWarnings(stats::cor(
        original[, column],
        transformed[, column],
        method = "spearman"
      ))
    }, numeric(1)))
  }
  correlations <- c(
    minimum_correlation(raw, scaled_raw),
    minimum_correlation(rank, scaled_rank),
    minimum_correlation(cost, scaled_cost),
    minimum_correlation(spatial, scaled_spatial),
    minimum_correlation(halfspace, scaled_halfspace)
  )
  homogeneity_errors <- c(
    max(abs(scaled_raw - scale * raw)),
    max(abs(scaled_rank - scale * rank)),
    max(abs(scaled_cost - scale * cost)),
    max(abs(scaled_spatial - spatial)),
    max(abs(scaled_halfspace - halfspace))
  )
  data.frame(
    check = c(
      "proposed_raw", "transport_rank_raw", "transport_cost",
      "metric_spatial_depth", "metric_halfspace_depth"
    ),
    scale = scale,
    expected_degree = c(1, 1, 1, 0, 0),
    max_homogeneity_error = homogeneity_errors,
    minimum_ordering_correlation = correlations,
    ordering_preserved = homogeneity_errors < 1e-10,
    stringsAsFactors = FALSE
  )
}

run_population_checks <- function(config, output_dir, workers = 1L) {
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  levels <- names(config$population)
  tasks <- do.call(c, lapply(levels, function(level_name) {
    lapply(seq_len(nrow(config$scenarios)), function(index) {
      list(level_name = level_name, scenario_index = index)
    })
  }))
  run_one <- function(task) {
    population_level_result(
      config$scenarios[task$scenario_index, , drop = FALSE],
      config,
      task$level_name,
      config$population[[task$level_name]]
    )
  }
  results <- run_parallel(tasks, run_one, workers)
  failures <- vapply(results, inherits, logical(1), what = "try-error")
  if (any(failures)) {
    stop("A population calculation failed: ",
         as.character(results[[which(failures)[[1L]]]]))
  }
  combine <- function(field) do.call(rbind, lapply(results, `[[`, field))
  summary <- combine("summary")
  curves <- combine("curves")
  tau <- combine("tau")
  modality <- combine("modality")
  depth_equivalence <- combine("depth_equivalence")
  scale <- scale_equivariance_check(config)
  assert_finite_frame(summary, "Population summary")
  assert_finite_frame(curves, "Population curves")
  assert_finite_frame(tau, "Population tau diagnostics")
  assert_finite_frame(depth_equivalence, "Population depth equivalence")
  assert_finite_frame(scale, "Scale-equivariance checks")
  write_csv_atomic(summary, file.path(output_dir, "population_summary.csv"))
  write_csv_atomic(curves, file.path(output_dir, "population_curves.csv"))
  write_csv_atomic(tau, file.path(output_dir, "population_tau_selection.csv"))
  write_csv_atomic(modality, file.path(output_dir, "modality_checks.csv"))
  write_csv_atomic(
    depth_equivalence,
    file.path(output_dir, "population_depth_equivalence.csv")
  )
  write_csv_atomic(scale, file.path(output_dir, "scale_equivariance.csv"))
  save_rds_atomic_comparison(
    list(
      config = config,
      config_signature = configuration_fingerprint(config),
      levels = names(config$population),
      expected_summary_rows = length(config$population) *
        nrow(config$scenarios) * length(comparison_method_levels),
      expected_curve_rows = nrow(curves),
      completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ),
    file.path(output_dir, "population_manifest.rds")
  )
  invisible(list(
    summary = summary,
    curves = curves,
    tau = tau,
    modality = modality,
    depth_equivalence = depth_equivalence,
    scale = scale
  ))
}
