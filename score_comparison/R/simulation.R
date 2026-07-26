comparison_score_matrix_from_profiles <- function(
    query,
    profiles,
    fit,
    selected_tau) {
  query <- as.numeric(query)
  profiles <- as.matrix(profiles)
  if (nrow(profiles) != length(query)) {
    stop("Query points and profile rows differ.")
  }
  answer <- cbind(
    proposed_scores_from_profiles(profiles, fit, selected_tau)[, 1L],
    transport_rank_scores_from_profiles(profiles, fit),
    transport_cost_centrality_from_profiles(profiles, fit),
    metric_spatial_depth(query, fit$training, fit$chunk_size),
    metric_halfspace_depth(query, fit$training, fit$chunk_size)
  )
  colnames(answer) <- comparison_method_levels
  answer
}

comparison_implementation_fingerprint <- function() {
  function_names <- c(
    "simulate_comparison_replicate", "comparison_score_matrix_from_profiles",
    "set_reproducible_rng", "comparison_seed", "numeric_fingerprint",
    "sorted_distance_profiles",
    "profile_integral_weights", "fit_profile_scores",
    "proposed_raw_scores_from_profiles", "transport_rank_raw_from_profiles",
    "mean_absolute_profile_cost", "metric_spatial_depth",
    "metric_spatial_depth_continuous_r1", "metric_halfspace_depth",
    "metric_halfspace_depth_continuous_r1", "manuscript_threshold",
    "select_tau_by_midpoints", "ensure_score_domain",
    "score_level_intervals", "empirical_depth_cells_r1",
    "empirical_depth_intervals_r1",
    "empirical_depth_score_from_counts", "merge_closed_intervals",
    "population_level_set_geometry",
    "oracle_hdr_geometry", "total_interval_length",
    "interval_symmetric_difference_length", "prediction_set_metrics",
    "score_curve_metrics", "draw_scenario", "scenario_cdf",
    "scenario_density", "scenario_true_modes",
    "scenario_true_antimode", "verify_scenario_modality"
  )
  bodies <- lapply(function_names, function(name) {
    paste(deparse(body(get(name, mode = "function"))), collapse = "\n")
  })
  names(bodies) <- function_names
  configuration_fingerprint(list(
    r_major_minor = paste(R.version$major, R.version$minor, sep = "."),
    function_bodies = bodies
  ))
}

simulate_comparison_replicate <- function(
    scenario,
    scenario_index,
    replicate,
    config) {
  seed <- comparison_seed(config, scenario_index, replicate)
  set_reproducible_rng(seed)
  training <- draw_scenario(config$n_train, scenario)
  calibration <- draw_scenario(config$n_calibration, scenario)
  testing <- draw_scenario(config$n_test, scenario)
  fingerprint <- numeric_fingerprint(training, calibration, testing)

  fit <- fit_profile_scores(training, config$tau_grid, config$chunk_size)
  calibration_profiles <- profiles_for_fit(calibration, fit)
  testing_profiles <- profiles_for_fit(testing, fit)
  tau_selection <- select_tau_by_midpoints(
    fit,
    calibration_profiles,
    alpha = config$alpha,
    penalty = config$tau_penalty
  )
  selected_tau <- tau_selection$tau

  calibration_scores <- comparison_score_matrix_from_profiles(
    calibration, calibration_profiles, fit, selected_tau
  )
  testing_scores <- comparison_score_matrix_from_profiles(
    testing, testing_profiles, fit, selected_tau
  )
  thresholds <- calibration_threshold_vector(
    calibration_scores,
    alpha = config$alpha,
    n_training = config$n_train
  )
  testing_membership <- sweep(testing_scores, 2L, thresholds, ">=")
  depth_methods <- c(
    "Metric spatial depth",
    "Metric halfspace (Tukey) depth"
  )
  depth_method_keys <- c(
    "Metric spatial depth" = "spatial",
    "Metric halfspace (Tukey) depth" = "halfspace"
  )
  exact_depth_intervals <- lapply(depth_methods, function(method) {
    empirical_depth_intervals_r1(
      reference = training,
      threshold = thresholds[[method]],
      method = depth_method_keys[[method]]
    )
  })
  names(exact_depth_intervals) <- depth_methods
  unbounded_methods <- depth_methods[vapply(
    exact_depth_intervals,
    function(intervals) is.infinite(total_interval_length(intervals)),
    logical(1)
  )]

  score_function <- function(points) {
    point_profiles <- profiles_for_fit(points, fit)
    comparison_score_matrix_from_profiles(
      points, point_profiles, fit, selected_tau
    )
  }
  domain <- ensure_score_domain(
    initial_grid = config$candidate_grid,
    score_function = score_function,
    thresholds = thresholds,
    outside_margin_points = config$outside_margin_points,
    expansion_fraction = config$domain_expansion_fraction,
    max_expansions = config$domain_max_expansions,
    unbounded_methods = unbounded_methods
  )
  grid <- domain$grid
  grid_scores <- domain$scores
  density <- scenario_density(grid, scenario)
  cdf_function <- function(point) scenario_cdf(point, scenario)
  oracle <- oracle_hdr_geometry(
    grid = grid,
    density = density,
    probability = 1 - config$alpha,
    cdf_function = cdf_function,
    outside_margin_points = config$outside_margin_points
  )
  modality <- verify_scenario_modality(scenario)

  target_points <- unique(c(
    0,
    modality$locations,
    if (is.finite(modality$antimode)) modality$antimode else numeric()
  ))
  target_scores <- score_function(target_points)
  target_index <- function(value) match(value, target_points)

  grid_depth_intervals <- lapply(depth_methods, function(method) {
    if (method %in% unbounded_methods) {
      return(data.frame(left = -Inf, right = Inf))
    }
    column <- match(method, comparison_method_levels)
    score_level_intervals(
      grid,
      grid_scores[, column],
      thresholds[[method]],
      outside_margin_points = config$outside_margin_points,
      enforce_margin = TRUE
    )
  })
  names(grid_depth_intervals) <- depth_methods

  intervals <- vector("list", length(comparison_method_levels))
  names(intervals) <- comparison_method_levels
  metric_rows <- lapply(seq_along(comparison_method_levels), function(column) {
    method <- comparison_method_levels[[column]]
    intervals[[method]] <<- if (method %in% depth_methods) {
      exact_depth_intervals[[method]]
    } else {
      score_level_intervals(
        grid,
        grid_scores[, column],
        thresholds[[column]],
        outside_margin_points = config$outside_margin_points,
        enforce_margin = TRUE
      )
    }
    mode_membership <- vapply(modality$locations, function(mode) {
      target_scores[target_index(mode), column] >= thresholds[[column]]
    }, logical(1))
    center_membership <-
      target_scores[target_index(0), column] >= thresholds[[column]]
    antimode_membership <- if (is.finite(modality$antimode)) {
      target_scores[target_index(modality$antimode), column] >=
        thresholds[[column]]
    } else {
      NA
    }
    prediction <- prediction_set_metrics(
      intervals = intervals[[method]],
      test_membership = testing_membership[, column],
      oracle_intervals = oracle$intervals,
      mode_membership = mode_membership,
      center_membership = center_membership,
      antimode_membership = antimode_membership
    )
    curve <- score_curve_metrics(
      grid,
      grid_scores[, column],
      density,
      modality$locations
    )
    method_tau <- if (identical(method, "Proposed C-score")) {
      selected_tau
    } else if (identical(method, "Transport rank")) {
      0
    } else {
      NA_real_
    }
    data.frame(
      scenario = scenario$id[[1L]],
      scenario_label = scenario$label[[1L]],
      regime = scenario$regime[[1L]],
      method = method,
      replicate = as.integer(replicate),
      seed = seed,
      dataset_fingerprint = fingerprint,
      tau = method_tau,
      score_threshold = thresholds[[column]],
      coverage = prediction$coverage,
      total_set_length = prediction$total_set_length,
      finite_replication_total_length =
        prediction$finite_replication_total_length,
      set_unbounded = prediction$set_unbounded,
      components = prediction$components,
      hdr_symmetric_difference_length =
        prediction$hdr_symmetric_difference_length,
      finite_replication_hdr_symmetric_difference_length =
        prediction$finite_replication_hdr_symmetric_difference_length,
      oracle_hdr_length = prediction$oracle_hdr_length,
      excess_length = prediction$excess_length,
      center_included = prediction$center_included,
      valley_excluded = prediction$valley_excluded,
      all_modes_included = prediction$all_modes_included,
      peak_inclusion_fraction = prediction$peak_inclusion_fraction,
      density_spearman = curve$density_spearman,
      center_loss = curve$center_loss,
      mode_minus_center = curve$mode_minus_center,
      stringsAsFactors = FALSE
    )
  })
  metrics <- do.call(rbind, metric_rows)

  initial_boundary_scores <- score_function(range(config$candidate_grid))
  depth_diagnostics <- do.call(rbind, lapply(depth_methods, function(method) {
    column <- match(method, comparison_method_levels)
    exact_length <- total_interval_length(exact_depth_intervals[[method]])
    grid_length <- total_interval_length(grid_depth_intervals[[method]])
    length_difference <- if (is.infinite(exact_length) &&
        is.infinite(grid_length)) {
      0
    } else {
      abs(exact_length - grid_length)
    }
    data.frame(
      scenario = scenario$id[[1L]],
      scenario_label = scenario$label[[1L]],
      replicate = as.integer(replicate),
      seed = seed,
      method = method,
      n_training = length(training),
      n_calibration = length(calibration),
      alpha = config$alpha,
      nonconformity_order_index =
        ceiling((length(calibration) + 1) * (1 - config$alpha)),
      conformity_order_index =
        floor((length(calibration) + 1) * config$alpha),
      score_orientation = "higher_is_more_conforming",
      acceptance_inequality = "score >= threshold",
      calibration_threshold = thresholds[[method]],
      calibration_minimum = min(calibration_scores[, column]),
      calibration_maximum = max(calibration_scores[, column]),
      calibration_tail_score_count =
        sum(calibration_scores[, column] == 0),
      depth_tail_limit_left = 0,
      depth_tail_limit_right = 0,
      negative_depth_tail_limit_left = 0,
      negative_depth_tail_limit_right = 0,
      initial_left = min(config$candidate_grid),
      initial_right = max(config$candidate_grid),
      initial_left_score = initial_boundary_scores[1L, column],
      initial_right_score = initial_boundary_scores[2L, column],
      evaluated_left = min(grid),
      evaluated_right = max(grid),
      evaluated_left_score = grid_scores[1L, column],
      evaluated_right_score = grid_scores[nrow(grid_scores), column],
      expansions = domain$expansions,
      exact_unbounded = as.numeric(is.infinite(exact_length)),
      exact_total_length = exact_length,
      exact_components = nrow(exact_depth_intervals[[method]]),
      grid_total_length = grid_length,
      grid_components = nrow(grid_depth_intervals[[method]]),
      grid_exact_absolute_length_difference = length_difference,
      grid_exact_component_equal =
        nrow(grid_depth_intervals[[method]]) ==
          nrow(exact_depth_intervals[[method]]),
      stringsAsFactors = FALSE
    )
  }))

  grid_profiles <- profiles_for_fit(grid, fit)
  c0_scores <- proposed_scores_from_profiles(grid_profiles, fit, tau = 0)[, 1L]
  c0_raw <- proposed_raw_scores_from_profiles(grid_profiles, fit, tau = 0)[, 1L]
  published_raw <- published_transport_rank_raw(grid, training)
  published_shift <- published_raw - c0_raw
  spatial_column <- match("Metric spatial depth", comparison_method_levels)
  halfspace_column <- match(
    "Metric halfspace (Tukey) depth", comparison_method_levels
  )
  spatial_grid <- grid_scores[, spatial_column]
  halfspace_grid <- grid_scores[, halfspace_column]
  grid_membership <- sweep(grid_scores, 2L, thresholds, ">=")
  spatial_exact_length <- total_interval_length(
    exact_depth_intervals[["Metric spatial depth"]]
  )
  halfspace_exact_length <- total_interval_length(
    exact_depth_intervals[["Metric halfspace (Tukey) depth"]]
  )
  spatial_halfspace_exact_length_difference <- if (
      is.infinite(spatial_exact_length) &&
        is.infinite(halfspace_exact_length)) {
    0
  } else {
    abs(spatial_exact_length - halfspace_exact_length)
  }
  checks <- data.frame(
    scenario = scenario$id[[1L]],
    replicate = as.integer(replicate),
    seed = seed,
    dataset_fingerprint = fingerprint,
    c0_rank_max_abs = max(abs(
      c0_scores -
        grid_scores[, match("Transport rank", comparison_method_levels)]
    )),
    published_raw_shift_range = diff(range(published_shift)),
    published_rank_correlation = suppressWarnings(stats::cor(
      published_raw, c0_raw, method = "spearman"
    )),
    spatial_halfspace_transform_max_abs = max(abs(
      spatial_grid - 4 * halfspace_grid * (1 - halfspace_grid)
    )),
    spatial_halfspace_rank_correlation = suppressWarnings(stats::cor(
      spatial_grid, halfspace_grid, method = "spearman"
    )),
    spatial_halfspace_grid_set_equal = identical(
      grid_membership[, spatial_column],
      grid_membership[, halfspace_column]
    ),
    spatial_halfspace_test_set_equal = identical(
      testing_membership[, spatial_column],
      testing_membership[, halfspace_column]
    ),
    spatial_halfspace_unbounded_equal =
      is.infinite(spatial_exact_length) ==
        is.infinite(halfspace_exact_length),
    spatial_halfspace_component_equal =
      nrow(exact_depth_intervals[["Metric spatial depth"]]) ==
        nrow(exact_depth_intervals[["Metric halfspace (Tukey) depth"]]),
    spatial_halfspace_exact_length_difference =
      spatial_halfspace_exact_length_difference,
    maximum_depth_grid_exact_length_difference = max(
      depth_diagnostics$grid_exact_absolute_length_difference
    ),
    all_depth_grid_exact_components_equal =
      all(depth_diagnostics$grid_exact_component_equal),
    n_training = length(training),
    n_calibration = length(calibration),
    n_test = length(testing),
    stringsAsFactors = FALSE
  )

  tau_rows <- data.frame(
    scenario = scenario$id[[1L]],
    replicate = as.integer(replicate),
    seed = seed,
    tau = fit$tau_grid,
    midpoint_proportion = tau_selection$proportion,
    criterion = tau_selection$criterion,
    selected = seq_along(fit$tau_grid) == tau_selection$index,
    midpoint_count = tau_selection$midpoint_count,
    stringsAsFactors = FALSE
  )

  domain_row <- data.frame(
    scenario = scenario$id[[1L]],
    replicate = as.integer(replicate),
    seed = seed,
    initial_left = min(config$candidate_grid),
    initial_right = max(config$candidate_grid),
    evaluated_left = min(grid),
    evaluated_right = max(grid),
    grid_step = validate_regular_grid(grid),
    expansions = domain$expansions,
    minimum_left_outside_points = domain$left_outside_points,
    minimum_right_outside_points = domain$right_outside_points,
    oracle_density_cutoff = oracle$density_cutoff,
    oracle_probability_content = oracle$probability_content,
    oracle_total_length = oracle$total_length,
    oracle_components = oracle$components,
    unbounded_method_count = length(unbounded_methods),
    unbounded_methods = paste(unbounded_methods, collapse = ";"),
    stringsAsFactors = FALSE
  )

  refinement <- data.frame()
  if (replicate %in% config$refinement_replicates) {
    factor <- as.integer(config$refinement_factor)
    fine_grid <- seq(
      from = min(grid),
      to = max(grid),
      length.out = (length(grid) - 1L) * factor + 1L
    )
    fine_scores <- score_function(fine_grid)
    coarse_geometry <- list(
      grid = grid,
      lengths = vapply(intervals, total_interval_length, numeric(1)),
      components = vapply(intervals, nrow, integer(1))
    )
    refinement <- refined_geometry_comparison(
      coarse_geometry = coarse_geometry,
      fine_grid = fine_grid,
      fine_scores = fine_scores,
      thresholds = thresholds,
      outside_margin_points = config$outside_margin_points * factor,
      unbounded_methods = unbounded_methods,
      exact_intervals = exact_depth_intervals
    )
    refinement$scenario <- scenario$id[[1L]]
    refinement$replicate <- as.integer(replicate)
    refinement$seed <- seed
    refinement$dataset_fingerprint <- fingerprint
    refinement <- refinement[c(
      "scenario", "replicate", "seed", "dataset_fingerprint",
      setdiff(names(refinement), c(
        "scenario", "replicate", "seed", "dataset_fingerprint"
      ))
    )]
  }

  assert_finite_frame(
    metrics,
    "Replication metrics",
    allow_infinite_columns = c(
      "total_set_length",
      "hdr_symmetric_difference_length",
      "excess_length"
    )
  )
  assert_finite_frame(checks, "Equivalence checks")
  assert_finite_frame(tau_rows, "Tau-selection diagnostics")
  assert_finite_frame(domain_row, "Domain diagnostics")
  assert_finite_frame(
    depth_diagnostics,
    "Exact depth-geometry diagnostics",
    allow_infinite_columns = c(
      "calibration_threshold", "exact_total_length", "grid_total_length"
    )
  )
  if (nrow(refinement)) {
    assert_finite_frame(
      refinement,
      "Grid-refinement diagnostics",
      allow_infinite_columns = c(
        "coarse_length", "fine_length"
      )
    )
  }
  probability_columns <- c(
    "coverage", "set_unbounded", "center_included", "valley_excluded",
    "all_modes_included", "peak_inclusion_fraction"
  )
  invalid_probability <- vapply(probability_columns, function(column) {
    values <- metrics[[column]]
    any(values < 0 | values > 1, na.rm = TRUE)
  }, logical(1))
  if (any(invalid_probability)) {
    stop("A replication produced an invalid probability.")
  }

  list(
    scenario = scenario$id[[1L]],
    replicate = as.integer(replicate),
    seed = seed,
    dataset_fingerprint = fingerprint,
    metrics = metrics,
    tau_selection = tau_rows,
    checks = checks,
    domain = domain_row,
    refinement = refinement,
    depth_diagnostics = depth_diagnostics
  )
}

run_comparison_simulation <- function(config, output_dir, workers = 1L, force = FALSE) {
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  cache_dir <- file.path(output_dir, "replicates")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  config_signature <- configuration_fingerprint(config)
  implementation_signature <- comparison_implementation_fingerprint()
  cache_signature <- configuration_fingerprint(list(
    config_signature = config_signature,
    implementation_signature = implementation_signature
  ))

  tasks <- do.call(c, lapply(seq_len(nrow(config$scenarios)), function(index) {
    lapply(seq_len(config$replicates), function(replicate) {
      list(scenario_index = index, replicate = replicate)
    })
  }))

  run_one <- function(task) {
    scenario <- config$scenarios[task$scenario_index, , drop = FALSE]
    path <- file.path(
      cache_dir,
      sprintf("%s_replicate_%04d.rds", scenario$id[[1L]], task$replicate)
    )
    if (file.exists(path) && !force) {
      saved <- tryCatch(readRDS(path), error = function(error) NULL)
      valid <- is.list(saved) &&
        identical(saved$cache_signature, cache_signature) &&
        is.list(saved$result) &&
        identical(saved$result$scenario, scenario$id[[1L]]) &&
        identical(saved$result$replicate, as.integer(task$replicate)) &&
        all(c(
          "metrics", "tau_selection", "checks", "domain", "refinement",
          "depth_diagnostics"
        ) %in% names(saved$result))
      if (valid) return(saved$result)
    }
    result <- simulate_comparison_replicate(
      scenario,
      task$scenario_index,
      task$replicate,
      config
    )
    save_rds_atomic_comparison(
      list(cache_signature = cache_signature, result = result),
      path
    )
    result
  }

  results <- run_parallel(tasks, run_one, workers)
  failures <- vapply(results, inherits, logical(1), what = "try-error")
  if (any(failures)) {
    stop(
      sum(failures), " replication(s) failed. First failure: ",
      as.character(results[[which(failures)[[1L]]]])
    )
  }

  metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
  tau <- do.call(rbind, lapply(results, `[[`, "tau_selection"))
  checks <- do.call(rbind, lapply(results, `[[`, "checks"))
  domains <- do.call(rbind, lapply(results, `[[`, "domain"))
  depth_diagnostics <- do.call(rbind, lapply(results, `[[`, "depth_diagnostics"))
  refinement_parts <- lapply(results, `[[`, "refinement")
  refinement_parts <- refinement_parts[vapply(
    refinement_parts, nrow, integer(1)
  ) > 0L]
  refinement <- if (length(refinement_parts)) {
    do.call(rbind, refinement_parts)
  } else {
    data.frame()
  }
  expected_replicates <- nrow(config$scenarios) * config$replicates
  if (length(results) != expected_replicates ||
      nrow(metrics) != expected_replicates * length(comparison_method_levels) ||
      nrow(domains) != expected_replicates ||
      nrow(depth_diagnostics) != expected_replicates * 2L) {
    stop("The simulation returned an inconsistent number of results.")
  }
  key <- paste(metrics$scenario, metrics$replicate, metrics$method, sep = ":")
  if (anyDuplicated(key)) stop("Duplicate scenario-replicate-method rows were produced.")
  if (!identical(sort(unique(metrics$method)), sort(comparison_method_levels))) {
    stop("The expected five literature methods are not all present.")
  }

  write_csv_atomic(metrics, file.path(output_dir, "raw_replicates.csv"))
  write_csv_atomic(tau, file.path(output_dir, "tau_selection.csv"))
  write_csv_atomic(checks, file.path(output_dir, "equivalence_checks.csv"))
  write_csv_atomic(domains, file.path(output_dir, "domain_diagnostics.csv"))
  write_csv_atomic(
    depth_diagnostics,
    file.path(output_dir, "depth_geometry_diagnostics.csv")
  )
  write_csv_atomic(refinement, file.path(output_dir, "grid_refinement.csv"))
  save_rds_atomic_comparison(
    list(
      config = config,
      config_signature = config_signature,
      implementation_signature = implementation_signature,
      cache_signature = cache_signature,
      r_version = R.version.string,
      package_versions = c(
        parallel = as.character(utils::packageVersion("parallel")),
        stats = as.character(utils::packageVersion("stats"))
      ),
      completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      successful_replicates = expected_replicates,
      result_rows = nrow(metrics),
      paired_methods_per_replicate = length(comparison_method_levels),
      exact_depth_diagnostic_rows = nrow(depth_diagnostics)
    ),
    file.path(output_dir, "simulation_manifest.rds")
  )
  invisible(list(
    metrics = metrics,
    tau = tau,
    checks = checks,
    domains = domains,
    refinement = refinement,
    depth_diagnostics = depth_diagnostics
  ))
}
