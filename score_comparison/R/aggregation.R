aggregate_metric_columns <- function(raw, metrics) {
  groups <- unique(raw[c("scenario", "scenario_label", "regime", "method")])
  groups <- groups[order(
    groups$scenario,
    match(groups$method, comparison_method_levels)
  ), ]
  rows <- lapply(seq_len(nrow(groups)), function(index) {
    group <- groups[index, , drop = FALSE]
    selected <- raw[
      raw$scenario == group$scenario &
        raw$method == group$method,
      ,
      drop = FALSE
    ]
    answer <- group
    answer$replicates <- nrow(selected)
    for (metric in metrics) {
      result <- summarize_with_mcse(selected[[metric]])
      answer[[metric]] <- unname(result[["estimate"]])
      answer[[paste0(metric, "_mcse")]] <- unname(result[["mcse"]])
    }
    answer
  })
  do.call(rbind, rows)
}

summarize_tau_selection <- function(tau_rows, config) {
  selected <- tau_rows[tau_rows$selected, , drop = FALSE]
  scenarios <- unique(selected$scenario)
  rows <- lapply(scenarios, function(scenario) {
    values <- selected$tau[selected$scenario == scenario]
    data.frame(
      scenario = scenario,
      replicates = length(values),
      tau_mean = mean(values),
      tau_mcse = stats::sd(values) / sqrt(length(values)),
      tau_median = stats::median(values),
      tau_q1 = as.numeric(stats::quantile(values, 0.25, type = 1L)),
      tau_q3 = as.numeric(stats::quantile(values, 0.75, type = 1L)),
      probability_tau_zero = mean(values == min(config$tau_grid)),
      probability_tau_max = mean(values == max(config$tau_grid)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_refinement_comparison <- function(population) {
  primary <- population[population$level == "primary", , drop = FALSE]
  refined <- population[population$level == "refined", , drop = FALSE]
  key <- c("scenario", "method")
  columns <- c(
    "tau", "density_spearman", "center_loss", "mode_minus_center",
    "optimum_to_mode", "optimum_location", "population_score_threshold",
    "population_probability_content", "population_total_length",
    "population_components", "population_peak_inclusion_fraction",
    "population_all_modes_included", "population_valley_excluded"
  )
  answer <- merge(
    primary[c(key, columns)],
    refined[c(key, columns)],
    by = key,
    suffixes = c("_primary", "_refined"),
    all = TRUE
  )
  for (column in columns) {
    answer[[paste0(column, "_absolute_change")]] <-
      abs(answer[[paste0(column, "_primary")]] -
            answer[[paste0(column, "_refined")]])
  }
  answer
}

summarize_grid_refinement <- function(refinement) {
  groups <- unique(refinement[c("scenario", "method")])
  groups <- groups[order(
    groups$scenario,
    match(groups$method, comparison_method_levels)
  ), ]
  rows <- lapply(seq_len(nrow(groups)), function(index) {
    group <- groups[index, , drop = FALSE]
    selected <- refinement[
      refinement$scenario == group$scenario &
        refinement$method == group$method,
      ,
      drop = FALSE
    ]
    data.frame(
      scenario = group$scenario,
      method = group$method,
      checked_replicates = nrow(selected),
      maximum_absolute_length_change =
        max(selected$absolute_length_change),
      mean_absolute_length_change =
        mean(selected$absolute_length_change),
      all_component_counts_stable =
        all(selected$component_count_stable),
      coarse_step = max(selected$coarse_step),
      fine_step = min(selected$fine_step),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_domains <- function(domains, config) {
  rows <- lapply(config$scenarios$id, function(scenario) {
    selected <- domains[domains$scenario == scenario, , drop = FALSE]
    data.frame(
      scenario = scenario,
      replicates = nrow(selected),
      maximum_expansions = max(selected$expansions),
      expanded_replicates = sum(selected$expansions > 0L),
      replications_with_unbounded_method =
        sum(selected$unbounded_method_count > 0L),
      minimum_left_outside_points =
        min(selected$minimum_left_outside_points),
      minimum_right_outside_points =
        min(selected$minimum_right_outside_points),
      minimum_evaluated_left = min(selected$evaluated_left),
      maximum_evaluated_right = max(selected$evaluated_right),
      maximum_oracle_content_error = max(abs(
        selected$oracle_probability_content - (1 - config$alpha)
      )),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_depth_geometry_diagnostics <- function(depth_diagnostics) {
  groups <- unique(depth_diagnostics[c("scenario", "scenario_label", "method")])
  groups <- groups[order(
    groups$scenario,
    match(groups$method, comparison_method_levels)
  ), , drop = FALSE]
  rows <- lapply(seq_len(nrow(groups)), function(index) {
    group <- groups[index, , drop = FALSE]
    selected <- depth_diagnostics[
      depth_diagnostics$scenario == group$scenario &
        depth_diagnostics$method == group$method,
      ,
      drop = FALSE
    ]
    data.frame(
      scenario = group$scenario,
      scenario_label = group$scenario_label,
      method = group$method,
      replications = nrow(selected),
      genuinely_unbounded_replications = sum(selected$exact_unbounded),
      unbounded_rate = mean(selected$exact_unbounded),
      maximum_grid_exact_length_difference =
        max(selected$grid_exact_absolute_length_difference),
      mean_grid_exact_length_difference =
        mean(selected$grid_exact_absolute_length_difference),
      all_grid_exact_components_equal =
        all(selected$grid_exact_component_equal),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

validate_raw_results <- function(
    raw,
    tau,
    checks,
    domains,
    depth_diagnostics,
    grid_refinement,
    population,
    population_tau,
    modality,
    population_depth_equivalence,
    scale,
    simulation_manifest,
    population_manifest,
    config) {
  expected <- nrow(config$scenarios) * config$replicates
  scenario_ids <- config$scenarios$id
  level_ids <- names(config$population)
  replicate_key <- unique(raw[c("scenario", "replicate")])
  method_key <- paste(raw$scenario, raw$replicate, raw$method, sep = ":")
  fingerprints <- stats::aggregate(
    dataset_fingerprint ~ scenario + replicate,
    raw,
    function(value) length(unique(value))
  )
  selected_tau <- tau[tau$selected, , drop = FALSE]
  expected_tau_rows <- expected * length(config$tau_grid)
  expected_replicate_keys <- as.vector(outer(
    scenario_ids,
    seq_len(config$replicates),
    paste,
    sep = ":"
  ))
  expected_method_keys <- as.vector(outer(
    expected_replicate_keys,
    comparison_method_levels,
    paste,
    sep = ":"
  ))
  depth_methods <- c(
    "Metric spatial depth",
    "Metric halfspace (Tukey) depth"
  )
  expected_depth_diagnostic_keys <- as.vector(outer(
    expected_replicate_keys,
    depth_methods,
    paste,
    sep = ":"
  ))
  observed_depth_diagnostic_keys <- paste(
    depth_diagnostics$scenario,
    depth_diagnostics$replicate,
    depth_diagnostics$method,
    sep = ":"
  )
  expected_tau_keys <- as.vector(outer(
    expected_replicate_keys,
    format(config$tau_grid, scientific = FALSE, trim = TRUE),
    paste,
    sep = ":"
  ))
  observed_tau_keys <- paste(
    tau$scenario,
    tau$replicate,
    format(tau$tau, scientific = FALSE, trim = TRUE),
    sep = ":"
  )
  expected_population_keys <- as.vector(outer(
    as.vector(outer(level_ids, scenario_ids, paste, sep = ":")),
    comparison_method_levels,
    paste,
    sep = ":"
  ))
  observed_population_keys <- paste(
    population$level, population$scenario, population$method, sep = ":"
  )
  expected_population_tau_keys <- as.vector(outer(
    as.vector(outer(level_ids, scenario_ids, paste, sep = ":")),
    format(config$tau_grid, scientific = FALSE, trim = TRUE),
    paste,
    sep = ":"
  ))
  observed_population_tau_keys <- paste(
    population_tau$level,
    population_tau$scenario,
    format(population_tau$tau, scientific = FALSE, trim = TRUE),
    sep = ":"
  )
  expected_level_scenario_keys <- as.vector(outer(
    level_ids, scenario_ids, paste, sep = ":"
  ))
  observed_modality_keys <- paste(
    modality$level, modality$scenario, sep = ":"
  )
  observed_depth_keys <- paste(
    population_depth_equivalence$level,
    population_depth_equivalence$scenario,
    sep = ":"
  )
  expected_refinement_keys <- as.vector(outer(
    as.vector(outer(
      scenario_ids,
      unique(config$refinement_replicates),
      paste,
      sep = ":"
    )),
    comparison_method_levels,
    paste,
    sep = ":"
  ))
  observed_refinement_keys <- paste(
    grid_refinement$scenario,
    grid_refinement$replicate,
    grid_refinement$method,
    sep = ":"
  )
  config_signature <- configuration_fingerprint(config)
  probability_columns <- c(
    "coverage", "set_unbounded", "center_included", "valley_excluded",
    "all_modes_included", "peak_inclusion_fraction"
  )
  probability_ranges_valid <- all(vapply(probability_columns, function(column) {
    values <- raw[[column]]
    all(values >= 0 & values <= 1, na.rm = TRUE)
  }, logical(1)))
  depth_pair_columns <- c(
    "scenario", "replicate", "set_unbounded", "components", "total_set_length"
  )
  spatial_rows <- raw[
    raw$method == "Metric spatial depth",
    depth_pair_columns,
    drop = FALSE
  ]
  halfspace_rows <- raw[
    raw$method == "Metric halfspace (Tukey) depth",
    depth_pair_columns,
    drop = FALSE
  ]
  paired_depth_rows <- merge(
    spatial_rows,
    halfspace_rows,
    by = c("scenario", "replicate"),
    suffixes = c("_spatial", "_halfspace"),
    all = TRUE,
    sort = FALSE
  )
  paired_depth_length_difference <- ifelse(
    is.infinite(paired_depth_rows$total_set_length_spatial) &
      is.infinite(paired_depth_rows$total_set_length_halfspace),
    0,
    abs(
      paired_depth_rows$total_set_length_spatial -
        paired_depth_rows$total_set_length_halfspace
    )
  )

  conditions <- data.frame(
    check = c(
      "configuration_manifests",
      "unique_replication_count",
      "unique_five_method_rows",
      "one_dataset_per_replication",
      "one_selected_tau_per_replication",
      "complete_tau_diagnostics",
      "complete_equivalence_diagnostics",
      "complete_domain_diagnostics",
      "complete_depth_geometry_diagnostics",
      "complete_grid_refinement",
      "complete_population_summary",
      "complete_population_tau_diagnostics",
      "complete_modality_diagnostics",
      "complete_population_depth_equivalence",
      "c0_rank_numeric_equivalence",
      "published_rank_constant_shift",
      "published_rank_order_equivalence",
      "sample_depth_transform_equivalence",
      "sample_depth_rank_equivalence",
      "sample_depth_region_equivalence",
      "sample_depth_boundedness_equivalence",
      "sample_depth_component_equivalence",
      "sample_depth_length_equivalence",
      "population_depth_transform_equivalence",
      "population_depth_rank_equivalence",
      "population_depth_region_equivalence",
      "sample_counts",
      "probability_ranges",
      "raw_metric_finiteness_and_unbounded_pattern",
      "domain_outside_margins",
      "oracle_probability_content",
      "grid_length_refinement_tolerance",
      "grid_component_refinement_stability",
      "exact_depth_geometry_audit",
      "population_probability_content",
      "population_scale_homogeneity",
      "population_scale_ordering"
    ),
    passed = c(
      identical(simulation_manifest$config_signature, config_signature) &&
        identical(population_manifest$config_signature, config_signature) &&
        as.integer(simulation_manifest$successful_replicates) == expected &&
        as.integer(simulation_manifest$paired_methods_per_replicate) ==
          length(comparison_method_levels) &&
        as.integer(population_manifest$expected_summary_rows) ==
          length(expected_population_keys),
      nrow(replicate_key) == expected &&
        setequal(
          paste(replicate_key$scenario, replicate_key$replicate, sep = ":"),
          expected_replicate_keys
        ),
      !anyDuplicated(method_key) &&
        setequal(method_key, expected_method_keys),
      all(fingerprints$dataset_fingerprint == 1L),
      nrow(selected_tau) == expected &&
        !anyDuplicated(paste(selected_tau$scenario, selected_tau$replicate)),
      nrow(tau) == expected_tau_rows &&
        !anyDuplicated(observed_tau_keys) &&
        setequal(observed_tau_keys, expected_tau_keys),
      nrow(checks) == expected &&
        !anyDuplicated(paste(checks$scenario, checks$replicate)) &&
        setequal(
          paste(checks$scenario, checks$replicate, sep = ":"),
          expected_replicate_keys
        ),
      nrow(domains) == expected &&
        !anyDuplicated(paste(domains$scenario, domains$replicate)) &&
        setequal(
          paste(domains$scenario, domains$replicate, sep = ":"),
          expected_replicate_keys
        ),
      nrow(depth_diagnostics) == 2L * expected &&
        !anyDuplicated(observed_depth_diagnostic_keys) &&
        setequal(
          observed_depth_diagnostic_keys,
          expected_depth_diagnostic_keys
        ) &&
        as.integer(simulation_manifest$exact_depth_diagnostic_rows) ==
          2L * expected,
      !anyDuplicated(observed_refinement_keys) &&
        setequal(observed_refinement_keys, expected_refinement_keys),
      length(observed_population_keys) == length(expected_population_keys) &&
        !anyDuplicated(observed_population_keys) &&
        setequal(observed_population_keys, expected_population_keys),
      length(observed_population_tau_keys) ==
        length(expected_population_tau_keys) &&
        !anyDuplicated(observed_population_tau_keys) &&
        setequal(observed_population_tau_keys, expected_population_tau_keys),
      length(observed_modality_keys) == length(expected_level_scenario_keys) &&
        !anyDuplicated(observed_modality_keys) &&
        setequal(observed_modality_keys, expected_level_scenario_keys),
      length(observed_depth_keys) == length(expected_level_scenario_keys) &&
        !anyDuplicated(observed_depth_keys) &&
        setequal(observed_depth_keys, expected_level_scenario_keys),
      max(checks$c0_rank_max_abs) < 1e-12,
      max(checks$published_raw_shift_range) < 1e-10,
      min(checks$published_rank_correlation) > 0.999,
      max(checks$spatial_halfspace_transform_max_abs) < 1e-12,
      min(checks$spatial_halfspace_rank_correlation) > 0.999999,
      all(checks$spatial_halfspace_grid_set_equal) &&
        all(checks$spatial_halfspace_test_set_equal),
      nrow(paired_depth_rows) == expected &&
        all(
          paired_depth_rows$set_unbounded_spatial ==
            paired_depth_rows$set_unbounded_halfspace
        ),
      nrow(paired_depth_rows) == expected &&
        all(
          paired_depth_rows$components_spatial ==
            paired_depth_rows$components_halfspace
        ),
      nrow(paired_depth_rows) == expected &&
        max(paired_depth_length_difference) <= 1e-12,
      max(population_depth_equivalence$transform_max_abs) < 1e-12,
      min(population_depth_equivalence$rank_correlation) > 0.9999,
      max(
        population_depth_equivalence$region_symmetric_difference_length
      ) <= config$length_tolerance &&
        all(population_depth_equivalence$component_count_equal) &&
        all(population_depth_equivalence$all_modes_membership_equal),
      all(checks$n_training == config$n_train) &&
        all(checks$n_calibration == config$n_calibration) &&
        all(checks$n_test == config$n_test),
      probability_ranges_valid,
      {
        allowed_infinite <- c(
          "total_set_length",
          "hdr_symmetric_difference_length",
          "excess_length"
        )
        numeric_columns <- names(raw)[vapply(raw, is.numeric, logical(1))]
        ordinary_columns <- setdiff(numeric_columns, allowed_infinite)
        ordinary_valid <- all(vapply(
          raw[ordinary_columns],
          function(x) all(is.finite(x) | is.na(x)),
          logical(1)
        ))
        infinite_pattern_valid <- all(vapply(
          allowed_infinite,
          function(column) {
            is.infinite(raw[[column]]) == (raw$set_unbounded == 1)
          },
          logical(nrow(raw))
        ))
        ordinary_valid && infinite_pattern_valid
      },
      all(domains$minimum_left_outside_points >=
            config$outside_margin_points) &&
        all(domains$minimum_right_outside_points >=
              config$outside_margin_points),
      max(abs(
        domains$oracle_probability_content - (1 - config$alpha)
      )) < 1e-6,
      max(grid_refinement$absolute_length_change) <=
        config$length_tolerance + 1e-10,
      all(grid_refinement$component_count_stable),
      all(depth_diagnostics$grid_exact_component_equal) &&
        max(depth_diagnostics$grid_exact_absolute_length_difference) <=
          2 * diff(config$candidate_grid)[[1L]] + 1e-10 &&
        all(
          depth_diagnostics$exact_unbounded ==
            as.numeric(depth_diagnostics$calibration_threshold <= 0)
        ),
      max(abs(
        population$population_probability_content - (1 - config$alpha)
      )) < 1e-6,
      all(scale$max_homogeneity_error < 1e-10),
      all(scale$ordering_preserved)
    ),
    detail = c(
      paste(
        simulation_manifest$config_signature,
        population_manifest$config_signature
      ),
      sprintf("%d of %d", nrow(replicate_key), expected),
      sprintf("%d method rows", nrow(raw)),
      sprintf(
        "maximum distinct fingerprints: %d",
        max(fingerprints$dataset_fingerprint)
      ),
      sprintf("%d selected rows", nrow(selected_tau)),
      sprintf("%d of %d", nrow(tau), expected_tau_rows),
      sprintf("%d of %d", nrow(checks), expected),
      sprintf("%d of %d", nrow(domains), expected),
      sprintf(
        "%d of %d",
        nrow(depth_diagnostics),
        2L * expected
      ),
      sprintf(
        "%d of %d",
        length(observed_refinement_keys),
        length(expected_refinement_keys)
      ),
      sprintf(
        "%d of %d",
        nrow(population),
        length(expected_population_keys)
      ),
      sprintf(
        "%d of %d",
        nrow(population_tau),
        length(expected_population_tau_keys)
      ),
      sprintf(
        "%d of %d",
        nrow(modality),
        length(expected_level_scenario_keys)
      ),
      sprintf(
        "%d of %d",
        nrow(population_depth_equivalence),
        length(expected_level_scenario_keys)
      ),
      format(max(checks$c0_rank_max_abs), scientific = TRUE),
      format(max(checks$published_raw_shift_range), scientific = TRUE),
      format(min(checks$published_rank_correlation), digits = 8L),
      format(
        max(checks$spatial_halfspace_transform_max_abs),
        scientific = TRUE
      ),
      format(
        min(checks$spatial_halfspace_rank_correlation),
        digits = 8L
      ),
      paste(
        all(checks$spatial_halfspace_grid_set_equal),
        all(checks$spatial_halfspace_test_set_equal)
      ),
      as.character(all(
        paired_depth_rows$set_unbounded_spatial ==
          paired_depth_rows$set_unbounded_halfspace
      )),
      as.character(all(
        paired_depth_rows$components_spatial ==
          paired_depth_rows$components_halfspace
      )),
      format(max(paired_depth_length_difference), scientific = TRUE),
      format(
        max(population_depth_equivalence$transform_max_abs),
        scientific = TRUE
      ),
      format(
        min(population_depth_equivalence$rank_correlation),
        digits = 8L
      ),
      format(
        max(
          population_depth_equivalence$region_symmetric_difference_length
        ),
        scientific = TRUE
      ),
      sprintf("%d/%d/%d", config$n_train, config$n_calibration, config$n_test),
      "all nonmissing reported probabilities in [0,1]",
      paste(
        "ordinary numeric metrics are finite;",
        "length, HDR, and excess are infinite exactly for unbounded sets"
      ),
      sprintf(
        "minimum margins %d/%d",
        min(domains$minimum_left_outside_points),
        min(domains$minimum_right_outside_points)
      ),
      format(
        max(abs(
          domains$oracle_probability_content - (1 - config$alpha)
        )),
        scientific = TRUE
      ),
      format(max(grid_refinement$absolute_length_change), scientific = TRUE),
      as.character(all(grid_refinement$component_count_stable)),
      paste(
        format(
          max(depth_diagnostics$grid_exact_absolute_length_difference),
          scientific = TRUE
        ),
        all(depth_diagnostics$grid_exact_component_equal)
      ),
      format(
        max(abs(
          population$population_probability_content - (1 - config$alpha)
        )),
        scientific = TRUE
      ),
      format(max(scale$max_homogeneity_error), scientific = TRUE),
      paste(scale$ordering_preserved, collapse = ",")
    ),
    stringsAsFactors = FALSE
  )
  if (any(!conditions$passed)) {
    stop(
      "Validation failed: ",
      paste(conditions$check[!conditions$passed], collapse = ", ")
    )
  }

  alternate_coverage <- tapply(
    raw$coverage,
    list(raw$scenario, raw$method),
    base::mean
  )
  alternate_hdr <- tapply(
    raw$hdr_symmetric_difference_length,
    list(raw$scenario, raw$method),
    base::mean
  )
  alternate_length <- tapply(
    raw$total_set_length,
    list(raw$scenario, raw$method),
    base::mean
  )
  alternate_peak <- tapply(
    raw$peak_inclusion_fraction,
    list(raw$scenario, raw$method),
    base::mean
  )
  list(
    checks = conditions,
    alternate_coverage = alternate_coverage,
    alternate_hdr = alternate_hdr,
    alternate_length = alternate_length,
    alternate_peak = alternate_peak
  )
}

aggregate_comparison_results <- function(config, input_dir, output_dir = input_dir) {
  required <- c(
    "raw_replicates.csv", "tau_selection.csv", "equivalence_checks.csv",
    "domain_diagnostics.csv", "depth_geometry_diagnostics.csv",
    "grid_refinement.csv",
    "population_summary.csv", "population_tau_selection.csv",
    "population_depth_equivalence.csv", "modality_checks.csv",
    "scale_equivariance.csv", "simulation_manifest.rds",
    "population_manifest.rds"
  )
  absent <- required[!file.exists(file.path(input_dir, required))]
  if (length(absent)) {
    stop("Required result files are missing: ", paste(absent, collapse = ", "))
  }
  read_result <- function(name) {
    utils::read.csv(file.path(input_dir, name), check.names = FALSE)
  }
  raw <- read_result("raw_replicates.csv")
  tau <- read_result("tau_selection.csv")
  checks <- read_result("equivalence_checks.csv")
  domains <- read_result("domain_diagnostics.csv")
  depth_diagnostics <- read_result("depth_geometry_diagnostics.csv")
  grid_refinement <- read_result("grid_refinement.csv")
  population <- read_result("population_summary.csv")
  population_tau <- read_result("population_tau_selection.csv")
  modality <- read_result("modality_checks.csv")
  population_depth_equivalence <-
    read_result("population_depth_equivalence.csv")
  scale <- read_result("scale_equivariance.csv")
  simulation_manifest <- readRDS(file.path(input_dir, "simulation_manifest.rds"))
  population_manifest <- readRDS(file.path(input_dir, "population_manifest.rds"))

  metric_names <- c(
    "coverage", "total_set_length", "finite_replication_total_length",
    "set_unbounded", "components", "hdr_symmetric_difference_length",
    "finite_replication_hdr_symmetric_difference_length",
    "oracle_hdr_length", "excess_length",
    "center_included", "valley_excluded", "all_modes_included",
    "peak_inclusion_fraction", "density_spearman", "center_loss",
    "mode_minus_center"
  )
  aggregate <- aggregate_metric_columns(raw, metric_names)
  tau_summary <- summarize_tau_selection(tau, config)
  population_refinement <- make_refinement_comparison(population)
  grid_refinement_summary <- summarize_grid_refinement(grid_refinement)
  domain_summary <- summarize_domains(domains, config)
  depth_geometry_summary <-
    summarize_depth_geometry_diagnostics(depth_diagnostics)
  validation <- validate_raw_results(
    raw,
    tau,
    checks,
    domains,
    depth_diagnostics,
    grid_refinement,
    population,
    population_tau,
    modality,
    population_depth_equivalence,
    scale,
    simulation_manifest,
    population_manifest,
    config
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv_atomic(aggregate, file.path(output_dir, "aggregate_summary.csv"))
  write_csv_atomic(tau_summary, file.path(output_dir, "tau_summary.csv"))
  write_csv_atomic(
    population_refinement,
    file.path(output_dir, "population_refinement.csv")
  )
  write_csv_atomic(
    grid_refinement_summary,
    file.path(output_dir, "grid_refinement_summary.csv")
  )
  write_csv_atomic(domain_summary, file.path(output_dir, "domain_summary.csv"))
  write_csv_atomic(
    depth_geometry_summary,
    file.path(output_dir, "depth_geometry_summary.csv")
  )
  write_csv_atomic(
    validation$checks,
    file.path(output_dir, "validation_checks.csv")
  )
  write_csv_atomic(
    as.data.frame(as.table(validation$alternate_coverage)),
    file.path(output_dir, "independent_coverage_recalculation.csv")
  )
  write_csv_atomic(
    as.data.frame(as.table(validation$alternate_hdr)),
    file.path(output_dir, "independent_hdr_recalculation.csv")
  )
  write_csv_atomic(
    as.data.frame(as.table(validation$alternate_length)),
    file.path(output_dir, "independent_total_length_recalculation.csv")
  )
  write_csv_atomic(
    as.data.frame(as.table(validation$alternate_peak)),
    file.path(output_dir, "independent_peak_recalculation.csv")
  )
  invisible(list(
    aggregate = aggregate,
    tau_summary = tau_summary,
    population_refinement = population_refinement,
    grid_refinement_summary = grid_refinement_summary,
    domain_summary = domain_summary,
    depth_geometry_summary = depth_geometry_summary,
    validation = validation
  ))
}
