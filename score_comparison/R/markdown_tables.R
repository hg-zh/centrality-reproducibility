format_number_mcse <- function(estimate, mcse, digits = 3L, multiplier = 1) {
  if (is.na(estimate)) return("\u2014")
  if (is.infinite(estimate)) return("$\\infty$ (MCSE undefined)")
  if (is.na(mcse)) {
    return(formatC(
      estimate * multiplier,
      digits = digits,
      format = "f"
    ))
  }
  paste0(
    formatC(estimate * multiplier, digits = digits, format = "f"),
    " (",
    formatC(mcse * multiplier, digits = digits, format = "f"),
    ")"
  )
}

format_length_summary <- function(
    estimate,
    mcse,
    finite_estimate,
    finite_mcse,
    digits = 3L) {
  if (is.infinite(estimate)) {
    return(paste0(
      "$\\infty$; finite-only ",
      format_number_mcse(finite_estimate, finite_mcse, digits = digits)
    ))
  }
  format_number_mcse(estimate, mcse, digits = digits)
}

format_tau <- function(value) {
  if (is.na(value)) "\u2014" else format(value, trim = TRUE, scientific = FALSE)
}

markdown_table <- function(headers, rows, align = NULL) {
  if (is.null(align)) align <- rep("---", length(headers))
  lines <- c(
    paste0("| ", paste(headers, collapse = " | "), " |"),
    paste0("| ", paste(align, collapse = " | "), " |")
  )
  body <- apply(rows, 1L, function(row) {
    paste0("| ", paste(row, collapse = " | "), " |")
  })
  paste(c(lines, body), collapse = "\n")
}

population_geometry_label <- function(row) {
  component_value <- as.numeric(row[["population_components"]])
  all_modes_value <- as.numeric(row[["population_all_modes_included"]])
  valley_value <- as.numeric(row[["population_valley_excluded"]])
  components <- paste0(
    as.integer(component_value),
    if (component_value == 1) " component" else " components"
  )
  peaks <- if (all_modes_value == 1) {
    "peak(s) in"
  } else {
    "peak(s) not all in"
  }
  valley <- if (is.na(valley_value)) {
    "no antimode"
  } else if (valley_value == 1) {
    "valley out"
  } else {
    "valley in"
  }
  paste(components, peaks, valley, sep = "; ")
}

make_results_markdown <- function(config, input_dir, output_path) {
  read_result <- function(name) {
    utils::read.csv(file.path(input_dir, name), check.names = FALSE)
  }
  aggregate <- read_result("aggregate_summary.csv")
  population <- read_result("population_summary.csv")
  tau_summary <- read_result("tau_summary.csv")
  grid_refinement <- read_result("grid_refinement_summary.csv")
  domains <- read_result("domain_summary.csv")
  depth_geometry <- read_result("depth_geometry_summary.csv")
  depth_equivalence <- read_result("population_depth_equivalence.csv")
  validation <- read_result("validation_checks.csv")

  primary <- population[population$level == "primary", , drop = FALSE]
  refined <- population[population$level == "refined", , drop = FALSE]
  scenario_order <- config$scenarios$id
  method_order <- comparison_method_levels
  expected_keys <- as.vector(outer(
    scenario_order, method_order, paste, sep = ":"
  ))
  primary_keys <- paste(primary$scenario, primary$method, sep = ":")
  refined_keys <- paste(refined$scenario, refined$method, sep = ":")
  if (anyDuplicated(primary_keys) || anyDuplicated(refined_keys) ||
      !setequal(primary_keys, expected_keys) ||
      !setequal(refined_keys, expected_keys)) {
    stop("Population summary is incomplete or duplicated.")
  }
  refined$tau_primary <- primary$tau[match(refined_keys, primary_keys)]
  order_rows <- function(frame) {
    frame[order(
      match(frame$scenario, scenario_order),
      match(frame$method, method_order)
    ), , drop = FALSE]
  }
  aggregate <- order_rows(aggregate)
  refined <- order_rows(refined)
  merged <- merge(
    aggregate,
    refined[c(
      "scenario", "method", "tau", "tau_primary", "mode_count",
      "mode_locations", "density_spearman", "center_loss",
      "mode_minus_center", "optimum_to_mode", "optimum_location",
      "population_score_threshold", "population_probability_content",
      "population_total_length", "population_components",
      "population_center_included", "population_peak_inclusion_fraction",
      "population_all_modes_included", "population_valley_excluded",
      "oracle_hdr_length", "oracle_hdr_components"
    )],
    by = c("scenario", "method"),
    suffixes = c("_finite", "_population"),
    sort = FALSE
  )
  merged <- order_rows(merged)
  merged$population_tau_display <- mapply(
    function(primary_value, refined_value) {
      if (is.na(refined_value)) return("\u2014")
      if (isTRUE(all.equal(primary_value, refined_value))) {
        return(format_tau(refined_value))
      }
      paste0(format_tau(primary_value), "\u2192", format_tau(refined_value))
    },
    merged$tau_primary,
    merged$tau,
    USE.NAMES = FALSE
  )
  merged$population_geometry <- apply(
    merged,
    1L,
    population_geometry_label
  )

  unimodal_ids <- config$scenarios$id[config$scenarios$regime == "unimodal"]
  unimodal <- merged[merged$scenario %in% unimodal_ids, , drop = FALSE]
  unimodal_rows <- cbind(
    unimodal$scenario_label,
    unimodal$method,
    unimodal$population_tau_display,
    unimodal$population_geometry,
    mapply(
      format_number_mcse,
      unimodal$coverage,
      unimodal$coverage_mcse,
      MoreArgs = list(digits = 2L, multiplier = 100)
    ),
    mapply(
      format_length_summary,
      unimodal$total_set_length,
      unimodal$total_set_length_mcse,
      unimodal$finite_replication_total_length,
      unimodal$finite_replication_total_length_mcse,
      MoreArgs = list(digits = 3L)
    ),
    mapply(
      format_number_mcse,
      unimodal$set_unbounded,
      unimodal$set_unbounded_mcse,
      MoreArgs = list(digits = 1L, multiplier = 100)
    ),
    mapply(
      format_number_mcse,
      unimodal$peak_inclusion_fraction,
      unimodal$peak_inclusion_fraction_mcse,
      MoreArgs = list(digits = 1L, multiplier = 100)
    ),
    mapply(
      format_length_summary,
      unimodal$hdr_symmetric_difference_length,
      unimodal$hdr_symmetric_difference_length_mcse,
      unimodal$finite_replication_hdr_symmetric_difference_length,
      unimodal$finite_replication_hdr_symmetric_difference_length_mcse,
      MoreArgs = list(digits = 3L)
    )
  )
  colnames(unimodal_rows) <- NULL

  mixture_ids <- config$scenarios$id[config$scenarios$regime != "unimodal"]
  mixtures <- merged[merged$scenario %in% mixture_ids, , drop = FALSE]
  mixture_labels <- paste0(
    sub("Gaussian mixture, ", "", mixtures$scenario_label),
    " (", mixtures$mode_count, " mode",
    ifelse(mixtures$mode_count == 1, "", "s"), ")"
  )
  mixture_rows <- cbind(
    mixture_labels,
    mixtures$method,
    mixtures$population_tau_display,
    mixtures$population_geometry,
    mapply(
      format_number_mcse,
      mixtures$coverage,
      mixtures$coverage_mcse,
      MoreArgs = list(digits = 2L, multiplier = 100)
    ),
    mapply(
      format_length_summary,
      mixtures$total_set_length,
      mixtures$total_set_length_mcse,
      mixtures$finite_replication_total_length,
      mixtures$finite_replication_total_length_mcse,
      MoreArgs = list(digits = 3L)
    ),
    mapply(
      format_number_mcse,
      mixtures$set_unbounded,
      mixtures$set_unbounded_mcse,
      MoreArgs = list(digits = 1L, multiplier = 100)
    ),
    mapply(
      format_number_mcse,
      mixtures$all_modes_included,
      mixtures$all_modes_included_mcse,
      MoreArgs = list(digits = 1L, multiplier = 100)
    ),
    mapply(
      format_number_mcse,
      mixtures$valley_excluded,
      mixtures$valley_excluded_mcse,
      MoreArgs = list(digits = 1L, multiplier = 100)
    ),
    mapply(
      format_number_mcse,
      mixtures$components,
      mixtures$components_mcse,
      MoreArgs = list(digits = 2L)
    ),
    mapply(
      format_length_summary,
      mixtures$hdr_symmetric_difference_length,
      mixtures$hdr_symmetric_difference_length_mcse,
      mixtures$finite_replication_hdr_symmetric_difference_length,
      mixtures$finite_replication_hdr_symmetric_difference_length_mcse,
      MoreArgs = list(digits = 3L)
    )
  )
  colnames(mixture_rows) <- NULL

  tau_summary <- merge(
    config$scenarios[c("id", "label")],
    tau_summary,
    by.x = "id",
    by.y = "scenario",
    sort = FALSE
  )
  tau_summary <- tau_summary[match(scenario_order, tau_summary$id), ]
  tau_rows <- cbind(
    tau_summary$label,
    paste0(
      vapply(tau_summary$tau_median, format_tau, character(1)),
      " [",
      vapply(tau_summary$tau_q1, format_tau, character(1)),
      ", ",
      vapply(tau_summary$tau_q3, format_tau, character(1)),
      "]"
    ),
    mapply(
      format_number_mcse,
      tau_summary$probability_tau_zero,
      sqrt(
        tau_summary$probability_tau_zero *
          (1 - tau_summary$probability_tau_zero) /
          tau_summary$replicates
      ),
      MoreArgs = list(digits = 1L, multiplier = 100)
    ),
    mapply(
      format_number_mcse,
      tau_summary$probability_tau_max,
      sqrt(
        tau_summary$probability_tau_max *
          (1 - tau_summary$probability_tau_max) /
          tau_summary$replicates
      ),
      MoreArgs = list(digits = 1L, multiplier = 100)
    )
  )
  colnames(tau_rows) <- NULL

  probability_percent <- format(
    100 * (1 - config$alpha), trim = TRUE, scientific = FALSE
  )
  tau_max_label <- paste0(
    "Pr($\\tau=", format_tau(max(config$tau_grid)), "$) %, estimate (MCSE)"
  )
  maximum_refinement_error <- max(
    grid_refinement$maximum_absolute_length_change
  )
  all_components_stable <- all(
    grid_refinement$all_component_counts_stable
  )
  minimum_left_margin <- min(domains$minimum_left_outside_points)
  minimum_right_margin <- min(domains$minimum_right_outside_points)
  maximum_population_depth_region_difference <- max(
    depth_equivalence$region_symmetric_difference_length
  )
  maximum_depth_grid_exact_length_difference <- max(
    depth_geometry$maximum_grid_exact_length_difference
  )

  text <- c(
    "# Five-method score-comparison results",
    "",
    paste0(
      "All scores are oriented so that larger means more central. ",
      "Parentheses give Monte Carlo standard errors (MCSEs). The ",
      "population calculations are deterministic and therefore have no MCSE."
    ),
    "",
    "## Symmetric unimodal distributions",
    "",
    markdown_table(
      c(
        "Distribution", "Method", "Population $\\tau$",
        "Population geometry", "Coverage %, estimate (MCSE)",
        "Total length, estimate (MCSE)",
        "Unbounded %, estimate (MCSE)",
        "Peak inclusion %, estimate (MCSE)",
        paste0(
          probability_percent,
          "% HDR symmetric-difference length, estimate (MCSE)"
        )
      ),
      unimodal_rows
    ),
    "",
    paste0(
      "Coverage is the Monte Carlo estimate of ",
      "$P\\{X_{\\mathrm{new}}\\in\\widehat C\\}$. Total length is ",
      "$\\lambda(\\widehat C)$, summed over every connected component. ",
      "Peak inclusion evaluates the exact population mode at zero. The HDR ",
      "column is $\\lambda(\\widehat C\\mathbin{\\triangle}H_{",
      probability_percent,
      "\\%})$; it has no denominator and is neither a probability-mass ",
      "capture rate nor a Jaccard index. If any replication is the whole ",
      "real line, the unconditional mean length and HDR error are infinite; ",
      "the table then also gives the finite-replication conditional mean and ",
      "the unbounded-set rate."
    ),
    "",
    "## Gaussian-mixture transition and bimodal distributions",
    "",
    markdown_table(
      c(
        "Separation", "Method", "Population $\\tau$",
        "Population geometry", "Coverage %, estimate (MCSE)",
        "Total length, estimate (MCSE)",
        "Unbounded %, estimate (MCSE)",
        "All true modes in %, estimate (MCSE)",
        "True valley excluded %, estimate (MCSE)",
        "Components, mean (MCSE)",
        paste0(
          probability_percent,
          "% HDR symmetric-difference length, estimate (MCSE)"
        )
      ),
      mixture_rows
    ),
    "",
    paste0(
      "Population geometry is the upper score level set whose probability ",
      "content is ", probability_percent,
      "%. “Peak(s) in” evaluates the true density mode locations directly; ",
      "“valley out” evaluates the true antimode at zero directly. “No ",
      "antimode” marks the one-mode transition setting. Finite-sample ",
      "depth lengths use exact empirical order-statistic cells; the other ",
      "three methods use linearly interpolated threshold crossings."
    ),
    "",
    "## Prespecified tuning-rule diagnostics",
    "",
    markdown_table(
      c(
        "Distribution", "Selected $\\tau$, median [IQR]",
        "Pr($\\tau=0$) %, estimate (MCSE)",
        tau_max_label
      ),
      tau_rows
    ),
    "",
    "## Numerical and theoretical checks",
    "",
    paste0(
      "The main length grid had spacing ",
      format(diff(config$candidate_grid)[[1L]], digits = 8L, trim = TRUE),
      "; selected replications were repeated at spacing ",
      format(
        diff(config$candidate_grid)[[1L]] / config$refinement_factor,
        digits = 8L,
        trim = TRUE
      ),
      ". The maximum observed absolute length change was ",
      formatC(maximum_refinement_error, digits = 6L, format = "fg"),
      ", and all component counts were stable: ",
      tolower(as.character(all_components_stable)),
      ". The smallest verified outside margins were ",
      minimum_left_margin, " grid points on the left and ",
      minimum_right_margin, " on the right; domains expand automatically ",
      "if either margin fails. Before replacing depth-grid interpolation, ",
      "the largest independent grid-versus-exact depth-length difference was ",
      formatC(
        maximum_depth_grid_exact_length_difference,
        digits = 6L,
        format = "fg"
      ),
      "; the reported depth lengths now use the exact calculation."
    ),
    "",
    paste0(
      "For continuous real-valued laws, metric spatial depth satisfies ",
      "$D_S=4D_H(1-D_H)$, so it and metric halfspace depth have the same ",
      "population ordering. Both depth definitions are retained, but they are ",
      "theoretically redundant rather than independent evidence. Their ",
      "largest numerical population-region symmetric difference was ",
      formatC(
        maximum_population_depth_region_difference,
        digits = 6L,
        format = "fg"
      ),
      "; all finite-sample grid and test memberships agreed in the ",
      "continuous simulations."
    ),
    "",
    paste0(
      "**Design note.** This ", config$profile, " run used ",
      config$replicates, " independent replications, ",
      config$n_train, " training, ", config$n_calibration,
      " calibration, and ", config$n_test,
      " independent test observations per replication. The nominal ",
      "coverage was ", probability_percent,
      "%. Every method used the same simulated samples, split, candidate ",
      "points, and test observations. The proposed score selected ",
      "$\\tau\\in\\{", paste(config$tau_grid, collapse = ","),
      "\\}$ by minimizing the manuscript midpoint proportion plus ",
      format(config$tau_penalty, scientific = FALSE),
      "$\\tau$. The calibration sample is reused for this tuning and final ",
      "threshold, so fixed-score split-conformal validity does not ",
      "automatically extend to the tuned C-score."
    ),
    "",
    paste0(
      "All ", sum(validation$passed), " saved-result validation checks passed."
    )
  )
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(text, output_path, useBytes = TRUE)
  invisible(output_path)
}
