comparison_method_levels <- c(
  "Proposed C-score",
  "Transport rank",
  "Transport cost",
  "Metric spatial depth",
  "Metric halfspace (Tukey) depth"
)

comparison_script_root <- function() {
  argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(argument)) {
    stop("This entry point must be run with Rscript.")
  }
  script <- normalizePath(sub("^--file=", "", argument[[1L]]), mustWork = TRUE)
  normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
}

parse_comparison_arguments <- function(
    args = commandArgs(trailingOnly = TRUE),
    default_config = "config/pilot_config.R",
    default_output = "results/pilot") {
  answer <- list(
    config = default_config,
    output_dir = default_output,
    workers = 1L,
    force = FALSE
  )
  for (argument in args) {
    if (grepl("^--config=", argument)) {
      answer$config <- sub("^--config=", "", argument)
    } else if (grepl("^--output-dir=", argument)) {
      answer$output_dir <- sub("^--output-dir=", "", argument)
    } else if (grepl("^--workers=", argument)) {
      answer$workers <- as.integer(sub("^--workers=", "", argument))
    } else if (identical(argument, "--force")) {
      answer$force <- TRUE
    } else {
      stop("Unknown command-line argument: ", argument)
    }
  }
  if (is.na(answer$workers) || answer$workers < 1L) {
    stop("--workers must be a positive integer.")
  }
  answer
}

resolve_from_root <- function(root, path, must_work = FALSE) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    normalizePath(path, mustWork = must_work)
  } else {
    normalizePath(file.path(root, path), mustWork = must_work)
  }
}

load_comparison_config <- function(root, path) {
  config_path <- resolve_from_root(root, path, must_work = TRUE)
  environment <- new.env(parent = baseenv())
  environment$.comparison_module_root <- root
  sys.source(config_path, envir = environment)
  if (!exists("comparison_config", envir = environment, inherits = FALSE)) {
    stop("The configuration must define comparison_config().")
  }
  config <- environment$comparison_config()
  validate_comparison_config(config)
  config
}

validate_comparison_config <- function(config) {
  required <- c(
    "profile", "replicates", "seed_base", "alpha", "n_train", "n_calibration",
    "n_test", "candidate_grid", "tau_grid", "tau_penalty", "scenarios",
    "population", "chunk_size", "outside_margin_points",
    "domain_expansion_fraction", "domain_max_expansions",
    "refinement_replicates", "refinement_factor", "length_tolerance"
  )
  absent <- setdiff(required, names(config))
  if (length(absent)) stop("Configuration fields are missing: ", paste(absent, collapse = ", "))
  if (config$replicates < 1L || config$n_train < 2L ||
      config$n_calibration < 2L || config$n_test < 1L) {
    stop("Replication and sample sizes must be positive.")
  }
  if (!is.numeric(config$alpha) || config$alpha <= 0 || config$alpha >= 1) {
    stop("alpha must lie strictly between zero and one.")
  }
  if (!is.numeric(config$candidate_grid) ||
      length(config$candidate_grid) < 3L ||
      is.unsorted(config$candidate_grid, strictly = TRUE)) {
    stop("candidate_grid must be a strictly increasing numeric vector.")
  }
  if (!is.numeric(config$tau_grid) || any(config$tau_grid < 0) ||
      is.unsorted(config$tau_grid, strictly = TRUE) ||
      !any(config$tau_grid == 0)) {
    stop("tau_grid must be strictly increasing, nonnegative, and contain zero.")
  }
  if (config$outside_margin_points < 1L ||
      config$domain_expansion_fraction <= 0 ||
      config$domain_max_expansions < 0L ||
      config$refinement_factor < 2L ||
      config$length_tolerance <= 0 ||
      any(config$refinement_replicates < 1L) ||
      any(config$refinement_replicates > config$replicates)) {
    stop("Domain-expansion and grid-refinement settings are invalid.")
  }
  if (!is.data.frame(config$scenarios) ||
      !all(c("id", "label", "family", "parameter", "regime") %in%
           names(config$scenarios)) ||
      anyDuplicated(config$scenarios$id)) {
    stop("scenarios must be a data frame with unique scenario identifiers.")
  }
  invisible(config)
}

write_csv_atomic <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(basename(path), ".temporary-"),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(object, temporary, row.names = FALSE, na = "")
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) stop("Could not complete output: ", path)
  invisible(path)
}

save_rds_atomic_comparison <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(basename(path), ".temporary-"),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) stop("Could not complete output: ", path)
  invisible(path)
}

comparison_seed <- function(config, scenario_index, replicate) {
  seed <- as.double(config$seed_base) +
    as.double(scenario_index) * 100000 +
    as.double(replicate)
  if (seed > .Machine$integer.max) stop("Configured seed exceeds the R integer range.")
  as.integer(seed)
}

numeric_fingerprint <- function(...) {
  values <- unlist(list(...), use.names = FALSE)
  if (!length(values) || any(!is.finite(values))) {
    stop("A dataset fingerprint requires finite numeric values.")
  }
  index <- seq_along(values)
  paste(
    formatC(sum(values), digits = 16L, format = "fg", flag = "#"),
    formatC(sum(values^2), digits = 16L, format = "fg", flag = "#"),
    formatC(sum(values * ((index %% 1009L) + 1L)), digits = 16L, format = "fg", flag = "#"),
    length(values),
    sep = ":"
  )
}

configuration_text <- function(config) {
  paste(capture.output(dput(config)), collapse = "\n")
}

configuration_fingerprint <- function(config) {
  raw <- charToRaw(configuration_text(config))
  modulus <- 2147483629
  state <- 104729
  for (value in as.integer(raw)) {
    state <- (state * 257 + value + 1) %% modulus
  }
  sprintf("%08x", as.integer(state))
}

run_parallel <- function(items, function_to_run, workers) {
  workers <- min(as.integer(workers), length(items))
  if (workers <= 1L || .Platform$OS.type == "windows") {
    return(lapply(items, function_to_run))
  }
  parallel::mclapply(
    items,
    function_to_run,
    mc.cores = workers,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
}

assert_finite_frame <- function(
    frame,
    context,
    allow_infinite_columns = character()) {
  numeric_columns <- vapply(frame, is.numeric, logical(1))
  numeric_columns[names(numeric_columns) %in% allow_infinite_columns] <- FALSE
  bad <- vapply(
    frame[numeric_columns],
    function(column) any(!is.finite(column) & !is.na(column)),
    logical(1)
  )
  if (any(bad)) {
    stop(context, " contains non-finite numeric values in: ",
         paste(names(bad)[bad], collapse = ", "))
  }
  invisible(frame)
}
