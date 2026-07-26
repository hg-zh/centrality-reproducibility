#!/usr/bin/env Rscript

script_argument <- grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)
if (!length(script_argument)) stop("This script must be run with Rscript.")
script_path <- normalizePath(
  sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE
)
module_root <- dirname(script_path)

available <- c("colon", "temperature", "mnist", "taxi")
arguments <- commandArgs(trailingOnly = TRUE)

if (any(arguments %in% c("-h", "--help"))) {
  cat(
    "Usage: Rscript --vanilla dataIllustrations/run_all.R ",
    "[--analyses=colon,temperature,mnist,taxi]\n",
    sep = ""
  )
  quit(status = 0L)
}

unknown_arguments <- arguments[!grepl("^--analyses=", arguments)]
if (length(unknown_arguments)) {
  stop("Unknown argument(s): ", paste(unknown_arguments, collapse = ", "))
}

analysis_argument <- arguments[grepl("^--analyses=", arguments)]
if (length(analysis_argument) > 1L) {
  stop("Specify --analyses at most once.")
}

analyses <- if (length(analysis_argument)) {
  value <- sub("^--analyses=", "", analysis_argument[[1L]])
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
} else {
  available
}

if (!length(analyses) || any(!nzchar(analyses))) {
  stop("--analyses must contain at least one analysis name.")
}
unknown_analyses <- setdiff(analyses, available)
if (length(unknown_analyses)) {
  stop(
    "Unknown analysis name(s): ", paste(unknown_analyses, collapse = ", "),
    ". Choose from: ", paste(available, collapse = ", "), "."
  )
}
analyses <- unique(analyses)

input_files <- list(
  colon = c(
    "data/colon/I2000.txt",
    "data/colon/names2.txt",
    "data/colon/tissues.txt"
  ),
  temperature = "data/temperature/tampa.csv",
  mnist = c(
    "data/mnist/train-images.idx3-ubyte",
    "data/mnist/train-labels.idx1-ubyte"
  ),
  taxi = "data/taxi/taxi2019gl.rds"
)
scripts <- c(
  colon = "colon.R",
  temperature = "temp.R",
  mnist = "mnist.R",
  taxi = "taxi.R"
)

missing_inputs <- unlist(lapply(analyses, function(analysis) {
  paths <- file.path(module_root, input_files[[analysis]])
  missing <- input_files[[analysis]][!file.exists(paths)]
  if (length(missing)) paste0(analysis, ": ", missing) else character()
}), use.names = FALSE)

if (length(missing_inputs)) {
  stop(
    "Required input files are missing:\n- ",
    paste(missing_inputs, collapse = "\n- "),
    "\nSee dataIllustrations/README.md for data access and placement."
  )
}

rscript <- file.path(R.home("bin"), "Rscript")
for (analysis in analyses) {
  message("Running ", analysis, " analysis...")
  status <- system2(
    rscript,
    c("--vanilla", shQuote(file.path(module_root, scripts[[analysis]])))
  )
  if (!identical(status, 0L)) {
    stop("The ", analysis, " analysis exited with status ", status, ".")
  }
}

message("Requested real-data analyses completed successfully.")
