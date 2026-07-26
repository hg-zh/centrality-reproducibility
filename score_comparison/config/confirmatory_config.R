comparison_config <- function() {
  environment <- new.env(parent = baseenv())
  sys.source(
    file.path(.comparison_module_root, "config", "paper_config.R"),
    envir = environment
  )
  config <- environment$comparison_config()
  config$profile <- "confirmatory"
  config$seed_base <- 25072026L
  config
}
