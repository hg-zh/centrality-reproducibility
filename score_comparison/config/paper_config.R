comparison_config <- function() {
  list(
    profile = "paper",
    replicates = 200L,
    seed_base = 24072026L,
    alpha = 0.05,
    n_train = 200L,
    n_calibration = 200L,
    n_test = 2000L,
    candidate_grid = seq(-10, 10, by = 0.01),
    tau_grid = seq(0, 20, by = 1),
    tau_penalty = 0.0001,
    chunk_size = 500L,
    outside_margin_points = 10L,
    domain_expansion_fraction = 0.25,
    domain_max_expansions = 8L,
    refinement_replicates = c(1L, 50L, 100L, 150L, 200L),
    refinement_factor = 2L,
    length_tolerance = 0.01,
    scenarios = data.frame(
      id = c(
        "normal", "laplace", "student_t3",
        "mixture_075", "mixture_125", "mixture_250"
      ),
      label = c(
        "Normal", "Laplace (variance 1)", "Student t3 (variance 1)",
        "Gaussian mixture, delta=0.75",
        "Gaussian mixture, delta=1.25",
        "Gaussian mixture, delta=2.50"
      ),
      family = c(
        "normal", "laplace", "student_t",
        "normal_mixture", "normal_mixture", "normal_mixture"
      ),
      parameter = c(NA, NA, 3, 0.75, 1.25, 2.50),
      regime = c(
        "unimodal", "unimodal", "unimodal",
        "transition", "bimodal", "bimodal"
      ),
      stringsAsFactors = FALSE
    ),
    population = list(
      primary = list(
        reference_n = 1201L,
        tuning_reference_n = 301L,
        grid_limits = c(-10, 10),
        grid_step = 0.01
      ),
      refined = list(
        reference_n = 1801L,
        tuning_reference_n = 401L,
        grid_limits = c(-10, 10),
        grid_step = 0.005
      )
    )
  )
}
