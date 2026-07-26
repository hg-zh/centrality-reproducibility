# Core implementation of centrality scores (C-scores)
#
# Random objects are supplied as a list. Supported object types are:
#   euclidean    numeric vectors with the Euclidean metric
#   network      equally sized matrices with the Frobenius metric
#   function     functions evaluated on a common grid (`optns$supp`)
#   measure      empirical one-dimensional distributions with W2 distance
#   compositional unit vectors with the spherical geodesic metric
#
# The main entry point is `tolReg()`. The historical function name is retained
# so that the analysis scripts and earlier project code remain compatible.

expit <- function(x) {
  stats::plogis(x)
}

logit <- function(x) {
  stats::qlogis(x)
}

.validate_options <- function(y, optns) {
  if (is.null(y) || !is.list(y) || length(y) < 2L) {
    stop("`y` must be a list containing at least two random objects.")
  }

  defaults <- list(
    type = "euclidean",
    alpha = 0.05,
    candidates = c(0, 1, 2, 3, 5, 10, 20),
    penalty = 1e-4,
    max_midpoints = 10000L,
    seed = 1L
  )
  for (name in names(defaults)) {
    if (is.null(optns[[name]])) optns[[name]] <- defaults[[name]]
  }

  supported <- c("euclidean", "network", "function", "measure",
                 "compositional", "correlation")
  if (!optns$type %in% supported) {
    stop("Unsupported object type: ", optns$type)
  }
  if (!is.numeric(optns$alpha) || length(optns$alpha) != 1L ||
      optns$alpha <= 0 || optns$alpha >= 1) {
    stop("`optns$alpha` must be a single number strictly between 0 and 1.")
  }
  if (!is.null(optns$tau) &&
      (!is.numeric(optns$tau) || length(optns$tau) != 1L ||
       !is.finite(optns$tau) || optns$tau < 0)) {
    stop("`optns$tau` must be a single nonnegative finite number.")
  }
  if (optns$type == "function") {
    if (is.null(optns$supp) || length(optns$supp) < 2L ||
        is.unsorted(optns$supp, strictly = TRUE)) {
      stop("Functional data require a strictly increasing `optns$supp` grid.")
    }
  }
  optns
}

.trapezoid_weights <- function(grid) {
  dx <- diff(grid)
  c(dx[1L] / 2, (dx[-1L] + dx[-length(dx)]) / 2,
    dx[length(dx)] / 2)
}

.empirical_w2_weights <- function(n) {
  if (n == 1L) return(1)
  c(0.5, rep(1, max(0, n - 2L)), 0.5) / n
}

.quantile_matrix <- function(y, target_size = max(lengths(y))) {
  sizes <- lengths(y)
  if (any(sizes < 1L)) stop("Empirical distributions cannot be empty.")
  probabilities <- ((seq_len(target_size) - 0.5) / target_size)
  t(vapply(
    y,
    function(values) {
      stats::quantile(
        values,
        probs = probabilities,
        type = 1,
        names = FALSE,
        na.rm = FALSE
      )
    },
    numeric(target_size)
  ))
}

.prepare_objects <- function(y, optns) {
  if (any(vapply(y, function(x) any(!is.finite(x)), logical(1)))) {
    stop("All random objects must contain finite numeric values.")
  }

  if (optns$type == "measure") {
    matrix_data <- .quantile_matrix(y)
    metric_weights <- .empirical_w2_weights(ncol(matrix_data))
  } else {
    object_lengths <- lengths(y)
    if (length(unique(object_lengths)) != 1L) {
      stop("All objects must have the same size for type `", optns$type, "`.")
    }
    matrix_data <- matrix(unlist(y, use.names = FALSE),
                          nrow = length(y), byrow = TRUE)
    metric_weights <- rep(1, ncol(matrix_data))

    if (optns$type == "function") {
      if (ncol(matrix_data) != length(optns$supp)) {
        stop("Each function must have one value per point in `optns$supp`.")
      }
      metric_weights <- .trapezoid_weights(optns$supp)
    }
  }

  list(matrix = matrix_data, weights = metric_weights)
}

getPairDist <- function(data, distFun = NULL, optns = list()) {
  data <- as.matrix(data)
  if (nrow(data) < 2L) stop("At least two rows are required.")

  if (optns$type == "compositional") {
    gram <- tcrossprod(data)
    gram <- pmin(1, pmax(-1, gram))
    return(acos(gram))
  }

  weights <- optns$metric_weights
  if (is.null(weights)) {
    if (optns$type == "function") {
      weights <- .trapezoid_weights(optns$supp)
    } else if (optns$type == "measure") {
      weights <- .empirical_w2_weights(ncol(data))
    } else {
      weights <- rep(1, ncol(data))
    }
  }
  if (length(weights) != ncol(data)) {
    stop("Metric weights do not match the number of object coordinates.")
  }

  as.matrix(stats::dist(sweep(data, 2L, sqrt(weights), `*`)))
}

distFun <- function(x, y, optns = list()) {
  if (length(x) != length(y)) stop("Objects have incompatible sizes.")
  if (optns$type == "compositional") {
    return(acos(pmin(1, pmax(-1, sum(x * y)))))
  }

  weights <- optns$metric_weights
  if (is.null(weights)) {
    if (optns$type == "function") {
      weights <- .trapezoid_weights(optns$supp)
    } else if (optns$type == "measure") {
      weights <- .empirical_w2_weights(length(x))
    } else {
      weights <- rep(1, length(x))
    }
  }
  sqrt(sum(weights * (x - y)^2))
}

.cross_distances <- function(reference, query, optns) {
  reference <- as.matrix(reference)
  query <- as.matrix(query)
  if (ncol(reference) != ncol(query)) {
    stop("Reference and query objects have incompatible sizes.")
  }

  if (optns$type == "compositional") {
    gram <- query %*% t(reference)
    return(acos(pmin(1, pmax(-1, gram))))
  }

  weights <- optns$metric_weights
  if (is.null(weights)) {
    if (optns$type == "function") {
      weights <- .trapezoid_weights(optns$supp)
    } else if (optns$type == "measure") {
      weights <- .empirical_w2_weights(ncol(reference))
    } else {
      weights <- rep(1, ncol(reference))
    }
  }
  reference_weighted <- sweep(reference, 2L, sqrt(weights), `*`)
  query_weighted <- sweep(query, 2L, sqrt(weights), `*`)
  squared <- outer(rowSums(query_weighted^2), rowSums(reference_weighted^2), `+`) -
    2 * tcrossprod(query_weighted, reference_weighted)
  squared[squared < 0] <- 0
  sqrt(squared)
}

transRank <- function(pd, optns) {
  pd <- as.matrix(pd)
  if (nrow(pd) != ncol(pd)) stop("`pd` must be a square distance matrix.")
  if (is.null(optns$tau)) stop("A value of `optns$tau` is required.")

  n <- nrow(pd)
  profiles <- t(apply(pd, 1L, sort.int, method = "quick"))
  weights <- exp(-optns$tau * seq_len(n) / n)
  profile_integrals <- as.vector(profiles %*% weights)
  raw_scores <- (sum(profile_integrals) - n * profile_integrals) / n^2
  expit(raw_scores)
}

brct <- function(y, optns) {
  y <- as.matrix(y)
  if (optns$type == "compositional") {
    center <- colMeans(y)
    norm_center <- sqrt(sum(center^2))
    if (norm_center == 0) stop("The spherical barycenter is not unique.")
    center / norm_center
  } else {
    colMeans(y)
  }
}

.scores_for_profiles <- function(reference_profiles, query_profiles, tau) {
  n <- nrow(reference_profiles)
  weights <- exp(-tau * seq_len(n) / n)
  reference_integrals <- as.vector(reference_profiles %*% weights)
  query_integrals <- as.vector(query_profiles %*% weights)
  expit((sum(reference_integrals) - n * query_integrals) / n^2)
}

.select_tau <- function(y_matrix, pair_dist, optns) {
  n <- nrow(y_matrix)
  pair_index <- utils::combn(seq_len(n), 2L)
  if (ncol(pair_index) > optns$max_midpoints) {
    set.seed(optns$seed)
    pair_index <- pair_index[, sample.int(ncol(pair_index),
                                         optns$max_midpoints), drop = FALSE]
  }

  midpoint_matrix <- vapply(
    seq_len(ncol(pair_index)),
    function(j) brct(y_matrix[pair_index[, j], , drop = FALSE], optns),
    numeric(ncol(y_matrix))
  )
  midpoint_matrix <- t(midpoint_matrix)

  midpoint_dist <- .cross_distances(y_matrix, midpoint_matrix, optns)

  reference_profiles <- t(apply(pair_dist, 1L, sort.int, method = "quick"))
  midpoint_profiles <- t(apply(midpoint_dist, 1L, sort.int, method = "quick"))
  criterion <- vapply(
    optns$candidates,
    function(tau) {
      reference_scores <- .scores_for_profiles(
        reference_profiles, reference_profiles, tau
      )
      midpoint_scores <- .scores_for_profiles(
        reference_profiles, midpoint_profiles, tau
      )
      cutoff <- as.numeric(stats::quantile(
        reference_scores, probs = optns$alpha, names = FALSE, type = 1
      ))
      mean(midpoint_scores >= cutoff) + optns$penalty * tau
    },
    numeric(1)
  )

  list(
    tau = optns$candidates[which.min(criterion)],
    criterion = stats::setNames(criterion, optns$candidates)
  )
}

tolReg <- function(y = NULL, optns = list()) {
  optns <- .validate_options(y, optns)
  prepared <- .prepare_objects(y, optns)
  optns$metric_weights <- prepared$weights
  pair_dist <- getPairDist(prepared$matrix, optns = optns)

  if (is.null(optns$tau)) {
    selected <- .select_tau(prepared$matrix, pair_dist, optns)
    optns$tau <- selected$tau
    optns$selection_criterion <- selected$criterion
  }

  scores <- transRank(pair_dist, optns)
  n <- length(y)
  central_indices <- order(scores, decreasing = TRUE)[
    seq_len(ceiling(n * (1 - optns$alpha)))
  ]

  structure(
    list(
      confSet = y[central_indices],
      confSetIndex = central_indices,
      size = min(scores[central_indices]),
      barycenter = brct(prepared$matrix, optns),
      transRank = scores,
      pairDist = pair_dist,
      distProf = t(apply(pair_dist, 1L, sort.int, method = "quick")),
      yM = prepared$matrix,
      optns = optns
    ),
    class = "cscore_fit"
  )
}

test <- function(fit, yTest = NULL) {
  if (!inherits(fit, "cscore_fit")) stop("`fit` must be returned by `tolReg()`.")
  if (is.null(yTest) || !is.list(yTest)) {
    stop("`yTest` must be a nonempty list.")
  }

  optns <- fit$optns
  if (optns$type == "measure") {
    test_matrix <- .quantile_matrix(yTest, target_size = ncol(fit$yM))
    prepared_test <- list(
      matrix = test_matrix,
      weights = .empirical_w2_weights(ncol(test_matrix))
    )
  } else {
    prepared_test <- .prepare_objects(yTest, optns)
  }
  if (ncol(prepared_test$matrix) != ncol(fit$yM)) {
    stop("Training and test objects have incompatible representations.")
  }

  test_dist <- .cross_distances(fit$yM, prepared_test$matrix, optns)
  test_profiles <- t(apply(test_dist, 1L, sort.int, method = "quick"))
  test_scores <- .scores_for_profiles(fit$distProf, test_profiles, optns$tau)
  cutoff <- as.numeric(stats::quantile(
    fit$transRank, probs = optns$alpha, names = FALSE, type = 1
  ))

  list(
    inout = test_scores >= cutoff,
    scores = test_scores,
    cutoff = cutoff,
    optns = optns
  )
}

# Classical MDS coordinates for a Euclidean feature representation. This is
# equivalent to applying `cmdscale()` to its Euclidean distance matrix, while
# avoiding formation and eigendecomposition of a second n-by-n matrix.
feature_mds <- function(features, k = 2L, weights = NULL) {
  features <- as.matrix(features)
  if (is.null(weights)) weights <- rep(1, ncol(features))
  weighted <- sweep(features, 2L, sqrt(weights), `*`)
  centered <- scale(weighted, center = TRUE, scale = FALSE)

  if (requireNamespace("RSpectra", quietly = TRUE) &&
      min(dim(centered)) > k + 1L) {
    decomposition <- RSpectra::svds(centered, k = k)
    coordinates <- sweep(decomposition$u, 2L, decomposition$d, `*`)
  } else {
    decomposition <- stats::prcomp(centered, center = FALSE,
                                   scale. = FALSE, rank. = k)
    coordinates <- decomposition$x[, seq_len(k), drop = FALSE]
  }
  colnames(coordinates) <- paste0("Coordinate ", seq_len(k))
  coordinates
}
