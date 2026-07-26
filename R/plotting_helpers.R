require_release_plot_packages <- function() {
  packages <- c("ggplot2", "cowplot", "KernSmooth", "scales")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required package(s): ", paste(missing, collapse = ", "))
  }
  invisible(packages)
}

distance_profile_values <- function(point, training, distance_grid, metric = "euclidean") {
  distances <- as.numeric(pairwise_distance(point, training, metric))
  vapply(distance_grid, function(radius) mean(distances <= radius), numeric(1))
}

profile_plot_theme <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 10)
    )
}

smooth_profile <- function(x, y, bandwidth) {
  result <- KernSmooth::locpoly(x, y, bandwidth = bandwidth, degree = 1)
  data.frame(x = result$x, y = result$y)
}

transport_curve <- function(
    source_cdf,
    target_quantile,
    distance_grid,
    probability_grid,
    first_bandwidth,
    second_bandwidth,
    cap_direction = c("none", "above_identity"),
    sign = c("none", "positive", "negative"),
    scale = 1) {
  cap_direction <- match.arg(cap_direction)
  sign <- match.arg(sign)
  transported <- stats::approx(
    probability_grid,
    target_quantile,
    xout = source_cdf,
    rule = 2
  )$y
  if (identical(cap_direction, "above_identity")) {
    transported[transported >= distance_grid] <-
      distance_grid[transported >= distance_grid]
  }
  first <- smooth_profile(distance_grid, transported, first_bandwidth)
  second <- smooth_profile(
    c(0, first$x),
    c(0, first$y),
    second_bandwidth
  )
  second$y[1L] <- 0
  second$y[nrow(second)] <- max(distance_grid)
  displacement <- scale * (second$y - second$x)
  if (identical(sign, "positive")) displacement[displacement < 0] <- 0
  if (identical(sign, "negative")) displacement[displacement > 0] <- 0
  data.frame(x = second$x, y = displacement)
}

reproduce_figure1 <- function(settings, output_directory) {
  cfg <- settings$figure1
  set_reproducible_rng(cfg$seed)
  training <- t(vapply(
    seq_len(cfg$n),
    function(i) MASS::mvrnorm(1L, cfg$mean, cfg$covariance),
    numeric(2)
  ))
  data <- data.frame(x = training[, 1L], y = training[, 2L])
  colors <- c(center = "red", middle = "purple", far = "blue")
  selected <- cfg$selected_points[c("center", "middle", "far"), , drop = FALSE]

  scatter_panels <- lapply(seq_len(nrow(selected)), function(i) {
    point <- data.frame(x = selected[i, 1L], y = selected[i, 2L])
    ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(size = 0.3) +
      ggplot2::geom_point(
        data = point, color = "gray", alpha = 0.5, size = 10
      ) +
      ggplot2::geom_point(
        data = point, color = "gray", alpha = 0.3, size = 20
      ) +
      ggplot2::geom_point(
        data = point, color = "gray", alpha = 0.2, size = 30
      ) +
      ggplot2::geom_point(data = point, color = unname(colors[i])) +
      ggplot2::coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
      profile_plot_theme()
  })

  profiles <- t(vapply(
    seq_len(nrow(selected)),
    function(i) distance_profile_values(
      selected[i, ], training, cfg$distance_grid, "euclidean"
    ),
    numeric(length(cfg$distance_grid))
  ))
  cdf_smooth <- lapply(seq_len(nrow(profiles)), function(i) {
    value <- smooth_profile(cfg$distance_grid, profiles[i, ], 0.8)
    value$y[value$y < 0] <- 0
    value
  })
  cdf_panels <- lapply(seq_along(cdf_smooth), function(i) {
    ggplot2::ggplot(cdf_smooth[[i]], ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_line(color = unname(colors[i])) +
      ggplot2::labs(x = "r", y = "") +
      profile_plot_theme()
  })

  density_smooth <- lapply(seq_len(nrow(profiles)), function(i) {
    increments <- c(0, diff(profiles[i, ]))
    value <- smooth_profile(cfg$distance_grid, increments, 0.8)
    value$y <- value$y / sum(value$y) * length(value$x) /
      diff(range(cfg$distance_grid))
    if (i == 3L) value$y[value$y < 0] <- 0
    value
  })
  density_panels <- lapply(seq_along(density_smooth), function(i) {
    ggplot2::ggplot(density_smooth[[i]], ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_line(color = unname(colors[i])) +
      ggplot2::labs(x = "r", y = "") +
      profile_plot_theme()
  })

  probability_grid <- make_probability_grid(
    cfg$probability_grid_size,
    include_endpoints = TRUE
  )
  adjusted_profiles <- profiles
  adjusted_profiles[, 1L] <- 0
  quantiles <- t(vapply(
    seq_len(nrow(adjusted_profiles)),
    function(i) {
      value <- stats::approx(
        adjusted_profiles[i, ],
        cfg$distance_grid,
        xout = probability_grid,
        rule = 2,
        ties = min
      )$y
      value[1L] <- 0
      value[length(value)] <- max(cfg$distance_grid)
      value
    },
    numeric(length(probability_grid))
  ))

  center_to_far <- transport_curve(
    adjusted_profiles[1L, ], quantiles[3L, ],
    cfg$distance_grid, probability_grid, 1.2, 0.5,
    sign = "positive"
  )
  center_to_middle <- transport_curve(
    adjusted_profiles[1L, ], quantiles[2L, ],
    cfg$distance_grid, probability_grid, 1.2, 0.5,
    sign = "positive", scale = 0.5
  )
  middle_to_far <- transport_curve(
    adjusted_profiles[2L, ], quantiles[3L, ],
    cfg$distance_grid, probability_grid, 0.5, 0.7
  )
  middle_to_center <- transport_curve(
    adjusted_profiles[2L, ], quantiles[1L, ],
    cfg$distance_grid, probability_grid, 0.5, 0.7,
    cap_direction = "above_identity", sign = "negative", scale = 2
  )
  far_to_center <- transport_curve(
    adjusted_profiles[3L, ], quantiles[1L, ],
    cfg$distance_grid, probability_grid, 0.5, 0.7,
    cap_direction = "above_identity", sign = "negative"
  )
  far_to_middle <- transport_curve(
    adjusted_profiles[3L, ], quantiles[2L, ],
    cfg$distance_grid, probability_grid, 0.5, 0.7,
    cap_direction = "above_identity", sign = "negative"
  )

  displacement_panel <- function(first, second, first_color, second_color) {
    ggplot2::ggplot(first, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 2, color = "gray") +
      ggplot2::geom_line(color = first_color) +
      ggplot2::geom_line(
        data = second,
        ggplot2::aes(x = x, y = y),
        color = second_color
      ) +
      ggplot2::labs(x = "r", y = " ") +
      profile_plot_theme()
  }
  transport_panels <- list(
    displacement_panel(center_to_far, center_to_middle, "blue", "purple"),
    displacement_panel(middle_to_far, middle_to_center, "blue", "red"),
    displacement_panel(far_to_center, far_to_middle, "red", "purple")
  )

  figure <- cowplot::plot_grid(
    plotlist = c(scatter_panels, cdf_panels, density_panels, transport_panels),
    nrow = 4
  )
  ggplot2::ggsave(
    file.path(output_directory, cfg$output),
    figure,
    width = 6,
    height = 7,
    device = grDevices::pdf
  )
  invisible(figure)
}

reproduce_figure2 <- function(settings, output_directory) {
  cfg <- settings$figure2
  set_reproducible_rng(cfg$seed)
  pair <- generate_gaussian_mixture_pair(
    cfg$n_train,
    cfg$mixing,
    cfg$means,
    cfg$covariances
  )
  search <- make_euclidean_search_grid(cfg$grid_n, cfg$grid_limits)
  probability_grid <- make_probability_grid(
    cfg$probability_grid_size,
    include_endpoints = TRUE
  )
  fit <- fit_centrality_score(
    pair$training,
    "euclidean",
    cfg$distance_grid,
    probability_grid
  )
  training_scores <- centrality_scores_from_quantiles(
    fit$training_quantiles, fit, cfg$tau
  )
  calibration_scores <- centrality_scores(pair$calibration, fit, cfg$tau)
  search_scores <- centrality_scores(search$observations, fit, cfg$tau)
  thresholds <- calibration_thresholds(
    calibration_scores,
    settings$alpha,
    cfg$n_train,
    settings$threshold$method,
    settings$threshold$quantile_type
  )
  membership <- sweep(search_scores, 2L, thresholds, ">=")
  training_data <- data.frame(
    x = pair$training[, 1L],
    y = pair$training[, 2L]
  )

  panels <- lapply(seq_along(cfg$tau), function(index) {
    scored <- transform(training_data, score = training_scores[, index])
    region <- search$data[membership[, index], , drop = FALSE]
    ggplot2::ggplot(scored, ggplot2::aes(x = x, y = y, color = score)) +
      ggplot2::geom_point(size = 0.4) +
      ggplot2::scale_color_gradient(low = "blue", high = "coral", name = NULL) +
      ggplot2::geom_point(
        data = region,
        ggplot2::aes(x = x, y = y),
        inherit.aes = FALSE,
        color = "gray",
        alpha = 0.03
      ) +
      ggplot2::scale_x_continuous(
        limits = cfg$grid_limits,
        expand = c(0, 0)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(-4.1, 5),
        expand = c(0, 0)
      ) +
      ggplot2::geom_hline(yintercept = 4.1, color = "black", linewidth = 1) +
      ggplot2::annotate(
        "rect",
        xmin = cfg$grid_limits[1L],
        xmax = cfg$grid_limits[2L],
        ymin = 4.1,
        ymax = 5,
        fill = "grey80"
      ) +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 4.78,
        label = paste0("tau == ", cfg$tau[index]),
        parse = TRUE,
        vjust = 1,
        fontface = "bold",
        size = 3
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "x", y = "y") +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(
          colour = "black", fill = NA, linewidth = 1
        ),
        panel.background = ggplot2::element_blank(),
        panel.spacing = grid::unit(0, "lines"),
        plot.margin = ggplot2::margin(5, 2, 5, 2, unit = "pt"),
        axis.title = ggplot2::element_text(size = 8),
        axis.text = ggplot2::element_text(size = 8),
        legend.text = ggplot2::element_text(size = 8),
        legend.key.width = grid::unit(0.3, "cm"),
        legend.key.height = grid::unit(0.6, "cm")
      )
  })
  legend <- cowplot::get_legend(panels[[1L]])
  panels <- lapply(panels, function(panel) {
    panel + ggplot2::theme(legend.position = "none")
  })
  figure <- cowplot::plot_grid(
    plotlist = c(panels, list(legend)),
    nrow = 1,
    rel_widths = c(rep(1, length(panels)), 0.3)
  )
  ggplot2::ggsave(
    file.path(output_directory, cfg$output),
    figure,
    width = 12.8 / 2.54,
    height = 17.8 / (3.2 * 2.54),
    device = grDevices::pdf
  )
  invisible(figure)
}

reproduce_figures_1_2 <- function(settings, project_root) {
  require_release_plot_packages()
  output_directory <- release_output_directory(
    project_root, "figures", settings$profile
  )
  reproduce_figure1(settings, output_directory)
  reproduce_figure2(settings, output_directory)
  invisible(c(
    file.path(output_directory, settings$figure1$output),
    file.path(output_directory, settings$figure2$output)
  ))
}

sample_size_labels <- function(values) {
  paste0("n=", values)
}

sample_colors <- function(values) {
  palette <- c("red", "blue", "green")
  labels <- sample_size_labels(values)
  stats::setNames(palette[seq_along(labels)], labels)
}

prepare_metric_summary <- function(result) {
  data <- result$metric_summary
  levels <- sample_size_labels(result$settings$sample_sizes)
  data$sample_label <- factor(
    sample_size_labels(data$sample_size),
    levels = levels
  )
  data
}

minimal_size_points <- function(data) {
  pieces <- split(data, data$sample_label)
  do.call(rbind, lapply(pieces, function(piece) {
    piece[which.min(piece$size), , drop = FALSE]
  }))
}

line_panel_theme <- function(base_size = 18) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        colour = "black", fill = NA, linewidth = 1
      ),
      panel.background = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

plot_figure3 <- function(result, output_directory) {
  data <- prepare_metric_summary(result)
  colors <- sample_colors(result$settings$sample_sizes)
  manuscript <- identical(result$profile, "manuscript")

  coverage <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = tau, y = coverage, color = sample_label)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(
      yintercept = 0.95, linetype = "dashed", color = "black"
    ) +
    ggplot2::labs(
      x = "tuning parameter", y = "coverage", color = "sample size"
    ) +
    ggplot2::scale_color_manual(values = colors) +
    line_panel_theme() +
    ggplot2::theme(legend.position = "none")
  if (manuscript) {
    coverage <- coverage +
      ggplot2::scale_x_continuous(
        limits = c(0, 100),
        expand = ggplot2::expansion(mult = 0.04)
      ) +
      ggplot2::coord_cartesian(ylim = c(0.925, 0.975))
  }

  minima <- minimal_size_points(data)
  size <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = tau, y = size, color = sample_label)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(
      data = minima,
      ggplot2::aes(x = tau, y = size),
      color = "black",
      size = 2,
      shape = 17
    ) +
    ggplot2::labs(
      x = "tuning parameter", y = "size", color = "sample size"
    ) +
    ggplot2::scale_color_manual(values = colors) +
    line_panel_theme()
  if (manuscript) {
    size <- size + ggplot2::coord_cartesian(ylim = c(2.5, 3))
  }
  legend <- cowplot::get_legend(
    size +
      ggplot2::theme(
        legend.position = "right",
        legend.text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 12)
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(keywidth = grid::unit(1.5, "cm")))
  )
  legend_canvas <- cowplot::ggdraw(legend)
  size <- size + ggplot2::theme(legend.position = "none")
  line_figure <- cowplot::plot_grid(
    coverage,
    size,
    legend_canvas,
    ncol = 3,
    rel_widths = c(1, 1, 0.35)
  )
  ggplot2::ggsave(
    file.path(output_directory, "mix_2d_1.pdf"),
    line_figure,
    width = 600 / 72,
    height = 262 / 72,
    device = grDevices::pdf
  )

  heat <- result$heatmap
  heat$mixing_label <- factor(
    paste0("a=", format(heat$mixing, trim = TRUE)),
    levels = paste0("a=", format(result$settings$heat_mixing, trim = TRUE))
  )
  heat$tau_label <- factor(
    paste0("tau=", format(heat$tau, trim = TRUE)),
    levels = paste0("tau=", format(result$settings$heat_tau, trim = TRUE))
  )
  heat_figure <- ggplot2::ggplot(
    heat,
    ggplot2::aes(x = x, y = y, fill = probability)
  ) +
    ggplot2::geom_tile(width = 1, height = 1) +
    ggplot2::scale_fill_gradient(
      low = "#FFF7BC",
      high = "#2C7FB8",
      limits = c(0, 1),
      oob = scales::squish,
      breaks = seq(0, 1, 0.25),
      labels = scales::percent
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(mixing_label),
      cols = ggplot2::vars(tau_label),
      switch = "y"
    ) +
    ggplot2::labs(fill = NULL, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(size = 18, face = "bold"),
      strip.text.x = ggplot2::element_text(size = 18, face = "bold"),
      legend.position = "right",
      legend.key.width = grid::unit(0.45, "cm"),
      legend.key.height = grid::unit(1.6, "cm"),
      legend.text = ggplot2::element_text(size = 12),
      panel.spacing = grid::unit(0.6, "lines"),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
  ggplot2::ggsave(
    file.path(output_directory, "mix_2d_heatmap_2x3.pdf"),
    heat_figure,
    width = 600 / 72,
    height = 350 / 72,
    device = grDevices::pdf
  )
  invisible(list(lines = line_figure, heatmap = heat_figure))
}

plot_figure4 <- function(sphere, wasserstein, output_directory) {
  sphere_data <- prepare_metric_summary(sphere)
  wasserstein_data <- prepare_metric_summary(wasserstein)
  colors <- sample_colors(sphere$settings$sample_sizes)
  manuscript <- identical(sphere$profile, "manuscript")

  make_coverage <- function(data, y_label) {
    plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = tau, y = coverage, color = sample_label)
    ) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_hline(
        yintercept = 0.95, linetype = "dashed", color = "black"
      ) +
      ggplot2::labs(x = NULL, y = y_label, color = NULL) +
      ggplot2::scale_color_manual(values = colors) +
      line_panel_theme() +
      ggplot2::theme(legend.position = "none")
    if (manuscript) {
      plot <- plot + ggplot2::coord_cartesian(
        xlim = c(0, 20), ylim = c(0.925, 0.975)
      )
    }
    plot
  }
  make_size <- function(data, y_label, sphere_panel = FALSE) {
    plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = tau, y = size, color = sample_label)
    ) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(
        data = minimal_size_points(data),
        ggplot2::aes(x = tau, y = size),
        color = "black",
        size = 2,
        shape = 17
      ) +
      ggplot2::labs(x = expression(tau), y = y_label, color = NULL) +
      ggplot2::scale_color_manual(values = colors) +
      line_panel_theme() +
      ggplot2::theme(legend.position = "none")
    if (manuscript && sphere_panel) {
      plot <- plot + ggplot2::coord_cartesian(ylim = c(2, 3))
    }
    plot
  }

  sphere_coverage <- make_coverage(sphere_data, "Coverage") +
    ggplot2::labs(title = "S2")
  wasserstein_coverage <- make_coverage(wasserstein_data, NULL) +
    ggplot2::labs(title = "W2")
  sphere_size <- make_size(sphere_data, "Size", TRUE)
  wasserstein_size <- make_size(wasserstein_data, NULL, FALSE)

  legend_plot <- ggplot2::ggplot(
    sphere_data,
    ggplot2::aes(x = tau, y = size, color = sample_label)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(color = NULL) +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = ggplot2::element_text(size = 16),
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        keywidth = grid::unit(2, "cm")
      )
    )
  legend <- cowplot::get_legend(legend_plot)
  legend_canvas <- cowplot::ggdraw(legend)
  panels <- cowplot::plot_grid(
    cowplot::plot_grid(
      sphere_coverage, wasserstein_coverage,
      nrow = 1, align = "hv", axis = "tblr"
    ),
    cowplot::plot_grid(
      sphere_size, wasserstein_size,
      nrow = 1, align = "hv", axis = "tblr"
    ),
    ncol = 1,
    rel_heights = c(1, 1)
  )
  line_figure <- cowplot::plot_grid(
    panels,
    legend_canvas,
    ncol = 1,
    rel_heights = c(1, 0.12)
  )
  ggplot2::ggsave(
    file.path(output_directory, "S2_W2_lines_2x2.pdf"),
    line_figure,
    width = 600 / 72,
    height = 600 / 72,
    device = grDevices::pdf
  )

  heat_scale <- ggplot2::scale_fill_gradient(
    low = "#FFF7BC",
    high = "#2C7FB8",
    limits = c(0, 1),
    oob = scales::squish,
    guide = "none"
  )
  heat_theme <- ggplot2::theme_minimal(base_size = 18) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 20),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 12),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
  sphere_heat <- ggplot2::ggplot(
    sphere$heatmap,
    ggplot2::aes(x = phi, y = theta, fill = probability)
  ) +
    ggplot2::geom_tile(width = 0.05, height = 0.05) +
    heat_scale +
    ggplot2::labs(x = expression(phi), y = expression(theta)) +
    heat_theme
  wasserstein_heat <- ggplot2::ggplot(
    wasserstein$heatmap,
    ggplot2::aes(x = mu, y = sigma, fill = probability)
  ) +
    ggplot2::geom_tile(width = 0.055, height = 0.014) +
    ggplot2::coord_cartesian(
      xlim = c(0, 5), ylim = c(0, 1), expand = FALSE
    ) +
    heat_scale +
    ggplot2::labs(x = expression(mu), y = expression(sigma)) +
    heat_theme
  legend_data <- data.frame(
    x = 1,
    y = seq(0.001, 0.999, length.out = 256L),
    probability = seq(0, 1, length.out = 256L)
  )
  legend_plot <- ggplot2::ggplot(
    legend_data,
    ggplot2::aes(x = x, y = y, fill = probability)
  ) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::scale_fill_gradient(
      low = "#FFF7BC", high = "#2C7FB8", guide = "none"
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      position = "right",
      breaks = c(0.001, 0.25, 0.5, 0.75, 0.999),
      labels = c("0%", "25%", "50%", "75%", "100%"),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y.left = ggplot2::element_blank(),
      axis.text.y.right = ggplot2::element_text(size = 12, color = "black"),
      axis.ticks.y.left = ggplot2::element_blank(),
      axis.ticks.y.right = ggplot2::element_line(color = "white"),
      axis.ticks.length.y.right = grid::unit(0.08, "cm"),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 20, 0, 0)
    )
  legend_panel <- cowplot::plot_grid(
    NULL,
    legend_plot,
    NULL,
    ncol = 1,
    rel_heights = c(0.37, 0.26, 0.37)
  )
  heat_column <- cowplot::plot_grid(
    sphere_heat + ggplot2::theme(legend.position = "none"),
    wasserstein_heat + ggplot2::theme(legend.position = "none"),
    ncol = 1,
    align = "v",
    axis = "lr"
  )
  heat_figure <- cowplot::plot_grid(
    heat_column,
    legend_panel,
    nrow = 1,
    rel_widths = c(1, 0.2)
  )
  ggplot2::ggsave(
    file.path(output_directory, "S2W2_heat_2_by_1.pdf"),
    heat_figure,
    width = 490 / 72,
    height = 730 / 72,
    device = grDevices::pdf
  )
  invisible(list(lines = line_figure, heatmap = heat_figure))
}

plot_all_simulations <- function(settings, project_root) {
  require_release_plot_packages()
  result_directory <- file.path(project_root, "results", settings$profile)
  paths <- file.path(
    result_directory,
    c(
      "figure3_euclidean.rds",
      "figure4_sphere.rds",
      "figure4_wasserstein.rds"
    )
  )
  if (any(!file.exists(paths))) {
    stop(
      "Clean simulation results are missing. Run ",
      "scripts/02_run_figures_3_4_simulations.R first."
    )
  }
  results <- lapply(paths, readRDS)
  scenarios <- c("euclidean", "sphere", "wasserstein")
  expected_analysis <- analysis_settings(settings)
  valid <- vapply(seq_along(results), function(index) {
    result <- results[[index]]
    scenario <- scenarios[index]
    identical(result$schema_version, 2L) &&
      identical(result$profile, settings$profile) &&
      identical(result$scenario, scenario) &&
      identical(result$settings, settings[[scenario]]) &&
      identical(result$analysis, expected_analysis) &&
      identical(result$rng, settings$rng)
  }, logical(1))
  if (any(!valid)) {
    stop(
      "A saved result does not match the requested scientific configuration. ",
      "Rerun scripts/02_run_figures_3_4_simulations.R."
    )
  }
  output_directory <- release_output_directory(
    project_root, "figures", settings$profile
  )
  plot_figure3(results[[1L]], output_directory)
  plot_figure4(results[[2L]], results[[3L]], output_directory)
  invisible(output_directory)
}
