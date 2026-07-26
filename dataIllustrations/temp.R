#!/usr/bin/env Rscript

# Tampa temperature application
# Reads daily station records, constructs one summer minimum-temperature
# distribution per year, computes C-scores, and recreates the paper figure.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  script_file <- sub("^--file=", "", script_arg[1L])
  script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
  dirname(normalizePath(script_file))
} else {
  normalizePath(getwd())
}
source(file.path(script_dir, "tolReg.R"))

required_packages <- c("ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install the following R packages before running this script: ",
       paste(missing_packages, collapse = ", "))
}

data_dir <- file.path(script_dir, "data", "temperature")
output_dir <- file.path(script_dir, "results", "temperature")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

temperature <- read.csv(file.path(data_dir, "tampa.csv"))
temperature$DATE <- as.Date(temperature$DATE)
if (!all(c("DATE", "TMIN") %in% names(temperature))) {
  stop("The Tampa data must contain DATE and TMIN columns.")
}

years <- 1939:2023
summer_distributions <- lapply(years, function(year) {
  start <- as.Date(sprintf("%d-06-01", year))
  end <- as.Date(sprintf("%d-09-01", year))
  temperature$TMIN[temperature$DATE >= start & temperature$DATE < end]
})
names(summer_distributions) <- years

# Two transcription errors in the source record are corrected as in the
# analysis used for the manuscript.
summer_distributions[["2019"]][62] <- 90
summer_distributions[["2019"]][c(26, 27)] <- 78

if (any(vapply(summer_distributions, anyNA, logical(1))) ||
    length(unique(lengths(summer_distributions))) != 1L) {
  stop("Each summer must contain the same number of nonmissing daily records.")
}
summer_distributions <- lapply(summer_distributions, sort)

fit <- tolReg(
  summer_distributions,
  optns = list(type = "measure")
)

score_table <- data.frame(
  year = years,
  c_score = fit$transRank,
  rank_from_periphery = rank(fit$transRank, ties.method = "first")
)
write.csv(score_table, file.path(output_dir, "temperature_c_scores.csv"),
          row.names = FALSE)

profile_distance <- getPairDist(
  fit$distProf,
  optns = list(type = "measure")
)
profile_coordinates <- stats::cmdscale(profile_distance)
object_coordinates <- stats::cmdscale(fit$pairDist)

density_list <- lapply(
  summer_distributions, stats::density, from = 50, to = 90, adjust = 1
)
density_frame <- data.frame(
  temperature = unlist(lapply(density_list, `[[`, "x"), use.names = FALSE),
  density = unlist(lapply(density_list, `[[`, "y"), use.names = FALSE),
  year = rep(years, each = length(density_list[[1L]]$x)),
  c_score = rep(fit$transRank, each = length(density_list[[1L]]$x))
)
most_peripheral <- order(fit$transRank)[1:2]

p_density <- ggplot2::ggplot(
  density_frame,
  ggplot2::aes(
    x = temperature, y = density, color = c_score,
    group = year
  )
) +
  ggplot2::geom_line(
    alpha = rep(
      ifelse(density_frame$year %in% years[most_peripheral], 1, 0.3)
    )
  ) +
  ggplot2::geom_label(
    data = data.frame(
      temperature = c(82.5, 80.8),
      density = c(0.075, 0.135),
      year = years[most_peripheral]
    ),
    ggplot2::aes(x = temperature, y = density, label = year),
    inherit.aes = FALSE, size = 3
  ) +
  ggplot2::scale_color_gradient(
    name = "Centrality score", low = "coral", high = "blue",
    breaks = c(0.495, 0.5)
  ) +
  ggplot2::labs(
    x = "Daily minimum temperature", y = "Density", tag = "(a)"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 20),
    legend.position = "top",
    legend.key.width = grid::unit(1, "cm")
  )

embedding_plot <- function(coordinates, tag) {
  ggplot2::ggplot(
    data.frame(
      x = coordinates[, 1L], y = coordinates[, 2L],
      year = years, c_score = fit$transRank
    ),
    ggplot2::aes(x = x, y = y, label = year, color = c_score)
  ) +
    ggplot2::geom_text() +
    ggplot2::scale_color_gradient(
      name = "Centrality score", low = "coral", high = "blue",
      breaks = c(0.495, 0.5)
    ) +
    ggplot2::labs(x = "Coordinate 1", y = "Coordinate 2", tag = tag) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 20),
      legend.position = "top",
      legend.key.width = grid::unit(1, "cm")
    )
}

p_profile <- embedding_plot(profile_coordinates, "(b)")
p_object <- embedding_plot(object_coordinates, "(c)")
figure <- p_density + p_profile + p_object + patchwork::guide_area() +
  patchwork::plot_layout(guides = "collect", widths = c(5, 5, 5, 3)) &
  ggplot2::theme(legend.position = "right")

grDevices::pdf(file.path(output_dir, "temperature.pdf"), width = 15, height = 4)
print(figure)
grDevices::dev.off()
message("Temperature analysis complete. Results: ", output_dir)
