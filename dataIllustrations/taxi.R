#!/usr/bin/env Rscript

# Manhattan taxi-network application
# Reads daily 2019 graph Laplacians, normalizes edge weights, computes C-scores,
# recreates profile MDS, and reproduces the DBSCAN cluster summaries.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  script_file <- sub("^--file=", "", script_arg[1L])
  script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
  dirname(normalizePath(script_file))
} else {
  normalizePath(getwd())
}
source(file.path(script_dir, "tolReg.R"))

required_packages <- c("ggplot2", "dbscan")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install the following R packages before running this script: ",
       paste(missing_packages, collapse = ", "))
}

data_dir <- file.path(script_dir, "data", "taxi")
output_dir <- file.path(script_dir, "results", "taxi")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

networks <- readRDS(file.path(data_dir, "taxi2019gl.rds"))
dates <- as.Date(names(networks))
if (length(networks) != 365L || anyNA(dates) ||
    any(vapply(networks, function(x) !identical(dim(x), c(13L, 13L)),
               logical(1)))) {
  stop("Unexpected structure in the 2019 taxi-network data.")
}

week_day_levels <- c(
  "Monday", "Tuesday", "Wednesday", "Thursday",
  "Friday", "Saturday", "Sunday"
)
week_day_number <- as.POSIXlt(dates, tz = "UTC")$wday
week_day_name <- c(
  "Sunday", "Monday", "Tuesday", "Wednesday",
  "Thursday", "Friday", "Saturday"
)[week_day_number + 1L]
week_day <- factor(week_day_name, levels = week_day_levels)
weekend <- week_day_number %in% c(0L, 6L)

normalized_networks <- lapply(networks, function(network) {
  scale_value <- max(abs(network[upper.tri(network)]))
  if (!is.finite(scale_value) || scale_value <= 0) {
    stop("Encountered a network with no positive off-diagonal scale.")
  }
  network / scale_value
})

fit <- tolReg(
  normalized_networks,
  optns = list(type = "network")
)

n_outliers <- floor(length(networks) * 0.05)
outlier_index <- order(fit$transRank)[seq_len(n_outliers)]
score_table <- data.frame(
  date = dates,
  day_of_week = week_day,
  weekend = weekend,
  c_score = fit$transRank,
  rank_from_periphery = rank(fit$transRank, ties.method = "first"),
  peripheral_5_percent = seq_along(networks) %in% outlier_index
)
write.csv(score_table, file.path(output_dir, "taxi_c_scores.csv"),
          row.names = FALSE)

profile_distance <- getPairDist(
  fit$distProf,
  optns = list(
    type = "measure",
    metric_weights = c(
      0.5, rep(1, ncol(fit$distProf) - 2L), 0.5
    ) / ncol(fit$distProf)
  )
)
profile_cluster <- dbscan::dbscan(
  stats::as.dist(profile_distance), eps = 0.08, minPts = 3
)$cluster
object_cluster <- dbscan::dbscan(
  stats::as.dist(fit$pairDist), eps = 0.5, minPts = 3
)$cluster

cluster_table <- function(cluster, method) {
  tab <- as.data.frame.matrix(table(cluster, week_day))
  tab$cluster <- rownames(tab)
  tab$method <- method
  rownames(tab) <- NULL
  tab[, c("method", "cluster",
          "Monday", "Tuesday", "Wednesday", "Thursday",
          "Friday", "Saturday", "Sunday")]
}
write.csv(
  rbind(
    cluster_table(profile_cluster, "Distance profiles"),
    cluster_table(object_cluster, "Original networks")
  ),
  file.path(output_dir, "taxi_dbscan_clusters.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    date = dates,
    day_of_week = week_day,
    profile_cluster = profile_cluster,
    object_cluster = object_cluster
  ),
  file.path(output_dir, "taxi_cluster_assignments.csv"),
  row.names = FALSE
)

profile_coordinates <- stats::cmdscale(profile_distance)
plot_frame <- data.frame(
  x = profile_coordinates[, 1L],
  y = profile_coordinates[, 2L],
  c_score = fit$transRank,
  weekend = factor(weekend, levels = c(FALSE, TRUE),
                   labels = c("Weekday", "Weekend"))
)

figure <- ggplot2::ggplot(
  plot_frame,
  ggplot2::aes(x = x, y = y, shape = weekend, color = c_score)
) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::scale_color_gradient(
    name = "Metric\ncentrality", low = "coral", high = "blue",
    breaks = c(0.495, 0.5),
    guide = ggplot2::guide_colorbar(order = 1, title.position = "left")
  ) +
  ggplot2::scale_shape_discrete(
    name = NULL, guide = ggplot2::guide_legend(order = 2)
  ) +
  ggplot2::labs(x = "Coordinate 1", y = "Coordinate 2") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 20),
    legend.position = "top",
    legend.key.width = grid::unit(1, "cm")
  )

grDevices::pdf(
  file.path(output_dir, "taxi_profile_mds.pdf"), width = 10, height = 6
)
print(figure)
grDevices::dev.off()
message("Taxi analysis complete. Results: ", output_dir)
