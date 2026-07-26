#!/usr/bin/env Rscript

# Colon gene-expression application
# Reads and preprocesses the Alon et al. colon data, constructs one optimal
# transport map per gene, computes C-scores, and recreates the paper figure.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  script_file <- sub("^--file=", "", script_arg[1L])
  script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
  dirname(normalizePath(script_file))
} else {
  normalizePath(getwd())
}

source(file.path(script_dir, "tolReg.R"))

required_packages <- c("ggplot2", "patchwork", "pracma")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install the following R packages before running this script: ",
       paste(missing_packages, collapse = ", "))
}

data_dir <- file.path(script_dir, "data", "colon")
output_dir <- file.path(script_dir, "results", "colon")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expression <- as.matrix(read.table(file.path(data_dir, "I2000.txt")))
tissue_id <- read.table(file.path(data_dir, "tissues.txt"))[[1L]]
annotations <- read.delim(
  file.path(data_dir, "names2.txt"),
  header = FALSE, sep = "\t", quote = "", fill = TRUE,
  comment.char = "", strip.white = TRUE
)
gene_id <- annotations[[2L]]

if (!identical(dim(expression), c(2000L, 62L)) ||
    length(tissue_id) != ncol(expression) ||
    length(gene_id) != nrow(expression)) {
  stop("Unexpected dimensions in the colon data files.")
}

# The sample labelled 12 is the normal-tissue outlier excluded in the paper.
keep <- tissue_id != 12
expression <- expression[, keep, drop = FALSE]
tissue_id <- tissue_id[keep]
normal <- expression[, tissue_id > 0, drop = FALSE]
tumor <- expression[, tissue_id %in% -tissue_id[tissue_id > 0],
                    drop = FALSE]

if (!identical(dim(normal), c(2000L, 21L)) ||
    !identical(dim(tumor), c(2000L, 21L))) {
  stop("Could not recover the 21 matched normal/tumor tissue pairs.")
}

support <- seq(0, 1, by = 0.01)
transport_map <- function(normal_values, tumor_values, gene_index) {
  origin <- sort(
    (normal_values - min(normal_values)) / diff(range(normal_values))
  )
  destination <- sort(
    (tumor_values - min(normal_values)) / diff(range(normal_values))
  )

  # This is the single tied knot handled in the original analysis.
  if (gene_index == 1358L) origin[9L] <- origin[8L] + 1e-5
  pracma::spinterp(origin, destination, xp = support)
}

transport_maps <- lapply(
  seq_len(nrow(expression)),
  function(i) transport_map(normal[i, ], tumor[i, ], i)
)
names(transport_maps) <- gene_id

# The analysis setting below is the one reported for this application.
fit <- tolReg(
  transport_maps,
  optns = list(type = "function", supp = support, tau = 0)
)

n_outliers <- floor(length(transport_maps) * 0.05)
outlier_index <- order(fit$transRank)[seq_len(n_outliers)]
is_outlier <- seq_along(transport_maps) %in% outlier_index

score_table <- data.frame(
  gene = gene_id,
  c_score = fit$transRank,
  rank_from_periphery = rank(fit$transRank, ties.method = "first"),
  outlier = is_outlier
)
write.csv(score_table, file.path(output_dir, "colon_c_scores.csv"),
          row.names = FALSE)

profile_distance <- getPairDist(
  fit$distProf,
  optns = list(type = "measure")
)
profile_coordinates <- stats::cmdscale(profile_distance)
object_coordinates <- stats::cmdscale(fit$pairDist)

map_frame <- data.frame(
  x = rep(support, length(transport_maps)),
  y = unlist(transport_maps, use.names = FALSE),
  gene_index = rep(seq_along(transport_maps), each = length(support)),
  central = rep(!is_outlier, each = length(support))
)

p_map <- ggplot2::ggplot(
  map_frame,
  ggplot2::aes(x = x, y = y, group = gene_index, color = central)
) +
  ggplot2::geom_line(alpha = 0.3) +
  ggplot2::scale_color_discrete(
    name = NULL, breaks = FALSE, labels = "Outlier"
  ) +
  ggplot2::labs(
    x = "Gene expression (normal)",
    y = "Gene expression (tumor)",
    tag = "(a)"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 20),
    legend.position = "top",
    legend.key.width = grid::unit(1, "cm")
  )

embedding_frame <- function(coordinates) {
  data.frame(
    x = coordinates[, 1L],
    y = coordinates[, 2L],
    c_score = fit$transRank,
    outlier = is_outlier
  )
}

p_profile <- ggplot2::ggplot(
  embedding_frame(profile_coordinates),
  ggplot2::aes(x = x, y = y, color = c_score, shape = outlier)
) +
  ggplot2::geom_point(size = 3, alpha = ifelse(is_outlier, 1, 0.4)) +
  ggplot2::scale_color_gradient(
    name = "Centrality\nscore", low = "coral", high = "blue",
    guide = ggplot2::guide_colorbar(
      direction = "vertical", title.position = "top",
      barheight = grid::unit(3, "cm")
    )
  ) +
  ggplot2::scale_shape_discrete(
    name = NULL, breaks = TRUE, labels = "Outlier"
  ) +
  ggplot2::labs(x = "Coordinate 1", y = "Coordinate 2", tag = "(b)") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 20),
    legend.position = "top",
    legend.key.width = grid::unit(1, "cm")
  )

p_object <- ggplot2::ggplot(
  embedding_frame(object_coordinates),
  ggplot2::aes(x = x, y = y, color = c_score, shape = outlier)
) +
  ggplot2::geom_point(size = 3, alpha = ifelse(is_outlier, 1, 0.4)) +
  ggplot2::scale_color_gradient(
    name = "Centrality\nscore", low = "coral", high = "blue",
    guide = ggplot2::guide_colorbar(
      direction = "vertical", title.position = "top",
      barheight = grid::unit(3, "cm")
    )
  ) +
  ggplot2::scale_shape_discrete(
    name = NULL, breaks = TRUE, labels = "Outlier"
  ) +
  ggplot2::labs(x = "Coordinate 1", y = "Coordinate 2", tag = "(c)") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 20),
    legend.position = "top",
    legend.key.width = grid::unit(1, "cm")
  )

highlight_genes <- c("T57619", "T48804", "T61609")
highlight_index <- match(highlight_genes, gene_id)
if (anyNA(highlight_index)) stop("The highlighted genes were not found.")
normal_density <- lapply(
  highlight_index, function(i) stats::density(normal[i, ], from = 0)
)
tumor_density <- lapply(
  highlight_index, function(i) stats::density(tumor[i, ], from = 0)
)
density_frame <- do.call(
  rbind,
  lapply(seq_along(highlight_genes), function(j) {
    rbind(
      data.frame(x = normal_density[[j]]$x, density = normal_density[[j]]$y,
                 tissue = "Normal", gene = highlight_genes[j]),
      data.frame(x = tumor_density[[j]]$x, density = tumor_density[[j]]$y,
                 tissue = "Tumor", gene = highlight_genes[j])
    )
  })
)

p_density <- ggplot2::ggplot(
  density_frame,
  ggplot2::aes(x = x, y = density, color = tissue)
) +
  ggplot2::geom_line() +
  ggplot2::scale_color_manual(
    name = NULL, values = c(Normal = "#619CFF", Tumor = "#F8766D")
  ) +
  ggplot2::facet_wrap(ggplot2::vars(gene), ncol = 3) +
  ggplot2::labs(x = "Gene expression", y = "Density", tag = "(d)") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 20),
    legend.key.width = grid::unit(1, "cm"),
    legend.position = "bottom"
  )

figure <- (
  p_map +
    (p_profile + p_object + patchwork::guide_area() +
       patchwork::plot_layout(guides = "collect", widths = c(5, 5, 2)) &
       ggplot2::theme(legend.position = "right")) +
    patchwork::plot_layout(widths = c(5, 12))
) / p_density

grDevices::pdf(file.path(output_dir, "colon.pdf"), width = 15, height = 9)
print(figure)
grDevices::dev.off()
message("Colon analysis complete. Results: ", output_dir)
