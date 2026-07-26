#!/usr/bin/env Rscript

# MNIST handwriting application
# Reads the original IDX files, samples 200 images per digit, computes C-scores
# within digit classes, and recreates the typical/atypical handwriting figure.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  script_file <- sub("^--file=", "", script_arg[1L])
  script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
  dirname(normalizePath(script_file))
} else {
  normalizePath(getwd())
}
source(file.path(script_dir, "tolReg.R"))

data_dir <- file.path(script_dir, "data", "mnist")
output_dir <- file.path(script_dir, "results", "mnist")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_idx_images <- function(path) {
  connection <- file(path, "rb")
  on.exit(close(connection))
  magic <- readBin(connection, "integer", n = 1L, size = 4L, endian = "big")
  n <- readBin(connection, "integer", n = 1L, size = 4L, endian = "big")
  n_row <- readBin(connection, "integer", n = 1L, size = 4L, endian = "big")
  n_col <- readBin(connection, "integer", n = 1L, size = 4L, endian = "big")
  if (magic != 2051L || n_row != 28L || n_col != 28L) {
    stop("Unexpected MNIST image-file header.")
  }
  pixels <- readBin(
    connection, "integer", n = n * n_row * n_col,
    size = 1L, signed = FALSE
  )
  matrix(pixels, nrow = n, byrow = TRUE)
}

read_idx_labels <- function(path) {
  connection <- file(path, "rb")
  on.exit(close(connection))
  magic <- readBin(connection, "integer", n = 1L, size = 4L, endian = "big")
  n <- readBin(connection, "integer", n = 1L, size = 4L, endian = "big")
  if (magic != 2049L) stop("Unexpected MNIST label-file header.")
  readBin(connection, "integer", n = n, size = 1L, signed = FALSE)
}

images <- read_idx_images(file.path(data_dir, "train-images.idx3-ubyte"))
labels <- read_idx_labels(file.path(data_dir, "train-labels.idx1-ubyte"))
if (nrow(images) != length(labels)) {
  stop("MNIST images and labels have incompatible lengths.")
}

# Match the orientation used in the manuscript figure.
images <- t(apply(images, 1L, function(x) {
  as.vector(matrix(x, nrow = 28L)[, 28:1])
}))

sample_size <- 200L
fits <- vector("list", 10L)
sampled_images <- vector("list", 10L)
sampled_indices <- vector("list", 10L)
score_tables <- vector("list", 10L)

for (digit in 0:9) {
  class_indices <- which(labels == digit)
  set.seed(1)
  selected <- sample(class_indices, sample_size)
  digit_images <- images[selected, , drop = FALSE] / 50
  digit_objects <- split(digit_images, row(digit_images))

  fit <- tolReg(
    digit_objects,
    optns = list(type = "euclidean")
  )

  fits[[digit + 1L]] <- fit
  sampled_images[[digit + 1L]] <- digit_images
  sampled_indices[[digit + 1L]] <- selected
  score_tables[[digit + 1L]] <- data.frame(
    digit = digit,
    original_training_index = selected,
    c_score = fit$transRank,
    rank_from_periphery = rank(fit$transRank, ties.method = "first")
  )
}

write.csv(
  do.call(rbind, score_tables),
  file.path(output_dir, "mnist_c_scores.csv"),
  row.names = FALSE
)

# Arrange digits 0--4 in the first two rows and 5--9 in the last two rows.
# Each digit contributes four representative and four atypical samples.
grDevices::png(
  file.path(output_dir, "mnist.png"),
  width = 8000, height = 1600, res = 200, bg = "black"
)
graphics::par(
  mfrow = c(4, 20), oma = rep(0, 4), mar = rep(0.08, 4),
  bg = "black", fg = "white"
)

draw_digit <- function(values) {
  graphics::image(
    matrix(values, ncol = 28L),
    col = grDevices::gray((0:255) / 255),
    axes = FALSE, useRaster = TRUE
  )
  graphics::box(col = "white")
}

for (digit in 0:4) {
  fit <- fits[[digit + 1L]]
  digit_images <- sampled_images[[digit + 1L]]
  ordering <- order(fit$transRank)
  representative <- rev(ordering)[1:4]
  atypical <- ordering[1:4]
  for (index in representative) draw_digit(digit_images[index, ])
}
for (digit in 0:4) {
  fit <- fits[[digit + 1L]]
  digit_images <- sampled_images[[digit + 1L]]
  for (index in order(fit$transRank)[1:4]) draw_digit(digit_images[index, ])
}
for (digit in 5:9) {
  fit <- fits[[digit + 1L]]
  digit_images <- sampled_images[[digit + 1L]]
  ordering <- order(fit$transRank)
  representative <- rev(ordering)[1:4]
  for (index in representative) draw_digit(digit_images[index, ])
}
for (digit in 5:9) {
  fit <- fits[[digit + 1L]]
  digit_images <- sampled_images[[digit + 1L]]
  for (index in order(fit$transRank)[1:4]) draw_digit(digit_images[index, ])
}
grDevices::dev.off()

message("MNIST analysis complete. Results: ", output_dir)
