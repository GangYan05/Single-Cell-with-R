# FILE: /single-cell-analysis/single-cell-analysis/scripts/utils.R
# This file contains utility functions for single-cell analysis.

load_data <- function(file_path) {
  # Function to load single-cell data from a specified file path
  data <- readRDS(file_path)
  return(data)
}

plot_qc_metrics <- function(qc_data) {
  # Function to plot quality control metrics
  library(ggplot2)
  ggplot(qc_data, aes(x = total_counts, y = percent_mito)) +
    geom_point() +
    theme_minimal() +
    labs(title = "Quality Control Metrics",
         x = "Total Counts",
         y = "Percent Mitochondrial")
}

save_results <- function(data, file_path) {
  # Function to save results to a specified file path
  saveRDS(data, file_path)
}

normalize_counts <- function(counts) {
  # Function to normalize single-cell counts
  library(scater)
  normalized_counts <- logNormCounts(counts)
  return(normalized_counts)
}