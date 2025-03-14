# FILE: /single-cell-analysis/single-cell-analysis/scripts/utils.R
# This file contains utility functions for single-cell analysis.

# -----------------------------------------------------------------------------
# 1. Install and Load Libraries 
# -----------------------------------------------------------------------------
# Function to install packages (Bioconductor and CRAN)
install_and_load_packages <- function(bioc_pkgs = NULL, cran_pkgs = NULL) {
    if (!is.null(bioc_pkgs)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        for (pkg in bioc_pkgs) {
          if (!requireNamespace(pkg, quietly = TRUE)) {
            BiocManager::install(pkg)
          } else {
            print(paste("package", pkg, "already installed"))
          }
          library(pkg, character.only = TRUE)
        }
    }

    if (!is.null(cran_pkgs)) {
      for (pkg in cran_pkgs) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            install.packages(pkg)
          } else {
            print(paste("package", pkg, "already installed"))
          }
         library(pkg, character.only = TRUE)
      }
    }
}

# -----------------------------------------------------------------------------
# 2. Data Loading Functions
# -----------------------------------------------------------------------------
load_data_from_tabular <- function(file_path, sparse = FALSE) {
  if(sparse) {
     message(paste("Loading data from text file:", file_path, "(sparse)"))
      sparse.mat <- readSparseCounts(file_path)[,-1]
      return(sparse.mat)
  } else {
     message(paste("Loading data from text file:", file_path, "(dense)"))
      mat <- as.matrix(read.delim(file_path, check.names = FALSE))
      return(mat)
  }
}

# grepl("ENSMUSG*", rownames(sce_count))

load_data_from_excel <- function(file_path) {
    message(paste("Loading data from Excel file:", file_path))
    gunzip(file_path, destname = gsub(".gz$","",file_path), remove=FALSE, overwrite=TRUE)
    all.counts <- read_excel(gsub(".gz$","",file_path))
    gene.names <- all.counts$ID
    all.counts <- as.matrix(all.counts[,-1])
    rownames(all.counts) <- gene.names
    return(all.counts)
}
load_data_from_cellranger <- function(file_path) {
    message(paste("Loading data from Cell Ranger output:", file_path))
    untar(file_path, exdir = file.path(dirname(file_path),"tenx-output"))
    sce <- read10xCounts(file.path(dirname(file_path),"tenx-output/raw_gene_bc_matrices/GRCh38/"))
    return(sce)
}

create_single_cell_object <- function(counts, metadata = NULL, spike_mat=NULL, length_data=NULL){
  message("Creating SingleCellExperiment Object")
    sce <- SingleCellExperiment(list(counts=counts))
    if(!is.null(spike_mat)) {
      assay(sce,"spike_counts") <- spike_mat
      message("  Added spike in counts as assay")
    }

    if(!is.null(metadata)){
      colData(sce) <- metadata
      message("  Added metadata to colData")
    }
    if(!is.null(length_data)) {
        rowData(sce)$Length <- length_data
        message("  Added gene lengths to rowData")
    }
   return(sce)
}

add_metadata_from_file <- function(sce, meta_path, data_file_name){
  message(paste("Adding metadata from file", meta_path))
  coldata <- read.delim(meta_path, check.names = F)
  coldata <- coldata[coldata[, "Derived Array Data File"] == data_file_name, ]
  coldata <- DataFrame(genotype=coldata[, "Characteristics[genotype]"], 
                       phenotype=coldata[, "Characteristics[phenotype]"], 
                       spike_in=coldata[, "Factor Value[spike-in addition]"], 
                       row.names = coldata[, "Source Name"])
    colData(sce) <- coldata
    return(sce)
}


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


# Define required 10X files (both compressed and uncompressed versions)
check_10x_files <- function(data_dir) {
    required_files <- list(
        compressed = c(
            "barcodes.tsv.gz",
            "features.tsv.gz",
            "matrix.mtx.gz"
        ),
        uncompressed = c(
            "barcodes.tsv",
            "features.tsv",
            "matrix.mtx"
        )
    )
    
    # Check for compressed files
    compressed_exists <- all(sapply(required_files$compressed, 
        function(f) file.exists(file.path(data_dir, f))))
    
    # Check for uncompressed files
    uncompressed_exists <- all(sapply(required_files$uncompressed, 
        function(f) file.exists(file.path(data_dir, f))))
    
    if(!compressed_exists && !uncompressed_exists) {
        message("Missing required files. Neither compressed nor uncompressed files found.")
        message("Expected either:")
        message(paste(" -", required_files$compressed, collapse = "\n"))
        message("OR")
        message(paste(" -", required_files$uncompressed, collapse = "\n"))
        return(FALSE)
    }
    
    # Print which format was found
    if(compressed_exists) {
        message("Found compressed 10X files")
    } else {
        message("Found uncompressed 10X files")
    }
    return(TRUE)
}
