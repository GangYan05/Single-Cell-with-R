# Single-Cell RNA-Seq Analysis Workflows in R

This repository contains comprehensive R scripts for performing single-cell RNA sequencing (scRNA-seq) analysis. Each script represents a complete, end-to-end workflow for a specific public dataset, demonstrating key steps from data loading to biological interpretation.

## Project Structure

- **`data/`**: This directory is intended to store the input data for the analysis scripts. The scripts expect specific subdirectories corresponding to the dataset being analyzed.
  - e.g., `data/GSE72056/` for the Melanoma dataset.
  - e.g., `data/raw_data/GSE176078/` for the TNBC dataset.

- **`scripts/`**: Contains the R scripts, each performing a full analysis workflow.

  - **`01_data_input`**: Reading count matrices from various formats (e.g., 10x Genomics, flat files).
  - **`02_quality_control`**: Calculating cell-level metrics (library size, feature count, mitochondrial/ribosomal content) and filtering out low-quality cells.
  - **`03_normalization`**: Correcting for library size differences using deconvolution-based size factors.
  - **`04_feature_selection`**:Identifying highly variable genes (HVGs) to focus on biological signal.
  - **`05_dimensionality_reduction`**: Performing PCA, UMAP, and t-SNE for data visualization and downstream analysis.
  - **`06_clustering`**: Grouping cells based on their expression profiles using graph-based methods.
  - **`07_marker_gene`**: Finding marker genes that define each cluster.
  - **`08_cell_annotation`**: Automatically assigning cell identities using reference datasets with `SingleR`.
  - **`09_batch_removal`**: Correcting for technical differences between samples processed in different batches using methods like MNN.
  - **`10_TNBC_cell_lines.R`**: An analysis pipeline for Triple-Negative Breast Cancer (TNBC) cell line data (GSE176078).
  - **`11_melanoma_analysis.R`**: An analysis pipeline for a melanoma patient dataset (GSE72056), including steps like batch effect investigation and sub-clustering.
  - **`00_utils.R` (optional)**: A place for utility functions that can be shared across different analysis scripts.
  
- **`results/`**: This directory is automatically created by the scripts to store all outputs, such as plots, RDS objects, and tables. Results are organized into subdirectories for each analysis.
  - **`results/melanoma/`**: Outputs from the melanoma analysis.
  - **`results/TNBC/`**: Outputs from the TNBC analysis.
  - **`results/qc_metrics/`**: A shared location for quality control metric files.

## Requirements

This project relies heavily on the R/Bioconductor ecosystem. Key packages include:

-   `SingleCellExperiment`
-   `scater`
-   `scran`
-   `DropletUtils` (for 10x data)
-   `celldex` & `SingleR` (for annotation)
-   `ggplot2`
-   `pheatmap`
-   `igraph`

## License

This project is licensed under the MIT License.
