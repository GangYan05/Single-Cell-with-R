# Single-Cell Analysis Project

This project is designed for the analysis of single-cell RNA sequencing data. It includes various scripts for quality control, normalization, feature selection, dimensionality reduction, and data subsetting. The results of the analyses are organized into separate directories for easy access and interpretation.

## Project Structure

- **data/raw_data**: Contains the raw single-cell data files for processing.
  
- **scripts/**: 
  - **01_quality_control.R**: Handles quality control of the single-cell data, including metrics calculation and filtering based on quality thresholds.
  - **02_normalization.R**: Performs normalization of the single-cell count data to account for differences in sequencing depth and other technical variations.
  - **03_feature_selection.R**: Selects relevant features (genes) for downstream analysis, focusing on those that are most informative for the biological questions being addressed.
  - **04_dimensionality_reduction.R**: Implements dimensionality reduction techniques such as PCA, t-SNE, and UMAP to visualize the high-dimensional single-cell data.
  - **05_subsetting_combining.R**: Subsets the single-cell data based on specific criteria and combines datasets if necessary for comparative analysis.
  - **utils.R**: Contains utility functions used across the other scripts, such as data loading, plotting functions, and other helper functions.

- **results/**:
  - **qc_metrics/**: Stores the output of the quality control metrics generated from the analysis.
  - **normalized_data/**: Contains the normalized single-cell data files after processing.
  - **feature_selection/**: Stores the results of the feature selection process, including selected genes and their metrics.
  - **dimensionality_reduction/**: Contains the results of the dimensionality reduction analyses, including plots and transformed data.
  - **subsets/**: Stores the results of any subsetting operations performed on the data.

## Instructions

1. Place your raw single-cell data files in the `data/raw_data` directory.
2. Run the scripts in the `scripts` directory in the following order:
   - `01_quality_control.R`
   - `02_normalization.R`
   - `03_feature_selection.R`
   - `04_dimensionality_reduction.R`
   - `05_clustering.R`
   - `06_marker_gene`
   - `07_cell_annotation`
3. Check the `results` directory for the output of each analysis step.

## Requirements

Ensure you have the necessary R packages installed to run the scripts. You may need packages such as `scater`, `uwot`, and `scRNAseq`.

## License

This project is licensed under the MIT License.# Single-Cell-with-R
