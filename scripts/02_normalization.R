library(scater)

# Load the single-cell experiment object
sce <- readRDS("path/to/your/single_cell_experiment.rds")

# Normalization of the count data
sce <- logNormCounts(sce)

# Save the normalized data
saveRDS(sce, file = "results/normalized_data/normalized_sce.rds")