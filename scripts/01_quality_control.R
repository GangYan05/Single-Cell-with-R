library(scater)
library(scRNAseq)

# Load the single-cell data
sce <- MacoskoRetinaData()

# Identify mitochondrial genes for quality control
is.mito <- grepl("^MT-", rownames(sce))

# Calculate quality control metrics
qcstats <- perCellQCMetrics(sce, subsets = list(Mito = is.mito))

# Filter cells based on quality metrics
filtered <- quickPerCellQC(qcstats, percent_subsets = "subsets_Mito_percent")

# Subset the SingleCellExperiment object to retain only high-quality cells
sce <- sce[, !filtered$discard]

# Save quality control metrics to results directory
saveRDS(qcstats, file = "results/qc_metrics/qc_stats.rds")