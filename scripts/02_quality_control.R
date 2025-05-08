# Calculate quality control metrics
qc_df <- perCellQCMetrics(sce, subsets = list(Mito = is_mito))
qc_df

# Identify low-quality cells with fixed thresholds (optional)
qc_size <- qc_df$sum < 1e5
qc_nexpression <- qc_df$detected < 5e3
qc_spike <- qc_df$detected > 1e4
qc_mito <- qc_df$subsets_Mito_percent > 10
todrop <- qc_size | qc_nexpression | qc_spike | qc_mito
DataFrame(LibSize=sum(qc_size), NExprs=sum(qc_nexpression), Spike=sum(qc_spike), Mito=sum(qc_mito), total=sum(todrop))

# Adaptively filter out low-quality cells
tofilter <- quickPerCellQC(qc_df, percent_subsets=c("altexps_ERCC_percent", "subsets_Mito_percent"))
colSums(as.matrix(tofilter))


sce_qc <- addPerCellQC(sce)
sce_qc
colData(sce_qc)
plotColData(sce_qc, x="phenotype", y="sum")

colData(sce)
plotColData(sce, x="total_counts", y="percent_mito", colour_by="discard") +
    geom_point() +
    theme_minimal() +
    labs(title = "Quality Control Metrics",
         x = "Total Counts",
         y = "Percent Mitochondrial")
filtered <- perCellQCFilters(qc_df, sub.fields = c("subsets_Mito_percent"))

# Subset the SingleCellExperiment object to retain only high-quality cells
sce <- sce[, !filtered$discard]

# Save quality control metrics to results directory
saveRDS(qc_df, file = "results/qc_metrics/qc_stats.rds")
saveRDS(sce, file = "results/SingleCellExperiment.rds")
