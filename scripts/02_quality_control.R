# Load sce object
sce <- readRDS("results/sce.rds")
is_mito <- readRDS("results/is_mito.rds")
coldata <- readRDS("results/coldata.rds")

# Calculate quality control metrics
qc_df <- perCellQCMetrics(sce, subsets=list(Mito = is_mito))
qc_df

# Identify low-quality cells with fixed thresholds (optional)
qc_size <- qc_df$sum < 1e5
qc_nexpression <- qc_df$detected < 5e3
qc_spike <- qc_df$detected > 1e4
qc_mito <- qc_df$subsets_Mito_percent > 10
todrop <- qc_size | qc_nexpression | qc_spike | qc_mito
DataFrame(LibSize=sum(qc_size), NExprs=sum(qc_nexpression), Spike=sum(qc_spike), Mito=sum(qc_mito), total=sum(todrop))

# Adaptively filter out low-quality cells
sce_qc <- addPerCellQC(sce, subsets=list(Mito = is_mito))
tofilter <- quickPerCellQC(colData(sce_qc), sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent"))
sce_qc$discard <- tofilter$discard

gridExtra::grid.arrange(
plotColData(sce_qc, x="phenotype", y="sum", colour_by="discard") + 
    scale_y_log10() + ggtitle("total counts"),
plotColData(sce_qc, x="phenotype", y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("total counts"),
plotColData(sce_qc, x="phenotype", y="altexps_ERCC_percent", colour_by="discard") +
    ggtitle("ERCC percent")
)

# Subset the SingleCellExperiment object to retain only high-quality cells
sce <- sce[, !sce_qc$discard]

# Check the low-quality cells by calculating the average count for each feature.
lost_cells <- calculateAverage(counts(sce)[, sce_qc$discard])
keep_cells <- calculateAverage(counts(sce)[, !sce_qc$discard])
logged <- cpm(cbind(lost_cells, keep_cells), log=TRUE)
logFC <- logged[, 1] - logged[, 2]
abundance <- rowMeans(logged)

plot(abundance, logFC, pch=20, cex=0.5, xlab="average abundance", ylab="log fold change")
points(abundance[is_mito], logFC[is_mito], pch=40, cex=0.5, col="darkblue")


calculateAver# Save quality control metrics to results directory
saveRDS(qc_df, file = "results/qc_metrics/qc_stats.rds")
saveRDS(sce, file = "results/SingleCellExperiment.rds")
