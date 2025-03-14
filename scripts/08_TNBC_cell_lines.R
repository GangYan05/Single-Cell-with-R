data_dir <- "data/raw_data/GSE176078"
files_status <- check_10x_files(data_dir)
# Load the single-cell data and metadata
if(files_status) {
    # Read the data using DropletUtils
    sce_raw <- DropletUtils::read10xCounts(
        samples = data_dir,
        col.names = TRUE,
        type = "sparse"
    )
    
    # Print basic information about the loaded data
    message("Data loaded successfully:")
    message(sprintf("Number of genes: %d", nrow(sce_raw)))
    message(sprintf("Number of cells: %d", ncol(sce_raw)))
}

# Define mitochondrial genes (for human data)
is_mito <- grepl("^MT-", rowData(sce_raw)$ID)
rowData(sce_raw)$ID[is_mito]

# Calculate and visualize mitochondrial content
mito_stats <- colSums(counts(sce_raw)[is_mito,]) / colSums(counts(sce_raw))
hist(mito_stats * 100, 
     breaks = 50, 
     main = "Distribution of Mitochondrial Content",
     xlab = "Mitochondrial Content (%)")

# Calculate ribosomal content
is_ribo <- grepl("^RP[LS]", rowData(sce_raw)$ID)
ribo_stats <- colSums(counts(sce_raw)[is_ribo,]) / colSums(counts(sce_raw))

# Visualize distribution
hist(ribo_stats * 100,
     breaks = 50,
     main = "Distribution of Ribosomal Content",
     xlab = "Ribosomal Content (%)")

# Print summary statistics
summary(ribo_stats * 100)

# Create data frame for plotting
plot_data <- data.frame(
    mito = mito_stats * 100,
    ribo = ribo_stats * 100
)

# Calculate correlation
cor_val <- cor(plot_data$mito, plot_data$ribo, method = "pearson")

# Create enhanced scatter plot with correlation line
library(ggplot2)
ggplot(plot_data, aes(x = mito, y = ribo)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "lm", color = "red") +
    labs(
        x = "Mitochondrial Content (%)",
        y = "Ribosomal Content (%)",
        title = "Mitochondrial vs Ribosomal Content",
        subtitle = sprintf("Correlation: %.3f", cor_val)
    ) +
    theme_minimal()

# Calculate comprehensive QC metrics
qc_metrics <- perCellQCMetrics(sce_raw,
    subsets=list(
        Mito=is_mito,
        Ribo=is_ribo
        )
)

# Add QC metrics to the SingleCellExperiment object
colData(sce_raw) <- cbind(colData(sce_raw), qc_metrics)

# Calculate QC thresholds using adaptive thresholds
qc_filters <- perCellQCFilters(qc_metrics,
    sub.fields=c("subsets_Mito_percent", "subsets_Ribo_percent"),
    nmads=3
)

# Library size distribution
hist(log10(qc_metrics$sum), 
     breaks=20, 
     main="Library size distribution",
     xlab="Log10 library size")

# Number of expressed genes
hist(log10(qc_metrics$detected), 
     breaks=20, 
     main="Number of expressed genes",
     xlab="Log10 number of genes")

# Mitochondrial proportion
hist(qc_metrics$subsets_Mito_percent,
     breaks=20,
     main="Mitochondrial proportion",
     xlab="Mitochondrial percent")

# Scatter plots
plot(qc_metrics$sum, 
     qc_metrics$subsets_Mito_percent,
     xlab="Library size",
     ylab="Mitochondrial percent",
     pch=16)

plot(qc_metrics$sum,
     qc_metrics$detected,
     xlab="Library size",
     ylab="Number of expressed genes",
     pch=16)


# Filter cells and summarize results
cells_removed <- data.frame(
    "Total cells" = ncol(sce_raw),
    "Low library size" = sum(qc_filters$low_lib_size),
    "High mito percent" = sum(qc_filters$high_sub_Mito_percent),
    "Low expressed genes" = sum(qc_filters$low_n_features),
    "Total removed" = sum(qc_filters$discard)
)

# Filter cells
sce_filtered <- sce_raw[, !qc_filters$discard]

# Save QC results
saveRDS(qc_metrics, "results/qc_metrics/qc_metrics.rds")
saveRDS(qc_filters, "results/qc_metrics/qc_filters.rds")
saveRDS(sce_filtered, "results/sce_filtered.rds")



# Need to modify the code below
# Load the filtered SingleCellExperiment object
sce <- readRDS("results/sce_filtered.rds")
# 1. Normalize the data
# Use deconvolution-based size factors
set.seed(1234)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
sce <- logNormCounts(sce)

# 2. Feature selection - identify highly variable genes
dec <- modelGeneVar(sce)
# Get top 2000 most variable genes
top_hvgs <- getTopHVGs(dec, n=2000)

# 3. Dimensionality reduction
# PCA
sce <- runPCA(sce, subset_row=top_hvgs)

# UMAP
sce <- runUMAP(sce, dimred="PCA", n_dimred=20)

# t-SNE
sce <- runTSNE(sce, dimred="PCA", n_dimred=20)

# 4. Clustering
# Graph-based clustering
g <- buildSNNGraph(sce, k=10, use.dimred="PCA")
clusters <- igraph::cluster_louvain(g)
colLabels(sce) <- factor(clusters$membership)

# 5. Visualization
# Plot UMAP with clusters
plotReducedDim(sce, "UMAP", colour_by="label")

# 6. Differential expression analysis
# Find markers for each cluster
markers <- findMarkers(sce, groups=colLabels(sce), 
                      test.type="wilcox", 
                      direction="up")

# 7. Cell type annotation
# Load reference data
ref <- celldex::HumanPrimaryCellAtlasData()

# Perform annotation
pred <- SingleR(test=sce, 
                ref=ref,
                labels=ref$label.main,
                de.method="wilcox")

# Add cell type labels to the object
colData(sce)$cell_type <- pred$labels

# 8. Quality visualization
# Plot number of cells per cluster
table(colLabels(sce))

# Plot UMAP with cell types
plotReducedDim(sce, "UMAP", colour_by="cell_type")

# 9. Save results
saveRDS(sce, file="results/tnbc_analyzed.rds")
saveRDS(markers, file="results/cluster_markers.rds")

# 10. Generate cluster marker heatmap
top_markers <- lapply(markers, function(x) {
    x <- as.data.frame(x)
    head(rownames(x)[x$FDR < 0.05], 10)
})

plotHeatmap(sce, features=unique(unlist(top_markers)),
           columns=order(colLabels(sce)),
           colour_columns_by="label",
           cluster_rows=TRUE,
           show_colnames=FALSE)










# Add metadata to the SingleCellExperiment object
sce <- add_metadata_from_file(sce, meta_fp, "counts_Calero_20160113.tsv")
# sce <- addPerCellQC(sce)
# sce <- addPerFeatureQC(sce)

# Add gene annotation
ah <- AnnotationHub()
ens_mm <- ah[["AH75036"]]
gene_info <- AnnotationDbi::select(ens_mm, keys = rownames(sce), keytype = "GENEID", columns = c("SYMBOL", "SEQNAME", "GENEBIOTYPE"))
rowData(sce)$SYMBOL <- gene_info$SYMBOL[match(rownames(sce), gene_info$GENEID)]
rowData(sce)$SEQNAME <- gene_info$SEQNAME[match(rownames(sce), gene_info$GENEID)]
rowData(sce)$GENEBIOTYPE <- gene_info$GENEBIOTYPE[match(rownames(sce), gene_info$GENEID)]

# Identify Mitochondrial genes using symbols
is_mito <- grepl("^mt-", rowData(sce)$SYMBOL)

# Calculate quality control metrics
qc_df <- perCellQCMetrics(sce, subsets = list(Mito = is_mito))

# Identify low-quality cells with fixed thresholds (optional)
qc_size <- qc_df$sum < 1e5
qc_nexpression <- qc_df$detected < 5e3
qc_spike <- qc_df$detected > 1e4
qc_mito <- qc_df$subsets_Mito_percent > 10
todrop <- qc_size | qc_nexpression | qc_spike | qc_mito
DataFrame(LibSize=sum(qc_size), NExprs=sum(qc_nexpression), Spike=sum(qc_spike), Mito=sum(qc_mito), total=sum(todrop))

# Adaptively filter out low-quality cells
filtered <- perCellQCFilters(qc_df ,sub.fields = c("subsets_Mito_percent"))

# Subset the SingleCellExperiment object to retain only high-quality cells
sce <- sce[, !filtered$discard]

# Save quality control metrics to results directory
saveRDS(qc_df, file = "results/qc_metrics/qc_stats.rds")
saveRDS(sce, file = "results/SingleCellExperiment.rds")
