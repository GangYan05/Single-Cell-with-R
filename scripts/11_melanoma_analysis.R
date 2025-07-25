# 11_melanoma_analysis.R
# This script performs a complete single-cell analysis workflow on a public melanoma dataset.
# The workflow includes:
# 1. Data Loading
# 2. Quality Control
# 3. Normalization
# 4. Feature Selection
# 5. Dimensionality Reduction
# 6. Clustering
# 7. Differential Expression Analysis
# 8. Cell Type Annotation
# 9. Visualization

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
# Ensure necessary packages are loaded. 
# You can use the function from 00_utils.R or load them manually.
library(SingleCellExperiment)
library(scRNAseq)
library(scater)
library(scran)
library(ggplot2)
library(celldex)
library(SingleR)
library(pheatmap)
library(igraph)
library(Matrix)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install()


# Create directories for results if they don't exist
dir.create("results/melanoma", showWarnings = FALSE, recursive = TRUE)
dir.create("results/qc_metrics/melanoma", showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# 1. Data Loading
# -----------------------------------------------------------------------------
# Load the Tirosh et al. 2016 melanoma dataset from the scRNAseq package
sce_raw <- read.table("data/GSE72056/GSE72056_melanoma_single_cell_revised_v2.txt", 
header=TRUE, sep="\t")
# The first three rows contain metadata, so we will skip them.
count_df <- sce_raw[-(1:3), ]
info <- sce_raw[1:3, ]
# Find the duplicate gene names
dup_genes <- sce_raw$Cell[duplicated(sce_raw$Cell)]
sce_raw[(sce_raw$Cell %in% dup_genes),1:5]
# Aggregate counts of duplicate genes by summing their counts.
count_matrix <- as.matrix(count_df[ , !(colnames(count_df) %in% "Cell") ])
agg_counts <- rowsum(count_matrix, group = count_df$Cell)
rownames(agg_counts) <- rownames(agg_counts)
# Load the metadata separately
sce_meta <- read.table("data/GSE72056/melanoma_cluster_assignment_portal.txt", 
header=2, sep="\t", )[-1,-1]

sce <- SingleCellExperiment(assays = list(counts=agg_counts), colData = sce_meta)
colData(sce)

# The data is already a SingleCellExperiment object.
# Gene identifiers are symbols in the rownames.
# Metadata is in colData.

# Print basic information about the loaded data
message("Data loaded successfully:")
message(sprintf("Number of genes: %d", nrow(sce)))
message(sprintf("Number of cells: %d", ncol(sce)))
message("Available metadata columns:")
print(colnames(colData(sce)))


# -----------------------------------------------------------------------------
# 2. Quality Control (QC)
# -----------------------------------------------------------------------------

# Identify mitochondrial genes or ribosomal genes.
is_mito <- grepl("^MT-|mt-|Mito", rownames(sce))
message(sprintf("Found %d mitochondrial genes.", sum(is_mito)))
is_ribo <- grepl("^RPL|^RPS", rownames(sce))
message(sprintf("Found %d ribosomal genes.", sum(is_ribo)))

# Calculate comprehensive QC metrics and add QC metrics to the SingleCellExperiment object
sce <- addPerCellQC(sce, subsets = list(
    Ribo= is_ribo))

# Determine QC thresholds adaptively using MADs (Median Absolute Deviations)
# We filter on low library size, low number of features, and high mitochondrial content.
qc_filters <- perCellQCFilters(qc_metrics,
    sub.fields=c("sum", "detected", "subsets_Ribo_percent"),
    nmads=3
)

# Summarize how many cells are filtered by each criterion
colSums(as.matrix(qc_filters))

# Add the final discard decision to the colData
colData(sce)$discard <- qc_filters$discard

# Visualize QC metrics
# Library size distribution
p1 <- plotColData(sce, y = "sum", colour_by = "discard") + 
    scale_y_log10() + ggtitle("Library Size")

# Number of expressed genes
p2 <- plotColData(sce, y = "detected", colour_by = "discard") + 
    scale_y_log10() + ggtitle("Detected Features")

# Eibosomal gene proportion
p3 <- plotColData(sce, y = "subsets_Ribo_percent", colour_by = "discard") + 
    ggtitle("Mitochondrial %")

# Arrange plots
gridExtra::grid.arrange(p1, p2, p3, ncol=1)


# Filter cells based on the discard flag
sce_filtered <- sce[, !sce$discard]
message(sprintf("Cells before filtering: %d", ncol(sce)))
message(sprintf("Cells after filtering: %d", ncol(sce_filtered)))
message(sprintf("Removed %d cells.", sum(sce$discard)))


# Save QC results
saveRDS(qc_metrics, "results/qc_metrics/melanoma/melanoma_qc_metrics.rds")  
saveRDS(qc_filters, "results/qc_metrics/melanoma/melanoma_qc_filters.rds")
saveRDS(sce_filtered, "results/melanoma/melanoma_sce_filtered.rds")


# -----------------------------------------------------------------------------
# 3. Normalization
# -----------------------------------------------------------------------------
# Use deconvolution-based size factors for more accurate normalization.
# This method pools cells with similar expression profiles to estimate size factors.
set.seed(1234)
clusters <- quickCluster(sce_filtered)
sce_filtered <- computeSumFactors(sce_filtered, clusters=clusters)

# Check the size factors
summary(sizeFactors(sce_filtered))

# Apply log-transformation to the counts
sce_filtered <- logNormCounts(sce_filtered)


# -----------------------------------------------------------------------------
# 4. Feature Selection
# -----------------------------------------------------------------------------
# Identify highly variable genes (HVGs) to focus on biologically meaningful variation.
# We model the variance of each gene and select those with high biological components.
dec <- modelGeneVar(sce_filtered)

# Get top 2000 most variable genes
top_hvgs <- getTopHVGs(dec, n=2000)

# Visualize the mean-variance trend
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
points(dec$mean[top_hvgs], dec$total[top_hvgs], col="red")
curve(metadata(dec)$trend(x), col="dodgerblue", add=TRUE)


# -----------------------------------------------------------------------------
# 5. Dimensionality Reduction
# -----------------------------------------------------------------------------
# Perform PCA on the HVGs to capture the main axes of variation.
set.seed(1234)
sce_filtered <- runPCA(sce_filtered, subset_row=top_hvgs)

# Run UMAP and t-SNE for visualization. We use the first 20 PCs.
sce_filtered <- runUMAP(sce_filtered, dimred="PCA", n_dimred=20)
sce_filtered <- runTSNE(sce_filtered, dimred="PCA", n_dimred=20)


# -----------------------------------------------------------------------------
# 6. Clustering
# -----------------------------------------------------------------------------
# Perform graph-based clustering on the PCA results.
# This builds a shared nearest-neighbor graph and identifies communities (clusters).
g <- buildSNNGraph(sce_filtered, k=10, use.dimred="PCA")
clusters_louvain <- igraph::cluster_louvain(g)
colLabels(sce_filtered) <- factor(clusters_louvain$membership)

# Visualize clusters on the UMAP plot
plotReducedDim(sce_filtered, "UMAP", colour_by="label", text_by="label") +
    ggtitle("UMAP colored by Louvain Clusters")


# -----------------------------------------------------------------------------
# 7. Differential Expression Analysis
# -----------------------------------------------------------------------------
# Find marker genes for each cluster to help identify cell types.
# We use a Wilcoxon rank-sum test to find genes upregulated in each cluster vs. others.
markers <- findMarkers(sce_filtered, groups=colLabels(sce_filtered), test.type="wilcox", direction="up")

# Example: View top markers for cluster 1
cluster1_markers <- markers[["1"]]
head(cluster1_markers)


# -----------------------------------------------------------------------------
# 8. Cell Type Annotation
# -----------------------------------------------------------------------------
# Use SingleR to automatically annotate clusters based on reference datasets.
# We use the Human Primary Cell Atlas as a reference.
ref <- celldex::HumanPrimaryCellAtlasData()

# Perform annotation. This may take a few minutes.
pred <- SingleR(test=sce_filtered, ref=ref, labels=ref$label.main, de.method="wilcox")

# Add cell type labels to the object
colData(sce_filtered)$cell_type_main <- pred$labels

# Visualize UMAP with predicted cell types
plotReducedDim(sce_filtered, "UMAP", colour_by="cell_type_main") +
    ggtitle("UMAP colored by Predicted Cell Type (SingleR)")

# Tabulate predicted cell types vs. clusters
table(Predicted=colData(sce_filtered)$cell_type_main, Cluster=colLabels(sce_filtered))


# -----------------------------------------------------------------------------
# 9. Visualization & Saving
# -----------------------------------------------------------------------------
# Generate a heatmap of the top 10 marker genes for each cluster.
top_markers <- lapply(markers, function(x) {
    # Filter for significant markers and take the top 10
    x_df <- as.data.frame(x)
    sig_markers <- rownames(x_df[x_df$FDR < 0.05, ])
    head(sig_markers, 10)
})

# Plot the heatmap
plotHeatmap(sce_filtered, 
            features=unique(unlist(top_markers)),
            columns=order(colLabels(sce_filtered)),
            colour_columns_by=c("label", "cell_type_main"),
            cluster_rows=TRUE,
            show_colnames=FALSE)

# Save final results
saveRDS(sce_filtered, file="results/melanoma/melanoma_analyzed.rds")
saveRDS(markers, file="results/melanoma/melanoma_cluster_markers.rds")

message("Melanoma analysis complete. Results saved in 'results/melanoma/'.")