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
library(gridExtra)
library(igraph)
library(Matrix)

# Create directories for results if they don't exist
dir.create("results/melanoma", showWarnings = FALSE, recursive = TRUE)
dir.create("results/qc_metrics/melanoma", showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# 1. Data Loading
# -----------------------------------------------------------------------------
# Load the Tirosh et al. 2016 melanoma dataset from local files.
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

# Extract patient ID from complex cell names (e.g., "72" from "Cy72_...", "80" from "cy80-...", "89" from "CY89NEG...").
# The patient ID is the numeric part immediately following a case-insensitive "Cy" prefix.
patient_ids <- sub("^cy(\\d+).*", "\\1", rownames(colData(sce)), ignore.case = TRUE)

# For cell names that do not match the "Cy" pattern (e.g., "monika_...", "SS2_..."),
# the substitution fails and returns the original string. We will label these as "unknown".
# We identify these by checking for any non-digit characters in the result.
patient_ids[grepl("[^0-9]", patient_ids)] <- "unknown"
sce$patient_id <- patient_ids
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

# Identify ribosomal genes.
is_ribo <- grepl("^RPL|^RPS", rownames(sce))
message(sprintf("Found %d ribosomal genes.", sum(is_ribo)))

# Calculate comprehensive QC metrics and add QC metrics to the SingleCellExperiment object
sce <- addPerCellQC(sce, subsets = list(
    Ribo= is_ribo))

# The QC metrics are in colData(sce)
qc_metrics <- colData(sce)

# Determine QC thresholds adaptively using MADs (Median Absolute Deviations)
# We filter on low library size, low number of features, and high ribosomal content.
qc_filters <- perCellQCFilters(qc_metrics,
    sub.fields=c("sum", "detected", "subsets_Ribo_percent"),
    nmads = 3
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

# Ribosomal gene proportion
p3 <- plotColData(sce, y = "subsets_Ribo_percent", colour_by = "discard") + 
    ggtitle("Ribosomal %")

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

qc_metrics <- readRDS("results/qc_metrics/melanoma/melanoma_qc_metrics.rds")
qc_filters <- readRDS("results/qc_metrics/melanoma/melanoma_qc_filters.rds")
sce_filtered <- readRDS("results/melanoma/melanoma_sce_filtered.rds")
# -----------------------------------------------------------------------------
# 3. Normalization
# -----------------------------------------------------------------------------
# Use deconvolution-based size factors for more accurate normalization.
set.seed(1234)
clusters <- quickCluster(sce_filtered)
summary(clusters)
sce_filtered <- computeSumFactors(sce_filtered, clusters=clusters)

# Check the size factors
summary(sizeFactors(sce_filtered))

# Apply log-transformation to the counts
sce_filtered <- logNormCounts(sce_filtered)

# -----------------------------------------------------------------------------
# 4. Feature Selection
# -----------------------------------------------------------------------------
# Identify highly variable genes (HVGs) to focus on biologically meaningful variation.
dec <- modelGeneVar(sce_filtered)

# Get top 2000 most variable genes
top_hvgs <- getTopHVGs(dec, n=2000)

# Visualize the mean-variance trend
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
points(dec[top_hvgs, "mean"], dec[top_hvgs, "total"], col="red")
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
plotReducedDim(sce_filtered, "UMAP", colour_by="CLUSTER", text_by="CLUSTER") +
    ggtitle("UMAP colored by Louvain Clusters")
# export the figure
ggsave("results/melanoma/melanoma_umap_clusters.png", width=10, height=10)

# Visualize clusters on the t-SNE plot
plotReducedDim(sce_filtered, "TSNE", colour_by="CLUSTER", text_by="CLUSTER") +
    ggtitle("tSNE colored by Louvain Clusters")
# export the figure
ggsave("results/melanoma/melanoma_tSNE_clusters.png", width=10, height=10)

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


# -----------------------------------------------------------------------------
# 10. Comparative Analysis Between Groups
# (Sub-clustering within the Malignant Population)
# -----------------------------------------------------------------------------
# After identifying broad cell types, we can perform a more focused analysis
# on a specific population to find more subtle heterogeneity. Here, we will
# take all cells identified as 'malignant' and re-cluster them.

# We will use the 'CLUSTER' column from the original metadata for grouping.
if ("CLUSTER" %in% colnames(colData(sce_filtered))) {
    
    message("Performing sub-clustering on the 'malignant' cell population...")
    
    # 1. Subset the data to malignant cells
    sce_malignant <- sce_filtered[, sce_filtered$CLUSTER == "malignant"]
    
    message(sprintf("Found %d malignant cells for sub-clustering.", ncol(sce_malignant)))
    
    # 2. Re-run feature selection to find HVGs specific to this subset
    # This identifies genes that drive variation *within* the malignant cells.
    dec_malig <- modelGeneVar(sce_malignant)
    top_hvgs_malig <- getTopHVGs(dec_malig, n=1500) # Using 1500 HVGs for the subset
    
    # 3. Re-run dimensionality reduction on the malignant-specific HVGs
    set.seed(5678)
    sce_malignant <- runPCA(sce_malignant, subset_row=top_hvgs_malig, name="PCA_malignant")
    sce_malignant <- runUMAP(sce_malignant, dimred="PCA_malignant", name="UMAP_malignant")
    
    # 4. Re-run clustering to find sub-clusters
    g_malig <- buildSNNGraph(sce_malignant, use.dimred="PCA_malignant", k=10)
    clusters_malig <- igraph::cluster_louvain(g_malig)
    sce_malignant$sub_cluster <- factor(clusters_malig$membership)
    
    # 5. Visualize the new sub-clusters
    p_subcluster <- plotReducedDim(sce_malignant, "UMAP_malignant", colour_by="sub_cluster", text_by="sub_cluster") +
        ggtitle("Sub-clusters within Malignant Cells")
    
    print(p_subcluster)
    ggsave("results/melanoma/melanoma_umap_malignant_subclusters.png", plot = p_subcluster, width=8, height=8)
    
    message("Sub-clustering of malignant cells complete.")
    
    # 6. Save the new SCE object with sub-clustering results
    saveRDS(sce_malignant, file="results/melanoma/melanoma_malignant_subclustered.rds")
    
} else {
    message("Column 'CLUSTER' not found in colData. Skipping sub-clustering analysis.")
}
