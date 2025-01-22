library(scater)
library(uwot)

# Load the single-cell experiment object
sce <- readRDS("path/to/your/sce_object.rds")

# Perform PCA
sce <- scater::runPCA(sce)
pca_results <- reducedDim(sce, "PCA")

# Perform t-SNE
sce <- scater::runTSNE(sce, perplexity = 30)
tsne_results <- reducedDim(sce, "TSNE")

# Perform UMAP
umap_results <- uwot::umap(t(logcounts(sce)), n_neighbors = 15)
reducedDim(sce, "UMAP_UWOT") <- umap_results

# Save results
saveRDS(pca_results, file = "results/dimensionality_reduction/pca_results.rds")
saveRDS(tsne_results, file = "results/dimensionality_reduction/tsne_results.rds")
saveRDS(umap_results, file = "results/dimensionality_reduction/umap_results.rds")

# Optionally, create plots for visualization
pdf("results/dimensionality_reduction/pca_plot.pdf")
plot(pca_results, main = "PCA Plot")
dev.off()

pdf("results/dimensionality_reduction/tsne_plot.pdf")
plot(tsne_results, main = "t-SNE Plot")
dev.off()

pdf("results/dimensionality_reduction/umap_plot.pdf")
plot(umap_results, main = "UMAP Plot")
dev.off()