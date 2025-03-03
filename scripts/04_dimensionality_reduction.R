# Load the single-cell experiment object
sce <- ZeiselBrainData()
sce <- aggregateAcrossFeatures(sce, id=sub("_loc[0-9]+$", "", rownames(sce)))

# Gene annotation
rowData(sce)$Ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")

# Quality control
stats <- perCellQCMetrics(sce, subsets=list(Mt=rowData(sce)$featureType=="mito"))
qc <- quickPerCellQC(stats, percent_subsets=c("altexps_ERCC_percent", "subsets_Mt_percent"))
sce <- sce[, !qc$discard]

# Normalization
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusters)
sce <- logNormCounts(sce)

# Feature selection
dec.sce <- modelGeneVarWithSpikes(sce, "ERCC")
top.hvgs <- getTopHVGs(dec.sce, prop=0.1)

# Perform PCA
sce <- scater::runPCA(sce, subset_row=top.hvgs)
dim(reducedDim(sce))

# Choose the number of PCs to use
percent.var <- attr(reducedDim(sce, "PCA"), "percentVar")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
chosen.pcs <- PCAtools::findElbowPoint(percent.var)
chosen.pcs

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