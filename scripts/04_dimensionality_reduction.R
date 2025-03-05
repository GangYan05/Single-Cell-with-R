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

# PCA plot
plotReducedDim(sce, dimred="PCA", color_by = "level1class")

# Perform t-SNE
sce <- scater::runTSNE(sce, perplexity = 30)
# t-SNE plot
plotReducedDim(sce, dimred="TSNE", color_by = "level1class")

# Perform UMAP
sce <- runUMAP(sce)
plotReducedDim(sce, dimred="UMAP", color_by = "level1class")
