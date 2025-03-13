
# Find marker genes for each cluster using pairwise comparisons
markers <- findMarkers(sce)

# Get the top 5 markers for cluster 7
cluster7 <- markers[[7]]
best5 <- cluster7[cluster7$Top <=5,]
lfc_cluster7 <- getMarkerEffects(best5)

# plot the expression of the top 5 markers of each pairwise comparison for cluster 7
pheatmap(lfc_cluster7, breaks=seq(-5, 5, length.out=101))

# Find cluster-specific markers
markers_up <- findMarkers(sce, pval.type="all", direction="up")
cluster7_up <- markers_up[[7]]
cluster7_up[1:10, 1:3]

# Less stringent criteria
markers_up <- findMarkers(sce, pval.type="some", direction="up")
cluster7_up <- markers_up[[7]]
cluster7_up[1:10, 1:3]

# Focusing on upregulated genes with a log-fold change of at least 1
markers_up <- findMarkers(sce, direction="up", lfc=1)
cluster7_up <- markers_up[[7]]
best10 <- cluster7_up[cluster7_up$Top <=10,]
lfc_cluster7_up <- getMarkerEffects(best10)
pheatmap(lfc_cluster7_up, breaks=seq(-5, 5, length.out=101))

# Alternative Wilcoxon rank-sum test
markers_wmw <- findMarkers(sce, test="wilcox", direction="up")
cluster7_wmw <- markers_wmw[[7]]
best10_wmw <- cluster7_wmw[cluster7_wmw$Top <=10,]
AUC_cluster7_wmw <- getMarkerEffects(best10_wmw, prefix="AUC")
pheatmap(AUC_cluster7_wmw, breaks=seq(0, 1, length.out=101))
