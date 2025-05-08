# Graph-based clustering
g <- buildSNNGraph(sce, k=10, use.dimred="PCA")
clust <- igraph::cluster_walktrap(g)$membership
table(clust)

# Plot the t-SNE plot with the clusters
colLabels(sce) <- factor(clust)
plotReducedDim(sce, "TSNE", colour_by = "label")

# K-means clustering
set.seed(100)
clust_kmean <- kmeans(reducedDim(sce, "PCA"), centers=10)
table(clust_kmean$cluster)
# Plot the t-SNE plot with the clusters
colLabels(sce) <- factor(clust_kmean$cluster)
plotReducedDim(sce, "TSNE", colour_by = "label")

# Hierarchical clustering
dist <- dist(reducedDim(sce, "PCA"))
tree <- hclust(dist, "ward.D2")
# Plot the dendrogram
dend <- as.dendrogram(tree)
plot(dend)
clust_hclust <- cutreeDynamic(tree, distM=as.matrix(dist), minClusterSize=10, deepSplit=1)
# PLot the t-SNE plot with the clusters
colLabels(sce) <- factor(clust_hclust)
plotReducedDim(sce, "TSNE", colour_by = "label")
