library(scater)

# Load the single-cell experiment object
sce <- readRDS("path/to/your/sce_object.rds")

# Feature selection based on variance
# Calculate the mean and variance for each gene
mean_var <- rowData(sce) %>%
  summarize(mean = rowMeans(logcounts(sce)), 
            variance = rowVars(logcounts(sce)))

# Select features with high variance
high_var_genes <- mean_var[mean_var$variance > quantile(mean_var$variance, 0.75), ]

# Subset the single-cell experiment object to keep only high variance genes
sce_filtered <- sce[rownames(sce) %in% rownames(high_var_genes), ]

# Save the filtered results
saveRDS(sce_filtered, file = "results/feature_selection/filtered_sce.rds")

# Optionally, visualize the selected features
plot(rowMeans(logcounts(sce)), rowVars(logcounts(sce)), 
     xlab = "Mean Expression", ylab = "Variance", 
     main = "Feature Selection: Mean vs Variance")
points(rowMeans(logcounts(sce_filtered)), rowVars(logcounts(sce_filtered)), col = "red")