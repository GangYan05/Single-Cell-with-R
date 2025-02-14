library(scater)

# Load the single-cell experiment object
sce_lib <- readRDS(file = "results/normalized_data/lib_normalized_sce.rds")
sce_deconv <- readRDS(file = "results/normalized_data/deconv_normalized_sce.rds")

sce_model <- modelGeneVar(sce_lib)
sce_fit <- metadata(sce_model)


plot(sce_fit$mean, sce_fit$var, xlab="Mean of log-expression",
    ylab="Variance of log-expression")
curve(sce_fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)


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