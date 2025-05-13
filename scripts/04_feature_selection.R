# Load the single-cell experiment object
sce_lib <- readRDS(file = "results/normalized_data/lib_normalized_sce.rds")
sce_deconv <- readRDS(file = "results/normalized_data/deconv_normalized_sce.rds")

# Quality per-gene variation
sce_model <- modelGeneVar(sce_lib)
sce_fit <- metadata(sce_model)


# Plot the mean-variance trend
plot(sce_fit$mean, sce_fit$var, xlab="Mean of log-expression",
    ylab="Variance of log-expression")
curve(sce_fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)

sce_model[order(sce_model$bio, decreasing=TRUE), ]

# Select features with high variance
sce_variable <- getTopHVGs(sce_model, fdr.threshold = 0.05)

# Save the filtered results
saveRDS(sce_variable, file = "results/feature_selection/sce_variable.rds")

# Optionally, visualize the selected features
plot(rowMeans(logcounts(sce)), rowVars(logcounts(sce)), 
     xlab = "Mean Expression", ylab = "Variance", 
     main = "Feature Selection: Mean vs Variance")
points(rowMeans(logcounts(sce_filtered)), rowVars(logcounts(sce_filtered)), col = "red") 