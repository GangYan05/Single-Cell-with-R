library(scater)

# Load the single-cell experiment object
sce <- readRDS("results/SingleCellExperiment.rds")

# Library size factor normalization
lib_sf <- librarySizeFactors(sce)
hist(log10(lib_sf), breaks = 10, col = "grey", xlab = "log10(Library size factor)", main = "Library size factor distribution")

# Deconvolution normalization
set.seed(100)
clust_sce <- quickCluster(sce, min.size=20)
summary(clust_sce)
deconv_sf <- calculateSumFactors(sce, cluster=clust_sce)
summary(sce_deconv)

# Plot the library size factors against the deconvolution size factors
plot(lib_sf, deconv_sf, xlab="Library size factor",
    ylab="Deconvolution size factor", log='xy', pch=16)
abline(a=0, b=1, col="red")

# Library size normalization and log transformation
sce_lib <- logNormCounts(sce)

# Deconvolution normalization and log transformation
sce_deconv <- computeSumFactors(sce, cluster=clust_sce)
sce_count <- logNormCounts(sce_deconv)

# Save the normalized data
saveRDS(sce_lib, file = "results/normalized_data/lib_normalized_sce.rds")
saveRDS(sce_deconv, file = "results/normalized_data/deconv_normalized_sce.rds")
