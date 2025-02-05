library(scater)

# Load the single-cell experiment object
sce <- readRDS("results/SingleCellExperiment.rds")

# Library size factor normalization
lib_sf <- librarySizeFactors(sce)
hist(log10(sce_libsize), breaks = 10, col = "grey", xlab = "log10(Library size factor)", main = "Library size factor distribution")

# Deconvolution normalization
set.seed(100)
clust_sce <- quickCluster(sce, min.size=20)
summary(clust_sce)
sce <- calculateSumFactors(sce, cluster=clust_sce)
summary(sce_deconv)


# Plot the library size factors against the deconvolution size factors
plot(lib_sf, deconv_sf, xlab="Library size factor",
    ylab="Deconvolution size factor", log='xy', pch=16)
abline(a=0, b=1, col="red")

# Scaling and log normalization of the count data
logNormCounts(sce)

# Save the normalized data
saveRDS(sce, file = "results/normalized_data/normalized_sce.rds")