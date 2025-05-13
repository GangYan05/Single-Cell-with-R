
all.sce <- list(
    pbmc3k=TENxPBMCData('pbmc3k'),
    pbmc4k=TENxPBMCData('pbmc4k'),
    pbmc8k=TENxPBMCData('pbmc8k')
)
all.sce.filtered <- list()

for (name in names(all.sce)) {
    cat("Processing:", name)
    # Get the current SingleCellExperiment object
    current.sce <- all.sce[[name]]
    # Identify mitochondrial genes using standard gene symbols
    symbol_col <- NULL
    if ("Symbol_TENx" %in% colnames(rowData(current.sce))) {
        symbol_col <- "Symbol_TENx"
    } else if ("Symbol" %in% colnames(rowData(current.sce))) {
        symbol_col <- "Symbol"
    } else {
        warning(paste("Could not find standard gene symbol column (Symbol or Symbol_TENx) in", name, "- Skipping mitochondrial QC."))
        is_mito <- rep(FALSE, nrow(current.sce)) # Create dummy vector
    }
    if (!is.null(symbol_col)) {
        # Using "^MT-" which is standard for human mitochondrial genes
        is_mito <- grepl("^MT-", rowData(current.sce)[[symbol_col]], ignore.case = TRUE) 
        cat("  Found", sum(is_mito), "mitochondrial genes using column:", symbol_col, "\n")
        if (sum(is_mito) == 0) {
            warning(paste("No mitochondrial genes found matching pattern '^MT-' in", name))
            }
    }
    # Calculate quality control metrics
    qc_subset <- list()
    if (sum(is_mito) > 0) {
        qc_subset$Mito <- is_mito
    }
    qc_df <- perCellQCMetrics(current.sce, subsets = qc_subset)
    # Apply adaptive filtering for low-quality cells
    quc_sub_fields <- c("subsets_Mito_percent")
    if ("subsets_Mito_percent" %in% colnames(qc_df)) {
    qc_sub_fields <- c("subsets_Mito_percent")
    cat("  Applying adaptive filtering on library size, features, and mitochondrial percentage.\n")
    } else {
    cat("  Applying adaptive filtering on library size and features only (no mito metrics found/calculated).\n")
    }
    filtered_results <- perCellQCFilters(
        qc_df, 
        sub.fields = qc_sub_fields 
    )
    # Filter the SingleCellExperiment object to retain only high-quality cells
    current.sce.filtered <- current.sce[, !filtered_results$discard]
    # Store the filtered object in the list
    all.sce.filtered[[name]] <- current.sce.filtered
}

# Find intersection of features across all datasets
cat("\nFinding intersection of features across all datasets...\n")
gene_lists <- lapply(all.sce.filtered, rownames)
# Find the common genes
intersect_genes <- Reduce(intersect, gene_lists)
# Check if any genes remain
if (length(intersect_genes) == 0) {
    stop("No common genes found across datasets.")
} else {
    cat("  Found", length(intersect_genes), "common genes.\n")
}

# Subset sce list to only include common genes
all.sce.intersected <- lapply(all.sce.filtered, function(sce) {
    sce[intersect_genes, ]
})
# Verify dimensions (rows should now be identical)
print("Dimensions of intersected SCE objects:")
print(lapply(all.sce.intersected, dim))

# Normalize and log-transform the filtered SingleCellExperiment objects
all.sce.intersected <- lapply(all.sce.intersected, logNormCounts)

# Variance-modeling for feature selection
all.dec <- lapply(all.sce.intersected, modelGeneVar)
all.hvgs <- lapply(all.dec, getTopHVGs, prop=0.1)

# Combine HVGs from all datasets
cat("Combining HVGs from all datasets...\n")
all.hvgs <- unique(unlist(all.hvgs))

# Add batch information to colData
cat("Adding batch information to colData...\n")
for (name in names(all.sce.intersected)) {
    if(is.null(colData(all.sce.intersected[[name]])$batch)) {
         colData(all.sce.intersected[[name]])$batch <- name
         cat("  Added batch label '", name, "'\n", sep="")
    } else {
         cat("  Batch label already exists for '", name, "'\n", sep="")
    }
}

# Run fastMNN
cat("Running fastMNN for batch correction...\n")
set.seed(10101010) # Set seed for reproducibility

sce.corrected <- fastMNN(
    all.sce.intersected, # Input list of SCEs with common PCA
    batch = "batch",    # Specify the batch variable
    subset.row = all.hvgs, # Specify the genes used for PCA
    d = 50,             # Number of PCs to use for MNN
    k = 20,             # Number of neighbors for MNN graph 
    correct.all = TRUE  # Get corrected PCs for all cells
    # assay.type = "logcounts" # Specify assay used for PCA/variance calculation
)
cat("  fastMNN complete.\n")
















