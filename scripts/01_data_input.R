# Define the packages to install
bioc_packages <- c("SingleCellExperiment", "scuttle", "scran", "scater", "uwot", 
                   "rtracklayer", "DropletUtils", "batchelor", "bluster", "ensembldb", 
                   "org.Mm.eg.db", "org.Hs.eg.db", "DropletTestFiles", "scRNAseq", "AnnotationHub",
                   "PCAtools", "celldex", "SingleR", "TENxPBMCData")
cran_packages <- c("uwot", "dynamicTreeCut", "dplyr", "pheatmap", "Seurat")

# install and load the packages
install_and_load_packages(bioc_pkgs = bioc_packages, cran_pkgs = cran_packages)

# Data path of E-MTAB-5522
raw_fp <- "data/raw_data/E-MTAB-5522/counts_Calero_20160113.tsv"
meta_fp <- "data/raw_data/E-MTAB-5522/E-MTAB-5522.sdrf.txt" 

# Data path of GSE176078
# R.utils::gzip("data/raw_data/GSE176078/barcodes.tsv")
# R.utils::gzip("data/raw_data/GSE176078/features.tsv")
# R.utils::gzip("data/raw_data/GSE176078/matrix.mtx")
# data_dir <- "data/raw_data/GSE176078"

# Load the single-cell data and metadata
count <- read.delim(raw_fp, header = TRUE, row.names = 1, check.names = FALSE)

# Separate feature lengths and count matrix
all_feature_lengths <- count[,1] # First column is lengths
names(all_feature_lengths) <- rownames(count)
mat <- as.matrix(count[,-1])



# Identify endogenous gene rows (e.g., those starting with "ENSMUSG")
is_endogenous <- grepl("^ENSMUSG", rownames(mat), ignore.case = TRUE)
endogenous_feature_names <- rownames(mat)[is_endogenous]

# Identify ERCC spike-in rows (e.g., those starting with "ERCC-")
is_ercc <- grepl("^ERCC-", rownames(mat), ignore.case = TRUE)
ercc_feature_names <- rownames(mat)[is_ercc]

# Create count matrices for main and alternative experiments
main_counts <- mat[endogenous_feature_names, ]
ercc_counts <- mat[ercc_feature_names, ]

# Create a SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = main_counts))

# Add gene lengths to the SingleCellExperiment object
rowData(sce)$Length <- all_feature_lengths[endogenous_feature_names]

# If ERCC features exist, create a SummarizedExperiment for them and add as altExp
if (length(ercc_feature_names) > 0) {
    alt_se_ercc <- SummarizedExperiment(assays = list(counts = ercc_counts))
    # Optionally, add rowData (like lengths) to the alternative experiment
    rowData(alt_se_ercc)$Length <- all_feature_lengths[ercc_feature_names]
    
    # Add the ERCC data as an alternative experiment NAMED "ERCC"
    altExp(sce, "ERCC") <- alt_se_ercc
    message(paste("Added ERCCs as an alternative experiment named 'ERCC' with", 
                  nrow(altExp(sce, "ERCC")), "features."))
    # You can check the names:
    # print(altExpNames(sce)) # Should output "ERCC"
} else {
    message("No ERCC features found to add as an alternative experiment.")
}

# Import the metadata
coldata <- read.delim(meta_fp, check.names=FALSE)
coldata <- coldata[coldata[,"Derived Array Data File"]=="counts_Calero_20160113.tsv",]
coldata <- DataFrame(genotype=coldata[, "Characteristics[genotype]"], 
                    phenotype=coldata[, "Characteristics[phenotype]"], 
                    spike_in=coldata[, "Factor Value[spike-in addition]"], 
                    row.names = coldata[, "Source Name"])

# Check if the metadata file has the same number of cells as the count matrix
if (ncol(sce) != nrow(coldata)) {
    stop("The number of cells in the metadata file does not match the count matrix.")
} else {
    message("The number of cells in the metadata file matches the count matrix.")
}

# Add metadata to the SingleCellExperiment object
colData(sce) <- coldata

# AnnotationHub object to get gene information
ah <- AnnotationHub()
ens_mm <- ah[["AH75036"]]
gene_info <- AnnotationDbi::select(ens_mm, keys = rownames(sce), keytype = "GENEID", columns = c("SYMBOL", "SEQNAME", "GENEBIOTYPE"))
rowData(sce)$SYMBOL <- gene_info$SYMBOL[match(rownames(sce), gene_info$GENEID)]
rowData(sce)$SEQNAME <- gene_info$SEQNAME[match(rownames(sce), gene_info$GENEID)]
rowData(sce)$GENEBIOTYPE <- gene_info$GENEBIOTYPE[match(rownames(sce), gene_info$GENEID)]

# Identify Mitochondrial genes using symbols
is_mito <- grepl("^mt-", rowData(sce)$SYMBOL)

# Save quality control metrics to results directory
saveRDS(sce, file = "results/sce.rds")
saveRDS(coldata, file = "results/coldata.rds")
saveRDS(is_mito, file = "results/is_mito.rds")



