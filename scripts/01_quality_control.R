# Define the packages to install
bioc_packages <- c("SingleCellExperiment", "scuttle", "scran", "scater", "uwot", 
                   "rtracklayer", "DropletUtils", "batchelor", "bluster", "ensembldb", 
                   "org.Mm.eg.db", "org.Hs.eg.db", "DropletTestFiles", "scRNAseq", "AnnotationHub",
                   "PCAtools", "celldex", "SingleR")
cran_packages <- c("uwot", "dynamicTreeCut", "dplyr", "pheatmap")

# install and load the packages
install_and_load_packages(bioc_pkgs = bioc_packages, cran_pkgs = cran_packages)

# Data path
raw_fp <- "data/raw_data/E-MTAB-5522/counts_Calero_20160113.tsv"
meta_fp <- "data/raw_data/E-MTAB-5522/E-MTAB-5522.sdrf.txt" 

# Load the single-cell data and metadata
sce_count <- load_data_from_tabular(raw_fp, sparse = TRUE)

# Create a SingleCellExperiment object
sce <- create_single_cell_object(sce_count)

# Add metadata to the SingleCellExperiment object
sce <- add_metadata_from_file(sce, meta_fp, "counts_Calero_20160113.tsv")
# sce <- addPerCellQC(sce)
# sce <- addPerFeatureQC(sce)

# Add gene annotation
ah <- AnnotationHub()
ens_mm <- ah[["AH75036"]]
gene_info <- AnnotationDbi::select(ens_mm, keys = rownames(sce), keytype = "GENEID", columns = c("SYMBOL", "SEQNAME", "GENEBIOTYPE"))
rowData(sce)$SYMBOL <- gene_info$SYMBOL[match(rownames(sce), gene_info$GENEID)]
rowData(sce)$SEQNAME <- gene_info$SEQNAME[match(rownames(sce), gene_info$GENEID)]
rowData(sce)$GENEBIOTYPE <- gene_info$GENEBIOTYPE[match(rownames(sce), gene_info$GENEID)]

# Identify Mitochondrial genes using symbols
is_mito <- grepl("^mt-", rowData(sce)$SYMBOL)

# Calculate quality control metrics
qc_df <- perCellQCMetrics(sce, subsets = list(Mito = is_mito))

# Identify low-quality cells with fixed thresholds (optional)
qc_size <- qc_df$sum < 1e5
qc_nexpression <- qc_df$detected < 5e3
qc_spike <- qc_df$detected > 1e4
qc_mito <- qc_df$subsets_Mito_percent > 10
todrop <- qc_size | qc_nexpression | qc_spike | qc_mito
DataFrame(LibSize=sum(qc_size), NExprs=sum(qc_nexpression), Spike=sum(qc_spike), Mito=sum(qc_mito), total=sum(todrop))

# Adaptively filter out low-quality cells
filtered <- perCellQCFilters(qc_df ,sub.fields = c("subsets_Mito_percent"))

# Subset the SingleCellExperiment object to retain only high-quality cells
sce <- sce[, !filtered$discard]

# Save quality control metrics to results directory
saveRDS(qc_df, file = "results/qc_metrics/qc_stats.rds")
saveRDS(sce, file = "results/SingleCellExperiment.rds")
