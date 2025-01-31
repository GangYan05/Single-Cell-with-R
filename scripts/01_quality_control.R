# Define the packages to install
bioc_packages <- c("SingleCellExperiment", "scuttle", "scran", "scater", "uwot", 
                   "rtracklayer", "DropletUtils", "batchelor", "bluster", "ensembldb", 
                   "org.Mm.eg.db", "org.Hs.eg.db", "DropletTestFiles", "scRNAseq", "AnnotationHub")
cran_packages <- c("uwot", "dynamicTreeCut")

# install and load the packages
install_and_load_packages(bioc_pkgs = bioc_packages, cran_pkgs = cran_packages)

# Data path
raw_fp <- "data/raw_data/E-MTAB-5522/counts_Calero_20160113.tsv"
meta_fp <- "data/raw_data/E-MTAB-5522/E-MTAB-5522.sdrf.txt" 

# Load the single-cell data and metadata
sce_count <- load_data_from_tabular(raw_fp, sparse = TRUE)


# Create a SingleCellExperiment object
sce <- create_single_cell_object(sce_count)
sce <- scuttle::logNormCounts(sce)

# Add metadata to the SingleCellExperiment object
sce <- add_metadata_from_file(sce, meta_fp, "counts_Calero_20160113.tsv")
sce <- addPerCellQC(sce)
sce <- addPerFeatureQC(sce)

# Add gene annotation
ah <- AnnotationHub()
ens_mm <- ah[["AH75036"]]
gene_info <- AnnotationDbi::select(ens_mm, keys = rownames(sce), keytype = "GENEID", columns = c("SYMBOL", "SEQNAME", "GENEBIOTYPE"))
rowData(sce)$SYMBOL <- gene_info$SYMBOL[match(rownames(sce), gene_info$GENEID)]
rowData(sce)$SEQNAME <- gene_info$SEQNAME[match(rownames(sce), gene_info$GENEID)]
rowData(sce)$GENEBIOTYPE <- gene_info$GENEBIOTYPE[match(rownames(sce), gene_info$GENEID)]


is_mito_symbol <- grepl("^mt-", rowData(sce)$SYMBOL)
is_mito_location <- rowData(sce)$SEQNAME == "MT"

# Identify Mitochondrial genes using symbols
is_mito_symbol <- grepl("^mt-", rowData(sce)$SYMBOL)
is_mito_location <- rowData(sce)$SEQNAME == "MT"
is_spike <- grepl("^ERCC-", rowData(sce)$SYMBOL)


gene_anno[gene_anno$type == "gene"]
# names(gene_info) <- gene_info$gene_id
# gene_revelant <- grep("gene_", colnames(mcols(gene_info)))
# mcols(gene_info) <- mcols(gene_info)[, gene_revelant]
# rowRanges(sce) <- gene_info[rownames(sce),]

first_5 <- sce[,1:5]
colData(first_5)

# Identify mitochondrial genes for quality control
is.mito <- grepl("^MT-", rownames(sce))

# Calculate quality control metrics
qcstats <- perCellQCMetrics(sce)

# Filter cells based on quality metrics
filtered <- quickPerCellQC(qcstats, percent_subsets = "subsets_Mito_percent")

# Subset the SingleCellExperiment object to retain only high-quality cells
sce <- sce[, !filtered$discard]

# Save quality control metrics to results directory
saveRDS(qcstats, file = "results/qc_metrics/qc_stats.rds")