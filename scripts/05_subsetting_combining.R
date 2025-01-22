sce <- scuttle::addPerCellQC(sce)

# Subsetting the first 10 cells
first.10 <- sce[, 1:10]
ncol(counts(first.10))
colData(first.10)

# Subsetting to include only protein-coding genes
coding.only <- sce[rowData(sce)$gene_biotype == "protein_coding", ]
counts(coding.only)

# Combining datasets
sce2 <- cbind(sce, sce)
ncol(counts(sce2))
sce2 <- rbind(sce, sce)
nrow(counts(sce2))