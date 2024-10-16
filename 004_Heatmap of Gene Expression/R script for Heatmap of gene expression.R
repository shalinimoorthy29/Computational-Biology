# Load required libraries
library(pheatmap)
library(biomaRt)

# 1. Load the FPKM data from the file
fpkm_data <- read.table("GSE46056_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", sep = "\t", header = TRUE, row.names = 1)

# 2. Log-transform the FPKM data (adding 1 to avoid log(0))
log_fpkm_data <- log2(fpkm_data + 1)

# 3. Use biomaRt to map Entrez Gene IDs to HGNC gene symbols
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract Entrez Gene IDs from the row names of the data
ncbi_gene_ids <- rownames(fpkm_data)

# Query biomaRt to get Gene Symbols using NCBI Gene IDs
gene_annotations <- getBM(attributes = c('entrezgene_id', 'hgnc_symbol'),
                          filters = 'entrezgene_id',
                          values = ncbi_gene_ids,
                          mart = ensembl)

# Debug: Print the first few rows of gene_annotations to check mappings
print(head(gene_annotations))

# 4. Merge FPKM data with gene_annotations
# Convert fpkm_data rownames (NCBI Gene IDs) into a column for merging
fpkm_data$entrezgene_id <- rownames(fpkm_data)

# Merge FPKM data with gene_annotations to get Gene Symbols
fpkm_data <- merge(fpkm_data, gene_annotations, by = "entrezgene_id", all.x = TRUE)

# Debug: Check if merge was successful and gene symbols are included
print(head(fpkm_data))

# 5. Handle missing gene symbols
# Remove rows with missing or empty gene symbols
fpkm_data <- fpkm_data[!is.na(fpkm_data$hgnc_symbol) & fpkm_data$hgnc_symbol != "", ]

# Ensure no numeric-only gene symbols are present and remove them
fpkm_data <- fpkm_data[!grepl("^[0-9]+$", fpkm_data$hgnc_symbol), ]

# 6. Replace row names with unique gene symbols
rownames(fpkm_data) <- make.unique(fpkm_data$hgnc_symbol)

# Remove the 'entrezgene_id' and 'hgnc_symbol' columns to clean up
fpkm_data <- fpkm_data[, !(colnames(fpkm_data) %in% c("entrezgene_id", "hgnc_symbol"))]

# 7. Subset the data for visualization (e.g., top 100 genes)
subset_data <- fpkm_data[1:100, ]

# Apply z-score normalization (scaling by rows, i.e., genes)
scaled_data <- t(scale(t(subset_data)))  # Transpose, scale, and transpose back

# Check for any rows with NA/NaN/Inf values and remove them
scaled_data <- scaled_data[complete.cases(scaled_data), ]

# 8. Generate heatmap with scaled (z-score) data
pheatmap(scaled_data, 
         cluster_rows = TRUE,    # Cluster genes (rows)
         cluster_cols = TRUE,    # Cluster samples (columns)
         show_rownames = TRUE,   # Show gene names (row names)
         show_colnames = TRUE,   # Show sample names (column names)
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Color scale from blue to red
         main = "Heatmap of Z-Score Normalized Gene Expression (FPKM with Gene Symbols)")
