# Load required libraries
library(pheatmap)
library(biomaRt)

# 1. Load the FPKM data from the file
fpkm_data <- read.table("GSE46056_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", sep = "\t", header = TRUE, row.names = 1)

# 2. Log-transform the FPKM data (adding 1 to avoid log(0))
log_fpkm_data <- log2(fpkm_data + 1)

# 3. Use biomaRt to map Entrez Gene IDs to HGNC gene symbols
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract Entrez Gene IDs from the row names of the log-transformed data
ncbi_gene_ids <- rownames(log_fpkm_data)

# Query biomaRt to get Gene Symbols using NCBI Gene IDs
gene_annotations <- getBM(attributes = c('entrezgene_id', 'hgnc_symbol'),
                          filters = 'entrezgene_id',
                          values = ncbi_gene_ids,
                          mart = ensembl)

# 4. Merge log-transformed FPKM data with gene_annotations
# Convert log_fpkm_data rownames (NCBI Gene IDs) into a column for merging
log_fpkm_data$entrezgene_id <- rownames(log_fpkm_data)

# Merge log_fpkm_data with gene_annotations to get Gene Symbols
log_fpkm_data <- merge(log_fpkm_data, gene_annotations, by = "entrezgene_id", all.x = TRUE)

# 5. Handle missing gene symbols
# Remove rows with missing or empty gene symbols
log_fpkm_data <- log_fpkm_data[!is.na(log_fpkm_data$hgnc_symbol) & log_fpkm_data$hgnc_symbol != "", ]

# Ensure no numeric-only gene symbols are present and remove them
log_fpkm_data <- log_fpkm_data[!grepl("^[0-9]+$", log_fpkm_data$hgnc_symbol), ]

# 6. Replace row names with unique gene symbols
rownames(log_fpkm_data) <- make.unique(log_fpkm_data$hgnc_symbol)

# Remove the 'entrezgene_id' and 'hgnc_symbol' columns to clean up
log_fpkm_data <- log_fpkm_data[, !(colnames(log_fpkm_data) %in% c("entrezgene_id", "hgnc_symbol"))]

# 7. Subset the log-transformed data for visualization (e.g., top 100 genes)
subset_data <- log_fpkm_data[1:100, ]

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
         main = "Heatmap of Z-Score Normalized Gene Expression (Log FPKM with Gene Symbols)")
