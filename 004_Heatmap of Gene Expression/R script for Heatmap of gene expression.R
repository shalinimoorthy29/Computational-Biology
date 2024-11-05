# R Script for heatmap creation
# GEO Accession ID: GSE46056
# Required file type: GSE46056_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz
# URL: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE46056

install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "circlize", "biomaRt"))

# Load required libraries
library(ComplexHeatmap)
library(circlize)
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
log_fpkm_data$entrezgene_id <- rownames(log_fpkm_data)
log_fpkm_data <- merge(log_fpkm_data, gene_annotations, by = "entrezgene_id", all.x = TRUE)

# 5. Remove rows with missing or empty gene symbols
log_fpkm_data <- log_fpkm_data[!is.na(log_fpkm_data$hgnc_symbol) & log_fpkm_data$hgnc_symbol != "", ]

# Ensure no numeric-only gene symbols are present and remove them
log_fpkm_data <- log_fpkm_data[!grepl("^[0-9]+$", log_fpkm_data$hgnc_symbol), ]

# 6. Replace row names with unique gene symbols
rownames(log_fpkm_data) <- make.unique(log_fpkm_data$hgnc_symbol)
log_fpkm_data <- log_fpkm_data[, !(colnames(log_fpkm_data) %in% c("entrezgene_id", "hgnc_symbol"))]

# 7. Select top 10% most variable genes based on standard deviation
gene_sd <- apply(log_fpkm_data, 1, sd)        # Calculate standard deviation for each gene
top_genes <- order(gene_sd, decreasing = TRUE) # Rank genes by variability
top_10_percent <- top_genes[1:(length(gene_sd) * 0.1)] # Select top 10%

# Subset data to include only the top 10% most variable genes
subset_data <- log_fpkm_data[top_10_percent, ]

# Apply z-score normalization
scaled_data <- t(scale(t(subset_data)))
scaled_data <- scaled_data[complete.cases(scaled_data), ]

# Load condition data
col_data <- read.csv("col_data.csv")
conditions <- col_data$Condition  # Extract the Condition column for annotation

# Define the condition annotation using the 'Condition' column
column_ha <- HeatmapAnnotation(
  condition = conditions,
  col = list(condition = c("PRMT4-Knockdown" = "red", "Control" = "blue")),
  annotation_legend_param = list(condition = list(title = "Condition"))
)

Heatmap(scaled_data,
        name = "Expression",
        col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
        top_annotation = column_ha,
        show_row_names = FALSE,      # Hide gene names if too many rows
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_side = "right",
        column_names_side = "bottom",
        heatmap_legend_param = list(title = "Scaled Normalized Expression"))
