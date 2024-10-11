# Load required libraries for data processing and visualization
library(pheatmap)

# Load the FPKM data from the file
fpkm_data <- read.table("GSE46056_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", sep = "\t", header = TRUE, row.names = 1)

# Log-transform the FPKM data (adding 1 to avoid log(0))
log_fpkm_data <- log2(fpkm_data + 1)

# Subset the data (e.g., top 50 genes)
subset_data <- log_fpkm_data[1:50, ]

# Apply z-score normalization (scaling by rows, i.e., genes)
scaled_data <- t(scale(t(subset_data)))  # Transpose, scale, and transpose back

# Check for any rows with NA/NaN/Inf values and remove them
scaled_data <- scaled_data[complete.cases(scaled_data), ]

# Generate heatmap with scaled (z-score) data
pheatmap(scaled_data, 
         cluster_rows = TRUE,    # Cluster genes (rows)
         cluster_cols = TRUE,    # Cluster samples (columns)
         show_rownames = TRUE,   # Show gene names (row names)
         show_colnames = TRUE,   # Show sample names (column names)
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Color scale from blue to red
         main = "Heatmap of Z-Score Normalized Gene Expression (FPKM)")
