# R Script to perform PCA analysis
# GEO Accession ID: GSE46056
# Required file type: GSE46056_raw_counts_GRCh38.p13_NCBI.tsv.gz
# URL: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE46056

install.packages("BiocManager")
BiocManager::install("DESeq2")  # DESeq2 for RNA-Seq analysis and PCA
BiocManager::install("ggplot2")  # For plotting

library(DESeq2)
library(ggplot2)

# 1) Load counts data (genes as rows and samples as columns).
counts_data <- read.csv('counts_data.csv', row.names = 1)

# 2) Load metadata (sample information, e.g., conditions like Knockdown/Control).
col_data <- read.csv('col_data.csv', row.names = 1)

# 3) Ensure column names of counts_data match the row names of col_data.
all(colnames(counts_data) == rownames(col_data))  # This should return TRUE

# 4) Create DESeq2 dataset.
dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                              colData = col_data, 
                              design = ~ Knockdown)  # Replace 'Knockdown' with your specific condition

# 5) Filter low count genes. 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 6) Perform variance stabilising transformation.
vsd <- vst(dds, blind = FALSE)

# 7) Generate PCA data.
pca_data <- plotPCA(vsd, intgroup = "Knockdown", returnData = TRUE)  
percent_var <- round(100 * attr(pca_data, "percentVar"))


# 8) Add sample names or IDs as labels.
pca_data$SampleID <- rownames(pca_data)  # If you have sample names in metadata, you can adjust this

# 9) Create PCA plot with sample labels.
# Assuming you don't have a column named 'Sample', use the rownames of pca_data to create it.
pca_data$Sample <- rownames(pca_data)  # Assign sample names or IDs from row names

# 10) Now plot with the correct sample labels.
ggplot(pca_data, aes(PC1, PC2, color = Knockdown)) +  # Remove 'label = Sample' from aes()
  geom_point(size = 3) +
  geom_text(aes(label = Sample), vjust = -0.5, hjust = 0.5, show.legend = FALSE) +  # Add labels without affecting legend
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_minimal() +
  labs(title = "PCA Plot with Sample Labels") +
  coord_cartesian(xlim = range(pca_data$PC1) + c(-2, 2),  # Extend X axis limits
                  ylim = range(pca_data$PC2) + c(-2, 2)) +  # Extend Y axis limits
  # Center the plot title
  theme(plot.title = element_text(hjust = 0.5))  # Center the title
