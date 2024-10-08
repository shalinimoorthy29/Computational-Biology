# R Script to perform differential gene expression analysis using DESeq2
# GEO Accession ID: GSE46056
# Required file type: GSE46056_raw_counts_GRCh38.p13_NCBI.tsv.gz
# URL: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE46056

# Prerequisite:
# Set working directory.
# Download file and extract from zip folder.
# Copy paste the data in an excel sheet and export as csv.
# Transfer file to working directory.
# Create another csv file with the GSM IDs, sample names and condition
# NOTE - GSM IDs which would be seen as column headers on the counts data
# should be in the same order, but as rows in the column data. 
# Create it using excel sheet, export as csv file.
# Transfer file to working directory

# -----------------------------------------------------------------------------
# Package description:

# DESeq2 is a package specifically designed for analyzing count data from
# high-throughput sequencing assays, particularly RNA-Seq data.

# The tidyverse is a collection of R packages designed for data science
# providing a cohesive framework for data manipulation, exploration, 
# and visualization.

# The readr package is part of the tidyverse and is specifically designed for 
# reading and writing data files efficiently.


# 1) Load required R packages
install.packages("BiocManager")
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(readr)

# 2) Read in counts data to R environment
counts_data <- read.csv('counts_data.csv', row.names = 1)
head(counts_data)

# 3) Read in column data (Sample information)
col_data <- read.csv('col_data.csv', row.names = 1)

# 4) Ensure the row names in col_data matches to 
# column names in counts_data as well as in the same order
all(colnames(counts_data) %in% rownames(col_data))
all(colnames(counts_data) == rownames(col_data))

# 5) Set up a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = col_data,
                       design = ~ Knockdown)

# 6) Remove rows with low gene counts (less than 10)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 7) Set the factor level to enable comparison between the conditions
dds$Knockdown <- relevel(dds$Knockdown, ref = "Untreated")

# NOTE: If there are technical replicates in the dataset, make sure to 
# collapse them. But for this dataset, it is not required.

# 8) Run DEseq and get result
dds <- DESeq(dds)
result <- results(dds)
summary(result)

# 9) Get significant result with an adjusted p-value threshold of 0.01
# instead of 0.1.
significant_result <- results(dds, alpha = 0.01)
summary(significant_result)

# Save significant results to CSV
write.csv(as.data.frame(significant_result), "significant_results.csv")

# 10) Plot MA plots for result and significant result
plotMA(significant_result, main="MA Plot")

# 11) Volcano Plot Preparation
# Convert the results to a data frame and calculate -log10(p-value)
result_df <- as.data.frame(result)
result_df <- result_df %>%
  mutate(log10_pvalue = -log10(pvalue)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1)

# 12) Generate Volcano Plot using ggplot2
volcano_plot <- ggplot(result_df, aes(x = log2FoldChange, y = log10_pvalue)) +
  geom_point(aes(color = significant), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(volcano_plot)