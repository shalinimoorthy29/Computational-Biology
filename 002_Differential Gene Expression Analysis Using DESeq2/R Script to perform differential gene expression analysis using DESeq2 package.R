# R Script to perform differential gene expression analysis using DESeq2
# GEO Accession ID: GSE46056
# Required file type: GSE46056_raw_counts_GRCh38.p13_NCBI.tsv.gz
# URL: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE46056

# 1) Install and load required R packages and libraries.
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("tidyverse")
BiocManager::install("ggplot2")
BiocManager::install("readr")
BiocManager::install("dplyr")
BiocManager::install("biomaRt")
BiocManager::install("ggrepel")

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(readr)
library(dplyr)
library(biomaRt)
library(ggrepel)

# 2) Read in counts data to R environment.
counts_data <- read.csv('counts_data.csv', row.names = 1)
head(counts_data)

# 3) Read in column data (Sample information).
col_data <- read.csv('col_data.csv', row.names = 1)

# 4) Ensure the row names in col_data matches to 
# column names in counts_data as well as in the same order.
all(colnames(counts_data) %in% rownames(col_data))
all(colnames(counts_data) == rownames(col_data))

# 5) Set up a DESeqDataSet object.
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = col_data,
                              design = ~ Knockdown)

# 6) Remove rows with low gene counts (less than 10).
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 7) Set the factor level to enable comparison between the conditions.
dds$Knockdown <- relevel(dds$Knockdown, ref = "Untreated")

# NOTE: If there are technical replicates in the dataset, make sure to 
# collapse them. But for this dataset, it is not required.

# 8) Run DEseq and get result.
dds <- DESeq(dds)
result <- results(dds)
summary(result)

# 9) Get significant result with an adjusted p-value threshold of 0.01
# instead of 0.1.
significant_result <- results(dds, alpha = 0.01)
summary(significant_result)

# Save significant results to CSV.
write.csv(as.data.frame(significant_result), "significant_results.csv")

# 10) Plot MA plots for result and significant result.
plotMA(significant_result, main="MA Plot")

# --------------------------------------------------------------------------------

# Load the data with ENTREZ IDs as row names
data <- read.csv("significant_results.csv", row.names = 1)

# Convert row names to a column for merging purposes, naming it as 'entrez_id'
data <- data %>%
  tibble::rownames_to_column(var = "entrez_id")

# Ensure that entrez_id is numeric for compatibility with biomaRt
data$entrez_id <- as.numeric(as.character(data$entrez_id))

# Use biomaRt to convert ENTREZ IDs to gene symbols
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(filters = "entrezgene_id",
               attributes = c("entrezgene_id", "hgnc_symbol"),
               values = unique(data$entrez_id),
               mart = mart)

# Merge gene symbols into the data frame based on entrez_id
data <- data %>%
  left_join(genes, by = c("entrez_id" = "entrezgene_id"))

# Remove rows with NA values in both hgnc_symbol and padj columns
data <- data %>% filter(!is.na(hgnc_symbol) & !is.na(padj))

# Add columns for plotting: -log10(p-value), group classification, and significance flag
plotdata <- data %>%
  mutate(log10_pvalue = -log10(pvalue),
         group = case_when(
           padj <= 0.01 & log2FoldChange > 1 ~ "Upregulated",
           padj <= 0.01 & log2FoldChange < -1 ~ "Downregulated",
           TRUE ~ "Not Significant"
         ))

# Filter significant results for labeling
plotdata_sig <- plotdata %>% filter(group != "Not Significant")

# Generate the volcano plot with centered title and axis labels
ggplot(plotdata, aes(x = log2FoldChange, y = log10_pvalue)) +
  geom_point(aes(color = group), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")) +
  ggrepel::geom_label_repel(data = plotdata_sig, aes(label = hgnc_symbol), size = 3, max.overlaps = 10) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value") +
  guides(color = guide_legend(title = "Differential Expression")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5)
  )