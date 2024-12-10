
# Differential Gene Expression Analysis Template Using DESeq2

## Overview
This template provides a step-by-step guide to perform differential gene expression analysis using the **DESeq2** package in R. Replace placeholders (indicated with `< >`) with your dataset-specific details.

---

## Prerequisites

1. **Set Your Working Directory**  
   Update your working directory to the folder containing your files:
   ```r
   setwd("<path_to_your_working_directory>")
   ```

2. **Prepare Your Data**:
   - **Counts Matrix CSV** (`<counts_data.csv>`): Rows as gene IDs, columns as sample IDs.
   - **Sample Metadata CSV** (`<col_data.csv>`): Rows as sample IDs, columns for condition/treatment groups.

---

## Step-by-Step Guide

### 1. Install and Load Required Packages
Install and load the necessary R packages:
```r
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
```

---

### 2. Read in Data
Load the counts and metadata files into R. Replace `<counts_data.csv>` and `<col_data.csv>` with your file names.
```r
counts_data <- read.csv('<counts_data.csv>', row.names = 1)
col_data <- read.csv('<col_data.csv>', row.names = 1)
```

---

### 3. Ensure Proper Alignment of Data
Ensure that column names in `counts_data` match row names in `col_data`:
```r
# Check if column names in counts_data match row names in col_data
all(colnames(counts_data) %in% rownames(col_data))
all(colnames(counts_data) == rownames(col_data))
```

---

### 4. Create DESeqDataSet Object
Replace `<condition_column>` with the name of the column in your metadata describing experimental conditions (e.g., "Treatment").
```r
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = col_data,
                              design = ~ <condition_column>)
```

---

### 5. Filter Low-Count Genes
Remove rows with counts below the threshold (e.g., `10`):
```r
keep <- rowSums(counts(dds)) >= <threshold_value>
dds <- dds[keep,]
```

---

### 6. Set Reference Level
Specify the control or untreated group as the reference:
```r
dds$<condition_column> <- relevel(dds$<condition_column>, ref = "<reference_condition>")
```

---

### 7. Run Differential Expression Analysis
Run DESeq2 to compute differential expression:
```r
dds <- DESeq(dds)
result <- results(dds)
summary(result)
```

---

### 8. Filter Significant Results
Filter results based on an adjusted p-value threshold (e.g., `0.01`):
```r
significant_result <- results(dds, alpha = <p_value_threshold>)
summary(significant_result)
```

---

### 9. Save Results to a CSV File
Save significant results to a CSV file:
```r
write.csv(as.data.frame(significant_result), "<output_filename.csv>")
```

---

## Visualisation

### 1. MA Plot
Create an MA plot:
```r
plotMA(significant_result, main = "MA Plot")
```

---

### 2. Volcano Plot
Create a volcano plot (replace `<entrez_id>` and `<hgnc_symbol>` with appropriate columns):
```r
# Prepare data for plotting
data <- read.csv("<output_filename.csv>", row.names = 1)
data <- data %>% tibble::rownames_to_column(var = "<entrez_id>")
data$<entrez_id> <- as.numeric(as.character(data$<entrez_id>))

# Add gene symbols using biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(filters = "entrezgene_id",
               attributes = c("entrezgene_id", "hgnc_symbol"),
               values = unique(data$<entrez_id>),
               mart = mart)

data <- data %>% left_join(genes, by = c("<entrez_id>" = "entrezgene_id"))

# Add significance flags
plotdata <- data %>% mutate(
  log10_pvalue = -log10(pvalue),
  group = case_when(
    padj <= <padj_threshold> & log2FoldChange > <log2fc_threshold> ~ "Upregulated",
    padj <= <padj_threshold> & log2FoldChange < -<log2fc_threshold> ~ "Downregulated",
    TRUE ~ "Not Significant"
  )
)

# Plot volcano
ggplot(plotdata, aes(x = log2FoldChange, y = log10_pvalue)) +
  geom_point(aes(color = group), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")) +
  ggrepel::geom_label_repel(data = filter(plotdata, group != "Not Significant"),
                            aes(label = hgnc_symbol), size = 3, max.overlaps = 10) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value")
```

---

## Notes
- Replace all placeholders (`< >`) with values specific to your dataset.
- Thresholds such as `p_value_threshold` and `log2fc_threshold` can be adjusted based on your analysis goals.
