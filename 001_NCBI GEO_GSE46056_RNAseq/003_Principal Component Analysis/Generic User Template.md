
# PCA Analysis Template Using DESeq2

## Overview
This template provides a step-by-step guide to perform Principal Component Analysis (PCA) on RNA-Seq data using the **DESeq2** package in R. Replace placeholders (indicated with `< >`) with your dataset-specific details.

---

## Prerequisites

1. **Set Your Working Directory**  
   Update your working directory to the folder containing your files:
   ```r
   setwd("<path_to_your_working_directory>")
   ```

2. **Prepare Your Data**:
   - **Counts Matrix CSV file** (`<counts_matrix_file.csv>`): Rows as gene IDs, columns as sample IDs.
   - **Sample Metadata CSV file** (`<sample_metadata_file.csv>`): Rows as sample IDs, columns for condition/treatment groups.

---

## Step-by-Step Guide

### 1. Install and Load Required Packages
Install and load the necessary R packages:
```r
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("ggplot2")

library(DESeq2)
library(ggplot2)
```

---

### 2. Load Data
Load the counts matrix and sample metadata files into R. Replace `<counts_matrix_file.csv>` and `<sample_metadata_file.csv>` with your file names:
```r
counts_data <- read.csv('<counts_matrix_file.csv>', row.names = 1)
col_data <- read.csv('<sample_metadata_file.csv>', row.names = 1)
```

---

### 3. Verify Data Consistency
Ensure that column names in `counts_data` match row names in `col_data`:
```r
# Check if column names in counts_data match row names in col_data
all(colnames(counts_data) == rownames(col_data))  # This should return TRUE
```

---

### 4. Create DESeqDataSet Object
Replace `<condition_column>` with the name of the column in your metadata describing experimental conditions (e.g., "Treatment"):
```r
dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                              colData = col_data, 
                              design = ~ <condition_column>)
```

---

### 5. Filter Low-Count Genes
Remove rows with counts below the threshold (e.g., `10`):
```r
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```

---

### 6. Perform Variance Stabilising Transformation
Transform the dataset for PCA:
```r
vsd <- vst(dds, blind = FALSE)
```

---

### 7. Generate PCA Data
Generate the PCA data and compute the percentage of variance for each principal component:
```r
pca_data <- plotPCA(vsd, intgroup = "<condition_column>", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
```

---

### 8. Add Sample Names or IDs
Assign sample names or IDs from row names to the PCA data:
```r
pca_data$Sample <- rownames(pca_data)
```

---

### 9. Create PCA Plot
Create a PCA plot with sample labels:
```r
ggplot(pca_data, aes(PC1, PC2, color = <condition_column>)) +  # Replace <condition_column> with your condition column
  geom_point(size = 3) +
  geom_text(aes(label = Sample), vjust = -0.5, hjust = 0.5, show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_minimal() +
  labs(title = "PCA Plot with Sample Labels") +
  coord_cartesian(xlim = range(pca_data$PC1) + c(-2, 2),
                  ylim = range(pca_data$PC2) + c(-2, 2)) +
  theme(plot.title = element_text(hjust = 0.5))
```

---

## Notes
- Replace all placeholders (`< >`) with values specific to your dataset.
- Ensure your counts matrix and sample metadata files are formatted correctly.
