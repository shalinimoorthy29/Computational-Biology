
# Heatmap Visualisation Template Using ComplexHeatmap

## Overview
This template provides a step-by-step guide to visualise gene expression heatmaps using **ComplexHeatmap** in R. Replace placeholders (indicated with `< >`) with your dataset-specific details.

---

## Prerequisites

1. **Set Your Working Directory**  
   Update your working directory to the folder containing your files:
   ```r
   setwd("<path_to_your_working_directory>")
   ```

2. **Prepare Your Data**:
   - **Normalised FPKM Counts File** (`<normalised_fpkm_file.tsv>`): Rows as gene IDs, columns as sample IDs.
   - **Sample Metadata CSV file** (`<sample_metadata_file.csv>`): Rows as sample IDs, columns for condition/treatment groups.

---

## Step-by-Step Guide

### 1. Install and Load Required Packages
Install and load the necessary R packages:
```r
install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "circlize", "biomaRt"))

library(ComplexHeatmap)
library(circlize)
library(biomaRt)
```

---

### 2. Load FPKM Data
Load the normalised FPKM data into R. Replace `<normalised_fpkm_file.tsv>` with your file name:
```r
fpkm_data <- read.table("<normalised_fpkm_file.tsv>", sep = "	", header = TRUE, row.names = 1)
```

---

### 3. Log-Transform the FPKM Data
Perform a log transformation (adding 1 to avoid `log(0)`):
```r
log_fpkm_data <- log2(fpkm_data + 1)
```

---

### 4. Map Gene IDs to HGNC Symbols
Use **biomaRt** to map Entrez Gene IDs to HGNC gene symbols:
```r
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ncbi_gene_ids <- rownames(log_fpkm_data)

gene_annotations <- getBM(attributes = c('entrezgene_id', 'hgnc_symbol'),
                          filters = 'entrezgene_id',
                          values = ncbi_gene_ids,
                          mart = ensembl)

log_fpkm_data$entrezgene_id <- rownames(log_fpkm_data)
log_fpkm_data <- merge(log_fpkm_data, gene_annotations, by = "entrezgene_id", all.x = TRUE)
log_fpkm_data <- log_fpkm_data[!is.na(log_fpkm_data$hgnc_symbol) & log_fpkm_data$hgnc_symbol != "", ]
log_fpkm_data <- log_fpkm_data[!grepl("^[0-9]+$", log_fpkm_data$hgnc_symbol), ]
rownames(log_fpkm_data) <- make.unique(log_fpkm_data$hgnc_symbol)
log_fpkm_data <- log_fpkm_data[, !(colnames(log_fpkm_data) %in% c("entrezgene_id", "hgnc_symbol"))]
```

---

### 5. Select Top 10% Most Variable Genes
Identify and subset the top 10% most variable genes based on standard deviation:
```r
gene_sd <- apply(log_fpkm_data, 1, sd)
top_genes <- order(gene_sd, decreasing = TRUE)
top_10_percent <- top_genes[1:(length(gene_sd) * 0.1)]

subset_data <- log_fpkm_data[top_10_percent, ]
```

---

### 6. Z-Score Normalisation
Apply z-score normalisation to the subset data:
```r
scaled_data <- t(scale(t(subset_data)))
scaled_data <- scaled_data[complete.cases(scaled_data), ]
```

---

### 7. Load Sample Metadata
Load the sample metadata file and extract conditions for annotation. Replace `<sample_metadata_file.csv>` with your file name:
```r
col_data <- read.csv("<sample_metadata_file.csv>")
conditions <- col_data$Condition
```

---

### 8. Define Condition Annotations
Create condition annotations for the heatmap using the `Condition` column in the metadata:
```r
column_ha <- HeatmapAnnotation(
  condition = conditions,
  col = list(condition = c("PRMT4-Knockdown" = "red", "Control" = "blue")),
  annotation_legend_param = list(condition = list(title = "Condition"))
)
```

---

### 9. Plot Heatmap
Generate the heatmap with annotations:
```r
Heatmap(scaled_data,
        name = "Expression",
        col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
        top_annotation = column_ha,
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_side = "right",
        column_names_side = "bottom",
        heatmap_legend_param = list(title = "Scaled Normalised Expression"))
```

---

## Notes
- Replace all placeholders (`< >`) with values specific to your dataset.
- Ensure your FPKM counts file and sample metadata files are formatted correctly.
