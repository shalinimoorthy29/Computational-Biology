
# Gene Ontology (GO) Enrichment and Gene Set Enrichment Analysis (GSEA) Template

## Overview
This template provides a step-by-step guide to perform Gene Ontology (GO) enrichment and Gene Set Enrichment Analysis (GSEA) using **clusterProfiler** and related R packages. Replace placeholders (indicated with `< >`) with your dataset-specific details.

---

## Prerequisites

1. **Set Your Working Directory**  
   Update your working directory to the folder containing your files:
   ```r
   setwd("<path_to_your_working_directory>")
   ```

2. **Prepare Your Data**:
   - **Significant Results CSV file** (`<deseq2_results_file.csv>`): Results from a previous DESeq2 analysis.
   - Ensure that the file includes `gene`, `padj`, and `log2FoldChange` columns.

---

## Step-by-Step Guide

### 1. Install and Load Required Packages
Install and load the necessary R packages:
```r
install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ggplot2", "enrichplot", "DOSE", "biomaRt", 
                       "tibble", "msigdbr", "dplyr", "patchwork"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(DOSE)
library(biomaRt)
library(tibble)
library(msigdbr)
library(dplyr)
library(patchwork)
```

---

### 2. Load and Prepare the Data
Load the significant results from the DESeq2 analysis into R. Replace `<deseq2_results_file.csv>` with your file name:
```r
data <- read.csv('<deseq2_results_file.csv>', row.names = 1)
head(data)
```

Map gene symbols using **biomaRt**:
```r
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ncbi_gene_ids <- rownames(data)
gene_annotations <- getBM(attributes = c('entrezgene_id', 'hgnc_symbol'),
                          filters = 'entrezgene_id',
                          values = ncbi_gene_ids,
                          mart = ensembl)
data$entrezgene_id <- rownames(data)
data <- merge(data, gene_annotations, by = "entrezgene_id", all.x = TRUE)
```

Filter and clean the data:
```r
data <- data[!is.na(data$padj) & !is.na(data$hgnc_symbol) & data$hgnc_symbol != "", ]
rownames(data) <- make.unique(data$hgnc_symbol)
```

---

### 3. Identify Significant and Background Genes
Identify significant genes and background genes:
```r
# Select genes with padj ≤ 0.01 and |log2FoldChange| ≥ 2
significant_genes <- data %>%
  filter(padj <= 0.01, abs(log2FoldChange) >= 2) %>%
  pull(gene)

# Select all genes with baseMean > 0 as background genes
background_genes <- data %>%
  filter(baseMean != 0) %>%
  pull(gene)
```

---

## Part 1: GO Enrichment Analysis

### Perform GO Enrichment Analysis
```r
ego <- enrichGO(
  gene = significant_genes,
  universe = background_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = TRUE
)
```

### Visualise GO Terms
Generate bar and dot plots for GO terms:
```r
barplot(ego, showCategory = 20) + ggtitle("Top 20 Enriched GO Terms - Bar Plot")
dotplot(ego, showCategory = 20) + ggtitle("Top 20 Enriched GO Terms - Dot Plot")
```

Create an enrichment map and a gene-concept network:
```r
emapplot(pairwise_termsim(ego), showCategory = 30) + ggtitle("Enrichment Map of GO Terms")
cnetplot(ego, showCategory = 10, foldChange = data$log2FoldChange) + ggtitle("Gene-Concept Network for GO Terms")
```

### Create Heatmaps for GO Terms
```r
p1 <- heatplot(ego, showCategory = 20) + ggtitle("Heatmap of GO Terms")
p2 <- heatplot(ego, foldChange = data$log2FoldChange, showCategory = 20) + ggtitle("Heatmap of GO Terms with Fold Changes")
```

---

## Part 2: Gene Set Enrichment Analysis (GSEA)

### Prepare Data for GSEA
Convert gene symbols to Entrez IDs and create a ranked gene list:
```r
significant_genes_map <- bitr(data$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_list <- data$log2FoldChange
names(gene_list) <- data$gene
```

Load hallmark gene sets and run GSEA:
```r
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)
em2 <- GSEA(gene_list, TERM2GENE = m_t2g)
```

Visualise GSEA results with dot and ridge plots:
```r
dotplot(em2, showCategory = 20) + ggtitle("Dot Plot of GSEA Enriched Pathways")
ridgeplot(em2, showCategory = 20) + ggtitle("Ridge Plot of GSEA Enriched Pathways")
```

---

## Part 3: Disease Ontology Enrichment Analysis

### Run Disease Ontology Enrichment Analysis
Perform Disease Ontology enrichment analysis with significant genes:
```r
dgn_enrich <- enrichDGN(
  gene = significant_genes_map$ENTREZID,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
```

Visualise disease terms with bar and dot plots:
```r
barplot(dgn_enrich, showCategory = 20) + ggtitle("Top Enriched Disease Terms - Bar Plot")
dotplot(dgn_enrich, showCategory = 20) + ggtitle("Top Enriched Disease Terms - Dot Plot")
emapplot(pairwise_termsim(dgn_enrich), showCategory = 20) + ggtitle("Enrichment Map of Disease Terms")
```

---

## Notes
- Replace all placeholders (`< >`) with values specific to your dataset.
- Ensure your DESeq2 results file contains the required columns (`gene`, `padj`, `log2FoldChange`).
