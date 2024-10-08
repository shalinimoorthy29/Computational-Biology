# Differential Gene Expression Analysis Using DESeq2

## Overview
This project demonstrates the analysis of RNA-Seq data to identify differentially expressed genes using the DESeq2 package in R. The dataset used for this analysis is from the GEO Accession ID: **GSE46056**.

## Prerequisites
- Set your working directory to the folder containing the necessary data files.
- Download the count data file `GSE46056_raw_counts_GRCh38.p13_NCBI.tsv.gz` from the following URL:  
  [GSE46056 Dataset](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE46056)
  
- Prepare two CSV files:
  1. `counts_data.csv`: Contains the raw counts matrix with gene IDs as rows and sample IDs (GSM IDs) as columns.
  2. `col_data.csv`: Contains sample information (GSM IDs, sample names, and experimental conditions). The GSM IDs must match the column names of the `counts_data.csv` file and be in the same order.

## Required R Packages
This analysis relies on the following R packages:
- **DESeq2**: For differential expression analysis of RNA-Seq data.
- **tidyverse**: A collection of R packages used for data manipulation and visualization.
- **ggplot2**: For creating visualizations such as MA and volcano plots.
- **readr**: For reading CSV files into R.

## Analysis Stages

### 1. Loading Data
The counts matrix and sample metadata are loaded into the R environment. The counts matrix contains raw gene expression data (with gene IDs as rows and sample IDs as columns), while the metadata contains information about the experimental conditions for each sample.

### 2. Data Consistency Checks
Before the analysis, the script checks that the sample names in the counts matrix match the corresponding entries in the sample metadata. This ensures that the samples are properly aligned for the analysis.

### 3. Creating a DESeqDataSet
A `DESeqDataSet` object is created using the raw counts and sample metadata. Genes with low expression levels (less than 10 counts) are filtered out to improve the reliability of the differential expression analysis. Additionally, the experimental condition being studied (e.g., "Knockdown" vs. "Untreated") is set as the reference for comparisons.

### 4. Differential Expression Analysis
The DESeq2 package is used to perform the differential gene expression analysis. This identifies genes whose expression levels significantly differ between the experimental conditions. The results include log2 fold changes, p-values, and adjusted p-values (to account for multiple testing).

### 5. Filtering Significant Genes
Significant differentially expressed genes are identified based on an adjusted p-value (padj) threshold of 0.01. The significant results are saved as a CSV file for further exploration or reporting.

### 6. Visualisation

#### MA Plot
An MA plot is generated to visualise the relationship between the mean expression of genes and the log2 fold change in expression. This plot helps highlight genes that have large changes in expression levels compared to their average expression.

#### Volcano Plot
A volcano plot is created to show the relationship between log2 fold changes and statistical significance (-log10 p-values) for each gene. Genes that are both statistically significant and have large fold changes are highlighted in red, making it easy to identify the most biologically relevant genes.

### 7. Summary of Results
The analysis successfully identifies differentially expressed genes between the experimental conditions. Genes with significant changes in expression (adjusted p-value < 0.01 and absolute log2 fold change > 1) are highlighted in the volcano plot, providing a clear visual representation of important gene expression changes.

## Conclusions
This project demonstrates a complete workflow for performing RNA-Seq differential gene expression analysis using DESeq2. It includes steps for data loading, quality control, differential expression analysis, and visualisation of results. The analysis reveals significant changes in gene expression between the experimental conditions, which are easily identifiable through the visualizations (MA and volcano plots).

## Output
- **significant_results.csv**: A CSV file containing the differentially expressed genes identified through the analysis.
- **MA Plot**: A visualisation showing the mean expression versus log2 fold change for all genes.
- **Volcano Plot**: A visualisation showing the significance and fold change of genes, highlighting those that are differentially expressed.

## License
This project is licensed under the MIT License - see the LICENSE file for details.
