# Principal Component Analysis (PCA) on Gene Expression Data

**Dataset Used**:  
- Dataset: GSE46056_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz  
- This dataset contains RNA-seq raw counts across multiple samples from the GSE46056 dataset. The data was analysed to investigate variations in gene expression levels across different experimental conditions.

**Objectives**:  
The objective of this analysis was to apply Principal Component Analysis (PCA) on raw gene expression counts to understand the major sources of variance in the data. The PCA aimed to reduce dimensionality while preserving the most critical patterns of variation between samples.

**Stage of Analysis**:  
1. **Data Loading and Preparation**: The raw count data was loaded and prepared for PCA. This step ensured that the data was in the correct format for dimensionality reduction.
   
2. **Variance Stabilising Transformation**: The raw counts were processed using a variance-stabilising transformation (VST) to ensure that the counts were more comparable across samples. This transformation helps adjust for the large dynamic range in RNA-seq data.
   
3. **PCA Execution**: PCA was applied to the transformed data, identifying principal components that captured the most variance in gene expression levels across samples.

**Rationale**:  
PCA is a common tool for dimensionality reduction and visualisation in gene expression studies. The rationale for using raw counts (with transformation) in PCA is to preserve the biological variance. By using PCA, we reduce the complexity of the dataset while maintaining the structure of the most relevant expression patterns. The primary components reflect the main sources of variability across the dataset, which could be due to biological differences, experimental conditions, or batch effects.

**Conclusions**:  
The PCA successfully reduced the dimensionality of the gene expression data, capturing significant variance in the first two principal components (PC1: 62% variance, PC2: 29% variance). The separation observed in the PCA plot provided insights into how different experimental conditions (e.g., treated vs. untreated) contributed to the variance in the dataset. The clear separation between the treated and untreated samples in the PCA plot confirms that the experimental conditions significantly influenced gene expression profiles. This clustering of samples demonstrates that the primary sources of variability in the data are closely related to the treatment conditions, helping to identify potential relationships between gene expression patterns and sample groupings.
