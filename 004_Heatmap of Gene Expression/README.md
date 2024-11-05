# Heatmap Analysis of Gene Expression Dat

**Dataset Used**:  
- **Dataset**: GSE46056_norm_counts_FPKM_GRCh38.p13_NCBI.tsv  
- This dataset contains RNA-seq FPKM (Fragments Per Kilobase of transcript per Million mapped reads) counts across multiple samples from the GSE46056 dataset. The data was analysed to investigate variations in gene expression levels across different experimental conditions.

**Objectives**:  
The objective of this analysis was to create a heatmap of gene expression data to visually explore the top 10% most variable genes. This heatmap aimed to highlight the differences in expression patterns between experimental conditions.

**Stage of Analysis**:  
1. **Data Loading and Preparation**: The FPKM data was loaded and prepared for analysis. This step ensured that the data was in the correct format for visualisation. FPKM counts provide a normalised measure of gene expression that accounts for sequencing depth and gene length. This makes them suitable for EDA as they allow for meaningful comparisons across genes and samples. By using FPKM counts, researchers can explore variance in gene expression while minimising the impact of outliers, leading to more robust interpretations.
   
2. **Log-Transformation of FPKM Data**: The raw FPKM counts were log-transformed (log2 transformation with a pseudocount of 1) to stabilise variance and make the expression data more comparable across samples.

3. **Gene Annotation**: BiomaRt was used to map Entrez Gene IDs to HGNC gene symbols, enhancing the interpretability of the heatmap.

4. **Data Filtering**: Genes with missing or invalid symbols were removed, ensuring only valid gene annotations were used.

5. **Selection of Most Variable Genes**: The analysis focused on the top 10% most variable genes based on standard deviation, which are often the most informative for exploratory data analysis. Focusing on the top 10% most variable genes allows researchers to highlight the genes that exhibit the greatest variability across samples. These genes are often the most biologically relevant and can reveal significant patterns in response to different experimental conditions. This approach reduces dimensionality, making it easier to visualise and interpret the data, while capturing the main sources of variability.

6. **Z-Score Normalisation**: The selected subset of genes was normalised using z-scores to allow for easier comparison of expression levels across samples.

7. **Heatmap Generation**: A heatmap was created to visualise the expression patterns of the top 10% most variable genes, with condition annotations for treated and control samples. Although Entrez IDs were mapped to genes, the plot does not display gene labels due to space restrictions. This approach allows for a clearer view of expression patterns across more genes.

**Rationale for Heatmap Visualisation**:  
Heatmaps are powerful tools for visualising high-dimensional data such as RNA-seq expression profiles. They provide an intuitive way to identify patterns, clusters, and relationships among samples. In the context of RNA-seq analysis, heatmaps allow researchers to visually assess how gene expression changes across different conditions, facilitating insights into biological processes and treatment effects.

**Conclusions**:  
The heatmap successfully illustrated the expression patterns of the most variable genes, providing insights into the differences in gene expression between conditions. The clustering of samples indicated how treatment conditions influenced gene expression profiles, confirming that biological variations were captured effectively in the visualisation.
