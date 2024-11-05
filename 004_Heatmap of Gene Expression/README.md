# Heatmap Analysis of Gene Expression Data

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
Conclusion
The heatmap displays the expression patterns of the top 10% most variable genes across different experimental conditions, specifically comparing the PRMT4-Knockdown and Control samples.

Distinct Clustering: The heatmap indicates a clear separation between the two conditions, with distinct clustering observed. This suggests that the PRMT4 knockdown significantly alters the expression profiles of the genes analysed, highlighting the biological impact of this treatment.

Expression Patterns: The colours in the heatmap represent z-score normalised expression levels, where red indicates higher expression and blue indicates lower expression. The differential expression of genes shows that certain genes are upregulated in the PRMT4-Knockdown samples compared to the Control, while others are downregulated.

Biological Relevance: The clustering and differential expression patterns observed in the heatmap may point to underlying biological mechanisms affected by PRMT4 knockdown. This could involve pathways related to gene regulation, cell signalling, or other processes that warrant further investigation.

Variability and Information: The use of the top 10% most variable genes allows for a focused analysis on those genes that exhibit the greatest variability, which are often the most informative in understanding the response to the knockdown treatment.

Overall, the heatmap serves as a valuable visual tool in identifying how gene expression is influenced by experimental conditions, supporting the hypothesis that PRMT4 plays a significant role in regulating gene expression patterns within the samples analysed. Further functional analyses could help elucidate the specific pathways involved and their relevance to the biological questions at hand.
