# Heatmap of Z-Score Normalized Gene Expression Data

**Dataset Used**:  
- Dataset: GSE46056_norm_counts_FPKM_GRCh38.p13_NCBI.csv (Gene Symbols)  
- This dataset contains normalized gene expression counts (FPKM) across multiple samples. The gene IDs were converted to Gene Symbols using the DAVID tool.

**Objectives**:  
The objective of this analysis was to create a heatmap of gene expression data to visualize patterns of expression across multiple samples. Z-score normalization was applied to the data to focus on relative changes in expression, making it easier to compare the variation across genes and samples.

**Stage of Analysis**:  
1. **Data Preprocessing**: The normalized counts were prepared by converting gene IDs to Gene Symbols using DAVID and saving the data as a CSV file. The processed file was loaded for further analysis.

2. **Log Transformation**: A log transformation was applied to the FPKM data to compress the range of values and make the differences in expression more interpretable. This is crucial when handling a dataset with a wide dynamic range in expression levels.

3. **Z-Score Normalization**: Z-score normalization was applied to the log-transformed data to standardize the expression values for each gene across samples. This ensures that the heatmap captures relative expression differences rather than absolute values.

4. **Heatmap Generation**: A heatmap was generated to visualize gene expression patterns. The data was clustered both by genes (rows) and samples (columns) to highlight relationships in the data.

**Rationale**:  
The heatmap is a powerful visualization tool that enables the comparison of expression patterns across genes and samples. By using FPKM data for the heatmap, we account for gene length and sequencing depth. Z-score normalization ensures that the focus is on relative changes in expression, which are often more biologically meaningful. Clustering both genes and samples helps identify co-expressed genes and sample groups with similar expression profiles.

**Conclusions**:  
The heatmap provided a clear visual representation of the gene expression patterns across the dataset. Clustering of genes and samples revealed distinct expression groups and potential biological insights. The use of z-score normalization allowed for a focus on relative differences in expression across the dataset, highlighting key patterns that may be related to experimental conditions or biological differences between samples.