# Heatmap of Z-Score Normalised Gene Expression Data (with Gene Symbol Mapping)

**Dataset Used**:  
- Dataset: `GSE46056_norm_counts_FPKM_GRCh38.p13_NCBI.tsv`  
- This dataset contains normalised gene expression counts (FPKM) across multiple samples, originally identified using NCBI Entrez Gene IDs.  

**Objectives**:  
The objective of this analysis was to create a heatmap of gene expression data to visualise patterns of expression across multiple samples. Z-score normalisation was applied to the data to focus on relative changes in expression, making it easier to compare the variation across genes and samples. The gene IDs were mapped to human-readable gene symbols using the **biomaRt** package for clearer biological interpretation.

**Stage of Analysis**:  
1. **Data Preprocessing**: The normalised counts were loaded for analysis. Entrez Gene IDs were mapped to HGNC gene symbols using the Ensembl database via the **biomaRt** package. Rows with missing or unmapped gene symbols were removed from the dataset.

2. **Log Transformation**: A log transformation was applied to the FPKM data (adding 1 to avoid log(0)) to compress the range of values and make the differences in expression more interpretable. This is crucial when handling a dataset with a wide dynamic range in expression levels.

3. **Z-Score Normalisation**: Z-score normalisation was applied to the log-transformed data to standardise the expression values for each gene across samples. This ensures that the heatmap captures relative expression differences rather than absolute values.

4. **Handling Missing/Unmapped Data**: After mapping Entrez Gene IDs to gene symbols, rows without valid gene symbols were removed. This ensures that only genes with corresponding human-readable gene symbols are displayed in the final heatmap. Duplicate gene symbols were handled using unique identifiers to avoid naming conflicts.

5. **Heatmap Generation**: A heatmap was generated to visualise gene expression patterns. The data was clustered both by genes (rows) and samples (columns) to highlight relationships in the data. The colour scale (blue-white-red) was used to represent the range of z-scores, from low to high expression.

**Rationale**:  
The heatmap is a powerful visualisation tool that enables the comparison of expression patterns across genes and samples. Mapping Entrez Gene IDs to gene symbols makes the data more interpretable and biologically relevant. Z-score normalisation ensures that the focus is on relative changes in expression, which are often more biologically meaningful. Clustering both genes and samples helps identify co-expressed genes and sample groups with similar expression profiles.

**Conclusions**:  
The heatmap provided a clear visual representation of the gene expression patterns across the dataset. The mapping of gene IDs to gene symbols improved the biological relevance of the data, and the removal of unmapped or invalid rows ensured data accuracy. Clustering of genes and samples revealed distinct expression groups and potential biological insights. The use of z-score normalisation allowed for a focus on relative differences in expression across the dataset, highlighting key patterns that may be biologically meaningful.
