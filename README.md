# Welcome to the Computational Biology Repository

This repository will showcase various sequencing data analysis techniques using R. Analysis will include, but is not limited to, Bulk RNA-Seq, ChIP-Seq, and Single-Cell RNA-Seq.

**Currently, I have developed R scripts for the following:**

- **Differential Gene Expression Analysis**: Identify genes with significant expression changes across conditions using DESeq2.
- **Principal Component Analysis (PCA)**: Visualise sample clustering to assess data quality and batch effects.
- **Heatmap Visualisations**: Generate informative heatmaps to display gene expression patterns, complete with gene symbol mapping.
- **Gene Ontology (GO) Enrichment and Gene Set Enrichment Analysis (GSEA)**: Discover enriched biological pathways and processes based on expression data.

## Packages:

- **BiocManager**: I’ll use BiocManager to manage Bioconductor package installations and dependencies, ensuring that required packages and libraries for analyses are up-to-date and compatible.

- **clusterProfiler**: This will facilitate functional enrichment analysis, enabling me to explore Gene Ontology (GO) terms, KEGG pathways, and other biological functions in clusters of genes or cells. It’s especially useful for interpreting single-cell clustering results in terms of biological relevance.

- **Seurat**: A comprehensive toolkit for single-cell RNA-Seq data analysis. I’ll use Seurat for data filtering, normalisation, dimensionality reduction, clustering, and differential expression. Seurat also supports visualisation techniques like t-SNE and UMAP, making it ideal for end-to-end single-cell workflows.

I will be continuously updating this repository, so please stay tuned for ongoing updates!
