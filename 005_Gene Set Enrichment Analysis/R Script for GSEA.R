# R Script for Gene ontology (GO) enrichment and gene set enrichment analysis (GSEA)
# GEO Accession ID: GSE46056
# Required file type: GSE46056_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz
# URL: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE46056

# significant_results.CSV obtained from previous DEseq2 analysis will be used in this script

# Install necessary packages
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ggplot2")
BiocManager::install("enrichplot")
BiocManager::install("DOSE")
BiocManager::install("biomaRt")
BiocManager::install("tibble")
BiocManager::install("msigdbr")
BiocManager::install("dplyr")
BiocManager::install("patchwork")
BiocManager::install("DOSE")

# Load required libraries
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
library(DOSE)

# Load and prepare the data
data <- read.csv('significant_results.csv', row.names = 1)
head(data)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ncbi_gene_ids <- rownames(data)
gene_annotations <- getBM(attributes = c('entrezgene_id', 'hgnc_symbol'),
                          filters = 'entrezgene_id',
                          values = ncbi_gene_ids,
                          mart = ensembl)
data$entrezgene_id <- rownames(data)
data <- merge(data, gene_annotations, by = "entrezgene_id", all.x = TRUE)
data <- data[!is.na(data$hgnc_symbol) & data$hgnc_symbol != "", ]
data <- data[!grepl("^[0-9]+$", data$hgnc_symbol), ]
rownames(data) <- make.unique(data$hgnc_symbol)
data <- data[, !(colnames(data) %in% c("entrezgene_id", "hgnc_symbol"))]
data <- data[!is.na(data$padj), ]

# Add gene symbols as a column
data <- data %>%
  rownames_to_column(var = "gene")

# Identify top genes and significant genes
# Arrange genes by padj and log2FoldChange, display the top 30
top_genes <- data %>%
  arrange(padj, desc(log2FoldChange)) %>%
  head(n = 30)

# Select genes that are significantly differentially expressed (e.g., padj ≤ 0.01 and |log2FoldChange| ≥ 2).
significant_genes <- data %>%
  filter(padj <= 0.01, abs(log2FoldChange) >= 2) %>%
  pull(gene)  # pull() to get a vector

# Select all genes with baseMean > 0 as background genes
background_genes <- data %>%
  filter(baseMean != 0) %>%
  pull(gene)  # pull() to get a vector

#-------------------------------------------------------------------------------

# PART 1 

# GO term enrichment analysis using Biological Process (BP) ontology

# GO Enrichment Analysis (ego)
ego <- enrichGO(
  gene = significant_genes,
  universe = background_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",      # Use Gene Symbols as identifiers
  ont = "BP",              # 'BP' for Biological Process
  pAdjustMethod = "BH",    # Benjamini-Hochberg correction for multiple testing
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Bar plot of GO Terms
barplot(ego, showCategory = 20) + 
  ggtitle("Top 20 Enriched GO Terms - Bar Plot")
ggsave("Top 20 Enriched GO Terms - Bar Plot.png", width = 16, height = 12, bg = "white")

# Dot plot of GO Terms
dotplot(ego, showCategory = 20) + 
  ggtitle("Top 20 Enriched GO Terms - Dot Plot")
ggsave("Top 20 Enriched GO Terms - Dot Plot.png", width = 16, height = 12, bg = "white")

# Calculate pairwise similarity for enriched GO terms
# Helps identify and visualise clusters of related GO terms, 
# enhancing biological interpretation by revealing functional relationships.
ego_sim <- pairwise_termsim(ego)

# Enrichment map of GO Terms
emapplot(ego_sim, showCategory = 30) +  
  ggtitle("Enrichment Map of GO Terms") +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),  # Increase margins
    text = element_text(size = 10)  # Adjust text size
  )
ggsave("Enrichment Map of GO Terms.png", width = 16, height = 12, bg = "white")

# Gene-Concept Network for GO Terms
cnetplot(ego, showCategory = 10, foldChange = data$log2FoldChange) + 
  ggtitle("Gene-Concept Network for GO Terms") +
  theme(
    plot.margin = unit(c(1,1,1,1), "cm"),
    text = element_text(size = 10)
  )
ggsave("Gene-Concept Network for GO Terms.png", width = 16, height = 12, bg = "white")

# Create geneList with log2FoldChange values, naming each value with the corresponding gene symbol
geneList <- data$log2FoldChange
names(geneList) <- data$gene

# Heatmap with and without fold change of GO Terms
p1 <- heatplot(ego, showCategory = 20) + 
  ggtitle("Heatmap of GO Terms") +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        plot.margin = unit(c(1, 1, 1, 1), "cm"))  # Center title and add margin

p2 <- heatplot(ego, foldChange = geneList, showCategory = 20) + 
  ggtitle("Heatmap of GO Terms with Fold Changes") +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        plot.margin = unit(c(1, 1, 1, 1), "cm"))  # Center title and add margin
# Arrange the two heatmaps side-by-side
combined_plot <- cowplot::plot_grid(p1, p2, ncol = 1, labels = LETTERS[1:2])
ggsave("Combined_Heatmaps of GO Terms.png", plot = combined_plot, width = 16, height = 12, bg = "white")

# Tree plot of GO Terms
# Generate the tree plot and assign it to a variable
tree_plot <- treeplot(ego_sim) + 
  ggtitle("Tree Plot of GO Enriched Terms")
ggsave("Tree Plot of GO Enriched Terms.png", width = 16, height = 12, bg = "white")

#-------------------------------------------------------------------------------

# PART 2

# Gene Set Enrichment Analysis (GSEA)

# Convert significant genes and background genes to Entrez IDs
significant_genes_map <- bitr(data$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
background_genes_map <- bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Prepare ranked list for GSEA based on log2FoldChange and p-value
res_df <- data %>% 
  mutate(signed_rank_stats = sign(log2FoldChange) * -log10(pvalue)) %>%
  left_join(significant_genes_map, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%  # Filter out rows with NA in ENTREZID
  arrange(desc(signed_rank_stats))

# Create ranked gene list for GSEA
gene_list <- res_df$signed_rank_stats
names(gene_list) <- res_df$ENTREZID

# Load MSigDB hallmark gene sets for Homo sapiens in TERM2GENE format
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

# Run GSEA
em2 <- GSEA(gene_list, TERM2GENE = m_t2g)

# Select top up- and down-regulated pathways for plotting
top_up_pathways <- em2@result %>%
  filter(p.adjust < 0.05, NES > 0) %>%
  arrange(desc(NES)) %>%
  head(10) %>%
  pull(Description)

top_down_pathways <- em2@result %>%
  filter(p.adjust < 0.05, NES < 0) %>%
  arrange(NES) %>%
  head(10) %>%
  pull(Description)

selected_pathways <- c(top_up_pathways, top_down_pathways)

# Create GSEA plots for selected pathways
plots <- list()
for (pathway in selected_pathways) {
  plot <- gseaplot(em2, geneSetID = pathway, by = "runningScore", title = pathway)
  plots[[pathway]] <- plot + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
}

# Define function to create pages with 4 plots each to prevent overcrowding of plots
create_plot_pages <- function(plots, plots_per_page = 4) {
  num_pages <- ceiling(length(plots) / plots_per_page)
  
  for (page in 1:num_pages) {
    start_index <- (page - 1) * plots_per_page + 1
    end_index <- min(page * plots_per_page, length(plots))
    
    final_plot <- wrap_plots(plots[start_index:end_index], ncol = 2) + 
      plot_layout(guides = "collect") +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    
    # Save each page of plots
    ggsave(paste0("GSEA_plots_page_", page, ".png"), plot = final_plot, width = 16, height = 12, bg = "white")
  }
}
# Generate and save GSEA plot pages
create_plot_pages(plots, plots_per_page = 4)

# Dot Plot of GSEA Enriched Pathways
dotplot(em2, showCategory = 20) +
  ggtitle("Dot Plot of GSEA Enriched Pathways")
ggsave("Dot Plot of GSEA Enriched Pathways.png", width = 16, height = 12, bg = "white")

# Ridge Plot of GSEA Enriched Pathways
ridgeplot(em2, showCategory = 20) + 
  ggtitle("Ridge Plot of GSEA Enriched Pathways") +
  theme(axis.text.y = element_text(size = 8))
ggsave("Ridge Plot of GSEA Enriched Pathways.png", width = 16, height = 12, bg = "white")

# ------------------------------------------------------------------------------

# PART 3

# Disease Ontology Enrichment Analysis

# Run Disease Ontology Enrichment Analysis (without OrgDb argument)
# Run DisGeNET Disease Enrichment Analysis with Entrez IDs
# Extract Entrez IDs from the significant_genes_map data frame
significant_genes_entrez <- significant_genes_map$ENTREZID

# Run DisGeNET Disease Enrichment Analysis with Entrez IDs
dgn_enrich <- enrichDGN(gene = significant_genes_entrez,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

# Calculate pairwise similarity for enriched disease terms
dgn_enrich_sim <- pairwise_termsim(dgn_enrich)

# Bar Plot for Disease Terms 
barplot(dgn_enrich, showCategory = 20) +
  ggtitle("Top Enriched Disease Terms - Bar Plot")
ggsave("Top_Enriched_Disease_Terms_Bar Plot.png", width = 16, height = 12, bg = "white")

# Dot Plot for Disease Terms
dotplot(dgn_enrich, showCategory = 20) +
  ggtitle("Top Enriched Disease Terms - Dot Plot")
ggsave("Top_Enriched_Disease_Terms_Dot Plot.png", width = 16, height = 12, bg = "white")

# Enrichment map plot (requires pairwise_termsim)
emapplot(dgn_enrich_sim, showCategory = 20) + 
  ggtitle("Enrichment Map of Disease Terms")
ggsave("Enrichment Map_Disease_Terms.png", width = 16, height = 12, bg = "white")


