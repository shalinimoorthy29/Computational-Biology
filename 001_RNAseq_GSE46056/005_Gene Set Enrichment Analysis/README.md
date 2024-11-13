# Gene Ontology (GO) and Gene Set Enrichment Analysis (GSEA) for RNAseq Data of PRMT4 Knockdown in CD34+ Cells

**Study Accession**: [GSE46056](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE46056)  
**Organism**: Homo sapiens  
**Summary**: This project explores the role of PRMT4, a type I arginine methyltransferase, in hematopoiesis and its potential oncogenic function in acute myelogenous leukemia (AML). The overexpression of PRMT4 in AML patient samples and its repressive action on myeloid differentiation in human stem/progenitor cells (HSPCs) suggest a key regulatory role. PRMT4 knockdown (KD) induces myeloid differentiation, hinting at its therapeutic targeting potential in AML.

## Project Overview

This repository contains an enrichment analysis pipeline applied to RNAseq data from the GSE46056 dataset. I used gene ontology (GO) enrichment, gene set enrichment analysis (GSEA), and disease ontology enrichment to uncover pathways and biological processes regulated by PRMT4 in CD34+ hematopoietic cells. The analyses performed aim to clarify the impact of PRMT4 knockdown on gene expression and potential downstream biological effects, particularly those relevant to myeloid differentiation and AML.

## Pipeline Description

The script performs the following steps:

1. **Data Loading and Preparation**:  
   I load the gene expression results from DESeq2 analysis and map gene symbols to Entrez IDs for enrichment analysis.

2. **GO Term Enrichment Analysis**:  
   GO enrichment analysis was conducted on differentially expressed genes (DEGs) to explore biological processes affected by PRMT4 knockdown. This included visualisations such as bar plots, dot plots, enrichment maps, and heatmaps.

3. **GSEA**:  
   Using GSEA, I assessed the enrichment of hallmark gene sets (from MSigDB) to identify pathways influenced by PRMT4 knockdown, focusing on pathways with high normalised enrichment scores (NES) for both upregulated and downregulated gene sets.

4. **Disease Ontology Enrichment Analysis**:  
   I performed disease ontology enrichment to relate DEGs to known disease associations, which can potentially shed light on diseases or conditions linked to PRMT4 dysregulation, particularly hematological malignancies.

## Interpretation of Results in the Context of GSE46056

### GO Term Enrichment Analysis
The GO enrichment results reveal biological processes associated with hematopoiesis, cell differentiation, and immune response pathways. Key findings from GO term enrichment include:
- **Myeloid Differentiation**: GO terms related to immune response and hematopoietic processes were enriched, supporting the study's observation that PRMT4 knockdown induces differentiation in myeloid progenitor cells.
- **Transcriptional Repression**: PRMT4’s role in repressing gene expression was indicated by GO terms associated with transcriptional regulation, which aligns with the study’s findings that PRMT4 represses miR-223 expression, thus blocking myeloid differentiation.

### Gene Set Enrichment Analysis (GSEA)
The hallmark pathways identified through GSEA provide insight into broader cellular processes affected by PRMT4 knockdown:
- **Myc Targets and G2M Checkpoint Pathways**: Enrichment in MYC target and cell cycle-related pathways supports PRMT4's role in maintaining a proliferative, undifferentiated state, as observed in AML.
- **Immune and Inflammatory Responses**: Pathways related to immune response, TNF-alpha signalling, and hypoxia were also enriched, suggesting that PRMT4 knockdown not only promotes differentiation but may also impact immune signalling pathways in hematopoietic cells.
- **Potential Therapeutic Implications**: Targeting PRMT4 may disrupt oncogenic pathways in AML, providing a basis for therapeutic intervention, as PRMT4 depletion decreases leukemia cell proliferation.

### Disease Ontology Enrichment
The disease enrichment analysis highlights connections between PRMT4-regulated genes and various hematological malignancies, particularly:
- **Hematologic Malignancies**: Enriched disease terms related to leukemia and other blood-related disorders align with the study’s finding that PRMT4 overexpression is common in AML and may contribute to leukemogenesis.
- **Immune and Autoimmune Diseases**: Diseases associated with immune dysregulation were also identified, which could provide additional context for understanding PRMT4’s role in the immune landscape of hematopoietic cells.

## Repository Contents

- **Script**: `GSEA_GO_DO_Enrichment.R` - R script for performing GO enrichment, GSEA, and disease enrichment analysis.
- **Figures**: Output plots illustrating enrichment results, including:
  - GO term bar and dot plots
  - GSEA enrichment plots for top pathways
  - Disease enrichment maps and network plots

## Summary

This project provides a comprehensive analysis of the transcriptional impact of PRMT4 knockdown on human CD34+ cells. The results underscore PRMT4's role in repressing myeloid differentiation and suggest that targeting PRMT4 could potentially disrupt oncogenic pathways in AML. This study contributes to the understanding of PRMT4 as a therapeutic target, with potential applications in treating hematological malignancies.

## How to Use

1. **Requirements**: Install required R packages (see script for package list).
2. **Data**: The script uses a CSV file (`significant_results.csv`) containing DESeq2 analysis results from PRMT4 KD vs. control.
3. **Running the Script**: Run the R script to reproduce the enrichment analyses and generate visual outputs.

## References

- Original Study: RNAseq of PRMT4KD in human cord blood derived CD34+ cells (GSE46056)
- ClusterProfiler: A universal enrichment tool for interpreting omics data
- MSigDB Hallmark Gene Sets: Broad Institute

## License

This project is licensed under the MIT License.

---

By exploring the biological processes, pathways, and diseases associated with PRMT4-regulated genes, this analysis provides further insight into PRMT4’s role in hematopoiesis and leukemia. The results suggest that PRMT4 is not only a regulator of myeloid differentiation but also a promising target for AML therapy.
