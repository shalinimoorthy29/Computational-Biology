# Gene Set and Disease Ontology Enrichment Analysis of PRMT4 Knockdown in Human CD34+ Cells**

## Dataset
[GSE46056](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE46056)  

## Objectives
- To investigate the transcriptional impact of PRMT4 knockdown** in human CD34+ cells, focusing on gene expression changes associated with myeloid differentiation.
- To identify enriched biological processes and pathways** influenced by PRMT4 knockdown, using Gene Ontology (GO) enrichment analysis to reveal processes such as immune response, inflammatory regulation, and ion homeostasis.
- To explore disease associations of differentially expressed genes** through Disease Ontology (DO) and DisGeNET enrichment analyses, highlighting potential connections to hematological and developmental disorders.
- To assess PRMT4’s potential as a therapeutic target** in acute myeloid leukemia (AML) and other conditions driven by epigenetic dysregulation, based on pathway and disease enrichment insights.

## Project Overview
This study explores the role of PRMT4, a type I arginine methyltransferase, in hematopoiesis and its potential oncogenic function in acute myelogenous leukemia (AML). The overexpression of PRMT4 in AML patient samples and its repressive action on myeloid differentiation in human stem/progenitor cells (HSPCs) suggest a key regulatory role. PRMT4 knockdown (KD) induces myeloid differentiation, hinting at its therapeutic targeting potential in AML.

This project contains an enrichment analysis pipeline applied to RNAseq data from the GSE46056 dataset. I used gene ontology (GO) enrichment, gene set enrichment analysis (GSEA), and disease ontology enrichment to uncover pathways and biological processes regulated by PRMT4 in CD34+ hematopoietic cells. The analyses performed aim to clarify the impact of PRMT4 knockdown on gene expression and potential downstream biological effects, particularly those relevant to myeloid differentiation and AML.

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
 
## Results

### GO Term Enrichment Analysis

**Top 20 Enriched GO Terms - Bar and Dot Plots**
- The bar and dot plots collectively reveal key biological processes enriched following PRMT4 knockdown, including immune response (e.g., chemotaxis, regulation of inflammatory response) and cellular signalling pathways (e.g., ion transport and homeostasis). These processes align with the observed myeloid differentiation in CD34+ cells.
- While the bar plot shows the absolute **counts** of genes per term, highlighting pathways with the most genes involved, the dot plot’s **GeneRatio** metric emphasises terms with a high proportion of relevant genes relative to the total. Together, these metrics provide both absolute and proportional insights into the enriched biological processes.

**Enrichment Map of GO Terms**
- The enrichment map clusters related GO terms based on their biological similarity, showing distinct groupings around immune response, ion homeostasis, and cell signalling pathways. This visual grouping highlights how PRMT4 knockdown influences interconnected processes in myeloid differentiation.
- Notable clusters, such as those involving ion transport and inflammatory response, underscore PRMT4’s multifaceted role in regulating both immune function and cellular homeostasis, which are crucial for hematopoietic cell behaviour.

**Gene-Concept Network for GO Terms**
- This network reveals the connections between enriched GO terms and associated genes, illustrating the complex interplay of PRMT4-regulated genes across pathways such as immune response, cytokine production, and ion homeostasis.
- The network structure highlights key genes linked to multiple GO terms, suggesting that PRMT4 knockdown impacts critical regulatory hubs within myeloid differentiation processes, particularly those related to immune modulation and cell signalling.

**Heatmap of GO Terms (A) and Heatmap of GO Terms with Fold Changes (B)**
- **General Heatmap (A)**: This plot shows the presence of significant genes across enriched GO terms related to immune response, ion transport, and inflammatory processes, emphasising the pathways most affected by PRMT4 knockdown.
- **Heatmap with Fold Changes (B)**: By incorporating fold change values, this plot highlights the direction and magnitude of gene regulation, offering insights into upregulated (red) and downregulated (blue) genes within the enriched pathways, which align with the expected outcomes of myeloid differentiation and immune activation.

**Tree Plot of GO Enriched Terms**
- **Clustered Immune and Inflammatory Pathways**: The tree plot shows closely related immune response and inflammatory pathways grouped together, reflecting PRMT4 knockdown's effect on activating immune defence and inflammatory responses.
- **Distinct Clusters of Ion Transport and Chemotaxis**: Separate clusters for ion homeostasis and chemotaxis highlight PRMT4's regulatory role in cellular signalling, transport, and movement, which are crucial in myeloid differentiation.

### GSEA

**GSEA Enrichment Plots - Summary**
The GSEA enrichment plots depict the running enrichment scores for various hallmark pathways, highlighting pathways that are upregulated or downregulated in response to PRMT4 knockdown in human CD34+ cells.
- **Top Enriched Pathways**: Key pathways with significant enrichment include **TNFA signaling via NFKB**, **hypoxia response**, **inflammatory response**, and **interferon gamma response**, suggesting an enhanced immune and inflammatory reaction due to PRMT4 knockdown.
- **Downregulated Pathways**: Pathways like **E2F targets**, **G2M checkpoint**, **MYC targets**, and **DNA repair** are among the downregulated sets, indicating a decrease in proliferation and cell cycle progression in the PRMT4 knockdown condition. This aligns with the goal of targeting PRMT4 to suppress oncogenic processes in AML.
These plots collectively reinforce PRMT4’s regulatory role in immune response pathways and cell proliferation, supporting its potential as a therapeutic target in myeloid malignancies.

**Dot Plot of GSEA Enriched Pathways**
- **Key Pathways**: This dot plot highlights significant pathways affected by PRMT4 knockdown, with the **MYC targets** and **E2F targets** pathways showing the highest enrichment ratios, indicating a strong influence of PRMT4 on cell cycle and proliferation mechanisms.
- **Immune and Inflammatory Response**: Pathways like **TNFA signaling via NFKB**, **interferon responses**, and **inflammatory response** are enriched, further suggesting that PRMT4 knockdown enhances immune-related pathways, which could contribute to differentiation and anti-proliferative effects in AML cells.
The dot size reflects the count of genes in each pathway, while the colour gradient indicates significance, with darker red signifying more significant enrichment.

**Ridge Plot of GSEA Enriched Pathways**
- **Distribution of Pathway Enrichment Scores**: This ridge plot visualises the distribution of enrichment scores across significant pathways, indicating the breadth and intensity of PRMT4 knockdown effects. Pathways with high peaks, like **TNFA signaling via NFKB** and **inflammatory response**, show strong and centralised enrichment, supporting the role of PRMT4 in immune modulation.
- **Cell Cycle and Proliferation Impact**: Pathways associated with cell cycle progression, such as **MYC targets** and **E2F targets**, also exhibit high enrichment, consistent with the hypothesis that PRMT4 knockdown may reduce proliferation in AML cells by promoting differentiation.
The colour gradient represents p-adjust values, with deeper shades of red indicating higher significance.

**Top Enriched Disease Terms - Bar Plot**
- **Disease Associations**: This bar plot reveals that PRMT4 knockdown in CD34+ cells is associated with enrichment in disease terms related to developmental abnormalities, such as **hypoplasia of corpus callosum** and **cerebral atrophy**. These associations may indicate underlying pathways linked to differentiation defects in hematopoietic cells.
- **Significance of Findings**: The p-adjusted values highlight a strong statistical association with specific disease phenotypes, suggesting that PRMT4 may play a role in pathways affecting both hematopoiesis and disease phenotypes in AML and related conditions.
The colour gradient reflects p-adjust values, with deeper red indicating higher statistical significance.

**Top Enriched Disease Terms - Dot Plot**
- **Disease Associations and Gene Ratio**: The dot plot shows enriched disease terms, highlighting the **gene ratio** for each term, which represents the proportion of genes associated with a particular disease term relative to the background. Higher gene ratios, such as for **Byzantine arch palate** and **low set ears**, suggest a stronger relevance of these terms in the context of PRMT4 knockdown.
- **Visual Representation of Significance**: The colour gradient in this plot represents the p-adjust values, similar to the bar plot, indicating statistical significance, while dot size represents the number of genes associated with each term, providing an additional layer of context.
This plot complements the bar plot by illustrating both the statistical and proportional relevance of each disease term.

**Enrichment Map of Disease Terms**
- **Clustering of Related Disease Terms**: The enrichment map visually clusters similar disease terms, revealing potential connections between conditions. For example, terms like **cerebral atrophy**, **cerebellar atrophy**, and **spasticity** are closely linked, indicating shared pathways or genetic associations affected by PRMT4 knockdown.
- **Significance and Association Strength**: The node size indicates the number of genes associated with each term, and colour gradient denotes the significance (p-adjust values), highlighting the most statistically enriched disease associations. This map provides an intuitive view of how PRMT4 knockdown may influence related disease pathways, especially in developmental and neurological conditions.

## Conclusions

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
  
## Summary

This project presents a detailed analysis of the transcriptional effects of PRMT4 knockdown on human CD34+ cells, highlighting PRMT4's multifaceted role in various biological and disease-related pathways. The results emphasise PRMT4's influence on immune response, inflammatory regulation, and cellular homeostasis, as well as its association with neurological and developmental diseases. These findings suggest that targeting PRMT4 may not only promote myeloid differentiation but also impact pathways relevant to oncogenic processes in AML and other hematological malignancies. This study enhances the understanding of PRMT4 as a therapeutic target, with broader implications for treatment strategies in cancer and potentially other disease contexts influenced by epigenetic regulation.

## License

This project is licensed under the MIT License.

---

By exploring the biological processes, pathways, and diseases associated with PRMT4-regulated genes, this analysis provides further insight into PRMT4’s role in hematopoiesis and leukemia. The results suggest that PRMT4 is not only a regulator of myeloid differentiation but also a promising target for AML therapy.
