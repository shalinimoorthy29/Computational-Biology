# Investigating Loss of Heterozygosity (LoH) in Glioblastoma

## Project Overview
Glioblastoma is a highly aggressive and recurrent brain tumour. Loss of Heterozygosity (LoH) refers to the loss of one allele at a specific locus, which can occur due to mechanisms such as mitotic recombination, chromosomal deletion, or uniparental disomy. LoH is significant in cancer as it often leads to genomic instability, inactivation of tumour suppressor genes, and therapeutic resistance.

This project investigates LoH events in glioblastoma using genomic data to better understand tumour evolution and progression.

## Objectives
1. Investigate evidence for loss of heterozygosity (LoH) in glioblastoma.
2. Identify genes and regions in the genome impacted by LoH and their frequency.
3. Determine whether LoH events are unique to recurrent glioblastoma.

## Datasets Used
1. `variants_gene_copy_number.csv`: Provides gene-level copy number variation (CNV) data.
2. `variants_titan_seg.csv`: Contains segment-level CNV data, enabling detection of regions with significant genomic alterations.
3. `ref_genes.csv`: Includes genomic coordinates for mapping LoH data to regions in the genome.

## Analysis Steps
### Data Preprocessing
- Checked for duplicates and missing values in all datasets.
- Unzipped compressed files and ensured proper loading of data.
- Conducted initial exploration of dataset structures.

### Investigating Evidence for LoH
- Summarised genomic alteration types using segment-level data.
- Categorised LoH events into NLOH (Copy Neutral), DLOH (Hemizygous Deletion), and ALOH (Amplified).
- Visualised counts and proportions of alteration types across all samples.

### Identifying Genes Impacted by LoH
- Mapped LoH segments to genes using genomic coordinates.
- Categorised LoH types for impacted genes.
- Summarised the top genes most frequently impacted by LoH and their breakdown by LoH type.

### Analysing Recurrent LoH Events
- Categorised samples into primary, recurrent, and metastatic based on barcodes.
- Identified unique LoH events in recurrent samples compared to primary tumours.
- Summarised recurrent LoH events by gene and sample type.

## Results
### 1. Evidence for LoH in Glioblastoma
- **NLOH** was identified as the most common LoH type, providing evidence of genomic instability in glioblastoma.

### 2. Genes and Regions Impacted by LoH
- **HLA-DRB1** was the most frequently impacted gene.
- Chromosome 1 exhibited the highest number of LoH events.

### 3. Unique LoH Events in Recurrent Tumours
- Unique LoH events in recurrent samples were predominantly **ALOH**.
- The **KDM5A** gene had the highest number of unique LoH events in recurrent samples.

## Conclusion
This analysis highlights the significance of LoH in glioblastoma progression:
- NLOH was prevalent, suggesting genomic instability across tumour samples.
- HLA-DRB1 and Chromosome 1 were highly impacted by LoH events.
- Unique recurrent LoH events were dominated by ALOH, particularly in genes like KDM5A.

Future work could expand on these findings by integrating transcriptomic data to further explore the functional impact of LoH events.
