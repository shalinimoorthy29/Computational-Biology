# Generating raw read counts per gene from BAM files using featureCounts.

### 1. Navigate to directory

```bash
cd ~/compbio
```

### 2. Create a folder in this directory called ensemble_gene_annotation

```bash
mkdir ensemble_gene_annotation
```

### 3. Download the gene annotation file from ENSEMBLE

Download Homo_sapiens.GRCh38.113.gtf.gz file from ENSEMBLE and save it into ensemble_gene_annotation folder.
(https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/)

### 4. Unzip the file

```bash
gunzip ~/compbio/ensemble_gene_annotation/Homo_sapiens.GRCh38.113.gtf.gz
```

### 5. Create a new folder called gene_counts so the raw counts can be saved here.

```bash
cd ~/compbio/GSE46056_RNAseq
mkdir ~/compbio/GSE46056_RNAseq/gene_counts
```

### 6. Perform featureCounts

```bash
featureCounts -T 4 -p -t exon -g gene_id \
-a ~/compbio/ensemble_gene_annotation/Homo_sapiens.GRCh38.113.gtf \
-o ~/compbio/GSE46056_RNAseq/gene_counts/gene_counts.txt \
~/compbio/GSE46056_RNAseq/bam_files/*.bam
```

### 7. Check first few lines of output text file

```bash
head -n 20 ~/compbio/GSE46056_RNAseq/gene_counts/gene_counts.txt
```

### 8. Clean up the gene_counts.txt file so only important columns are kept

```bash
cut -f1,7- ~/compbio/GSE46056_RNAseq/gene_counts/gene_counts.txt > ~/compbio/GSE46056_RNAseq/gene_counts/clean_gene_counts.txt
```

### 9. Check first few lines of cleaned text file

```bash
head ~/compbio/GSE46056_RNAseq/gene_counts/clean_gene_counts.txt
```

### 10. Convert to csv file

```bash
cat ~/compbio/GSE46056_RNAseq/gene_counts/clean_gene_counts.txt | tr '\t' ',' > ~/compbio/GSE46056_RNAseq/gene_counts/clean_gene_counts.csv
```