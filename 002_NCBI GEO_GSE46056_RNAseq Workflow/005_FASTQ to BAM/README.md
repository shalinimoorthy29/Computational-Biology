# Direct Conversion of FASTQ Files to BAM Files Using HISAT2 and SAMtools

This README describes the process of aligning FASTQ files to a reference genome and directly converting the output to BAM files to save disk space and processing time.

---

### 1: Navigate to the Directory for Alignment

```bash
cd ~/compbio/GSE46056_RNAseq
```

### 2: Create a Directory for BAM Files

```bash
mkdir ~/compbio/GSE46056_RNAseq/bam_files
```

### 3: Test Alignment for One FASTQ Pair
(The -p 4 option specifies the use of 4 threads for the alignment. Adjust this based on available CPU cores)

```bash
hisat2 -p 4 -x ~/compbio/hisat2_index/grch38/genome \
-1 ~/compbio/GSE46056_RNAseq/fastq_files/SRR827457_1.fastq \
-2 ~/compbio/GSE46056_RNAseq/fastq_files/SRR827457_2.fastq | \
samtools view -bS - > ~/compbio/GSE46056_RNAseq/bam_files/SRR827457.bam
```


Remove the Test BAM File (Optional)

```bash
rm ~/compbio/GSE46056_RNAseq/bam_files/SRR827457.bam
```

### 4: Automate the Alignment for All FASTQ Pairs
Use a loop to process all paired FASTQ files in the test_fastq_files directory and directly convert them to BAM files:

```bash
for file in ~/compbio/GSE46056_RNAseq/fastq_files/*_1.fastq; do
    base=$(basename $file _1.fastq)
    hisat2 -p 4 -x ~/compbio/hisat2_index/grch38/genome \
    -1 ~/compbio/GSE46056_RNAseq/fastq_files/${base}_1.fastq \
    -2 ~/compbio/GSE46056_RNAseq/fastq_files/${base}_2.fastq | \
    samtools view -bS - > ~/compbio/GSE46056_RNAseq/bam_files/${base}.bam
done
```

Monitor Alignment Output (will look like this):

```
27640449 reads; of these:
  27640449 (100.00%) were paired; of these:
    2087594 (7.55%) aligned concordantly 0 times
    23931404 (86.58%) aligned concordantly exactly 1 time
    1621451 (5.87%) aligned concordantly >1 times
    ----
    2087594 pairs aligned concordantly 0 times; of these:
      163721 (7.84%) aligned discordantly 1 time
    ----
    1923873 pairs aligned 0 times concordantly or discordantly; of these:
      3847746 mates make up the pairs; of these:
        2303600 (59.87%) aligned 0 times
        1302803 (33.86%) aligned exactly 1 time
        241343 (6.27%) aligned >1 times
95.83% overall alignment rate
```

### 5: BAM File Sorting

```bash
for bam in ~/compbio/GSE46056_RNAseq/bam_files/*.bam; do
    samtools sort -o ${bam%.bam}.sorted.bam $bam
    mv ${bam%.bam}.sorted.bam $bam
done
```

### 6: BAM File Indexing

```bash
for bam in ~/compbio/GSE46056_RNAseq/bam_files/*.bam; do
    samtools index $bam
done
```

### 7: Validate BAM File Headers

```bash
samtools view -H ~/compbio/GSE46056_RNAseq/bam_files/SRR827457.bam
```

Output will look like this:

```
@HD     VN:1.0  SO:coordinate
@SQ     SN:1    LN:248956422
@SQ     SN:10   LN:133797422
@SQ     SN:11   LN:135086622
...
@PG     ID:hisat2       PN:hisat2       VN:2.2.1        CL:"/usr/bin/hisat2-align-s --wrapper basic-0 -p 4 ..."
@PG     ID:samtools     PN:samtools     PP:hisat2       VN:1.19.2       CL:samtools view -bS -
...
```
