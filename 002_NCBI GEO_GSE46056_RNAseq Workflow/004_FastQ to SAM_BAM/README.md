# Generating BAM Files from FASTQ Files Using HISAT2

## Step 1: Create HISAT2 Index Directory and Download Genome Index

```bash
mkdir -p ~/compbio/hisat2_index
cd ~/compbio/hisat2_index
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
```

## Step 2: Extract the Downloaded Genome Index

```bash
tar -xzvf grch38_genome.tar.gz
```

## Step 3: Create Directories for SAM and BAM Files

```bash
mkdir -p ~/compbio/sratoolkit/test_sam_files
mkdir -p ~/compbio/sratoolkit/test_bam_files
```

## Step 4: Install HISAT2 and SAMtools

```bash
sudo apt update
sudo apt install hisat2 samtools
hisat2 --version
samtools --version
```

## Step 5: Align Reads and Generate BAM Files Directly

```bash
for file in ~/compbio/sratoolkit/test_fastq_files/*_1.fastq; do
    base=$(basename $file _1.fastq)
    hisat2 -p 4 -x ~/compbio/hisat2_index/grch38/genome \
    -1 ~/compbio/sratoolkit/test_fastq_files/${base}_1.fastq \
    -2 ~/compbio/sratoolkit/test_fastq_files/${base}_2.fastq | \
    samtools view -bS - > ~/compbio/sratoolkit/test_bam_files/${base}.bam
done
```

## Step 6: Sort and Index BAM Files

```bash
for bam in ~/compbio/sratoolkit/test_bam_files/*.bam; do
    samtools sort -o ${bam%.bam}.sorted.bam $bam
    mv ${bam%.bam}.sorted.bam $bam
done
for bam in ~/compbio/sratoolkit/test_bam_files/*.bam; do
    samtools index $bam
done
```
 