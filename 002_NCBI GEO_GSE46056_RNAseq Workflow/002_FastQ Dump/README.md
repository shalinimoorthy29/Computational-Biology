# Conversion of SRA Files to FASTQ Files

1. Check version of fastq-dump tool inside sra toolkit
   ```bash
   fastq-dump --version
   ```
   
Output should look like: fastq-dump : 3.1.1

2. Navigate to the working directory
   ```bash
   cd ~/compbio/GSE46056_RNAseq
   ```

3. Create a directory to store FASTQ files
   ```bash
   mkdir ~/compbio/GSE46056_RNAseq/fastq_files
   ```

4. Convert .sra files to FASTQ using fastq-dump
   ```bash
    for file in ~/compbio/GSE46056_RNAseq/sra_files/*.sra; do
        fastq-dump --split-files --outdir ~/compbio/GSE46056_RNAseq/fastq_files $file
    done
   ```

Output should look like:

Read 27640449 spots for /home/shali/compbio/GSE46056_RNAseq/sra_files/SRR827457.sra

Written 27640449 spots for /home/shali/compbio/GSE46056_RNAseq/sra_files/SRR827457.sra
