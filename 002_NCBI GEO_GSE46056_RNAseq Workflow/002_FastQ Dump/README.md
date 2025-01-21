# Conversion of SRA Files to FASTQ Files

1. Check version of fastq-dump tool inside sra toolkit
   ```bash
   fastq-dump --version
   ```
Output should look like: fastq-dump : 3.1.1

2. Navigate to the SRA Toolkit directory
   ```bash
   cd ~/compbio/sratoolkit

3. Create a directory to store .sra files
   ```bash
   mkdir ~/compbio/sratoolkit/test_sra_files

4. Copy .sra files to the new directory
   ```bash
   cp ~/compbio/sratoolkit/SRR827457/*.sra ~/compbio/sratoolkit/test_sra_files/
   cp ~/compbio/sratoolkit/SRR827459/*.sra ~/compbio/sratoolkit/test_sra_files/
   cp ~/compbio/sratoolkit/SRR827461/*.sra ~/compbio/sratoolkit/test_sra_files/
   cp ~/compbio/sratoolkit/SRR827474/*.sra ~/compbio/sratoolkit/test_sra_files/
   cp ~/compbio/sratoolkit/SRR827476/*.sra ~/compbio/sratoolkit/test_sra_files/
   cp ~/compbio/sratoolkit/SRR827478/*.sra ~/compbio/sratoolkit/test_sra_files/

5. List all .sra files in the directory
   ```bash
   ls ~/compbio/sratoolkit/test_sra_files
   ```

Output should look like: SRR827457.sra  SRR827459.sra  SRR827461.sra  SRR827474.sra  SRR827476.sra  SRR827478.sra ...

6. Create a directory to store FASTQ files
   ```bash
   mkdir ~/compbio/sratoolkit/test_fastq_files

7. Convert .sra files to FASTQ using fastq-dump
   ```bash
    for file in ~/compbio/sratoolkit/test_sra_files/*.sra; do
        fastq-dump --split-files --outdir ~/compbio/sratoolkit/test_fastq_files $file
    done
   ```

Output should look like:

Read 27640449 spots for /home/shali/compbio/sratoolkit/test_sra_files/SRR827457.sra

Written 27640449 spots for /home/shali/compbio/sratoolkit/test_sra_files/SRR827457.sra
