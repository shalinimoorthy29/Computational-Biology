# SRA Toolkit Installation and Retrieval of .sra Files
---

## SRA Toolkit Installation and Setup

### 1. Navigate to main working directory and create all the necessary sub-folders:

   ```bash
   cd /home/shali/compbio
   mkdir -p /home/shali/compbio/GSE46056_RNAseq
   cd /home/shali/compbio/GSE46056_RNAseq
   mkdir ~/compbio/GSE46056_RNAseq/sra_files
   mkdir ~/compbio/sratookit
   ```

### 2. Download the SRA Toolkit from [here](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit).
   ![SRA Toolkit Installation Screenshot](sra%20toolkit%20installation.png "SRA Toolkit Installation")

### 3. Copy the downloaded SRA Toolkit tarball to the working directory: 
 
   ```bash
   cp /mnt/c/Users/shali/Downloads/sratoolkit.3.1.1-ubuntu64.tar.gz ~/compbio/sratoolkit
   ```

### 4. Extract the tarball:

   ```bash
   tar -xvzf sratoolkit.3.1.1-ubuntu64.tar.gz
   ```

### 5. Add the toolkit binaries to the PATH:

   ```bash
   export PATH=$PATH:~/compbio/GSE46056_RNAseq/sratoolkit.3.1.1-ubuntu64/bin
   ```

### 6. Create a text file with a list of SRR accession numbers:

   ```bash
   nano ~/compbio/GSE46056_RNAseq/srr_accessions.txt
   ```

Example shown below:

SRR827457

SRR827462

SRR827467

... until the final ID

One SRR accession IDs have been entered into the nano terminal, press Ctrl+O to save the file and Ctrl+X to close the nano terminal. 

### 7. Use prefetch to download the .sra files for each accession number:

   ```bash
   cat ~/compbio/GSE46056_RNAseq/srr_accessions.txt | xargs -I {} prefetch --output-directory ~/compbio/GSE46056_RNAseq/sra_files {}
   ```

### 8. SRA files (.sra) are saved into individual folders and better to be moved to parent directory:

   ```bash
   cd ~/compbio/GSE46056_RNAseq/sra_files
   find . -type f -name "*.sra" -exec mv {} ./ \;
   find . -type d -empty -delete
   ```

This will bring them out of their folders and following that, the empty folders will be deleted. 

### 9. List all .sra files in the directory

   ```bash
   ls ~/compbio/GSE46056_RNAseq/sra_files

Output should look like: SRR827457.sra  SRR827462.sra  SRR827467.sra  SRR827472.sra  SRR827477.sra  SRR827482.sra ...


