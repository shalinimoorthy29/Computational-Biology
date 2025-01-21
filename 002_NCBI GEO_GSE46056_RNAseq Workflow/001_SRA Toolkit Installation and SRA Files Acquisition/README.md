# SRA Toolkit Installation and Retrieval of .sra Files

## SRA Toolkit Installation and Setup

1. Navigate to your desired directory:
   ```bash
   cd /home/shali/compbio
   mkdir -p /home/shali/compbio/sratoolkit
   cd /home/shali/compbio/sratoolkit

2. Download the SRA Toolkit from [here](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit).
   ![SRA Toolkit Installation Screenshot](sra%20toolkit%20installation.png "SRA Toolkit Installation")

3. Copy the downloaded SRA Toolkit tarball to the working directory:  
   ```bash
   cp /mnt/c/Users/shali/Downloads/sratoolkit.3.1.1-ubuntu64.tar.gz .
   
4. Extract the tarball:
   ```bash
   tar -xvzf sratoolkit.3.1.1-ubuntu64.tar.gz

5. Add the toolkit binaries to the PATH:
   ```bash
   export PATH=$PATH:~/compbio/sratoolkit/sratoolkit.3.1.1-ubuntu64/bin

6. Check the version of fastq-dump (or any other tool from the SRA Toolkit)
   ```bash
   fastq-dump --version

7. Create a text file with a list of SRR accession numbers:
   ```bash
   nano ~/compbio/sratoolkit/srr_accessions.txt
   ```
Example shown below:

SRR827457

SRR827458

SRR827459

... until the final ID

One SRR accession IDs have been entered into the nano terminal, press Ctrl+O to save the file and Ctrl+X to close the nano terminal. 

8. Use prefetch to download the .sra files for each accession number:
   ```bash
   cat ~/compbio/sratoolkit/srr_accessions.txt | xargs -n 1 prefetch

9. Check the downloaded files:
   ```bash
   ls
