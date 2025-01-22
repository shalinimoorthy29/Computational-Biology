# Generating FASTQC Reports From FASTQ Files

1. Check if FastQC is installed
   ```bash
   fastqc --version
   ```
If it is not installed, output will look like: Command 'fastqc' not found, but can be installed with: sudo apt install fastqc

2. Install FastQC if not already installed
   ```bash
   sudo apt update
   sudo apt install fastqc
   ```

3. Verify installation
   ```bash
   fastqc --version
   ```
Output will look like: FastQC v0.12.1


4. Navigate to the working directory
   ```bash
   cd ~/compbio/sratoolkit
   ```

5. Create a directory to store FastQC reports
   ```bash
   mkdir ~/compbio/sratoolkit/test_fastqc_files
   ```

6. Run FastQC on all FASTQ files in the specified directory
   ```bash
   fastqc ~/compbio/sratoolkit/test_fastq_files/*.fastq --outdir ~/compbio/sratoolkit/test_fastqc_files/
   ```

7. Check the FastQC reports in the output directory
   ```bash
   ls ~/compbio/sratoolkit/test_fastqc_files
