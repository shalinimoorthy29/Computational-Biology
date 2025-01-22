# Generating FASTQC and MULTIQC Reports From FASTQ Files

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


8. Check if MultiQC is installed
   ```bash
   multiqc --version
   ```

If it is not installed, output will look like: Command 'multiqc' not found, but can be installed with: sudo apt install multiqc


9. Install MultiQC if not already installed
   ```bash
   sudo apt update
   sudo apt install multiqc
   ```

10. Verify installation
   ```bash
   multiqc --version
   ```
Output will look like: MultiQC, version 1.18

11. Navigate to the working directory
   ```bash
   cd ~/compbio/sratoolkit
   ```

12. Create a directory to store MultiQC outputs
   ```bash
   mkdir ~/compbio/sratoolkit/test_multiqc_file
   ```

13. Run MultiQC on the FastQC output directory
   ```bash
   multiqc ~/compbio/sratoolkit/test_fastqc_files -o ~/compbio/sratoolkit/test_multiqc_file/
   ```

14. Check the generated MultiQC report
   ```bash
   ls ~/compbio/sratoolkit/test_multiqc_file
   ```
