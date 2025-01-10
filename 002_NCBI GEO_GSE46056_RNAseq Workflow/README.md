# RNAseq Workflow: A Bioinformatics Pipeline Using Linux and R

## Project Overview
This project demonstrates a bioinformatics pipeline for RNAseq data analysis. The workflow is divided into two stages:
1. **Data Acquisition and Pre-processing**: Conducted in Linux using tools like the SRA Toolkit, FastQC, and HISAT2.
2. **Analysis and Visualisation**: Performed in R using tools such as DESeq2 and ggplot2 for differential expression analysis and plotting.

The pipeline is designed for a beginner-friendly setup in WSL (Windows Subsystem for Linux) with Ubuntu.

---

## Steps Completed So Far

### **1. Data Acquisition and Pre-processing**

#### 1.1 Setting Up the Environment
- **Tools Installed**:
  - **SRA Toolkit** for downloading RNAseq data from the Sequence Read Archive (SRA).
  - `wget` for general file downloads.
- **SRA Toolkit Setup**:
  1. Downloaded the **SRA Toolkit** for Ubuntu Linux 64-bit architecture (`sratoolkit.3.1.1-ubuntu64.tar`) from the official [NCBI GitHub page](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit).
  2. Extracted the `.tar` file using the following command:
     ```bash
     tar -xvzf sratoolkit.3.1.1-ubuntu64.tar.gz
     ```
  3. Added the `bin` directory of the extracted SRA Toolkit to the system's `PATH`, enabling terminal access to its tools:
     ```bash
     export PATH=$PATH:~/compbio/sratoolkit/sratoolkit.3.1.1-ubuntu64/bin
     ```
     - **Explanation for Beginners**: The `PATH` variable tells the Linux system where to look for executable programs. By adding the `bin` directory of the SRA Toolkit to the `PATH`, I can run its tools (e.g., `prefetch`, `fastq-dump`) from anywhere in the terminal without needing to specify the full path.

#### 1.2 Data Acquisition
- Downloaded raw `.sra` files using the **SRA Toolkit's `prefetch` utility**.
- Managed an SRR accession list with 30 IDs to batch download the required data.
- Verified successful downloads and corrected incomplete `.sra` files.

#### 1.3 Preparing for Conversion
- Used the `fastq-dump` utility from the SRA Toolkit to convert `.sra` files into paired-end FASTQ files.
- Configured the process to generate `_1.fastq` and `_2.fastq` for forward and reverse reads.

---

### **2. Analysis and Visualisation**

#### Planned Steps
1. **Differential Expression Analysis**:
   - Perform statistical analysis using **DESeq2** or **edgeR** to identify differentially expressed genes.
2. **Visualisation**:
   - Create exploratory plots in R, including:
     - **PCA plots** to assess data variability.
     - **Heatmaps** to visualise gene expression patterns.
     - **Volcano plots** to highlight significant genes.

---

## Technologies Used
- **Linux** (Ubuntu via WSL) for data acquisition and pre-processing.
- **SRA Toolkit, FastQC, Trimmomatic, HISAT2** for processing raw sequencing data (to be implemented).
- **R** for statistical analysis and visualisation.

---

## Project Goals
1. Develop a reproducible pipeline for RNAseq analysis.
2. Gain proficiency with Linux command-line tools in bioinformatics.
3. Document and share the workflow for future reference and collaboration.
