# Direct Conversion of FASTQ Files to BAM Files Using HISAT2 and SAMtools

This README describes the process of aligning FASTQ files to a reference genome and directly converting the output to BAM files to save disk space and processing time.

---

## Step 1: Navigate to the Directory for Alignment

```bash
cd ~/compbio/sratoolkit
```

## Step 2: Create a Directory for BAM Files

```bash
mkdir ~/compbio/sratoolkit/test_bam_files
```

## Step 3: Test Alignment for One FASTQ Pair
(The -p 4 option specifies the use of 4 threads for the alignment. Adjust this based on available CPU cores)

```bash
hisat2 -p 4 -x ~/compbio/hisat2_index/grch38/genome \
-1 ~/compbio/sratoolkit/test_fastq_files/SRR827457_1.fastq \
-2 ~/compbio/sratoolkit/test_fastq_files/SRR827457_2.fastq | \
samtools view -bS - > ~/compbio/sratoolkit/test_bam_files/SRR827457.bam
```


Remove the Test BAM File (Optional)

```bash
rm ~/compbio/sratoolkit/test_bam_files/SRR827457.bam
```

## Step 4: Automate the Alignment for All FASTQ Pairs
Use a loop to process all paired FASTQ files in the test_fastq_files directory and directly convert them to BAM files:

```bash
for file in ~/compbio/sratoolkit/test_fastq_files/*_1.fastq; do
    base=$(basename $file _1.fastq)
    hisat2 -p 4 -x ~/compbio/hisat2_index/grch38/genome \
    -1 ~/compbio/sratoolkit/test_fastq_files/${base}_1.fastq \
    -2 ~/compbio/sratoolkit/test_fastq_files/${base}_2.fastq | \
    samtools view -bS - > ~/compbio/sratoolkit/test_bam_files/${base}.bam
done
```

Monitor Alignment Output (will look like this):

```
27640449 reads; of these:
  27640449 (100.00%) were paired; of these:
    2087594 (7.55%) aligned concordantly 0 times
    23931404 (86.58%) aligned concordantly exactly 1 time
    1621451 (5.87%) aligned concordantly >1 times
95.83% overall alignment rate
```

## Step 5: BAM File Sorting

```bash
for bam in ~/compbio/sratoolkit/test_bam_files/*.bam; do
    samtools sort -o ${bam%.bam}.sorted.bam $bam
    mv ${bam%.bam}.sorted.bam $bam
done
```

## Step 6: BAM File Indexing

```bash
for bam in ~/compbio/sratoolkit/test_bam_files/*.bam; do
    samtools index $bam
done
```

## Step 7: Validate BAM File Headers

```bash
samtools view -H ~/compbio/sratoolkit/test_bam_files/SRR827457.bam
```

Output will look like this:

@HD     VN:1.0  SO:coordinate
@SQ     SN:1    LN:248956422
@SQ     SN:10   LN:133797422
@SQ     SN:11   LN:135086622
...
@PG     ID:hisat2       PN:hisat2       VN:2.2.1        CL:"/usr/bin/hisat2-align-s --wrapper basic-0 -p 4 ..."
@PG     ID:samtools     PN:samtools     PP:hisat2       VN:1.19.2       CL:samtools view -bS -
...



