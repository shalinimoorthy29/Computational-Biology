# -----------------------------------------------------------------------------------------------
# Investigating Loss of Heterozygosity (LoH) in Glioblastoma
# -----------------------------------------------------------------------------------------------

# Download all required datasets to working directory:
# variants_gene_copy_number.csv.gz (Needs extracting)
# variants_titan_seg.csv
# ref_genes.csv

# Set working directory (will be different for another user of this script)
setwd("C:/Users/shali/Documents/L&D/GitHub Projects/Computational Biology/003_LoH_Analysis_in_Glioblastoma")

# Start Timer for Debugging
start_time <- Sys.time()

# Debugging: Environment Check
cat("\n--- DEBUGGING START ---\n")
cat("Current Working Directory: ", getwd(), "\n")

# Install required CRAN and Bioconductor packages
cran_packages <- c("R.utils", "tidyverse", "data.table", "ggplot2", "RColorBrewer")
install.packages(setdiff(cran_packages, installed.packages()[, "Package"]))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioconductor_packages <- c("GenomicRanges")
BiocManager::install(setdiff(bioconductor_packages, installed.packages()[, "Package"]))

# Load libraries
library(R.utils)         # For file utilities like gunzip
library(tidyverse)       # Data manipulation and visualisation
library(data.table)      # Efficient data handling
library(GenomicRanges)   # Genomic range operations
library(ggplot2)         # Visualisation
library(RColorBrewer)    # Colour palettes

# Debugging: File Check
input_files <- c("variants_gene_copy_number.csv.gz", "variants_titan_seg.csv", "ref_genes.csv")
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop("Missing Input Files: ", paste(missing_files, collapse = ", "))
} else {
  cat("All Input Files Found.\n")
}

# Unzip the file (this will keep the original .gz file)
if (!file.exists("variants_gene_copy_number.csv")) {
  gunzip("variants_gene_copy_number.csv.gz", remove = FALSE)
  cat("File 'variants_gene_copy_number.csv.gz' successfully unzipped.\n")
} else {
  cat("File 'variants_gene_copy_number.csv' already exists.\n")
}

# -----------------------------------------------------------------------------------------------
# Load all datasets required for this project
# -----------------------------------------------------------------------------------------------
cat("\nLoading Datasets...\n")
titan_seg <- read.csv("variants_titan_seg.csv", header = TRUE)
gene_copy_number <- read.csv("variants_gene_copy_number.csv", header = TRUE)
ref_genes <- read.csv("ref_genes.csv", header = TRUE)
cat("\nDatasets loaded successfully.\n")


# -----------------------------------------------------------------------------------------------
# Inspect titan_seg dataset
# -----------------------------------------------------------------------------------------------
cat("\nTITAN Segment Data Structure:\n")
str(titan_seg)

# Check for duplicates in titan_seg
duplicates <- titan_seg %>% filter(duplicated(.))
cat("\nNumber of Duplicate Rows in TITAN Segment Data: ", nrow(duplicates), "\n")

# Summarise the unique barcodes
unique_barcodes <- titan_seg %>% distinct(pair_barcode) %>% nrow()
cat("\nNumber of Unique Barcodes: ", unique_barcodes, "\n")

# Check for missing values in titan_call
total_na_titan_call <- titan_seg %>% filter(is.na(titan_call)) %>% nrow()
cat("\nNumber of Missing Values in titan_call Column: ", total_na_titan_call, "\n")

# Inspect unique values in titan_call
unique_titan_calls <- titan_seg %>% distinct(titan_call)
cat("\nUnique Values in titan_call Column:\n")
print(unique_titan_calls)

# -----------------------------------------------------------------------------------------------
# Inspect gene_copy_number dataset
# -----------------------------------------------------------------------------------------------
cat("\nGene Copy Number Data Structure:\n")
str(gene_copy_number)

# Check for duplicates in gene_copy_number
duplicates_gcn <- gene_copy_number %>% filter(duplicated(.))
cat("\nNumber of Duplicate Rows in Gene Copy Number Data: ", nrow(duplicates_gcn), "\n")

# Check for missing values in hlvl_call
total_na_hlvl_call <- gene_copy_number %>% filter(is.na(hlvl_call)) %>% nrow()
cat("\nNumber of Missing Values in hlvl_call Column: ", total_na_hlvl_call, "\n")

# Inspect unique values in hlvl_call
unique_hlvl_calls <- gene_copy_number %>% distinct(hlvl_call)
cat("\nUnique Values in hlvl_call Column:\n")
print(unique_hlvl_calls)


# -----------------------------------------------------------------------------------------------
# QUESTION 1: What evidence is there for loss of heterozygosity (LoH) in glioblastoma?
# -----------------------------------------------------------------------------------------------

# General Distribution of Genomic Alteration Events
general_distribution <- titan_seg %>%
  group_by(titan_call) %>%
  summarise(
    count = n(),
    frequency = (count / nrow(titan_seg)) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(count))

cat("\nGeneral Distribution of Genomic Alteration Events:\n")
print(general_distribution)

# Bar Plot: Genomic Alteration Types Across All Samples (COUNT)
count_plot <- ggplot(general_distribution, aes(x = reorder(titan_call, -count), y = count, fill = titan_call)) +
  geom_bar(stat = "identity", fill = "darkgrey") +
  labs(title = "Genomic Alteration Types Across All Samples", x = "Alteration Type", y = "Count") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),        # Axis tick labels: size 12, black
    axis.text.x = element_text(angle = 45, hjust = 1),           # Rotate x-axis labels
    axis.title = element_text(size = 14, color = "black"),       # Axis titles: size 14, black
    axis.title.x = element_text(margin = margin(t = 10)),        # Increase spacing above x-axis title
    axis.title.y = element_text(margin = margin(r = 10)),        # Increase spacing to the right of y-axis title
    plot.title = element_text(size = 16, color = "black",        # Plot title: size 16, black
                              hjust = 0.5,                       # Centre-align the title
                              margin = margin(b = 10)),          # Add spacing below the title
  )
print(count_plot)

# Save plot 
ggsave(filename = "Genomic Alteration Types Across All Samples_COUNT.png", 
       plot = count_plot, 
       width = 8, 
       height = 6, 
       dpi = 300,
       bg = "white")

print("Plot saved as Genomic Alteration Types Across All Samples_COUNT.png")

# Pie Chart: Genomic Alteration Types (FREQUENCY)
pie_chart_data <- general_distribution %>%
  mutate(cumulative_angle = cumsum(frequency),
         mid_angle = cumulative_angle - (frequency / 2),
         x_pos = 1.5 * sin(mid_angle * pi / 50),  # Adjusting position for label
         y_pos = 1.5 * cos(mid_angle * pi / 50)) # Adjusting position for label

frequency_pie <- ggplot(general_distribution, aes(x = "", y = frequency, fill = titan_call)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  labs(
    title = "Proportion of Genomic Alteration Types",
    fill = "Alteration Type"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 10)),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  ) +
  geom_text(aes(label = paste0(round(frequency, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 4)
print(frequency_pie)

# Save plot 
ggsave(filename = "Proportion of Genomic Alteration Types.png", 
       plot = frequency_pie, 
       width = 8, 
       height = 6, 
       dpi = 300,
       bg = "white")

print("Plot saved as Proportion of Genomic Alteration Types.png")


# Group by aliquot_barcode and hlvl_call to count occurrences
hlvl_call_summary <- gene_copy_number %>%
  group_by(hlvl_call) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(count))

# Bar plot: Distribution of hlvl_call across all barcodes
hlvl_call_plot <- ggplot(hlvl_call_summary, aes(x = factor(hlvl_call), y = count, fill = factor(hlvl_call))) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  labs(
    title = "Distribution of High-Level Copy Number Calls",
    x = "Copy Number Call",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),        # Axis tick labels: size 12, black
    axis.text.x = element_text(angle = 45, hjust = 1),           # Rotate x-axis labels
    axis.title = element_text(size = 14, color = "black"),       # Axis titles: size 14, black
    axis.title.x = element_text(margin = margin(t = 10)),        # Increase spacing above x-axis title
    axis.title.y = element_text(margin = margin(r = 10)),        # Increase spacing to the right of y-axis title
    plot.title = element_text(size = 16, color = "black",        # Plot title: size 16, black
                              hjust = 0.5,                       # Centre-align the title
                              margin = margin(b = 10)),          # Add spacing below the title
  ) +
  scale_y_continuous(
    limits = c(0, 30000000),      # Set y-axis range from 0 to 1500
    breaks = seq(0, 30000000, 5000000),  # Add breaks at intervals of 500
    labels = scales::comma      # Use normal number formatting with commas
  )
print(hlvl_call_plot)

# Save plot 
ggsave(filename = "Distribution of High-Level Copy Number Calls.png", 
       plot = hlvl_call_plot, 
       width = 8, 
       height = 6, 
       dpi = 300,
       bg = "white")

print("Plot saved as Distribution of High-Level Copy Number Calls.png")

# Filter for genes with high-level deletions (LoH)
loh_genes <- gene_copy_number %>%
  filter(hlvl_call == -2) %>%
  group_by(gene_symbol) %>%
  summarise(
    loh_count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(loh_count))

# Bar plot: Top 20 genes with the most high-level deletions
loh_plot <- ggplot(loh_genes[1:20, ], aes(x = reorder(gene_symbol, -loh_count), y = loh_count)) +
  geom_bar(stat = "identity", fill = "darkred") +
  labs(
    title = "Top 20 Genes with High-Level Deletions",
    x = "Gene Symbol",
    y = "Number of High-Level Deletions"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),        # Axis tick labels: size 12, black
    axis.text.x = element_text(angle = 45, hjust = 1),           # Rotate x-axis labels
    axis.title = element_text(size = 14, color = "black"),       # Axis titles: size 14, black
    axis.title.x = element_text(margin = margin(t = 10)),        # Increase spacing above x-axis title
    axis.title.y = element_text(margin = margin(r = 10)),        # Increase spacing to the right of y-axis title
    plot.title = element_text(size = 16, color = "black",        # Plot title: size 16, black
                              hjust = 0.5,                       # Centre-align the title
                              margin = margin(b = 10)),          # Add spacing below the title
  )

print(loh_plot)

# Save plot 
ggsave(filename = "Top 20 Genes with High-Level Deletions.png", 
       plot = loh_plot, 
       width = 8, 
       height = 6, 
       dpi = 300,
       bg = "white")

print("Plot saved as Top 20 Genes with High-Level Deletions.png")

# -----------------------------------------------------------------------------------------------
# QUESTION 2: Where in the genome does LoH occur in glioblastoma (which genes are impacted) and
# how frequently?
# -----------------------------------------------------------------------------------------------

# Create GRanges for genes
genes_gr <- GRanges(seqnames = ref_genes$chrom,
                    ranges = IRanges(start = ref_genes$start, end = ref_genes$end),
                    mcols = ref_genes$gene_symbol)

# Filter for LoH events in titan_seg
loh_segments_q2 <- titan_seg %>%
  filter(titan_call %in% c("DLOH", "NLOH", "ALOH")) %>%
  select(chrom, start, end, titan_call)

# Create GRanges for LoH segments
loh_gr_q2 <- GRanges(seqnames = loh_segments_q2$chrom,
                     ranges = IRanges(start = loh_segments_q2$start, end = loh_segments_q2$end),
                     mcols = loh_segments_q2$titan_call)

# Find overlaps
overlap_q2 <- findOverlaps(loh_gr_q2, genes_gr)

# Extract overlapping genes and their LoH types
loh_genes_q2 <- data.frame(
  gene_symbol = mcols(genes_gr)$mcols[subjectHits(overlap_q2)],
  titan_call = mcols(loh_gr_q2)$mcols[queryHits(overlap_q2)]
)

# Summarise impacted genes
loh_gene_summary <- loh_genes_q2 %>%
  group_by(gene_symbol, titan_call) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))

cat("\nSummary of Genes Impacted by LoH:\n")
print(loh_gene_summary)

# Extract the top 10 genes with highest LoH counts
top_genes <- loh_gene_summary %>%
  group_by(gene_symbol) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  arrange(desc(total_count)) %>%
  head(10)

cat("\nTop 10 Genes Most Frequently Impacted by LoH:\n")
print(top_genes)

# Bar plot: Top 10 Genes Most Frequently Impacted by LoH
top_genes_plot <- ggplot(top_genes, aes(x = reorder(gene_symbol, -total_count), y = total_count)) +
  geom_bar(stat = "identity", fill = "darkgrey") +
  labs(
    title = "Top 10 Genes Most Frequently Impacted by LoH",
    x = "Gene Symbol",
    y = "Total LoH Events"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),        # Axis tick labels: size 12, black
    axis.text.x = element_text(angle = 45, hjust = 1),           # Rotate x-axis labels
    axis.title = element_text(size = 14, color = "black"),       # Axis titles: size 14, black
    axis.title.x = element_text(margin = margin(t = 10)),        # Increase spacing above x-axis title
    axis.title.y = element_text(margin = margin(r = 10)),        # Increase spacing to the right of y-axis title
    plot.title = element_text(size = 16, color = "black",        # Plot title: size 16, black
                              hjust = 0.5,                       # Centre-align the title
                              margin = margin(b = 10)),          # Add spacing below the title
    legend.position = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 1500),      # Set y-axis range from 0 to 1500
    breaks = seq(0, 1500, 500),  # Add breaks at intervals of 500
    labels = scales::comma      # Use normal number formatting with commas
  )
print(top_genes_plot)

# Save plot 
ggsave(filename = "Top 10 Genes Most Frequently Impacted by LoH.png", 
       plot = top_genes_plot, 
       width = 8, 
       height = 6, 
       dpi = 300,
       bg = "white")

print("Plot saved as Top 10 Genes Most Frequently Impacted by LoH.png")


# Stacked bar plot: Breakdown of LoH Types for Top 10 Genes
loh_gene_details_top10 <- loh_gene_summary %>%
  filter(gene_symbol %in% top_genes$gene_symbol)

loh_gene_details_top10_plot <- ggplot(loh_gene_details_top10, aes(x = reorder(gene_symbol, -count), y = count, fill = titan_call)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Breakdown of LoH Types for Top 10 Genes",
    x = "Gene Symbol",
    y = "LoH Event Count",
    fill = "LoH Type"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),        # Axis tick labels: size 12, black
    axis.text.x = element_text(angle = 45, hjust = 1),           # Rotate x-axis labels
    axis.title = element_text(size = 14, color = "black"),       # Axis titles: size 14, black
    axis.title.x = element_text(margin = margin(t = 10)),        # Increase spacing above x-axis title
    axis.title.y = element_text(margin = margin(r = 10)),        # Increase spacing to the right of y-axis title
    plot.title = element_text(size = 16, color = "black",        # Plot title: size 16, black
                              hjust = 0.5,                       # Centre-align the title
                              margin = margin(b = 10)),          # Add spacing below the title
  ) +
  scale_y_continuous(
    limits = c(0, 1500),      # Set y-axis range from 0 to 1500
    breaks = seq(0, 1500, 500),  # Add breaks at intervals of 500
    labels = scales::comma      # Use normal number formatting with commas
  )
print(loh_gene_details_top10_plot)

# Save plot 
ggsave(filename = "Breakdown of LoH Types for Top 10 Genes.png", 
       plot = loh_gene_details_top10_plot, 
       width = 8, 
       height = 6, 
       dpi = 300,
       bg = "white")

print("Plot saved as Breakdown of LoH Types for Top 10 Genes.png")

# Bar plot: Chromosome-wide LoH Events
loh_chrom_summary <- loh_segments_q2 %>%
  group_by(chrom, titan_call) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))

loh_chrom_summary_plot <- ggplot(loh_chrom_summary, aes(x = factor(chrom), y = count, fill = titan_call)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Chromosome-Wide Distribution of LoH Events",
    x = "Chromosome",
    y = "LoH Event Count",
    fill = "LoH Type"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),        # Axis tick labels: size 12, black
    axis.text.x = element_text(angle = 45, hjust = 1),           # Rotate x-axis labels
    axis.title = element_text(size = 14, color = "black"),       # Axis titles: size 14, black
    axis.title.x = element_text(margin = margin(t = 10)),        # Increase spacing above x-axis title
    axis.title.y = element_text(margin = margin(r = 10)),        # Increase spacing to the right of y-axis title
    plot.title = element_text(size = 16, color = "black",        # Plot title: size 16, black
                              hjust = 0.5,                       # Centre-align the title
                              margin = margin(b = 10)),  
  ) +
  scale_y_continuous(
    limits = c(0, 50000),      # Set y-axis range from 0 to 50,000
    breaks = seq(0, 50000, 10000),  # Add breaks at intervals of 10,000
    labels = scales::comma      # Use normal number formatting with commas
  )
print(loh_chrom_summary_plot)

# Save plot 
ggsave(filename = "Chromosome-Wide Distribution of LoH Events.png", 
       plot = loh_chrom_summary_plot, 
       width = 8, 
       height = 6, 
       dpi = 300,
       bg = "white")

print("Plot saved as Chromosome-Wide Distribution of LoH Events.png")

# -----------------------------------------------------------------------------------------------
# QUESTION 3: Are there LoH events unique to recurrent tumour samples?
# -----------------------------------------------------------------------------------------------

# Filter for LoH events and categorise samples
loh_segments_q3 <- titan_seg %>%
  filter(titan_call %in% c("DLOH", "NLOH", "ALOH")) %>%
  select(pair_barcode, chrom, start, end, titan_call) %>%
  mutate(sample_type = case_when(
    str_detect(pair_barcode, "-TP-") ~ "Primary",
    str_detect(pair_barcode, "-R1-") ~ "Recurrence 1",
    str_detect(pair_barcode, "-R2-") ~ "Recurrence 2",
    str_detect(pair_barcode, "-R3-") ~ "Recurrence 3",
    str_detect(pair_barcode, "-R4-") ~ "Recurrence 4",
    str_detect(pair_barcode, "-R5-") ~ "Recurrence 5",
    str_detect(pair_barcode, "-M1-") ~ "Metastatic",
    TRUE ~ "Other"
  ))

# Check for unknown barcodes
unknown_barcodes <- loh_segments_q3 %>%
  filter(sample_type == "Other") %>%
  distinct(pair_barcode)
cat("\nBarcodes Categorised as Unknown:\n")
print(unknown_barcodes)

# Create GRanges for LoH segments
loh_gr_q3 <- GRanges(seqnames = as.character(loh_segments_q3$chrom),
                     ranges = IRanges(start = loh_segments_q3$start, end = loh_segments_q3$end),
                     pair_barcode = loh_segments_q3$pair_barcode,
                     titan_call = loh_segments_q3$titan_call,
                     sample_type = loh_segments_q3$sample_type)

# Ensure chromosome names match between loh_gr_q3 and genes_gr
seqlevelsStyle(loh_gr_q3) <- seqlevelsStyle(genes_gr)

# Ensure genes_gr only contains chromosomes present in loh_gr_q3
matching_chromosomes <- intersect(seqlevels(loh_gr_q3), seqlevels(genes_gr))
genes_gr <- keepSeqlevels(genes_gr, matching_chromosomes, pruning.mode = "coarse")


# Find overlaps
overlap_q3 <- findOverlaps(loh_gr_q3, genes_gr)

# Extract LoH genes and their sample types
loh_genes_q3 <- data.frame(
  gene_symbol = mcols(genes_gr)$mcols[subjectHits(overlap_q3)],
  titan_call = mcols(loh_gr_q3)$titan_call[queryHits(overlap_q3)],
  sample_type = mcols(loh_gr_q3)$sample_type[queryHits(overlap_q3)],
  pair_barcode = mcols(loh_gr_q3)$pair_barcode[queryHits(overlap_q3)] 
)

# Verify extracted data
cat("\nExtracted LoH Genes Data:\n")
print(head(loh_genes_q3))

# Identify LoH events unique to recurrence samples
loh_unique_recurrence <- loh_genes_q3 %>%
  filter(sample_type != "Primary") %>%
  anti_join(
    loh_genes_q3 %>% filter(sample_type == "Primary"),
    by = c("gene_symbol", "titan_call")
  )

cat("\nLoH Events Unique to Recurrence Samples:\n")
print(head(loh_unique_recurrence))

# Summarise unique recurrent LoH events by LoH type
loh_recurrence_summary <- loh_unique_recurrence %>%
  group_by(titan_call) %>%
  summarise(
    total_genes = n_distinct(gene_symbol),
    total_events = n(),
    .groups = "drop"
  )

cat("\nSummary of LoH Events Unique to Recurrence Samples:\n")
print(loh_recurrence_summary)

# Bar plot: Total LoH Events by Sample Type and LoH Type
loh_summary_q3 <- loh_genes_q3 %>%
  group_by(sample_type, titan_call) %>%
  summarise(count = n(), .groups = "drop")

loh_summary_q3_plot <- ggplot(loh_summary_q3, aes(x = titan_call, y = count, fill = sample_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Distribution of LoH Events by Sample Type and LoH Type",
    x = "LoH Type",
    y = "Total Event Count",
    fill = "Sample Type"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),        # Axis tick labels: size 12, black
    axis.text.x = element_text(angle = 45, hjust = 1),           # Rotate x-axis labels
    axis.title = element_text(size = 14, color = "black"),       # Axis titles: size 14, black
    axis.title.x = element_text(margin = margin(t = 10)),        # Increase spacing above x-axis title
    axis.title.y = element_text(margin = margin(r = 10)),        # Increase spacing to the right of y-axis title
    plot.title = element_text(size = 16, color = "black",        # Plot title: size 16, black
                              hjust = 0.5,                       # Centre-align the title
                              margin = margin(b = 10)),
  ) +
  scale_y_continuous(
    limits = c(0, 1200000),      # Set y-axis range from 0 to 1,200,000
    breaks = seq(0, 1200000, 200000),  # Add breaks at intervals of 200,000
    labels = scales::comma      # Use normal number formatting with commas
  )
print(loh_summary_q3_plot)

# Save plot 
ggsave(filename = "Distribution of LoH Events by Sample Type and LoH Type.png", 
       plot = loh_summary_q3_plot, 
       width = 8, 
       height = 6, 
       dpi = 300,
       bg = "white")

print("Plot saved as Distribution of LoH Events by Sample Type and LoH Type.png")

# Count the number of unique samples (pair_barcodes) in each sample type
sample_count_summary <- loh_segments_q3 %>%
  group_by(sample_type) %>%
  summarise(total_samples = n_distinct(pair_barcode), .groups = "drop") %>%
  arrange(desc(total_samples))

cat("\nNumber of Samples per Sample Type:\n")
print(sample_count_summary)

# Bar plot : Visualise sample numbers
sample_count_summary_plot <- ggplot(sample_count_summary, aes(x = reorder(sample_type, -total_samples), y = total_samples, fill = sample_type)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Number of Samples per Sample Type",
    x = "Sample Type",
    y = "Total Number of Samples",
    fill = "Sample Type"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),        # Axis tick labels: size 12, black
    axis.text.x = element_text(angle = 45, hjust = 1),           # Rotate x-axis labels
    axis.title = element_text(size = 14, color = "black"),       # Axis titles: size 14, black
    axis.title.x = element_text(margin = margin(t = 10)),        # Increase spacing above x-axis title
    axis.title.y = element_text(margin = margin(r = 10)),        # Increase spacing to the right of y-axis title
    plot.title = element_text(size = 16, color = "black",        # Plot title: size 16, black
                              hjust = 0.5,                       # Centre-align the title
                              margin = margin(b = 10)),          # Add spacing below the title
  )
print(sample_count_summary_plot)

# Save plot 
ggsave(filename = "Number of Samples per Sample Type.png", 
       plot = sample_count_summary_plot, 
       width = 8, 
       height = 6, 
       dpi = 300,
       bg = "white")

print("Plot saved as Number of Samples per Sample Type.png")


# Bar plot: Proportion of LoH Events Unique to Recurrence
# Summarise the unique LoH events by type
loh_unique_recurrence_summary <- loh_unique_recurrence %>%
  group_by(titan_call) %>%
  summarise(unique_count = n(), .groups = "drop")

# Add a column for percentage
loh_unique_recurrence_summary <- loh_unique_recurrence_summary %>%
  mutate(percentage = (unique_count / sum(unique_count)) * 100)

# Use a better colour palette
palette <- brewer.pal(n = 8, "Pastel2")[1:nrow(loh_unique_recurrence_summary)]

# Create a more appealing pie chart
loh_unique_recurrence_pie <- ggplot(loh_unique_recurrence_summary, aes(x = "", y = unique_count, fill = titan_call)) +
  geom_bar(width = 1, stat = "identity", color = NA) + # Remove the inner borders
  coord_polar("y", start = 0) +
  labs(
    title = "Proportion of Unique LoH Events in Recurrent Tumours",
    fill = "LoH Type"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 10)),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(), # Remove gridlines
    panel.background = element_blank(), # Remove panel background
    plot.background = element_blank() # Remove plot background
  ) +
  geom_text(aes(label = paste0(round(unique_count / sum(unique_count) * 100, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 4)
print(loh_unique_recurrence_pie)

# Save plot 
ggsave(filename = "Proportion of Unique LoH Events in Recurrent Tumours.png", 
       plot = loh_unique_recurrence_pie, 
       width = 8, 
       height = 6, 
       dpi = 300,
       bg = "white")

print("Plot saved as Proportion of Unique LoH Events in Recurrent Tumours.png")

# Bar Plot: Number of unique events for each gene, grouped by sample_type.
gene_count <- loh_unique_recurrence %>%
  group_by(gene_symbol, sample_type) %>%
  summarise(count = n(), .groups = "drop")

gene_count_plot <- ggplot(gene_count, aes(x = reorder(gene_symbol, -count), y = count, fill = sample_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Unique LoH Events per Gene by Sample Type",
    x = "Gene Symbol",
    y = "Number of Unique Events",
    fill = "Sample Type"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, color = "black"),        # Axis tick labels: size 12, black
    axis.text.x = element_text(angle = 45, hjust = 1),           # Rotate x-axis labels
    axis.title = element_text(size = 14, color = "black"),       # Axis titles: size 14, black
    axis.title.x = element_text(margin = margin(t = 10)),        # Increase spacing above x-axis title
    axis.title.y = element_text(margin = margin(r = 10)),        # Increase spacing to the right of y-axis title
    plot.title = element_text(size = 16, color = "black",        # Plot title: size 16, black
                              hjust = 0.5,                       # Centre-align the title
                              margin = margin(b = 10)),          # Add spacing below the title
    
  ) +
  scale_y_continuous(
    limits = c(0, 14),      # Set y-axis range from 0 to 14
    breaks = seq(0, 14, 2),  # Add breaks at intervals of 2
    labels = scales::comma      # Use normal number formatting with commas
  )
print(gene_count_plot)

# Save plot 
ggsave(filename = "Unique LoH Events per Gene by Sample Type.png", 
       plot = gene_count_plot, 
       width = 8, 
       height = 6, 
       dpi = 300,
       bg = "white")

print("Plot saved as Unique LoH Events per Gene by Sample Type.png")

# -----------------------------------------------------------------------------------------------
# Debugging: Finish and Runtime Tracking
# -----------------------------------------------------------------------------------------------
end_time <- Sys.time()
cat("\n--- DEBUGGING FINISH ---\n")
cat("Script Runtime: ", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

# Save workspace to revisit results later
save.image(file = "workspace_debug.RData")
cat("Workspace Saved to 'workspace_debug.RData'.\n")

# Display memory usage
cat("\nMemory Usage:\n")
print(memory.size())

# -----------------------------------------------------------------------------------------------
# END OF SCRIPT
# -----------------------------------------------------------------------------------------------

