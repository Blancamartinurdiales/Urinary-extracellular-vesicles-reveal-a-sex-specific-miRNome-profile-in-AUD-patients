# Urinary extracellular vesicles reveal a sex-specific miRNome profile in AUD patients

## Overview

*This repository contains the code and supporting files used for the analyses presented in the study:*

> Urinary extracellular vesicles reveal a sex-specific miRNome profile in AUD patients

The project investigates sex-specific microRNA expression profiles in urinary extracellular vesicles from patients with Alcohol Use Disorder (AUD).
This repository contains a bioinformatics workflow for the processing and analysis of miRNA-seq data, with a focus on extracellular vesicle (EV) samples. The pipeline covers quality control, adapter removal, read alignment, quantification, differential expression analysis, target annotation, and functional enrichment.

--- 

## Workflow overview

Raw FASTQ │ ▼ Quality Control FastQC + MultiQC │ ▼ Adapter & Quality Trimming Cutadapt │ ▼ Post-trimming QC FastQC + MultiQC │ ▼ Lane merging │ ▼ Read Mapping Bowtie │ ├──► Human genome │ └──► Mature miRNAs (miRBase) │ ▼ SAM → BAM Samtools │ ▼ Alignment statistics Samtools flagstat │ ▼ miRNA count matrix │ ▼ Normalization & Differential Expression DESeq2 │ ▼ miRNA target annotation TarBase + miRTarBase │ ▼ Functional enrichment GO + KEGG + mdgsa

---

## Pipeline Architecture

1. Quality Control (QC) — FastQC + MultiQC

The first step is to assess the quality of the sequencing reads.

> Because miRNAs are very short (~18–24 nt), adapter sequences can represent a substantial proportion of each read. Quality control is therefore particularly important for identifying adapter contamination, low-quality bases, and abnormal read-length distributions.

 - Tools
  - **FastQC** - generates individual quality reports for each sample.
  - **MultiQC** - aggregates FastQC reports into a single summary report.
 - The following metrics are inspected:
  - Per-base sequence quality
  - Read length distribution
  - GC content
  - Adapter contamination
  - Overrepresented sequences
  - Overall quality across samples and sequencing lanes

---
 
2. Adapter Trimming — Cutadapt
 
miRNA-seq libraries frequently contain a large proportion of adapter sequence because mature miRNAs are much shorter than the sequencing reads.

Cutadapt is used to:

 - Remove sequencing adapters.
 - Filter low-quality reads.
 - Remove reads below a minimum length.
 - Set an appropriate maximum read length.
 
Filtered reads are stored in a *filter/* directory with the suffix:

*_filtered.fastq.gz*

The process is automated using a *cutadapt.sh* script and submitted to the cluster queue with:

*sbatch cutadapt.sh*

 ```text
 ADAPTER="TGGAATTCTCGGGTGCCAAGG"
QUALITY=30
MIN_LENGTH=17
MAX_LENGTH=35 
 
 ```
 
---
 
3. Post-trimming QC

FastQC and MultiQC are run again after trimming to verify:

 - Adapter removal.
 - Improvement in read quality.
 - Expected read-length distribution.
 - Consistency between samples and lanes.


---
 
4. Lane Merging

A single biological sample may be sequenced across multiple lanes. Although each lane represents a separate technical sequencing unit, all lanes belonging to the same sample should ultimately be treated as a single biological sample.

Example:

 ```text

cat sample_L001_R1_filtered.fastq \
    sample_L002_R1_filtered.fastq \
    sample_L003_R1_filtered.fastq \
    sample_L004_R1_filtered.fastq \
    > sample_filtered_combined.fastq
    
 ```    

The resulting combined FASTQ file is used for downstream alignment and quantification.

---
 
5. Mapping — Bowtie

miRNA reads are very short and therefore require alignment strategies suitable for short sequences.

Bowtie 1 is used for sensitive alignment of short miRNA reads, with the number of allowed mismatches explicitly controlled.

Depending on the analysis strategy, mapping can be performed in two stages:

 - **Human genome (GRCh38/hg38)**
Used to characterize and filter reads that may originate from other genomic regions or RNA species.
 - **Mature human miRNAs from miRBase**
Used to identify reads corresponding to known mature miRNA sequences.

For miRBase-based mapping, mature human miRNAs (hsa) are extracted from the mature.fa reference.

```text

bowtie-build human_mirna_mature.fa miRNA_index

```

 ```text
 for FILE in "$DATADIR"/*.fastq.gz; do
    BASENAME=$(basename "$FILE" _filtered.fastq.gz)
    bowtie2 --no-unal -p 40 -L 6 -i S,0,0.5 --ignore-quals --norc --score-min L,-1,-0.6 -D 20 \
    -x "$INDEX" \
    -U "$FILE" \
    -S "$OUTDIR/${BASENAME}.sam" \
    2> "$OUTDIR/log.${BASENAME}.txt"
done

 
 ```
Alignment jobs are automated using shell scripts and submitted to the cluster using sbatch.

---
 
6. SSAM/BAM Processing and Aligment QC

SAM files are converted to BAM format to reduce storage requirements and facilitate downstream processing.

```text
samtools view -b sample_aligned.sam > sample_aligned.bam
samtools sort sample_aligned.bam -o sample_sorted.bam
samtools index sample_sorted.bam
```

Alignment statistics are generated using:

```text
samtools flagstat sample_sorted.bam
```

The main metrics evaluated are:

 - Total number of reads.
 - Number of mapped reads.
 - Mapping percentage.
 - Mapping quality.
 - Ambiguous or multiple alignments.

A summary table can be generated across all samples:

| Sample | Total_Reads | Mapped_Reads | Mapping_Percentage |
|--------|-------------|--------------|--------------------|
| sample_1| 100000| 98000| 98.00 | 
| sample_2| 95000| 91000| 95.79 | 

This provides an overview of alignment performance and helps identify potential outlier samples.

---
 
6. miRNA Quantification and Count Matrix

Following alignment, reads are assigned to individual miRNAs to generate a count matrix.

Expected structure:

| sample_1 | sample_2 | sample_3| 
|----------|----------|---------|
| hsa-miR-1 | 120 | 150 | 98 | 
| hsa-miR-2 | 45 | 32 | 51 | 
| hsa-miR-3 | 800 | 920 | 760 | 


 - Rows: miRNAs
 - Columns: samples
 - Values: read counts

This count matrix is used as input for downstream statistical analysis with **DESeq2**.

---
 
7. Normalization and Differential Expression — DESeq2

Differential miRNA expression is analyzed in R using DESeq2.

The analysis includes:

 - Importing the miRNA count matrix.
 - Defining sample metadata.
 - Filtering low-count miRNAs.
 - Normalizing sequencing depth.
 - Fitting the differential expression model.
 - Performing condition-specific comparisons.
 - Exporting results as CSV files.
 - Experimental design

The analysis considers:

 - Condition: ALC vs CTROL
 - Sex: male vs female

The main model is:

> design = ~ condition + sex

This allows the effect of condition to be analyzed while adjusting for sex.

Differential expression comparisons

Three main comparisons are performed:

 - ALC_male vs CTROL_male
 - ALC_female vs CTROL_female
 - ALC vs CTROL

The main output statistics include:

 - log2 Fold Change (log2FC) — magnitude and direction of expression change.
 - p-value — statistical significance.
 - Adjusted p-value (padj) — significance corrected for multiple testing.

Results are exported as CSV files for downstream visualization and interpretation.

---
 
8. miRNA Target Annotation

Differentially expressed miRNAs are used to identify their potential target genes.

Curated miRNA–mRNA interaction databases are used:

 - TarBase
 - miRTarBase

These databases provide experimentally supported miRNA–target interactions, allowing differential miRNA expression to be linked to potential changes in gene regulation.

---
 
9. Functional Enrichment Analysis

The identified target genes are analyzed to determine which biological processes and signaling pathways may be associated with the differential miRNAs.

Tools and databases

 - Gene Ontology (GO)

  - Biological Processes (BP)
  - Molecular Functions (MF)
  - Cellular Components (CC)

 - KEGG
 - Reactome

---

## Repository structure

```text
.
├── 01Functions/      # Custom R functions
├── 02Scripts/        # Analysis scripts
├── Data/             # Input datasets
├── metadata/         # Sample and study metadata
├── Results/          # Generated results, figures and tables
├── miRNAuEVsAUD.Rproj
├── README.md
└── .gitignore

```

---



## Reproducibility

Analyses were performed in R.

-Usage.

 -Clone the repository.

 -Open miRNAuEVsAUD.Rproj in RStudio.
 
 -Run the scripts in 02Scripts/ following the analysis workflow.
 
 -Results will be generated in the Results/ directory.

---

## Data availability

The repository contains the data required to reproduce the analyses presented in the manuscript. Any restrictions on data sharing should comply with the corresponding ethical approvals and institutional regulations.

---

## Author

*Blanca Martín Urdiales*

