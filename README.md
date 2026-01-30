# CHIP-seq-analysis-Pipeline

A comprehensive, production-ready Nextflow pipeline for ChIP-seq (Chromatin Immunoprecipitation Sequencing) data analysis - from raw FASTQ (suitable for both singel read and pari read) files to annotated peaks and quality reports.

## Table of Contents
- <ins>Overview<ins>
- <ins>What is Chip-seq?<ins>
- <ins>Pipeline Workflow<ins>
- <ins>Installation<ins>


## Overview
This pipelie automates the complete CHIP-seq analysis workflow. It processes raw sequencing reads through quality control, alignment, peak calling and annotation.

## Key Features
  **Automated Quality Control:**  FastQC + MultiQC reports \
  **Flexible Execution:** Docker, Singularity, Conda, or native installation \
  **HPC Compatible:** SLURM and local execution \
  **Reproducible:** Containerized tools with version control \
  **Production Ready:** Error handling, retry logic, and detailed logging \
  **Comprehensive Output:** BAM files, peaks, bigWigs, and annotations

  ## What is ChIP-seq?
  **ChIP-seq (Chromatin Immunoprecipitation Sequencing)** is a technique used to identify where proteins (like transcription factor or histones) bind to DNA across the entire genome.

  **The Biology**
  1. **Crosslink:** Proteins are crosslinked to DNA in living cells
  2. **Fragment:** DNA is broken into small pieces
  3. **Immunoprecipitate:** Antibodies capture DNA bound to your protein on interest
  4. **Sequence:** The captured DNA fragements are sequenced
  5. **Analyze:** Computational analysis  identifies enriched regions (peaks)

## What you Get
- **Binding Sites:** Where protein binds in the genome
- **Peak Regions:** Genomic coordinates of enriched areas
- **Gene Associations:** Which genes are near binding sites
- **Motifs:** DNA sequence patterns at binding sites

  ## Pipeline Workflow
```
  RAW FASTQ FILES
      ↓
┌─────────────────────────────────────┐
│  STEP 1: Quality Control (FastQC)  │
└─────────────────────────────────────┘
      ↓
┌─────────────────────────────────────┐
│  STEP 2: Adapter Trimming           │
│  (Trim Galore)                      │
└─────────────────────────────────────┘
      ↓
┌─────────────────────────────────────┐
│  STEP 3: Genome Alignment (BWA)     │
└─────────────────────────────────────┘
      ↓
┌─────────────────────────────────────┐
│  STEP 4: Remove PCR Duplicates      │
│  (Picard MarkDuplicates)            │
└─────────────────────────────────────┘
      ↓
┌─────────────────────────────────────┐
│  STEP 5: Quality Metrics            │
│  (PhantomPeakQualTools)             │
└─────────────────────────────────────┘
      ↓
┌─────────────────────────────────────┐
│  STEP 6: Create bigWig Files        │
│  (deepTools bamCoverage)            │
└─────────────────────────────────────┘
      ↓
┌─────────────────────────────────────┐
│  STEP 7: Peak Calling (MACS2)       │
└─────────────────────────────────────┘
      ↓
┌─────────────────────────────────────┐
│  STEP 8: Annotate Peaks (HOMER)     │
└─────────────────────────────────────┘
      ↓
┌─────────────────────────────────────┐
│  STEP 9: Summary Report (MultiQC)   │
└─────────────────────────────────────┘
      ↓
RESULTS: Peaks, Annotations, QC Reports

```

## Detailed Steps Explanations

**Step 1: Quality Control (FastQC)**
 **What:** Checks the quality of your raw sequencing data
 **Why:** Identifies problems like low quality bases, adapter contamination, or failed sequencing
 **Output:** HTML reports with quality metrics

**Step 2: Adapter Trimming (Trim Galore)**

 **What:** Removes adapter sequences and low-quality bases
 **Why:** Adapters are artificial sequences added during library prep that aren't part of your genome
 **Output:** Cleaned FASTQ files

**Step 3: Alignment (BWA)**

 **What:** Maps your reads to the reference genome
 **Why:** Determines where each DNA fragment came from in the genome
 **Output:** BAM files (aligned reads)

**Step 4: Remove Duplicates (Picard)**

 **What:** Identifies and removes PCR duplicates
 **Why:** PCR creates artificial copies during library prep; removing them prevents false enrichment
 **Output:** Deduplicated BAM files

**Step 5: Quality Assessment (PhantomPeakQualTools)**

 **What:** Calculates ChIP-seq specific quality metrics
 **Why:** Assesses if your ChIP worked (NSC, RSC scores, fragment length)
 **Output:** Quality metrics and fragment length estimates

**Step 6: Coverage Tracks (bamCoverage)**

 **What:** Creates normalized genome-wide coverage files
 **Why:** Allows visualization in genome browsers (IGV, UCSC)
 **Output:** bigWig files

**Step 7: Peak Calling (MACS2)**

 **What:** Identifies regions significantly enriched in ChIP vs Input
 **Why:** These are your binding sites!
 **Output:** Peak coordinates (BED/narrowPeak format)

**Step 8: Peak Annotation (HOMER)**

 **What:** Associates peaks with nearby genes and genomic features
 **Why:** Tells you which genes your protein might regulate
 **Output:** Annotated peak tables

**Step 9: Summary Report (MultiQC)**

 **What:** Aggregates all QC metrics into one interactive report
 **Why:** Provides overview of entire experiment quality
 **Output:** HTML report with all metrics

 # Installation

 ### Prerequisites
 - *Nextflow* (≥ 21.10.0)
 - *Java* (≥ 8)
 - One of: *Docker*, *Singularity*, or *Conda*

   
   ### Install Nextflow
   
```
# Download Nextflow
curl -s https://get.nextflow.io | bash

# Make it executable
chmod +x nextflow

# Move to a location in your PATH
sudo mv nextflow /usr/local/bin/

# Test installation
nextflow -version

```
### Clone This Repository

```
git clone https://github.com/yourusername/chipseq-pipeline.git
cd chipseq-pipeline

```

### *Quick Start*
Create a CSV file (samples.csv) with your samples(single end or pair end):

```
chip_sample,chip_fastq_r1,chip_fastq_r2,input_sample,input_fastq_r1,input_fastq_r2
MOLM-14_DMS01_5_BRD4,data/fastq/MOLM-14_DMS01_5_BRD4.fastq.gz,,MOLM-14_DMS01_5_Input,data/fastq/MOLM-14_DMS01_5_Input.fastq.gz,
MOLM-14_DMS02_6_BRD4,data/fastq/MOLM-14_DMS02_6_BRD4.fastq.gz,,MOLM-14_DMS02_6_Input,data/fastq/MOLM-14_DMS02_6_Input.fastq.gz,

```
**Column Descriptions:**
chip_sample: Unique name for your ChIP sample\
chip_fastq_r1: Path to R1 FASTQ file (required)\
chip_fastq_r2: Path to R2 FASTQ file (leave empty for single-end)\
input_sample: Name for the matched Input/Control sample\
input_fastq_r1: Path to Input R1 FASTQ\
input_fastq_r2: Path to Input R2 FASTQ (leave empty for single-end)

### Download Refrence Genome
```
mkdir -p genome

# Download human genome (hg38)
wget -O genome/hg38.fa.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip genome/hg38.fa.gz

```
# Run the pipeline
## With Docker
```
nextflow run main.nf \
  --input samples.csv \
  --genome hg38 \
  --genome_fasta genome/hg38.fa \
  --outdir results \
  -profile docker`

```
## With Singularity (for HPC)

```
nextflow run main.nf \
  --input samples.csv \
  --genome hg38 \
  --genome_fasta genome/hg38.fa \
  --outdir results \
  -profile singularity

```
# Detaild Usage
Currently supported genomes: hg38 - Human (GRCh38), hg19 - Human (GRCh37), mm10 - Mouse (GRCm38), mm9 - Mouse (GRCm37)

# Peak Calling Options
## For Transcription Factors *(narrow peaks):*
```
nextflow run main.nf --input samples.csv --broad_peaks false

```
## For Histone Marks *(braod peaks):*
```
nextflow run main.nf --input samples.csv --broad_peaks true

```

# Advance Parameters
```
nextflow run main.nf \
  --input samples.csv \
  --genome hg38 \
  --genome_fasta genome/hg38.fa \
  --outdir results \
  --macs_qvalue 0.01 \              # Stricter peak calling (default: 0.05)
  --broad_peaks true \               # For histone marks
  --skip_fastqc false \              # Include FastQC
  --skip_trimming false \            # Include trimming
  -profile singularity \
  -resume                            # Resume from last successful step

```

# Output Files
The pipeline creates the following directory structure:

```
results/
├── 01_fastqc/                    # Raw read quality reports
│   ├── sample1_R1_fastqc.html
│   └── sample1_R1_fastqc.zip
│
├── 02_trimmed/                   # Trimmed reads
│   ├── sample1_trimmed.fq.gz
│   └── sample1_trimming_report.txt
│
├── 03_alignment/                 # Aligned reads
│   ├── sample1.sorted.bam
│   ├── sample1.sorted.bam.bai
│   └── sample1.flagstat.txt
│
├── 04_dedup/                     # Deduplicated reads
│   ├── sample1.dedup.bam
│   ├── sample1.dedup.bam.bai
│   ├── sample1.dup_metrics.txt
│   └── sample1.dedup.flagstat.txt
│
├── 05_qc/                        # ChIP-seq quality metrics
│   ├── sample1.spp.out
│   └── sample1.spp.pdf
│
├── 06_bigwig/                    # Genome browser tracks
│   └── sample1.bw
│
├── 07_peaks/                     # Called peaks
│   ├── sample1_peaks.narrowPeak
│   ├── sample1_peaks.xls
│   ├── sample1_summits.bed
│   └── sample1_model.r
│
├── 08_annotation/                # Annotated peaks
│   └── sample1.annotated.txt
│
├── 09_multiqc/                   # Summary report
│   ├── multiqc_report.html
│   └── multiqc_data/
│
└── pipeline_info/                # Pipeline execution info
    ├── timeline.html
    ├── report.html
    ├── trace.txt
    └── dag.svg

```










