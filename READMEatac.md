# ATAC-seq Single-Cell TE Quantification Pipeline

Brief guide for processing ATAC-seq data from raw FASTQ to single-cell BAM with cell barcodes and multi-mapper support for TE quantification.

## Pipeline Overview

```
FASTQ Download → CellRanger ATAC (QC + CB assignment) → STAR (multi-mapper alignment) → BAM-to-BAM CB Transfer → scTE/scTEATAC
```

## Required Scripts (in order)

### 1. Download ATAC-seq Data
**Script:** `scripts/downloadATACseq_prefetch.sh`
- Downloads all 32 ATAC-seq samples from SRA
- Uses two-step method: `prefetch` + `fasterq-dump --include-technical`
- Creates 4 files per sample: `_1.fastq.gz` (R1), `_2.fastq.gz` (R2), `_3.fastq.gz` (barcode), `_4.fastq.gz` (index)
- Output: `~/sc/sra_downloads/ATAC-seq/SRR14514*.fastq.gz`

**Critical:** Must use `--include-technical` flag to get barcode file (_3.fastq.gz)

### 2. Run CellRanger ATAC (Cell QC + Barcode Assignment)
**Script:** `scripts/process_cellranger_continuous.sh`
- Continuously monitors for completed downloads and processes samples
- Runs CellRanger ATAC 2.2.0 for cell calling and barcode assignment
- Output: `~/sc/cellranger_atac_output/SRR*/outs/possorted_bam.bam` (has CB:Z: tags)

**Purpose:** CellRanger performs cell QC and assigns barcodes, but only keeps "best hit" per read (no multi-mapper support for TEs)

### 3. STAR Multi-mapper Alignment
**Script:** `scripts/align_atacseq_TEs.sh`
- Aligns with STAR using `--outFilterMultimapNmax 100` for TE quantification
- Keeps up to 100 alignments per read (essential for repetitive TEs)
- Output: `~/sc/atacseq_aligned/SRR*/Aligned.sortedByCoord.out.bam` (~320M reads, 9.4GB)

**Purpose:** STAR provides multi-mapper support that CellRanger lacks

### 4. BAM-to-BAM CB Transfer
**Script:** `scripts/transfer_cb_tags.sh`
- Transfers cell barcode tags (CB:Z:) from CellRanger BAM to STAR multi-mapper BAM
- Matches reads by read ID (98.7% match rate)
- Output: `~/sc/atacseq_aligned/SRR*/star_multimapper_final.bam` (9.0GB)

**Helper:** `scripts/add_cb_tags.py` - Python script that performs read ID matching

**Result:** Final BAM has both:
- ✓ Multi-mapper support (from STAR)
- ✓ Cell barcode tags (from CellRanger)
- ✓ Ready for scTE/scTEATAC single-cell TE quantification

## Quick Start

```bash
# 1. Download data (or use existing downloads)
./scripts/downloadATACseq_prefetch.sh

# 2. Start continuous CellRanger processing (runs in background)
cd ~/sc && nohup ./scripts/process_cellranger_continuous.sh > logs_atacseq/continuous_main.log 2>&1 &

# 3. For a specific sample, run STAR alignment
./scripts/align_atacseq_TEs.sh SRR14514130

# 4. Transfer CB tags to create final BAM
./scripts/transfer_cb_tags.sh SRR14514130

# Output: ~/sc/atacseq_aligned/SRR14514130/star_multimapper_final.bam
```

## Key Files

- **Genome:** `~/sc/genome/GRCh38.primary_assembly.genome.fa`
- **STAR Index:** `~/sc/star_index/` (built with sjdbOverhang 49 for 50bp reads)
- **TE Annotation:** `~/sc/annotations/hg38_rmsk_tetranscripts.gtf` (TEtranscripts format)
- **Gene Annotation:** `~/sc/annotations/gencode.v45.primary_assembly.annotation.gtf`

## Performance

- **SRR14514130 Stats:**
  - Input: 183M read pairs (732M total reads with technical)
  - STAR output: 320M aligned reads
  - CB transfer: 316M reads matched (98.7%)
  - Final BAM: 9.0GB

- **Hardware:** AWS EC2 with 16,000 IOPS / 1,000 MB/s EBS volume
- **Time:** ~3 hours per sample for CellRanger, ~1 hour for STAR + CB transfer

## Next Steps

Run scTE or scTEATAC on `star_multimapper_final.bam` for single-cell TE quantification.
