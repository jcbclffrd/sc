# Single-Nucleus RNA-seq TE Analysis in Alzheimer's Disease

Analysis of transposable element (TE) expression in Alzheimer's disease using single-nucleus RNA-seq data from **Morabito et al. (2021)**.

---

## Project Overview

**Study**: Single-nuclei epigenomic and transcriptomic landscape in Alzheimer's disease  
**Paper**: Morabito et al. (2021) Nature Genetics, https://doi.org/10.1038/s41588-021-00894-z  
**GEO**: GSE175952  
**BioProject**: PRJNA729525

This project reprocesses the raw snRNA-seq data with a focus on detecting and quantifying transposable element expression at single-nucleus resolution.

---

## Dataset

### Samples Processed

| Assay Type | Samples | Downloaded | Aligned | Status |
|------------|---------|------------|---------|--------|
| **10x snRNA-seq** | 152 | 150 | 150 | ✅ Complete |
| **Bulk RNA-seq** | 191 | 154 | 154 | ✅ Complete |
| **ATAC-seq** | 32 | - | - | ⏸️ Not processed |

### Raw Data
- **10x snRNA-seq**: 150 samples, 208bp paired-end (Illumina NovaSeq 6000)
- **Chemistry**: 10x Chromium Single Cell 3' v3
- **Cell barcodes**: 16bp + 12bp UMI
- **Total raw cells**: 574,298 nuclei (before QC filtering)
- **Expected after QC**: ~192,000 high-quality nuclei (~33% pass rate)

---

## Reference Annotations

### Human Genome (hg38 / GRCh38)
- **File**: `genome/hg38.fa`
- **Size**: 3.0 GB
- **Source**: UCSC Genome Browser

### Gene Annotations (GENCODE v45)
- **File**: `annotations/gencode.v45.primary_assembly.annotation.gtf`
- **Release**: September 2023 (Ensembl 111)
- **Lines**: 3,428,060 annotation entries
- **Genes**: 63,241 total genes
  - 20,073 protein-coding genes
  - 19,370 lncRNAs
  - 10,145 processed pseudogenes
  - 2,217 misc_RNA
  - 1,910 snRNA
  - 1,879 miRNA
  - Others: snoRNA, rRNA, TEC, etc.
- **Features**: Exons, CDS, UTRs, transcripts, start/stop codons
- **Provider**: GENCODE / Ensembl
- **URL**: https://www.gencodegenes.org/

### Repeat Element Annotations (RepeatMasker)
- **File**: `annotations/hg38_rmsk.gtf`
- **Lines**: 5,683,690 repeat instances
- **Coverage**: ~45-50% of human genome
- **Source**: UCSC RepeatMasker track

**Repeat Element Classes:**
- **SINE**: 1,854,435 instances (includes ~1.1M Alu elements)
- **LINE**: 1,581,845 instances (L1, L2, L3 families)
- **LTR**: 747,112 instances (endogenous retroviruses)
- **DNA transposons**: 503,288 instances (hAT-Charlie, Tigger, MER, etc.)
- **Simple repeats**: 724,562 instances (microsatellites, tandem repeats)
- **Low complexity**: 106,053 instances (homopolymers)
- **Satellites**: 9,133 instances (centromeric/telomeric)
- **SVA**: 5,974 instances (SINE-VNTR-Alu composite)
- **Other**: snRNA, rRNA, tRNA repeats

**Note**: RepeatMasker includes ALL repetitive elements, not just transposable elements. For TE-specific analysis, we focus on SINE, LINE, LTR, DNA, and SVA classes.

### Barcode Whitelist
- **File**: `annotations/3M-february-2018.txt`
- **Barcodes**: 3,686,400 valid 16bp cell barcodes
- **Source**: 10x Genomics CellRanger 10.0.0
- **Chemistry**: 10x Chromium v3

---

## Pipeline

### 1. Alignment (STARsolo)
- **Tool**: STAR 2.7.11b with `--soloType CB_UMI_Simple`
- **Parameters**: 
  - Multi-mapping: 100 (increased from default 10 for TE detection)
  - Barcode read length: 100bp (28bp barcode/UMI region)
- **Output**: 436GB BAM files (150 samples)

### 2. TE Quantification (scTE)
- **Tool**: scTE 1.0 (single-cell TE quantification)
- **Index**: hg38_gencode45.exclusive.idx (611MB, 4.8M features)
- **Threads**: 10 (reduced from 20 to avoid memory issues)
- **Output**: 7.8GB .h5ad files (150 samples)

### 3. Data Aggregation
- **Merged dataset**: 574,298 cells × 62,651 features (genes + TEs)
- **Total UMI counts**: 2.3 billion
- **Features detected**: 56,316 (90% of all features)
- **File**: `analysis/merged_samples.h5ad` (7.5GB)

---

## Results Summary

### TE Expression
- **TE features detected**: 4,763
- **Total TE counts**: 1.44 billion UMIs (62.7% of all reads)
- **Top TE families**: Alu (dominant), L2, MIR, L1

### Top 15 Expressed TEs
1. AluJb - 69.6M counts
2. AluY - 68.8M counts
3. AluSx1 - 68.1M counts
4. AluSx - 65.5M counts
5. AluSz - 59.5M counts
6. L2a - 48.1M counts
7. AluJr - 42.9M counts
8. AluJo - 38.3M counts
9. AluSq2 - 35.6M counts
10. MIRb - 33.3M counts
11. L2c - 32.0M counts
12. AluSp - 31.6M counts
13. AluSz6 - 27.4M counts
14. MIR - 26.9M counts
15. AluSg - 21.3M counts

### Top Expressed Genes (non-TEs)
1. MALAT1 (lncRNA) - 73.2M counts
2. TALAM1 - 65.4M counts
3. NEAT1 (lncRNA) - 3.4M counts
4. PLP1 (oligodendrocyte marker) - 2.7M counts
5. MEG3 (imprinted gene) - 2.7M counts

---

## Directory Structure

```
sc/
├── sra_downloads/          # 1.2TB FASTQ files (150 10x + 154 bulk samples)
├── starsolo_aligned/       # 436GB BAM files (150 10x samples)
├── star_bulk_aligned/      # 424GB BAM files (154 bulk samples)
├── scTE_output/            # 7.8GB .h5ad files (150 10x samples)
├── analysis/               # Merged datasets and results
│   ├── merged_samples.h5ad (7.5GB)
│   └── combinedCells.csv (3.4MB)
├── genome/                 # hg38 reference genome
├── annotations/            # Gene, TE, and barcode annotations
├── star_index/             # STAR genome index (28GB)
├── scripts/                # Analysis scripts
├── docs/                   # Documentation
└── data/                   # Sample lists and metadata
```

---

## Storage Requirements

- **Total project size**: ~2.1TB
- **FASTQ files**: 1.2TB
- **BAM files**: 860GB (436GB + 424GB)
- **Analysis files**: 15GB
- **Reference data**: 33GB

**Current disk usage**: 3.1TB / 3.7TB (88%), 462GB available

---

## Key Scripts

### Data Processing
- `scripts/download_10x_snRNA_samples.sh` - Download 10x samples from SRA
- `scripts/run_starsolo_alignment.sh` - STARsolo alignment pipeline
- `scripts/run_star_bulk_alignment.sh` - Bulk RNA-seq alignment
- `scripts/build_scte_index.sh` - Build scTE exclusive index
- `scripts/run_scte_quantification.sh` - scTE TE quantification
- `scripts/merge_and_aggregate_samples.py` - Merge all samples

### Monitoring
- `scripts/check_10x_download.sh` - Monitor download progress
- `scripts/check_alignment.sh` - Monitor 10x alignment
- `scripts/check_bulk_alignment.sh` - Monitor bulk alignment
- `scripts/monitor_scte.sh` - Monitor scTE progress

---

## Issues Resolved

1. **STAR Index Version Mismatch** - Rebuilt index with STAR 2.7.11b
2. **Barcode Length Issue** - Added `--soloBarcodeReadLength 100`
3. **scTE Path Handling Bug** - Used absolute paths and subshell execution
4. **Memory Constraints** - Reduced scTE threads from 20 to 10

---

## Next Steps

1. ✅ Download GEO sample metadata (GSE175952)
2. ✅ Map SRR → GSM → Condition (AD/Control)
3. ⏳ Apply QC filters (200-10,000 genes, <10% MT reads)
4. ⏳ Cell type clustering and annotation
5. ⏳ Differential TE expression analysis (AD vs Control)
6. ⏳ Visualization and interpretation

---

## Technical Notes

### Understanding "Samples" and Cell Counts

**What is a "Sample"?**
- **1 Sample = 1 SRR accession = 1 patient = 1 brain tissue piece**
- Each sample consists of:
  - 2 FASTQ files (R1: 28bp barcode+UMI, R2: 100bp RNA-seq read)
  - 1 10x Chromium capture (one microfluidic chip run)
  - ~3,000-5,000 cells from that capture

**Dataset Summary:**
- **Total available**: 152 samples (152 different patients)
- **We processed**: 150 samples = 150 patients × 2 FASTQ files = 300 FASTQ files
- **Morabito analyzed**: 18 samples (18 patients from the 152 available)
  - Sample-17, 19, 22, 27, 33, 37, 43, 45, 46, 47, 50, 52, 58, 66, 82, 90, 96, 100

**Cell Count Analysis:**

| Dataset | Samples | Total Cells | Cells/Sample | Mean UMI/Cell |
|---------|---------|-------------|--------------|---------------|
| **Our data** | 150 | 574,298 | 3,829 | 3,985 |
| **Morabito (GEO: GSE174367)** | 18 | 61,770 | 3,432 | 9,681 |

**Key Finding**: The cells per sample are very similar (3,829 vs 3,432)! The 9.3x difference in total cell count is primarily due to **sample selection**, not filtering stringency.

**Why the Discrepancy?**
1. **Sample selection (18/150 = 12%)** ← Main factor! They used only 18 of 152 available samples
2. **Slightly stricter filtering** (~10% fewer cells per sample)
3. **Higher quality threshold** (their mean UMI: 9,681 vs ours: 3,985)
   - They likely applied higher UMI cutoffs during CellRanger processing
   - Our data includes more low-quality cells from STARsolo's liberal filtering

**Paper's Reported Numbers:**
- **191,890 total nuclei** = ~130,000 ATAC-seq + ~61,000 snRNA-seq (not all RNA!)
- Their snRNA-seq subset: 61,770 cells from 18 samples
- Cell type distribution: ODC (60%), EX (10%), INH (10%), ASC (8%), MG (7%), OPC (5%), PER.END (1%)

### Quality Control Strategy

**CellRanger vs STARsolo Filtering:**
- **CellRanger** (Morabito's method): Conservative cell calling, higher quality thresholds
- **STARsolo** (our method): More liberal, includes lower-quality cells
- **scTE**: Reads from BAM directly, bypasses STARsolo's `filtered/` matrix
  - Result: We get ~38% more cells per sample than STARsolo's filtered count
  - Example: SRR14513977 has 3,345 filtered cells but scTE kept 4,629 cells

**Standard QC Filters** (to match paper):
- 200-10,000 genes per cell
- <10% mitochondrial reads
- Remove doublets and dying cells
- UMI threshold: ~7,800+ UMIs to match their median of 6,374

**Our Current Status:**
- **Raw cells**: 574,298 (all samples, minimal filtering)
- **After STARsolo filter**: ~500,000 cells (using `filtered/` barcodes)
- **After UMI threshold (≥7,847)**: ~71,685 cells
- **Target**: Match their 18-sample subset quality for comparison

Expected outcome: **We could match their cell quality by applying stricter UMI thresholds**

### TE Detection Strategy
- **Multi-mapping reads**: Set to 100 (vs default 10) to capture reads mapping to multiple TE copies
- **Exclusive index**: scTE builds index that separates overlapping genes/TEs
- **UMI deduplication**: Counts unique molecules, not raw reads
- **Result**: 62.7% of UMIs map to TEs (expected for repetitive genome regions)

---

## Citation

If you use this analysis, please cite the original paper:

**Morabito S, Miyoshi E, Michael N, et al.** Single-nucleus chromatin accessibility and transcriptomic characterization of Alzheimer's disease. *Nature Genetics*. 2021;53(8):1143-1155. doi:10.1038/s41588-021-00894-z

---

## Software Versions

- STAR: 2.7.11b
- scTE: 1.0
- scanpy: 1.10.0
- Python: 3.12
- samtools: 1.18
- SRA Toolkit: 3.0.3

---

## Contact

Project repository: `/home/jacobc/sc`  
Date: December 11, 2025
