# PRJNA729525 10x snRNA-seq TE Analysis - Project Status

**Date**: December 11, 2025  
**Project**: Single-nucleus TE expression analysis in Alzheimer's disease  
**Paper**: Morabito et al. (2021) Nature Genetics

---

## Dataset Classification

We identified **3 different assay types** in PRJNA729525:

| Assay Type | Sample Count | Read Length | Use Case |
|------------|--------------|-------------|----------|
| **10x snRNA-seq** | **152** | **208bp** | **PRIMARY - TE analysis** ✅ |
| Bulk RNA-seq | 191 | 200bp | Future comparison |
| ATAC-seq | 32 | varies | Chromatin accessibility |

**Total**: 375 samples

---

## Current Status

### ✅ Completed (Steps 1-7)

1. **Directory structure** created and organized
   - scripts/, docs/, data/, analysis/ directories
   - Git repository initialized with proper .gitignore

2. **Reference files** copied from hcaTE:
   - hg38.fa genome (3.0GB)
   - gencode.v45 annotations (1.5GB)
   - hg38_rmsk.gtf (TE annotations, 1.1GB)
   - 10x v3 whitelist: 3M-february-2018.txt (3.7M barcodes, 60MB)

3. **Sample identification and download**: 
   - Downloaded 150/152 10x snRNA-seq samples (1.2TB FASTQ files)
   - 2 samples missing from SRA

4. **STAR index** rebuilt with version 2.7.11b (28GB)
   - Fixed compatibility issue with STAR 2.7.4a index

5. **scTE index** built successfully
   - hg38_gencode45.exclusive.idx (611MB, 4.8M features)
   - Completed: Dec 9, 02:24:45

6. **STARsolo alignment** - ALL 150 SAMPLES COMPLETE ✅
   - Fixed barcode length issue (--soloBarcodeReadLength 100)
   - Multi-mapping set to 100 (vs default 10) for comprehensive TE detection
   - Output: 436GB BAM files in starsolo_aligned/
   - Completed: Dec 9, 00:09:18

7. **scTE quantification** - ALL 150 SAMPLES COMPLETE ✅
   - Generated .h5ad files with gene and TE counts per cell
   - Used 10 threads (reduced from 20 to avoid memory issues)
   - Output: 7.8GB of .h5ad files in scTE_output/
   - Completed: Dec 9, 17:37:52

8. **Data aggregation** - COMPLETE ✅
   - Merged all 150 samples into unified dataset
   - 574,298 cells × 62,651 features (genes + TEs)
   - Output: merged_samples.h5ad (7.5GB), combinedCells.csv (3.4MB)
   - Completed: Dec 10, 04:55

9. **Bulk RNA-seq alignment** - COMPLETE ✅
   - Downloaded 154/191 bulk RNA-seq samples
   - All 154 samples aligned with STAR (424GB BAM files)
   - Multi-mapping set to 100 for TE detection
   - Output: star_bulk_aligned/
   - Completed: Dec 11, 14:03:36

### 🔄 In Progress
- **Differential TE Analysis**
  - Need to add sample metadata (AD vs Control labels)
  - Need normalization (TPM/CPM) for differential expression
  - Quality control and filtering
  - Cell type clustering and annotation
  - Differential TE expression analysis

### 📋 To Do
1. Add sample metadata from GEO (GSE175952)
2. Implement normalization (consider gene length bias)
3. Quality control filtering (cells and features)
4. Differential TE expression analysis (AD vs Control)
5. Visualization and interpretation

---

## Key Files

**Sample Lists:**
- `data/snRNA_seq_10x_samples.txt` - 152 10x snRNA-seq samples (PRIMARY)
- `data/bulk_rnaseq_samples.txt` - 191 bulk RNA-seq samples
- `data/sra_metadata.csv` - Full metadata for all 375 samples

**Scripts (in scripts/ directory):**
- `download_10x_snRNA_samples.sh` - Download 10x samples
- `check_10x_download.sh` - Monitor download progress
- `run_starsolo_alignment.sh` - STARsolo alignment pipeline
- `run_star_bulk_alignment.sh` - Bulk RNA-seq STAR alignment
- `check_bulk_alignment.sh` - Monitor bulk alignment progress
- `build_scte_index.sh` - Build scTE exclusive index
- `run_scte_quantification.sh` - scTE TE quantification
- `merge_and_aggregate_samples.py` - Merge all samples and create aggregated counts
- `check_alignment.sh` - Monitor alignment progress
- `monitor_scte.sh` - Monitor scTE progress

**Reference Data:**
- `genome/hg38.fa` - Human genome reference (3.0GB)
- `annotations/gencode.v45.primary_assembly.annotation.gtf` - Gene annotations (1.5GB)
- `annotations/hg38_rmsk.gtf` - TE annotations (1.1GB)
- `annotations/3M-february-2018.txt` - 10x v3 whitelist (3.7M barcodes, 60MB)
- `star_index/` - STAR genome index (28GB, built with STAR 2.7.11b)
- `annotations/hg38_gencode45.exclusive.idx` - scTE index (611MB)

**Data Outputs:**
- `sra_downloads/` - 1.2TB FASTQ files (150 10x + 154 bulk samples)
- `starsolo_aligned/` - 436GB BAM files (150 10x samples)
- `star_bulk_aligned/` - 424GB BAM files (154 bulk samples)
- `scTE_output/` - 7.8GB .h5ad files (150 10x samples)
- `analysis/merged_samples.h5ad` - Unified snRNA-seq dataset (7.5GB, 574,298 cells × 62,651 features)
- `analysis/combinedCells.csv` - Aggregated count table (3.4MB)

---

## 10x Chromium v3 Chemistry Details

From the paper:
- **Platform**: 10x Chromium Single Cell 3' v3
- **Sequencing**: Illumina NovaSeq 6000, 100bp paired-end
- **Read structure**: 
  - 28bp: 16bp cell barcode + 12bp UMI
  - 91bp: cDNA sequence
  - Total: 208bp in SRA

**Required for STARsolo:**
- Barcode whitelist: `3M-february-2018.txt` (10x v3)
- Cell barcode: 16bp (start: 1, length: 16)
- UMI: 12bp (start: 17, length: 12)

---

## Storage Status

**Actual Usage:**
- **FASTQ files**: 1.2TB (sra_downloads/)
- **10x BAM files**: 436GB (starsolo_aligned/)
- **Bulk BAM files**: 424GB (star_bulk_aligned/)
- **scTE output**: 7.8GB (scTE_output/)
- **Analysis files**: 7.5GB (analysis/)
- **Reference data**: ~33GB (genome + annotations + indices)
- **Total project**: ~2.1TB

**Disk Space:**
- Total: 3.7TB
- Used: 3.1TB (88%)
- Available: 462GB

---

## Data Summary

**Merged Dataset:**
- **Total cells**: 574,298
- **Total features**: 62,651 (genes + TEs)
- **Total UMI counts**: 2,288,375,201 (~2.3 billion)
- **Features expressed**: 56,316 (90%)
- **Data format**: UMI-deduplicated molecular counts (not raw reads)

**Top Expressed Features:**
1. MALAT1 (lncRNA) - 73M counts
2. AluJb (TE) - 69M counts
3. AluY (TE) - 68M counts
4. AluSx1 (TE) - 68M counts
5. AluSx (TE) - 65M counts

---

## Issues Resolved

1. **STAR Index Version Mismatch** (Dec 8)
   - Problem: Index built with STAR 2.7.4a, running 2.7.11b
   - Solution: Rebuilt index with matching version

2. **Barcode Length Issue** (Dec 8)
   - Problem: 100bp barcode reads vs expected 28bp
   - Solution: Added --soloBarcodeReadLength 100

3. **scTE Path Handling Bug** (Dec 9)
   - Problem: scTE duplicates output path in internal file paths
   - Solution: cd into output directory and use absolute paths

4. **scTE Memory Issues** (Dec 9)
   - Problem: 20 threads caused 115GB/120GB memory usage
   - Solution: Reduced to 10 threads

---

## Next Steps for Analysis

1. **Add sample metadata**
   - Download GEO metadata (GSE175952)
   - Identify AD vs Control samples
   - Add condition labels to merged_samples.h5ad

2. **Normalization**
   - Implement TPM/CPM normalization
   - Consider gene length bias (longer transcripts capture more UMIs)
   - Standard scanpy workflow: normalize_total() + log1p()

3. **Quality Control**
   - Filter low-quality cells (min genes, max mitochondrial %)
   - Filter low-count features
   - Generate QC plots

4. **Differential Expression**
   - Cell type clustering and annotation
   - Compare TE expression: AD vs Control
   - Focus on specific TE families (Alu, LINE, etc.)

5. **Visualization**
   - UMAP colored by condition and cell type
   - Volcano plots for differential TEs
   - Heatmaps of top TEs

---

## Notes

- All 150 samples have identical 62,651 features in same order (can be directly concatenated)
- UMI counts are molecular counts, not raw read counts
- Multimapping set to 100 (vs default 10) for comprehensive TE detection
- TE detection is working: Alu elements among top expressed features
- Gene length normalization needed before differential expression analysis
- Git repository tracks only scripts/docs (~14MB), excludes data (1.7TB)
