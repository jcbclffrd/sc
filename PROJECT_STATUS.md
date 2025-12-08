# PRJNA729525 10x snRNA-seq TE Analysis - Project Status

**Date**: December 5, 2025  
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

### ✅ Completed
1. **Directory structure** created
2. **Reference files** copied from hcaTE:
   - hg38.fa genome
   - gencode.v45 annotations
   - hg38_rmsk.gtf (TE annotations)
   - 10x v2 whitelist (737K-august-2016.txt)
   - Need v3 whitelist for these samples

3. **Sample identification**: Filtered 152 10x snRNA-seq samples
4. **scTE installed** in hcaTE virtual environment

### 🔄 In Progress
- **Downloading 152 10x snRNA-seq samples** (208bp reads)
  - Monitor: `./check_10x_download.sh`
  - Log: `download_10x.log`
  - Individual logs: `logs_10x/`

### 📋 To Do
1. Complete 10x sample downloads (~152 samples)
2. Link STAR index
3. Build scTE exclusive index
4. Run STARsolo alignment (10x v3 chemistry)
5. Run scTE quantification
6. Differential TE analysis

---

## Key Files

**Sample Lists:**
- `snRNA_seq_10x_samples.txt` - 152 10x snRNA-seq samples (PRIMARY)
- `sra_metadata.csv` - Full metadata for all 375 samples

**Scripts:**
- `download_10x_snRNA_samples.sh` - Download 10x samples
- `check_10x_download.sh` - Monitor download progress

**Reference Data:**
- `genome/hg38.fa` - Human genome reference
- `annotations/gencode.v45.primary_assembly.annotation.gtf` - Gene annotations
- `annotations/hg38_rmsk.gtf` - TE annotations
- `annotations/737K-august-2016.txt` - 10x v2 whitelist (need v3!)

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

## Storage Estimates

- **10x samples**: 152 × ~3.5GB = **~530GB**
- **STAR alignments**: ~1-1.5TB (BAM files)
- **scTE output**: ~50-100GB (sparse matrices)
- **Total needed**: ~2TB

**Current usage**: 567GB (includes some bulk RNA-seq samples to clean up later)

---

## Next Steps When Downloads Complete

1. **Verify a downloaded sample has correct structure**
   ```bash
   zcat sra_downloads/SRR14513977/SRR14513977_1.fastq.gz | head -8
   # Should see 28bp reads with barcodes/UMIs
   ```

2. **Get 10x v3 whitelist**
   ```bash
   cp /path/to/3M-february-2018.txt annotations/
   ```

3. **Link STAR index**
   ```bash
   ln -s /home/jacobc/hcaTE/star_index star_index
   ```

4. **Build scTE index**
   ```bash
   source /home/jacobc/hcaTE/.venv/bin/activate
   scTE_build -g annotations/gencode.v45.primary_assembly.annotation.gtf \
              -r annotations/hg38_rmsk.gtf \
              -o annotations/hg38.exclusive \
              -build exclusive
   ```

---

## Notes

- Bulk RNA-seq (191 samples) can be analyzed later for comparison
- Some samples failed in initial download (154/375) but those were mostly bulk RNA-seq
- The 10x data will give us single-nucleus resolution for TE expression analysis
