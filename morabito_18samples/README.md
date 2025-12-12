# Morabito 18-Sample Subset Analysis

This folder contains the analysis of the **18 samples** that Morabito et al. (2021) used in their published study, extracted from our complete 150-sample dataset.

## Overview

**Goal**: Compare our scTE TE quantification results with Morabito's published analysis using the exact same samples and cell barcodes.

**Samples**: 18 patients (Sample-17, 19, 22, 27, 33, 37, 43, 45, 46, 47, 50, 52, 58, 66, 82, 90, 96, 100)
- **AD samples**: 12 patients
- **Control samples**: 6 patients

**Cells**: 61,770 high-quality nuclei (filtered by CellRanger/STARsolo)

### Why Only 18 of 150+ Samples?

Morabito et al. processed **152 total samples** but published only **18 samples** (GSE174367). The selection was likely based on:

1. **Quality Control**: Not all samples passed stringent QC metrics
   - Cell yield (samples with < 1,000 cells excluded)
   - Sequencing depth and gene detection rate
   - Doublet rate and mitochondrial content
   
2. **Batch Effects**: Selected samples from consistent sequencing batches to minimize technical variation

3. **Balanced Design**: Matched samples for:
   - Age distribution between AD and controls
   - Sex ratio (male/female balance)
   - Post-mortem interval (PMI)
   - Tissue quality and dissection consistency

4. **Statistical Power**: 18 high-quality samples with ~60,000 cells provides excellent power without adding noise from lower-quality samples

### Empty Droplets and Cell Filtering

**Raw 10x Chromium data** captures ~100,000 droplets per sample, but most are empty:
- **~98.7% of droplets are empty or low-quality** (< 100 UMIs)
- Only **~1-2% contain real single cells** (> 1,000 UMIs)
- For 18 samples: ~6.7M raw droplets → ~60K high-quality cells

**Processing differences:**
- **STARsolo/CellRanger**: Generates both `raw/` (all droplets) and `filtered/` (high-quality cells only) outputs
- **scTE**: Originally processed `raw/` outputs from STARsolo BAM files, including all ~500K droplets across samples
- **This analysis**: Filtered scTE results to match Morabito's exact cell barcodes from the `filtered/` outputs (61,770 cells)

The massive reduction (500K → 60K cells) is **expected and correct** - it represents removing empty droplets and low-quality cells to keep only real, high-quality single nuclei.

## Files

### Reference Data (from GEO: GSE174367)
- `GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5` - CellRanger filtered gene expression (261 MB)
- `GSE174367_snRNA-seq_cell_meta.csv.gz` - Cell metadata with annotations (425 KB)
- `sample_mapping.csv` - Mapping of Sample-X to SRR IDs

### Analysis Scripts
- `01_extract_18samples.py` - Extract 18 samples from our scTE data
- `02_filter_to_morabito_cells.py` - Filter to Morabito's exact cell barcodes
- `03_merge_and_qc.py` - Merge samples and perform QC
- `04_differential_analysis.py` - Differential TE/gene expression (AD vs Control)

### Output Data
- `scte_18samples/` - scTE .h5ad files for 18 samples only
- `merged_18samples.h5ad` - Combined dataset (genes + TEs, 61,770 cells)
- `merged_genes_only.h5ad` - Gene expression only (for comparison with Morabito)
- `differential_results/` - DE analysis results (TEs and genes)

## Workflow

### 1. Map Sample Names to SRR IDs
```bash
python 01_map_samples.py
```
Creates `sample_mapping.csv` with Sample-X → SRR mapping.

### 2. Extract 18 Samples
```bash
python 02_extract_18samples.py
```
Copies scTE output for the 18 samples to `scte_18samples/`.

### 3. Filter to Morabito's Cells
```bash
python 03_filter_to_morabito_cells.py
```
Uses Morabito's cell barcodes to filter our data to the exact same cells.

### 4. Merge and QC
```bash
python 04_merge_and_qc.py
```
Merges all 18 samples, adds metadata, performs QC.

### 5. Differential Analysis
```bash
python 05_differential_analysis.py
```
Performs DE analysis for:
- TEs: AD vs Control
- Genes: AD vs Control
- Cell-type specific analysis

## Analysis Goals

1. **Compare TE expression** between our scTE quantification and Morabito's gene-only data
2. **Identify AD-associated TEs** that are differentially expressed
3. **Correlate TE and gene expression** to understand regulatory relationships
4. **Cell-type specific TE analysis** (ODC, EX, INH, ASC, MG, OPC, PER.END)

## Expected Results

- **Differential TEs**: TEs upregulated/downregulated in AD
- **Differential genes**: Genes upregulated/downregulated in AD (compare with paper)
- **TE-gene correlations**: Co-expressed TEs and genes
- **Cell type specificity**: Which TEs are cell-type specific?

## Notes

- Morabito's data only contains **gene expression** (no TEs)
- Our scTE data contains **genes + TEs** from the same samples
- We can directly compare gene expression, and uniquely analyze TE expression
- Cell type labels are from Morabito's clustering/annotation

## Directory Structure

```
morabito_18samples/
├── README.md
├── 01_map_samples.py
├── 02_extract_18samples.py
├── 03_filter_to_morabito_cells.py
├── 04_merge_and_qc.py
├── 05_differential_analysis.py
├── sample_mapping.csv
├── scte_18samples/
│   ├── SRR14513977/
│   ├── SRR14513978/
│   └── ... (18 samples)
├── merged_18samples.h5ad
├── merged_genes_only.h5ad
└── differential_results/
    ├── TEs_AD_vs_Control.csv
    ├── Genes_AD_vs_Control.csv
    └── figures/
```
