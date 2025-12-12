# SoloTE Analysis: Morabito 18 Samples

Alternative TE quantification using SoloTE for comparison with scTE results.

## Overview

This folder contains an alternative analysis of the same 18 samples using **SoloTE** instead of scTE. SoloTE annotates aligned BAM files with TE information, providing a different quantification approach.

### Key Differences: scTE vs SoloTE

| Feature | scTE | SoloTE |
|---------|------|--------|
| **Input** | FASTQ files | BAM files (aligned) |
| **Alignment** | Built-in STAR alignment | Uses existing alignments |
| **Multi-mappers** | Custom handling | Tags reads in BAM |
| **Output** | Count matrix (.h5ad) | Annotated BAM + matrix |
| **TE annotation** | During alignment | Post-alignment |
| **Speed** | Slower (aligns reads) | Faster (uses BAM) |

## Dependencies

### Required Software

1. **Samtools** v1.16+ ✅ (installed: v1.19.2)
2. **BEDTools** v2.29.2+ ❌ (needs installation)
3. **R** v4+ ✅ (installed: v4.3.3)
4. **Python 3.9.5+** ✅ (installed: v3.12.3)

### Python Packages

- `pysam`
- `pandas` v1.5.0+

### Installation

```bash
# Install BEDTools (required)
sudo apt install bedtools

# OR install via conda
conda install -c bioconda bedtools

# Install Python dependencies
pip install pysam pandas
```

## Setup

### 1. Prepare TE Annotation

SoloTE needs RepeatMasker annotation in BED format:

```bash
cd SoloTE/
python SoloTE_RepeatMasker_to_BED.py -g hg38
```

This creates: `hg38_repeatmasker.bed`

### 2. Input Files

We'll use the BAM files from STARsolo alignment:

```
/home/jacobc/sc/starsolo_aligned/SRR14513984/Aligned.sortedByCoord.out.bam
/home/jacobc/sc/starsolo_aligned/SRR14514000/Aligned.sortedByCoord.out.bam
... (18 samples total)
```

### 3. Sample List

The 18 samples (from `../sample_mapping.csv`):
- Sample-100 → SRR14513984 (Control)
- Sample-17 → SRR14514000 (AD)
- Sample-19 → SRR14514008 (AD)
- ... (15 more)

## Running SoloTE

### Single Sample Example

```bash
python SoloTE/SoloTE_pipeline.py \
  --threads 8 \
  --bam /home/jacobc/sc/starsolo_aligned/SRR14513984/Aligned.sortedByCoord.out.bam \
  --teannotation SoloTE/hg38_repeatmasker.bed \
  --outputprefix SRR14513984 \
  --outputdir soloTE_output/SRR14513984
```

### Batch Processing (All 18 Samples)

```bash
bash run_soloTE_batch.sh
```

## Output Structure

For each sample, SoloTE creates:

```
soloTE_output/SRR14513984/
├── SRR14513984_annotated.bam          # BAM with TE tags
├── SRR14513984_annotated.bam.bai      # Index
├── SRR14513984_TE_counts.mtx          # Count matrix
├── SRR14513984_barcodes.tsv           # Cell barcodes
└── SRR14513984_features.tsv           # TE features
```

## Analysis Pipeline

1. **Run SoloTE** on 18 BAM files → TE-annotated BAMs + count matrices
2. **Merge matrices** → Single AnnData object (like scTE)
3. **Add metadata** → Cell type, diagnosis, etc. (same as scTE analysis)
4. **Differential analysis** → AD vs Control (genes + TEs)
5. **Compare with scTE** → Which TEs are consistent? Which differ?

## Comparison Analysis

After running both pipelines, we can compare:

### Questions to Answer

1. **Concordance**: Do scTE and SoloTE identify the same TEs?
2. **Sensitivity**: Which method detects more TEs?
3. **Reproducibility**: Are the AD-associated TEs consistent?
4. **Quantification**: Do expression levels correlate?

### Expected Files

```
soloTE_analysis/
├── SoloTE/                           # Cloned repository
├── soloTE_output/                    # Per-sample outputs (18 dirs)
├── merged_soloTE_18samples.h5ad      # Merged data
├── comparison/
│   ├── scTE_vs_soloTE_correlation.csv
│   ├── method_comparison.png
│   └── venn_diagrams.png
└── differential_results_soloTE/      # DE results from SoloTE
```

## Notes

### BAM File Compatibility

- SoloTE works with **any aligned BAM file**
- Our BAMs are from STAR (cell barcodes in CB tag)
- SoloTE will annotate reads with TE information using read name tags

### Cell Filtering

- SoloTE processes all cells in the BAM
- We'll filter to the same 60,328 cells as scTE analysis
- Use barcodes from: `../morabito_reference/GSE174367_snRNA-seq_cell_meta.csv.gz`

### TE Reference

- SoloTE uses RepeatMasker (same as scTE)
- Should have similar TE families detected
- Differences will be in quantification approach

## Timeline

**Estimated runtime per sample:** ~1-2 hours (depending on BAM size)
**Total for 18 samples:** ~18-36 hours (can parallelize)

## Pipeline Status

Check progress at any time:
```bash
bash monitor_progress.sh
```

## Full Analysis Pipeline

Once SoloTE completes, run the full analysis automatically:

```bash
bash run_full_pipeline.sh
```

This will:
1. ✅ Wait for SoloTE to finish (if still running)
2. Merge all 18 SoloTE outputs → `merged_soloTE_18samples.h5ad`
3. Run differential expression (AD vs Control) → `differential_results_soloTE/`
4. Compare scTE vs SoloTE → `comparison/` directory

## Comparison Output Files

The comparison generates CSV files with **joined data using TE names as keys**:

### 00_summary_statistics.csv
Overall comparison statistics:
- Number of TEs detected by each method
- Significant TEs in AD vs Control
- Upregulated/downregulated TEs
- Method overlap/concordance

### 01_TE_detection_comparison.csv
Columns:
- `TE_name`: TE subfamily name (key)
- `detected_scTE`: Boolean (found by scTE?)
- `detected_soloTE`: Boolean (found by SoloTE?)  
- `detection_method`: Both / scTE_only / SoloTE_only

### 02_expression_comparison.csv
For TEs detected by both methods:
- `TE_name`: TE subfamily name (key)
- `scTE_mean_expression`: Average counts across cells (scTE)
- `soloTE_mean_expression`: Average counts across cells (SoloTE)
- `log2_scTE`: Log2 transformed expression
- `log2_soloTE`: Log2 transformed expression
- `fold_difference`: Ratio of expression levels
- `pearson_r`: Correlation coefficient
- `spearman_r`: Rank correlation

### 03_differential_expression_comparison.csv
**Main comparison file with all DE statistics:**

Columns:
- `TE_name`: TE subfamily name (key)
- `scTE_logFC`: Log fold change AD vs Control (scTE)
- `scTE_pval`: P-value (scTE)
- `scTE_padj`: Adjusted p-value (scTE)
- `soloTE_logFC`: Log fold change AD vs Control (SoloTE)
- `soloTE_pval`: P-value (SoloTE)
- `soloTE_padj`: Adjusted p-value (SoloTE)
- `TE_family`: LINE, SINE, LTR, etc.
- `TE_class`: Repeat class
- `scTE_significant`: Boolean (padj < 0.05)
- `soloTE_significant`: Boolean (padj < 0.05)
- `scTE_direction`: Up / Down / None
- `soloTE_direction`: Up / Down / None
- `agreement`: Concordance category
  - `Both_significant_same_direction` ← **High confidence TEs**
  - `Both_significant_opposite_direction` ← Discordant
  - `Only_scTE_significant`
  - `Only_soloTE_significant`
  - `Neither_significant`
  - `Missing_data`
- `logFC_pearson_r`: Correlation of fold changes
- `logFC_spearman_r`: Rank correlation

## References

- **SoloTE paper**: [Nature Communications Biology (2022)](https://www.nature.com/articles/s42003-022-04020-5)
- **GitHub**: https://github.com/bvaldebenitom/SoloTE
- **scTE comparison**: See `../RESULTS_SUMMARY.md` for scTE findings

---

*Created: December 12, 2025*  
*Purpose: Compare TE quantification methods (scTE vs SoloTE)*  
*Status: SoloTE batch processing started 2025-12-12 04:58*
