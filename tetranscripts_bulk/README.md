# TEtranscripts Bulk RNA-seq Analysis

## Overview
Quantification and differential expression analysis of transposable elements (TEs) in bulk RNA-seq data from 18 Alzheimer's disease patients (11 AD, 7 Control) using TEtranscripts/TEcount.

## Analysis Date
December 26, 2024

## Sample Information
- **Total Samples**: 18
- **AD Samples**: 11 (Sample-17, 19, 22, 27, 33, 37, 43, 45, 46, 47, 50)
- **Control Samples**: 7 (Sample-52, 58, 66, 82, 90, 96, 100)
- **Patient Mapping**: See `../data/bulk_rnaseq_18patients.csv`

## Files

### Count Matrices
- **combined_counts_matrix.tsv** (3.7MB)
  - All genes (63,241) and TEs (1,330)
  - 18 samples × 64,571 features
  
- **combined_counts_matrix_genes_only.tsv** (3.6MB)
  - ENSEMBL genes only
  - 18 samples × 63,241 genes

- **combined_counts_matrix_TEs_only.tsv** (117KB)
  - Transposable elements only
  - 18 samples × 1,330 TE families

### Differential Expression Results
- **differential_TE_results_AD_vs_Control.csv** (182KB)
  - DESeq2 analysis results for all 1,145 tested TEs
  - Columns: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, TE, TE_family, TE_class, significance
  - **Key Finding**: 0 FDR-significant TEs, 20 nominally significant (p < 0.05)

### Visualizations
- **volcano_plot_TE_AD_vs_Control.pdf**
  - Log2 fold change vs -log10(p-value)
  
- **MA_plot_TE_AD_vs_Control.pdf**
  - Mean expression vs log2 fold change

- **PCA_plot_TE_expression.pdf**
  - Principal component analysis of TE expression

### Individual Sample Counts (not in git)
- `SRR*.cntTable` - TEcount output for each sample (18 files, 1.4MB each)

## Methods

### TE Quantification (TEcount)
```bash
TEcount --BAM <sample.bam> \
        --GTF gencode.v45.primary_assembly.annotation.gtf \
        --TE GRCh38_GENCODE_rmsk_TE.gtf \
        --mode multi \
        --stranded reverse \
        --sortByPos \
        --format BAM
```

### Differential Expression (DESeq2)
- Design: `~ diagnosis`
- Filtering: TEs expressed in ≥25% samples (≥10 reads)
- Tested: 1,145 of 1,330 TEs
- Multiple testing: Benjamini-Hochberg FDR

## Key Results

### Statistical Summary
- **TEs tested**: 1,145
- **FDR < 0.05**: 0
- **Nominal p < 0.05**: 20

### Top TE Families by Expression
1. L2 LINEs (L2a, L2c, L2b): >13.7M total reads
2. MIR SINEs: >7.6M total reads
3. Alu SINEs: Diverse subfamilies with substantial expression

### Nominal Trends (p < 0.05, uncorrected)
**Upregulated in AD** (11 TEs):
- LTR retrotransposons (ERVs): HERVE-int, LTR22A, LTR47B
- Alu SINEs: AluYe6, AluYj4

**Downregulated in AD** (9 TEs):
- DNA transposons (hAT-Charlie family): MER105, Charlie1, Charlie26a, Charlie9

## Interpretation

The lack of FDR-significant TEs likely reflects:
1. **Limited sample size** (n=18) - underpowered for multiple testing correction
2. **Bulk tissue heterogeneity** - cell-type-specific signals diluted
3. **High TE variability** - inherent noise in TE expression

**Recommendation**: Cell-type-specific TE analysis using matched snRNA-seq data to identify TE dysregulation in specific brain cell populations.

## Scripts

Located in `../scripts/`:
- `run_tetranscripts_bulk.sh` - Parallel TEcount quantification (6 concurrent jobs)
- `complete_missing_bulk_samples.sh` - Download/align/quantify missing samples
- `combine_tetranscripts_counts.py` - Merge individual count tables into matrices
- `differential_TE_analysis.R` - DESeq2 differential expression analysis
- `install_r_packages.R` - Install required R packages (DESeq2, ggplot2, etc.)

## Dependencies

### Software
- TEtranscripts v2.2.3 (Python venv: `../tetranscripts_env`)
- DESeq2 (R/Bioconductor)
- STAR aligner (for missing samples)
- Python: pandas, numpy
- R: ggplot2, pheatmap, dplyr, tidyr

### Annotations
- **Gene GTF**: gencode.v45.primary_assembly.annotation.gtf
- **TE GTF**: GRCh38_GENCODE_rmsk_TE.gtf (official TEtranscripts annotation, 706MB)

## Related Analyses
- **snRNA-seq**: `../morabito_18samples/` - Single-cell RNA-seq from same 18 patients
- **ATAC-seq**: Multi-omics chromatin accessibility data (to be analyzed)

## Citation
If using these results:
- TEtranscripts: Jin et al. (2015) Bioinformatics
- DESeq2: Love et al. (2014) Genome Biology
- Morabito et al. dataset: GSE174367

## Contact
Generated as part of multi-omics Alzheimer's disease analysis (scWSL project)
