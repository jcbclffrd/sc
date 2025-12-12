# Morabito 18-Sample Analysis Setup - COMPLETE ✓

**Date**: December 12, 2025  
**Location**: `/home/jacobc/sc/morabito_18samples/`

## What Was Created

### 1. Analysis Directory Structure
```
morabito_18samples/
├── README.md                          # Full documentation
├── QUICKSTART.sh                      # Quick reference
├── run_pipeline.sh                    # Master pipeline script
├── 01_map_samples.py                  # Map Sample-X → SRR
├── 02_extract_18samples.py            # Extract scTE data
├── 03_filter_to_morabito_cells.py     # Filter to exact cells
├── 04_merge_and_qc.py                 # Merge + metadata
├── 05_differential_analysis.py        # DE analysis (AD vs Control)
├── sample_mapping.csv                 # ✓ Sample-X → SRR mapping (created)
├── all_snrnaseq_samples.csv           # All snRNA-seq samples
└── scte_18samples/                    # ✓ Extracted scTE data (980 MB)
    ├── SRR14513984/                   # Sample-100 (Control)
    ├── SRR14514000/                   # Sample-17 (AD)
    ├── SRR14514008/                   # Sample-19 (AD)
    └── ... (18 samples total)
```

### 2. Sample Mapping Complete ✓

**Morabito's 18 samples mapped to SRR IDs:**

| Sample ID  | SRR ID      | GSM ID     | Diagnosis |
|------------|-------------|------------|-----------|
| Sample-100 | SRR14513984 | GSM5292838 | Control   |
| Sample-17  | SRR14514000 | GSM5292840 | AD        |
| Sample-19  | SRR14514008 | GSM5292841 | AD        |
| Sample-22  | SRR14514016 | GSM5292842 | AD        |
| Sample-27  | SRR14514024 | GSM5292843 | AD        |
| Sample-33  | SRR14514032 | GSM5292844 | AD        |
| Sample-37  | SRR14514040 | GSM5292845 | AD        |
| Sample-43  | SRR14514048 | GSM5292846 | AD        |
| Sample-45  | SRR14514056 | GSM5292847 | AD        |
| Sample-46  | SRR14514064 | GSM5292848 | AD        |
| Sample-47  | SRR14514072 | GSM5292849 | AD        |
| Sample-50  | SRR14514080 | GSM5292850 | AD        |
| Sample-52  | SRR14514088 | GSM5292851 | Control   |
| Sample-58  | SRR14514096 | GSM5292852 | Control   |
| Sample-66  | SRR14514104 | GSM5292853 | Control   |
| Sample-82  | SRR14514112 | GSM5292854 | Control   |
| Sample-90  | SRR14514120 | GSM5292855 | Control   |
| Sample-96  | SRR14514128 | GSM5292856 | Control   |

**Summary**: 12 AD samples, 6 Control samples

### 3. scTE Data Extracted ✓

**Total**: 18 samples, 980 MB  
**Source**: `/home/jacobc/sc/scTE_output/`  
**Destination**: `scte_18samples/`

All samples successfully copied!

## Next Steps - Ready to Run

### Option 1: Run Complete Pipeline
```bash
cd /home/jacobc/sc/morabito_18samples
./run_pipeline.sh
```

This will execute all remaining steps:
- ✓ Step 1: Sample mapping (DONE)
- ✓ Step 2: Extract samples (DONE)
- ⏳ Step 3: Filter to Morabito's cell barcodes
- ⏳ Step 4: Merge samples and add metadata
- ⏳ Step 5: Differential analysis (AD vs Control)

### Option 2: Run Steps Individually
```bash
# Step 3: Filter to Morabito's exact 61,770 cells
python3 03_filter_to_morabito_cells.py

# Step 4: Merge all 18 samples with metadata
python3 04_merge_and_qc.py

# Step 5: Differential TE/gene expression
python3 05_differential_analysis.py
```

## Expected Outputs

After running the complete pipeline, you will have:

### Data Files
- `scte_18samples_filtered/` - scTE data filtered to Morabito's exact cells
- `merged_18samples.h5ad` - Combined dataset (genes + TEs, ~61,770 cells)
- `merged_18samples_genes.h5ad` - Genes only (for comparison)

### Analysis Results
- `differential_results/genes_AD_vs_Control_all_cells.csv` - DE genes
- `differential_results/TEs_AD_vs_Control_all_cells.csv` - DE TEs
- `differential_results/<celltype>_AD_vs_Control.csv` - Cell-type specific
- `differential_results/celltype_summary.csv` - Summary table
- `differential_results/figures/volcano_plots.png` - Visualizations

## Key Analysis Goals

1. **Compare with Morabito**: Use exact same cells and samples
2. **TE Expression**: Identify AD-associated transposable elements
3. **Gene Expression**: Validate our gene quantification vs CellRanger
4. **Cell-Type Specificity**: TE/gene changes by cell type
5. **TE-Gene Correlations**: Co-expression patterns

## Reference Data (Already Downloaded)

From GEO GSE174367:
- `../analysis/morabito_reference/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5`
- `../analysis/morabito_reference/GSE174367_snRNA-seq_cell_meta.csv.gz`

Contains Morabito's:
- 61,770 filtered cells
- Cell type annotations (ODC, EX, INH, ASC, MG, OPC, PER.END)
- Diagnosis labels (AD, Control)
- Clinical metadata (age, sex, PMI, etc.)

## Notes

- **scTE data**: Contains BOTH genes and TEs (unlike Morabito's gene-only data)
- **Cell barcodes**: Will be matched exactly to Morabito's 61,770 cells
- **Filtering**: Step 3 filters our data to match their exact cell barcodes
- **DE analysis**: Will identify novel AD-associated TEs not in original paper

---

**Status**: ✓ Setup complete, ready to run pipeline  
**Next command**: `./run_pipeline.sh`
