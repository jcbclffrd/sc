# File Inventory: morabito_18samples/

Complete directory structure and file descriptions for the 18-sample Alzheimer's disease TE analysis.

---

## Directory Structure

```
morabito_18samples/
├── README.md                          # Main documentation
├── FILE_INVENTORY.md                  # This file - complete file inventory
├── RESULTS_SUMMARY.md                 # Comprehensive results and findings
├── ANALYSIS_SETUP_COMPLETE.md         # Setup and troubleshooting log
├── QUICKSTART.sh                      # Quick reference commands
├── run_pipeline.sh                    # Master pipeline execution script
│
├── 01_map_samples.py                  # Step 1: Map Sample-X to SRR IDs
├── 02_extract_18samples.py            # Step 2: Copy scTE .h5ad files
├── 03_filter_to_morabito_cells.py     # Step 3: Filter to exact cell barcodes
├── 04_merge_and_qc.py                 # Step 4: Merge samples + add metadata
├── 05_differential_analysis.py        # Step 5: DE analysis (AD vs Control)
│
├── sample_mapping.csv                 # Sample-X to SRR mapping (18 samples)
├── all_snrnaseq_samples.csv           # All 152 available snRNA-seq samples
├── te_names_rmsk.csv                  # TE annotations from RepeatMasker
│
├── scte_18samples/                    # Extracted scTE outputs (980 MB)
│   ├── SRR14513984/
│   │   └── SRR14513984.h5ad          # Sample-100 (Control, 4,520 cells)
│   ├── SRR14514000/
│   │   └── SRR14514000.h5ad          # Sample-17 (AD, 3,730 cells)
│   └── ... (16 more samples)
│
├── scte_18samples_filtered/           # Filtered to Morabito's barcodes
│   ├── SRR14513984/
│   │   └── SRR14513984.h5ad          # Filtered version
│   └── ... (18 samples)
│
├── merged_18samples.h5ad              # Merged dataset: 60,328 cells × 62,651 features
├── merged_18samples_genes.h5ad        # Genes only: 60,328 cells × 61,667 genes
│
└── differential_results/              # All differential expression results
    ├── genes_AD_vs_Control_all_cells.csv        # Overall gene DE (61,667 genes)
    ├── TEs_AD_vs_Control_all_cells.csv          # Overall TE DE (984 TEs)
    ├── ODC_AD_vs_Control.csv                    # Oligodendrocytes
    ├── EX_AD_vs_Control.csv                     # Excitatory neurons
    ├── INH_AD_vs_Control.csv                    # Inhibitory neurons
    ├── ASC_AD_vs_Control.csv                    # Astrocytes
    ├── OPC_AD_vs_Control.csv                    # Oligodendrocyte precursors
    ├── MG_AD_vs_Control.csv                     # Microglia
    ├── PER.END_AD_vs_Control.csv                # Pericytes/Endothelial
    ├── celltype_summary.csv                     # Cell-type DE summary table
    ├── TE_summary.csv                           # TE significance counts
    └── figures/
        └── volcano_plots.png                     # Volcano plots for visualization
```

---

## File Descriptions

### 📚 Documentation Files

| File | Size | Purpose |
|------|------|---------|
| **README.md** | ~20 KB | Main documentation: pipeline overview, sample info, usage instructions |
| **FILE_INVENTORY.md** | ~5 KB | This file - complete directory structure and file descriptions |
| **RESULTS_SUMMARY.md** | ~15 KB | Comprehensive results: top TEs, biological interpretation, comparisons |
| **ANALYSIS_SETUP_COMPLETE.md** | ~8 KB | Setup completion log, troubleshooting notes, pipeline validation |
| **QUICKSTART.sh** | ~2 KB | Quick reference commands for running the pipeline |

### 🐍 Python Scripts (Pipeline Steps)

| Script | Lines | Purpose | Input | Output |
|--------|-------|---------|-------|--------|
| **01_map_samples.py** | ~80 | Map Morabito's Sample-X names to SRR IDs | `../analysis/morabito_reference/` metadata | `sample_mapping.csv`, `all_snrnaseq_samples.csv` |
| **02_extract_18samples.py** | ~65 | Copy 18 scTE .h5ad files to working directory | `../scTE_output/`, `sample_mapping.csv` | `scte_18samples/` directory |
| **03_filter_to_morabito_cells.py** | ~90 | Filter to Morabito's exact cell barcodes | `scte_18samples/`, Morabito metadata | `scte_18samples_filtered/` |
| **04_merge_and_qc.py** | ~180 | Merge samples, add metadata, annotate TEs | `scte_18samples_filtered/`, RepeatMasker GTF | `merged_18samples.h5ad`, `merged_18samples_genes.h5ad` |
| **05_differential_analysis.py** | ~250 | Differential expression analysis (AD vs Control) | `merged_18samples.h5ad` | `differential_results/` directory |

### 📊 Data Files

#### Mapping and Metadata

| File | Rows | Purpose |
|------|------|---------|
| **sample_mapping.csv** | 18 | Maps Sample-X names to SRR IDs for the 18 samples used by Morabito |
| **all_snrnaseq_samples.csv** | 152 | Complete list of all available snRNA-seq samples in PRJNA729525 |
| **te_names_rmsk.csv** | 15,631 | TE annotations extracted from RepeatMasker hg38 GTF |

#### Main Data Objects

| File | Size | Shape | Description |
|------|------|-------|-------------|
| **merged_18samples.h5ad** | 926 MB | 60,328 cells × 62,651 features | Complete merged dataset with genes (61,667) + TEs (984) |
| **merged_18samples_genes.h5ad** | 697 MB | 60,328 cells × 61,667 genes | Genes only (for comparison with Morabito's original analysis) |

**Cell Metadata in .h5ad files:**
- `sample_id`: Sample-X identifier (e.g., Sample-100)
- `srr`: SRA run ID (e.g., SRR14513984)
- `diagnosis`: AD or Control
- `cell_type`: ODC, EX, INH, ASC, MG, OPC, PER.END
- `batch`: Batch number
- `age`, `sex`, `pmi`: Clinical metadata
- `tangle_stage`, `plaque_stage`: Neuropathology
- `n_counts`, `n_genes`: QC metrics
- `pct_counts_tes`: Percentage of UMIs from TEs

**Feature Metadata (.var):**
- `feature_type`: 'gene' or 'TE'

### 📁 Directories

#### scte_18samples/ (980 MB total)

Contains raw scTE output for 18 samples:
- **18 subdirectories** (one per SRR ID)
- Each contains: `{SRR}.h5ad` (scTE counts matrix)
- **Not filtered** - contains all cells from original scTE run
- Feature space: 62,651 features (genes + TEs)

#### scte_18samples_filtered/ (~900 MB total)

Filtered versions matching Morabito's exact cells:
- **18 subdirectories** (one per SRR ID)
- Each contains: `{SRR}.h5ad` (filtered to Morabito's barcodes)
- **60,328 cells total** across all 18 samples
- Used as input for merging step

#### differential_results/

Complete differential expression results:

**Overall Results:**
- `genes_AD_vs_Control_all_cells.csv` (61,667 rows) - All genes tested
- `TEs_AD_vs_Control_all_cells.csv` (984 rows) - All TEs tested

**Cell-Type Specific Results (7 files):**
- Each file contains genes + TEs for one cell type
- Format: `{CELLTYPE}_AD_vs_Control.csv`
- Cell types: ODC, EX, INH, ASC, OPC, MG, PER.END

**Summary Files:**
- `celltype_summary.csv` - Count of significant genes/TEs per cell type
- `TE_summary.csv` - Overall TE significance statistics

**Result File Columns:**
- `names`: Gene/TE name
- `scores`: Test statistic (z-score)
- `logfoldchanges`: Log2 fold change (AD vs Control)
- `pvals`: Raw p-values
- `pvals_adj`: Adjusted p-values (Benjamini-Hochberg FDR)
- `feature_type`: 'gene' or 'TE'

#### differential_results/figures/

Visualization outputs:
- `volcano_plots.png` - Multi-panel volcano plots (overall + cell types)

---

## Key Statistics

### Dataset Overview
- **Total cells:** 60,328 (37,875 AD, 22,453 Control)
- **Features:** 62,651 total (61,667 genes + 984 TEs)
- **Samples:** 18 (12 AD, 6 Control)
- **Cell types:** 7 major cell types
- **Mean UMIs/cell:** 4,698 (median: 3,199)
- **TE proportion:** ~64% of total UMI counts

### Cell Type Distribution
| Cell Type | Count | Percentage | Description |
|-----------|-------|------------|-------------|
| ODC | 36,920 | 61.2% | Oligodendrocytes (myelinating) |
| EX | 6,022 | 10.0% | Excitatory neurons |
| INH | 5,846 | 9.7% | Inhibitory neurons |
| ASC | 4,538 | 7.5% | Astrocytes |
| MG | 3,851 | 6.4% | Microglia |
| OPC | 2,714 | 4.5% | Oligodendrocyte precursors |
| PER.END | 437 | 0.7% | Pericytes/Endothelial |

### Differential Expression Summary
| Category | Count | Percentage |
|----------|-------|------------|
| TEs tested | 984 | 100% |
| TEs significant (padj < 0.05) | 374 | 38.0% |
| TEs upregulated in AD | 318 | 32.3% |
| TEs downregulated in AD | 56 | 5.7% |
| Genes tested | 61,667 | 100% |
| Genes significant (padj < 0.05) | ~15,000 | ~24% |

---

## File Dependencies

```
Pipeline Flow:

01_map_samples.py
  ↓ (creates sample_mapping.csv)
02_extract_18samples.py
  ↓ (creates scte_18samples/)
03_filter_to_morabito_cells.py
  ↓ (creates scte_18samples_filtered/)
04_merge_and_qc.py
  ↓ (creates merged_18samples.h5ad)
05_differential_analysis.py
  ↓ (creates differential_results/)
```

**External Dependencies:**
- `../scTE_output/` - Full scTE results for all 150 samples
- `../analysis/morabito_reference/` - Morabito's metadata and filtered matrix
- `/home/jacobc/sc/annotations/hg38_rmsk.gtf` - RepeatMasker TE annotations

---

## Disk Usage

| Item | Size | Count |
|------|------|-------|
| Python scripts | ~50 KB | 5 files |
| Documentation | ~50 KB | 5 files |
| CSV files | ~5 MB | 3 files |
| scte_18samples/ | 980 MB | 18 samples |
| scte_18samples_filtered/ | ~900 MB | 18 samples |
| merged_18samples.h5ad | 926 MB | 1 file |
| merged_18samples_genes.h5ad | 697 MB | 1 file |
| differential_results/ | ~150 MB | 15 files |
| **Total** | **~3.7 GB** | |

---

## Access Patterns

### Quick Analysis
```bash
# View top TEs
head -30 differential_results/TEs_AD_vs_Control_all_cells.csv

# Check cell type results
head differential_results/EX_AD_vs_Control.csv

# Summary statistics
cat differential_results/TE_summary.csv
```

### Load in Python
```python
import scanpy as sc
import pandas as pd

# Load merged data
adata = sc.read_h5ad('merged_18samples.h5ad')

# Load TE results
tes = pd.read_csv('differential_results/TEs_AD_vs_Control_all_cells.csv')
tes_sig = tes[tes['pvals_adj'] < 0.05]
```

### Re-run Pipeline
```bash
# Full pipeline
bash run_pipeline.sh

# Individual steps
python3 01_map_samples.py
python3 02_extract_18samples.py
python3 03_filter_to_morabito_cells.py
python3 04_merge_and_qc.py
python3 05_differential_analysis.py
```

---

## Version Control Notes

**Generated:** December 12, 2025  
**Python Environment:** `~/hcaTE/.venv/` (Python 3.12, scanpy 1.10.0)  
**Reference Genome:** hg38 (GRCh38)  
**TE Annotation:** RepeatMasker (15,537 TE subfamilies)

**Key Software Versions:**
- scanpy: 1.10.0
- anndata: 0.10.x
- pandas: 2.x
- numpy: 1.26.x
- scipy: 1.x
- matplotlib: 3.x

---

## Related Files (Outside This Directory)

```
../scTE_output/                        # Full scTE results (150 samples)
../analysis/morabito_reference/       # Morabito's published data
  ├── GSE174367_filtered_feature_bc_matrix.h5
  └── GSE174367_snRNA-seq_cell_meta.csv.gz
../annotations/
  ├── hg38_rmsk.gtf                   # RepeatMasker annotations
  └── gencode.v45.primary_assembly.annotation.gtf
```

---

## Notes

1. **Barcode Format:** Barcodes use format `{barcode}-{sample_number}` in Morabito's data but scTE outputs use clean barcodes. We create composite key: `{barcode_clean}_{sample_id}` for matching.

2. **Feature Type Annotation:** TEs are identified by matching feature names against RepeatMasker GTF. The `feature_type` column is added during the merge step (04_merge_and_qc.py).

3. **Cell Matching:** We match 60,328 cells vs Morabito's 61,770 (98% match rate). Small discrepancy due to strict barcode filtering.

4. **TE Quantification:** scTE uses both uniquely and multi-mapping reads for TE quantification, while CellRanger only uses unique reads (and doesn't quantify TEs at all).

---

*For detailed usage instructions, see README.md*  
*For results interpretation, see RESULTS_SUMMARY.md*
