# Pseudobulk Analysis Fix Instructions

## Problem Statement

The current implementation in `morabito_18samples/05_differential_analysis.py` contains a **critical pseudoreplication error** that treats individual cells as independent replicates when performing differential expression analysis. This approach violates the assumption of statistical independence and has resulted in **374 false positive significant transposable elements (TEs)**.

### Root Cause
- Individual cells from the same patient are not independent observations
- Current analysis: ~thousands of cells treated as independent replicates
- Correct approach: Aggregate cells by patient → 18 biological replicates

## Required Solution: Pseudobulk Aggregation

### Overview
Transform single-cell data into pseudobulk data by aggregating raw counts at the patient level before running differential expression analysis with DESeq2.

---

## Implementation Instructions

### Step 1: Create Pseudobulk Aggregation Function

Add the following function to aggregate counts by patient:

```python
def create_pseudobulk_counts(adata, group_by='patient_id', condition_col='condition'):
    """
    Aggregate single-cell counts to pseudobulk by summing raw counts per patient.
    
    Parameters:
    -----------
    adata : AnnData
        Single-cell AnnData object with raw counts
    group_by : str
        Column in adata.obs containing patient/sample identifiers (default: 'patient_id')
    condition_col : str
        Column in adata.obs containing condition labels (e.g., 'AD', 'Control')
    
    Returns:
    --------
    pseudobulk_df : pd.DataFrame
        DataFrame with genes/TEs as rows and patients as columns
    metadata_df : pd.DataFrame
        DataFrame with patient metadata (patient_id, condition, etc.)
    """
    import pandas as pd
    import numpy as np
    
    # Ensure we're working with raw counts
    if 'counts' in adata.layers:
        count_matrix = adata.layers['counts']
    else:
        count_matrix = adata.X
    
    # Convert to dense if sparse
    if hasattr(count_matrix, 'toarray'):
        count_matrix = count_matrix.toarray()
    
    # Create DataFrame for easier manipulation
    counts_df = pd.DataFrame(
        count_matrix,
        columns=adata.var_names,
        index=adata.obs.index
    )
    
    # Add patient and condition info
    counts_df[group_by] = adata.obs[group_by].values
    counts_df[condition_col] = adata.obs[condition_col].values
    
    # Aggregate by summing counts per patient
    pseudobulk_df = counts_df.groupby(group_by).sum()
    
    # Create metadata DataFrame
    metadata_df = adata.obs[[group_by, condition_col]].drop_duplicates()
    metadata_df = metadata_df.set_index(group_by)
    metadata_df = metadata_df.loc[pseudobulk_df.index]  # Ensure order matches
    
    print(f"Pseudobulk aggregation complete:")
    print(f"  - Original cells: {adata.n_obs}")
    print(f"  - Pseudobulk samples: {pseudobulk_df.shape[0]}")
    print(f"  - Features: {pseudobulk_df.shape[1]}")
    print(f"  - Conditions: {metadata_df[condition_col].value_counts().to_dict()}")
    
    return pseudobulk_df.T, metadata_df  # Transpose: genes as rows, patients as columns
```

### Step 2: Modify DESeq2 Analysis Pipeline

Replace the current cell-level DESeq2 analysis with pseudobulk analysis:

```python
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# Activate pandas conversion
pandas2ri.activate()

def run_deseq2_pseudobulk(pseudobulk_counts, metadata, condition_col='condition', 
                          ref_level='Control', test_level='AD'):
    """
    Run DESeq2 on pseudobulk data.
    
    Parameters:
    -----------
    pseudobulk_counts : pd.DataFrame
        Gene/TE counts matrix (features × patients)
    metadata : pd.DataFrame
        Patient metadata with condition information
    condition_col : str
        Column name for condition in metadata
    ref_level : str
        Reference condition level
    test_level : str
        Test condition level
        
    Returns:
    --------
    results_df : pd.DataFrame
        DESeq2 results with log2FC, p-values, adjusted p-values
    """
    
    # Import R packages
    deseq2 = importr('DESeq2')
    
    # Ensure counts are integers
    pseudobulk_counts = pseudobulk_counts.astype(int)
    
    # Create DESeq2 dataset
    ro.globalenv['counts_matrix'] = pseudobulk_counts
    ro.globalenv['metadata'] = metadata
    
    ro.r(f'''
    # Ensure factor levels are correct
    metadata${condition_col} <- factor(metadata${condition_col}, 
                                       levels = c("{ref_level}", "{test_level}"))
    
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(
        countData = counts_matrix,
        colData = metadata,
        design = ~ {condition_col}
    )
    
    # Filter low count genes (optional but recommended)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    # Run DESeq2
    dds <- DESeq(dds)
    
    # Extract results
    res <- results(dds, contrast = c("{condition_col}", "{test_level}", "{ref_level}"))
    res <- as.data.frame(res)
    ''')
    
    # Get results back to Python
    results_df = ro.r('res')
    results_df = pandas2ri.rpy2py(results_df)
    
    return results_df
```

### Step 3: Update Main Analysis Script

Modify `05_differential_analysis.py` to use pseudobulk approach:

```python
# Load single-cell data
adata = sc.read_h5ad('path_to_your_data.h5ad')

# CRITICAL: Ensure you have patient identifiers in adata.obs
# Check column names - might be 'patient', 'donor_id', 'sample_id', etc.
print("Available metadata columns:", adata.obs.columns.tolist())

# Create pseudobulk counts
pseudobulk_counts, patient_metadata = create_pseudobulk_counts(
    adata, 
    group_by='patient_id',  # ADJUST THIS TO YOUR ACTUAL COLUMN NAME
    condition_col='condition'  # ADJUST THIS TO YOUR ACTUAL COLUMN NAME
)

# Run DESeq2 on pseudobulk data
deseq2_results = run_deseq2_pseudobulk(
    pseudobulk_counts,
    patient_metadata,
    condition_col='condition',
    ref_level='Control',
    test_level='AD'
)

# Filter for significant results
sig_threshold_pval = 0.05
sig_threshold_lfc = 1.0  # log2 fold change

significant_TEs = deseq2_results[
    (deseq2_results['padj'] < sig_threshold_pval) & 
    (abs(deseq2_results['log2FoldChange']) > sig_threshold_lfc)
]

print(f"\nSignificant TEs with pseudobulk approach: {len(significant_TEs)}")
```

---

## Validation Steps

### 1. Sample Size Verification
```python
# Before aggregation
print(f"Number of cells: {adata.n_obs}")
print(f"Cells per patient: {adata.obs.groupby('patient_id').size()}")

# After aggregation
print(f"Number of pseudobulk samples: {pseudobulk_counts.shape[1]}")
assert pseudobulk_counts.shape[1] == 18, "Should have exactly 18 patient-level samples"
```

### 2. Count Aggregation Verification
```python
# Verify counts are summed correctly
patient_example = 'patient_001'  # Use actual patient ID
cells_from_patient = adata[adata.obs['patient_id'] == patient_example]
gene_example = 'TE_element_1'  # Use actual TE name

# Sum from single cells
manual_sum = cells_from_patient.layers['counts'][:, adata.var_names == gene_example].sum()
# From pseudobulk
pseudobulk_value = pseudobulk_counts.loc[gene_example, patient_example]

print(f"Manual sum: {manual_sum}")
print(f"Pseudobulk value: {pseudobulk_value}")
assert manual_sum == pseudobulk_value, "Counts should match!"
```

### 3. Results Comparison
```python
# Compare with previous (incorrect) analysis
print("\n=== Results Comparison ===")
print(f"Original (cell-level) significant TEs: 374")
print(f"Pseudobulk (patient-level) significant TEs: {len(significant_TEs)}")
print(f"Reduction: {374 - len(significant_TEs)} false positives eliminated")
print(f"Percentage reduction: {((374 - len(significant_TEs)) / 374 * 100):.1f}%")
```

### 4. Comparison with Bulk RNA-seq
```python
# Load bulk RNA-seq results (if available)
bulk_results = pd.read_csv('bulk_rnaseq_results.csv')  # Adjust path

# Find overlap
pseudobulk_sig_TEs = set(significant_TEs.index)
bulk_sig_TEs = set(bulk_results[bulk_results['padj'] < 0.05].index)

overlap = pseudobulk_sig_TEs & bulk_sig_TEs
print(f"\nOverlap with bulk RNA-seq:")
print(f"  - Pseudobulk significant TEs: {len(pseudobulk_sig_TEs)}")
print(f"  - Bulk significant TEs: {len(bulk_sig_TEs)}")
print(f"  - Overlap: {len(overlap)}")
print(f"  - Concordance: {len(overlap) / len(pseudobulk_sig_TEs) * 100:.1f}%")
```

### 5. Statistical Power Check
```python
# Ensure adequate depth per pseudobulk sample
total_counts_per_patient = pseudobulk_counts.sum(axis=0)
print("\nDepth per pseudobulk sample:")
print(total_counts_per_patient.describe())
print(f"Minimum depth: {total_counts_per_patient.min():,.0f}")
print(f"Median depth: {total_counts_per_patient.median():,.0f}")

# Should have sufficient counts (typically >1M per sample)
assert total_counts_per_patient.min() > 100000, "Warning: Low count depth detected"
```

---

## Expected Outcomes

### Before Fix (Current State)
- **Replicates:** ~Thousands of cells
- **Significant TEs:** 374
- **Problem:** Pseudoreplication → inflated statistical power → false positives

### After Fix (Expected State)
- **Replicates:** 18 patients
- **Significant TEs:** Expected 10-50 (substantial reduction)
- **Benefit:** Proper statistical inference, results match bulk RNA-seq

### Success Criteria
✅ Pseudobulk samples = 18 (one per patient)  
✅ Significant TEs reduced by >80% (from 374 to <75)  
✅ High concordance (>60%) with bulk RNA-seq results  
✅ No warnings from DESeq2 about sample size or dispersion  

---

## Important Notes

1. **Raw Counts Required:** Always use raw counts (not normalized) for DESeq2
2. **Patient Metadata:** Verify patient ID column name in your data
3. **Condition Labels:** Ensure AD/Control labels are correct and consistent
4. **R Dependencies:** Requires DESeq2 installed in R: `BiocManager::install("DESeq2")`
5. **Memory:** Pseudobulk aggregation significantly reduces memory requirements

## Troubleshooting

### Issue: "Cannot find column 'patient_id'"
**Solution:** Check actual column name with `print(adata.obs.columns)` and update `group_by` parameter

### Issue: "DESeq2 convergence failures"
**Solution:** Filter lowly expressed genes more stringently (e.g., `rowSums(counts(dds)) >= 20`)

### Issue: "All TEs significant"
**Solution:** You may still have pseudoreplication. Verify `pseudobulk_counts.shape[1] == 18`

### Issue: "No significant TEs"
**Solution:** May indicate true lack of signal or too stringent filtering. Check:
- Count depth per sample
- Dispersion estimates
- Relax FDR threshold slightly (e.g., 0.1)

---

## Additional Resources

- DESeq2 documentation: https://bioconductor.org/packages/DESeq2
- Pseudobulk best practices: Squair et al. (2021) Nature Communications
- Single-cell DE analysis: Crowell et al. (2020) Nature Communications

---

**Last Updated:** 2025-12-27  
**Status:** REQUIRES IMMEDIATE IMPLEMENTATION  
**Priority:** HIGH - Affects validity of all TE differential expression results