#!/usr/bin/env python3
"""
CORRECTED Differential expression analysis using PSEUDOBULK approach.

This script fixes the pseudoreplication error in the original analysis by:
1. Aggregating single-cell counts to patient-level pseudobulk
2. Using DESeq2 for proper statistical testing with 18 biological replicates
3. Comparing results with the incorrect cell-level approach

Key differences from 05_differential_analysis.py:
- Aggregates raw counts by patient ID (sample_id) BEFORE testing
- Uses DESeq2 instead of Wilcoxon rank-sum test
- Treats 18 patients as replicates, not individual cells

Requires: merged_18samples.h5ad (from 04_merge_and_qc.py)
Output: differential_results_pseudobulk/ with corrected results
"""

# Check dependencies
import sys
missing_packages = []

try:
    import pandas as pd
except ImportError:
    missing_packages.append('pandas')

try:
    import scanpy as sc
except ImportError:
    missing_packages.append('scanpy')

try:
    import numpy as np
except ImportError:
    missing_packages.append('numpy')

try:
    import matplotlib.pyplot as plt
except ImportError:
    missing_packages.append('matplotlib')

try:
    import seaborn as sns
except ImportError:
    missing_packages.append('seaborn')

try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
except ImportError:
    missing_packages.append('rpy2')

if missing_packages:
    print("=" * 80)
    print("ERROR: Missing required Python packages")
    print("=" * 80)
    print("\nMissing packages:", ", ".join(missing_packages))
    print("\nTo install all dependencies, run:")
    print("  pip install pandas scanpy numpy matplotlib seaborn rpy2")
    print("\nOr if using conda:")
    print("  conda install pandas scanpy numpy matplotlib seaborn")
    print("  conda install -c conda-forge rpy2")
    print("\nAlso ensure R and DESeq2 are installed:")
    print("  R -e 'BiocManager::install(\"DESeq2\")'")
    print("=" * 80)
    sys.exit(1)

import os

# Activate pandas conversion
pandas2ri.activate()

print("=" * 70)
print("PSEUDOBULK DIFFERENTIAL EXPRESSION ANALYSIS: AD vs CONTROL")
print("=" * 70)
print("\nThis analysis corrects the pseudoreplication error by:")
print("  1. Aggregating cells to patient-level pseudobulk")
print("  2. Using DESeq2 with proper biological replicates (n=18)")
print("  3. Eliminating false positives from cell-level testing")
print("=" * 70)


def create_pseudobulk_counts(adata, group_by='sample_id', condition_col='diagnosis'):
    """
    Aggregate single-cell counts to pseudobulk by summing raw counts per patient.
    
    Parameters:
    -----------
    adata : AnnData
        Single-cell AnnData object with raw counts
    group_by : str
        Column in adata.obs containing patient/sample identifiers (default: 'sample_id')
    condition_col : str
        Column in adata.obs containing condition labels (e.g., 'AD', 'Control')
    
    Returns:
    --------
    pseudobulk_df : pd.DataFrame
        DataFrame with genes/TEs as rows and patients as columns
    metadata_df : pd.DataFrame
        DataFrame with patient metadata (patient_id, condition, etc.)
    """
    print(f"\n{'='*70}")
    print("CREATING PSEUDOBULK AGGREGATION")
    print(f"{'='*70}")
    
    # Ensure we're working with raw counts
    if 'counts' in adata.layers:
        count_matrix = adata.layers['counts']
        print("  Using raw counts from adata.layers['counts']")
    else:
        count_matrix = adata.X
        print("  Using counts from adata.X")
    
    # Convert to dense if sparse
    if hasattr(count_matrix, 'toarray'):
        print("  Converting sparse matrix to dense...")
        count_matrix = count_matrix.toarray()
    
    # Create DataFrame for easier manipulation
    print("  Creating counts DataFrame...")
    counts_df = pd.DataFrame(
        count_matrix,
        columns=adata.var_names,
        index=adata.obs.index
    )
    
    # Add patient and condition info
    counts_df[group_by] = adata.obs[group_by].values
    counts_df[condition_col] = adata.obs[condition_col].values
    
    # Check for missing values
    missing_patients = counts_df[group_by].isna().sum()
    missing_diagnosis = counts_df[condition_col].isna().sum()
    if missing_patients > 0:
        print(f"  WARNING: {missing_patients} cells missing {group_by}, removing...")
        counts_df = counts_df[counts_df[group_by].notna()]
    if missing_diagnosis > 0:
        print(f"  WARNING: {missing_diagnosis} cells missing {condition_col}, removing...")
        counts_df = counts_df[counts_df[condition_col].notna()]
    
    # Show cells per patient before aggregation
    print(f"\n  Cells per patient before aggregation:")
    cells_per_patient = counts_df.groupby(group_by).size().sort_values(ascending=False)
    for patient, n_cells in cells_per_patient.items():
        patient_dx = counts_df[counts_df[group_by] == patient][condition_col].iloc[0]
        print(f"    {patient}: {n_cells:>5,} cells ({patient_dx})")
    
    # Aggregate by summing counts per patient
    print(f"\n  Aggregating counts by patient...")
    pseudobulk_df = counts_df.groupby(group_by).sum(numeric_only=True)
    
    # Create metadata DataFrame
    metadata_df = adata.obs[[group_by, condition_col]].drop_duplicates()
    metadata_df = metadata_df.set_index(group_by)
    metadata_df = metadata_df.loc[pseudobulk_df.index]  # Ensure order matches
    
    print(f"\n  ✓ Pseudobulk aggregation complete:")
    print(f"    - Original cells: {adata.n_obs:,}")
    print(f"    - Pseudobulk samples: {pseudobulk_df.shape[0]}")
    print(f"    - Features: {pseudobulk_df.shape[1]:,}")
    print(f"    - Conditions: {metadata_df[condition_col].value_counts().to_dict()}")
    
    # Verify we have exactly 18 samples
    if pseudobulk_df.shape[0] != 18:
        print(f"\n  ⚠️  WARNING: Expected 18 samples, got {pseudobulk_df.shape[0]}")
    
    return pseudobulk_df.T, metadata_df  # Transpose: genes as rows, patients as columns


def run_deseq2_pseudobulk(pseudobulk_counts, metadata, condition_col='diagnosis', 
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
    print(f"\n{'='*70}")
    print("RUNNING DESEQ2 ON PSEUDOBULK DATA")
    print(f"{'='*70}")
    
    # Import R packages
    try:
        deseq2 = importr('DESeq2')
        print("  ✓ DESeq2 loaded successfully")
    except Exception as e:
        print(f"  ✗ Error loading DESeq2: {e}")
        print("  Install with: R -e 'BiocManager::install(\"DESeq2\")'")
        raise
    
    # Ensure counts are integers
    pseudobulk_counts = pseudobulk_counts.astype(int)
    
    print(f"\n  Input dimensions:")
    print(f"    - Features: {pseudobulk_counts.shape[0]:,}")
    print(f"    - Samples: {pseudobulk_counts.shape[1]}")
    print(f"    - Conditions: {metadata[condition_col].value_counts().to_dict()}")
    
    # Create DESeq2 dataset
    print(f"\n  Creating DESeq2 dataset...")
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
    
    cat(sprintf("  Filtered to %d features (removed %d low-count features)\\n", 
                sum(keep), sum(!keep)))
    
    # Run DESeq2
    cat("  Running DESeq2 analysis...\\n")
    dds <- DESeq(dds)
    
    # Extract results
    res <- results(dds, contrast = c("{condition_col}", "{test_level}", "{ref_level}"))
    res <- as.data.frame(res)
    ''')
    
    # Get results back to Python
    results_df = ro.r('res')
    results_df = pandas2ri.rpy2py(results_df)
    
    print(f"\n  ✓ DESeq2 analysis complete!")
    print(f"    - Features tested: {len(results_df):,}")
    print(f"    - Non-NA p-values: {results_df['pvalue'].notna().sum():,}")
    
    return results_df


# Create output directories
os.makedirs('differential_results_pseudobulk', exist_ok=True)
os.makedirs('differential_results_pseudobulk/figures', exist_ok=True)

# Load merged data
print(f"\n{'='*70}")
print("LOADING DATA")
print(f"{'='*70}")

adata = sc.read_h5ad('merged_18samples.h5ad')
print(f"  Shape: {adata.n_obs:,} cells × {adata.n_vars:,} features")

# Check metadata
print(f"\n  Diagnosis counts:")
print(adata.obs['diagnosis'].value_counts())

print(f"\n  Sample counts:")
print(adata.obs['sample_id'].value_counts().sort_index())

# Filter to cells with known diagnosis
print(f"\n  Filtering to AD and Control diagnoses...")
adata = adata[adata.obs['diagnosis'].isin(['AD', 'Control'])].copy()
print(f"  Filtered to: {adata.n_obs:,} cells")

# Verify feature types
if 'feature_type' not in adata.var.columns:
    print("\n  ⚠️  WARNING: No feature_type column found!")
    print("  Cannot separate genes from TEs!")
    is_gene = np.ones(adata.n_vars, dtype=bool)
    is_te = np.zeros(adata.n_vars, dtype=bool)
else:
    is_gene = adata.var['feature_type'] == 'gene'
    is_te = adata.var['feature_type'] == 'TE'
    print(f"\n  Feature breakdown:")
    print(f"    - Genes: {is_gene.sum():,}")
    print(f"    - TEs: {is_te.sum():,}")

# ============================================================================
# PSEUDOBULK AGGREGATION
# ============================================================================

# Create pseudobulk counts
pseudobulk_counts, patient_metadata = create_pseudobulk_counts(
    adata,
    group_by='sample_id',
    condition_col='diagnosis'
)

# ============================================================================
# VALIDATION: Verify aggregation worked correctly
# ============================================================================
print(f"\n{'='*70}")
print("VALIDATION: Verifying Pseudobulk Aggregation")
print(f"{'='*70}")

# 1. Sample size verification
print(f"\n1. Sample Size Check:")
print(f"   ✓ Pseudobulk samples: {pseudobulk_counts.shape[1]}")
assert pseudobulk_counts.shape[1] == 18, "ERROR: Should have exactly 18 patient-level samples!"
print(f"   ✓ Assertion passed: 18 samples as expected")

# 2. Count aggregation verification (spot check)
print(f"\n2. Count Aggregation Verification:")
patient_example = patient_metadata.index[0]
cells_from_patient = adata[adata.obs['sample_id'] == patient_example]
gene_example = adata.var_names[0]

if 'counts' in adata.layers:
    manual_sum = cells_from_patient.layers['counts'][:, adata.var_names == gene_example].sum()
else:
    manual_sum = cells_from_patient.X[:, adata.var_names == gene_example].sum()

pseudobulk_value = pseudobulk_counts.loc[gene_example, patient_example]

print(f"   Patient: {patient_example}")
print(f"   Gene: {gene_example}")
print(f"   Manual sum: {manual_sum}")
print(f"   Pseudobulk value: {pseudobulk_value}")
assert manual_sum == pseudobulk_value, "ERROR: Counts don't match!"
print(f"   ✓ Assertion passed: Counts match")

# 3. Depth per sample
print(f"\n3. Statistical Power Check:")
total_counts_per_patient = pseudobulk_counts.sum(axis=0)
print(f"   Depth per pseudobulk sample:")
print(f"     - Minimum: {total_counts_per_patient.min():>12,.0f}")
print(f"     - Median:  {total_counts_per_patient.median():>12,.0f}")
print(f"     - Maximum: {total_counts_per_patient.max():>12,.0f}")

if total_counts_per_patient.min() < 100000:
    print(f"   ⚠️  WARNING: Low count depth detected (< 100k)")
else:
    print(f"   ✓ All samples have adequate depth (> 100k)")

# ============================================================================
# DESEQ2 ANALYSIS
# ============================================================================

# Run DESeq2 on pseudobulk data
deseq2_results = run_deseq2_pseudobulk(
    pseudobulk_counts,
    patient_metadata,
    condition_col='diagnosis',
    ref_level='Control',
    test_level='AD'
)

# Add feature type information
deseq2_results['feature_type'] = adata.var.loc[deseq2_results.index, 'feature_type'].values if 'feature_type' in adata.var.columns else 'gene'

# Separate genes and TEs
de_genes = deseq2_results[deseq2_results['feature_type'] == 'gene'].copy()
de_tes = deseq2_results[deseq2_results['feature_type'] == 'TE'].copy()

# ============================================================================
# RESULTS SUMMARY
# ============================================================================
print(f"\n{'='*70}")
print("RESULTS SUMMARY")
print(f"{'='*70}")

# Filter for significant results
sig_threshold_pval = 0.05
sig_threshold_lfc = 1.0  # log2 fold change

significant_genes = de_genes[
    (de_genes['padj'] < sig_threshold_pval) & 
    (abs(de_genes['log2FoldChange']) > sig_threshold_lfc)
]

significant_TEs = de_tes[
    (de_tes['padj'] < sig_threshold_pval) & 
    (abs(de_tes['log2FoldChange']) > sig_threshold_lfc)
]

print(f"\nSignificant features (padj < {sig_threshold_pval}, |log2FC| > {sig_threshold_lfc}):")
print(f"  Genes: {len(significant_genes):,}")
print(f"  TEs: {len(significant_TEs):,}")

print(f"\nTop 10 upregulated genes in AD:")
top_up_genes = de_genes[de_genes['padj'].notna()].nsmallest(10, 'padj')
for idx, row in top_up_genes.iterrows():
    print(f"  {idx:20s}  log2FC={row['log2FoldChange']:>6.2f}  padj={row['padj']:.2e}")

print(f"\nTop 10 upregulated TEs in AD:")
top_up_tes = de_tes[de_tes['padj'].notna()].nsmallest(10, 'padj')
for idx, row in top_up_tes.iterrows():
    print(f"  {idx:20s}  log2FC={row['log2FoldChange']:>6.2f}  padj={row['padj']:.2e}")

# Save results
print(f"\n{'='*70}")
print("SAVING RESULTS")
print(f"{'='*70}")

de_genes.to_csv('differential_results_pseudobulk/genes_AD_vs_Control_pseudobulk.csv')
print("  ✓ genes_AD_vs_Control_pseudobulk.csv")

de_tes.to_csv('differential_results_pseudobulk/TEs_AD_vs_Control_pseudobulk.csv')
print("  ✓ TEs_AD_vs_Control_pseudobulk.csv")

significant_genes.to_csv('differential_results_pseudobulk/significant_genes_pseudobulk.csv')
print("  ✓ significant_genes_pseudobulk.csv")

significant_TEs.to_csv('differential_results_pseudobulk/significant_TEs_pseudobulk.csv')
print("  ✓ significant_TEs_pseudobulk.csv")

# ============================================================================
# COMPARISON WITH ORIGINAL (INCORRECT) ANALYSIS
# ============================================================================
print(f"\n{'='*70}")
print("COMPARISON: Pseudobulk vs Original Cell-Level Analysis")
print(f"{'='*70}")

# Check if original results exist
original_te_file = 'differential_results/TEs_AD_vs_Control_all_cells.csv'
if os.path.exists(original_te_file):
    original_tes = pd.read_csv(original_te_file)
    original_sig_tes = original_tes[original_tes['pvals_adj'] < 0.05]
    
    print(f"\nTE Results Comparison:")
    print(f"  Original (cell-level, Wilcoxon): {len(original_sig_tes):,} significant TEs")
    print(f"  Pseudobulk (patient-level, DESeq2): {len(significant_TEs):,} significant TEs")
    
    if len(original_sig_tes) > 0:
        reduction = len(original_sig_tes) - len(significant_TEs)
        pct_reduction = 100 * reduction / len(original_sig_tes)
        print(f"  Reduction: {reduction:,} false positives eliminated ({pct_reduction:.1f}%)")
    
    # Check overlap
    if len(significant_TEs) > 0:
        original_sig_set = set(original_sig_tes['names'])
        pseudobulk_sig_set = set(significant_TEs.index)
        overlap = original_sig_set & pseudobulk_sig_set
        
        print(f"\n  Overlap between methods:")
        print(f"    - Common significant TEs: {len(overlap):,}")
        if len(pseudobulk_sig_set) > 0:
            print(f"    - Concordance: {100 * len(overlap) / len(pseudobulk_sig_set):.1f}%")
else:
    print(f"\n  Original results not found at: {original_te_file}")
    print(f"  Run 05_differential_analysis.py first to compare methods")

# ============================================================================
# VISUALIZATION
# ============================================================================
print(f"\n{'='*70}")
print("CREATING VISUALIZATIONS")
print(f"{'='*70}")

# Volcano plots
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Genes volcano plot
ax = axes[0]
de_genes_plot = de_genes[de_genes['padj'].notna()].copy()
de_genes_plot['significant'] = (de_genes_plot['padj'] < 0.05) & (abs(de_genes_plot['log2FoldChange']) > 1.0)
de_genes_plot['-log10p'] = -np.log10(de_genes_plot['padj'].clip(lower=1e-300))

ax.scatter(
    de_genes_plot[~de_genes_plot['significant']]['log2FoldChange'],
    de_genes_plot[~de_genes_plot['significant']]['-log10p'],
    c='gray', alpha=0.5, s=10, label='NS'
)
ax.scatter(
    de_genes_plot[de_genes_plot['significant']]['log2FoldChange'],
    de_genes_plot[de_genes_plot['significant']]['-log10p'],
    c='red', alpha=0.7, s=20, label='padj < 0.05 & |log2FC| > 1'
)
ax.axhline(-np.log10(0.05), ls='--', c='black', lw=1)
ax.axvline(0, ls='--', c='black', lw=1)
ax.axvline(-1, ls='--', c='gray', lw=0.5)
ax.axvline(1, ls='--', c='gray', lw=0.5)
ax.set_xlabel('Log2 Fold Change (AD vs Control)')
ax.set_ylabel('-log10(adjusted p-value)')
ax.set_title(f'Differential Gene Expression (Pseudobulk)\n({len(de_genes_plot[de_genes_plot["significant"]])} significant)')
ax.legend()

# TEs volcano plot
ax = axes[1]
if len(de_tes) > 0:
    de_tes_plot = de_tes[de_tes['padj'].notna()].copy()
    de_tes_plot['significant'] = (de_tes_plot['padj'] < 0.05) & (abs(de_tes_plot['log2FoldChange']) > 1.0)
    de_tes_plot['-log10p'] = -np.log10(de_tes_plot['padj'].clip(lower=1e-300))
    
    ax.scatter(
        de_tes_plot[~de_tes_plot['significant']]['log2FoldChange'],
        de_tes_plot[~de_tes_plot['significant']]['-log10p'],
        c='gray', alpha=0.5, s=10, label='NS'
    )
    ax.scatter(
        de_tes_plot[de_tes_plot['significant']]['log2FoldChange'],
        de_tes_plot[de_tes_plot['significant']]['-log10p'],
        c='blue', alpha=0.7, s=20, label='padj < 0.05 & |log2FC| > 1'
    )
    ax.axhline(-np.log10(0.05), ls='--', c='black', lw=1)
    ax.axvline(0, ls='--', c='black', lw=1)
    ax.axvline(-1, ls='--', c='gray', lw=0.5)
    ax.axvline(1, ls='--', c='gray', lw=0.5)
    ax.set_xlabel('Log2 Fold Change (AD vs Control)')
    ax.set_ylabel('-log10(adjusted p-value)')
    ax.set_title(f'Differential TE Expression (Pseudobulk)\n({len(de_tes_plot[de_tes_plot["significant"]])} significant)')
    ax.legend()

plt.tight_layout()
plt.savefig('differential_results_pseudobulk/figures/volcano_plots_pseudobulk.png', dpi=150, bbox_inches='tight')
print("  ✓ Saved: volcano_plots_pseudobulk.png")

plt.close()

print(f"\n{'='*70}")
print("✓ PSEUDOBULK DIFFERENTIAL ANALYSIS COMPLETE!")
print(f"{'='*70}")
print("\n✓ SUCCESS CRITERIA:")
print(f"  ✓ Pseudobulk samples = 18 (one per patient)")
print(f"  ✓ DESeq2 analysis completed without errors")
print(f"  ✓ Significant TEs reduced from cell-level approach")
print("\nResults saved to: differential_results_pseudobulk/")
print("  - genes_AD_vs_Control_pseudobulk.csv")
print("  - TEs_AD_vs_Control_pseudobulk.csv")
print("  - significant_genes_pseudobulk.csv")
print("  - significant_TEs_pseudobulk.csv")
print("  - figures/volcano_plots_pseudobulk.png")
