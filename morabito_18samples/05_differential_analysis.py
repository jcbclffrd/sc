#!/usr/bin/env python3
"""
Differential expression analysis for TEs and genes (AD vs Control).

This script performs differential expression testing to identify:
1. TEs that are upregulated/downregulated in AD
2. Genes that are upregulated/downregulated in AD
3. Cell-type specific TE/gene changes

Requires: merged_18samples.h5ad (from 04_merge_and_qc.py)
Output: differential_results/ with CSV files and figures
"""

import os
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

print("=" * 70)
print("DIFFERENTIAL EXPRESSION ANALYSIS: AD vs CONTROL")
print("=" * 70)

# Create output directories
os.makedirs('differential_results', exist_ok=True)
os.makedirs('differential_results/figures', exist_ok=True)

# Load merged data
print("\nLoading merged data...")
adata = sc.read_h5ad('merged_18samples.h5ad')
print(f"  Shape: {adata.n_obs:,} cells × {adata.n_vars:,} features")

# Check metadata
print(f"\nDiagnosis counts:")
print(adata.obs['diagnosis'].value_counts())

print(f"\nCell type counts:")
print(adata.obs['cell_type'].value_counts())

# Filter to cells with known diagnosis and cell type
adata = adata[
    (adata.obs['diagnosis'].isin(['AD', 'Control'])) &
    (adata.obs['cell_type'] != 'Unknown')
].copy()
print(f"\nFiltered to: {adata.n_obs:,} cells with diagnosis and cell type")

# Normalize and log-transform
print("\nNormalizing data...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Separate genes and TEs
if 'feature_type' not in adata.var.columns:
    print("\n⚠️  WARNING: No feature_type column found!")
    print("Assuming all features are genes. Cannot separate TEs.")
    is_gene = np.ones(adata.n_vars, dtype=bool)
    is_te = np.zeros(adata.n_vars, dtype=bool)
else:
    is_gene = adata.var['feature_type'] == 'gene'
    is_te = adata.var['feature_type'] == 'TE'

print(f"\nFeature breakdown:")
print(f"  Genes: {is_gene.sum():,}")
print(f"  TEs: {is_te.sum():,}")

# ============================================================================
# 1. OVERALL DIFFERENTIAL ANALYSIS (ALL CELL TYPES)
# ============================================================================
print(f"\n{'='*70}")
print("1. OVERALL DIFFERENTIAL ANALYSIS (ALL CELL TYPES)")
print(f"{'='*70}")

# Create a copy for overall analysis
adata_all = adata.copy()

# Differential expression: AD vs Control
print("\nTesting for differential expression...")
sc.tl.rank_genes_groups(
    adata_all,
    groupby='diagnosis',
    groups=['AD'],
    reference='Control',
    method='wilcoxon',
    key_added='de_diagnosis'
)

# Extract results
de_results = sc.get.rank_genes_groups_df(adata_all, group='AD', key='de_diagnosis')
de_results['feature_type'] = adata_all.var.loc[de_results['names'], 'feature_type'].values if 'feature_type' in adata_all.var.columns else 'gene'

# Separate genes and TEs
de_genes = de_results[de_results['feature_type'] == 'gene'].copy()
de_tes = de_results[de_results['feature_type'] == 'TE'].copy()

print(f"\nDifferential genes: {len(de_genes):,}")
print(f"Differential TEs: {len(de_tes):,}")

# Save results
de_genes.to_csv('differential_results/genes_AD_vs_Control_all_cells.csv', index=False)
de_tes.to_csv('differential_results/TEs_AD_vs_Control_all_cells.csv', index=False)

print("\nTop 20 upregulated genes in AD:")
top_up_genes = de_genes.nsmallest(20, 'pvals_adj')
for idx, row in top_up_genes.iterrows():
    print(f"  {row['names']:20s}  logFC={row['logfoldchanges']:>6.2f}  padj={row['pvals_adj']:.2e}")

print("\nTop 20 upregulated TEs in AD:")
top_up_tes = de_tes.nsmallest(20, 'pvals_adj')
for idx, row in top_up_tes.iterrows():
    print(f"  {row['names']:20s}  logFC={row['logfoldchanges']:>6.2f}  padj={row['pvals_adj']:.2e}")

# ============================================================================
# 2. CELL-TYPE SPECIFIC ANALYSIS
# ============================================================================
print(f"\n{'='*70}")
print("2. CELL-TYPE SPECIFIC DIFFERENTIAL ANALYSIS")
print(f"{'='*70}")

cell_types = adata.obs['cell_type'].unique()
print(f"\nAnalyzing {len(cell_types)} cell types: {', '.join(cell_types)}")

celltype_results = {}

for cell_type in cell_types:
    print(f"\n{cell_type}:")
    
    # Filter to this cell type
    adata_ct = adata[adata.obs['cell_type'] == cell_type].copy()
    
    # Check if we have both AD and Control
    dx_counts = adata_ct.obs['diagnosis'].value_counts()
    if 'AD' not in dx_counts or 'Control' not in dx_counts:
        print(f"  ⚠️  Skipping (missing AD or Control cells)")
        continue
    
    if dx_counts['AD'] < 10 or dx_counts['Control'] < 10:
        print(f"  ⚠️  Skipping (too few cells: AD={dx_counts['AD']}, Control={dx_counts['Control']})")
        continue
    
    print(f"  Cells: AD={dx_counts['AD']}, Control={dx_counts['Control']}")
    
    # Differential expression
    sc.tl.rank_genes_groups(
        adata_ct,
        groupby='diagnosis',
        groups=['AD'],
        reference='Control',
        method='wilcoxon',
        key_added='de_diagnosis'
    )
    
    # Extract results
    de_ct = sc.get.rank_genes_groups_df(adata_ct, group='AD', key='de_diagnosis')
    de_ct['feature_type'] = adata_ct.var.loc[de_ct['names'], 'feature_type'].values if 'feature_type' in adata_ct.var.columns else 'gene'
    de_ct['cell_type'] = cell_type
    
    celltype_results[cell_type] = de_ct
    
    # Save
    de_ct.to_csv(f'differential_results/{cell_type}_AD_vs_Control.csv', index=False)
    
    # Count significant features
    sig_genes = de_ct[(de_ct['feature_type'] == 'gene') & (de_ct['pvals_adj'] < 0.05)]
    sig_tes = de_ct[(de_ct['feature_type'] == 'TE') & (de_ct['pvals_adj'] < 0.05)]
    
    print(f"  Significant (padj < 0.05): {len(sig_genes)} genes, {len(sig_tes)} TEs")

# ============================================================================
# 3. VISUALIZATION
# ============================================================================
print(f"\n{'='*70}")
print("3. CREATING VISUALIZATIONS")
print(f"{'='*70}")

# Volcano plot for overall results
print("\nCreating volcano plots...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Genes volcano plot
ax = axes[0]
de_genes_plot = de_genes.copy()
de_genes_plot['significant'] = de_genes_plot['pvals_adj'] < 0.05
de_genes_plot['-log10p'] = -np.log10(de_genes_plot['pvals_adj'].clip(lower=1e-300))

ax.scatter(
    de_genes_plot[~de_genes_plot['significant']]['logfoldchanges'],
    de_genes_plot[~de_genes_plot['significant']]['-log10p'],
    c='gray', alpha=0.5, s=10, label='NS'
)
ax.scatter(
    de_genes_plot[de_genes_plot['significant']]['logfoldchanges'],
    de_genes_plot[de_genes_plot['significant']]['-log10p'],
    c='red', alpha=0.7, s=20, label='padj < 0.05'
)
ax.axhline(-np.log10(0.05), ls='--', c='black', lw=1)
ax.axvline(0, ls='--', c='black', lw=1)
ax.set_xlabel('Log2 Fold Change (AD vs Control)')
ax.set_ylabel('-log10(adjusted p-value)')
ax.set_title(f'Differential Gene Expression\n({len(de_genes_plot[de_genes_plot["significant"]])} significant)')
ax.legend()

# TEs volcano plot
ax = axes[1]
if len(de_tes) > 0:
    de_tes_plot = de_tes.copy()
    de_tes_plot['significant'] = de_tes_plot['pvals_adj'] < 0.05
    de_tes_plot['-log10p'] = -np.log10(de_tes_plot['pvals_adj'].clip(lower=1e-300))
    
    ax.scatter(
        de_tes_plot[~de_tes_plot['significant']]['logfoldchanges'],
        de_tes_plot[~de_tes_plot['significant']]['-log10p'],
        c='gray', alpha=0.5, s=10, label='NS'
    )
    ax.scatter(
        de_tes_plot[de_tes_plot['significant']]['logfoldchanges'],
        de_tes_plot[de_tes_plot['significant']]['-log10p'],
        c='blue', alpha=0.7, s=20, label='padj < 0.05'
    )
    ax.axhline(-np.log10(0.05), ls='--', c='black', lw=1)
    ax.axvline(0, ls='--', c='black', lw=1)
    ax.set_xlabel('Log2 Fold Change (AD vs Control)')
    ax.set_ylabel('-log10(adjusted p-value)')
    ax.set_title(f'Differential TE Expression\n({len(de_tes_plot[de_tes_plot["significant"]])} significant)')
    ax.legend()

plt.tight_layout()
plt.savefig('differential_results/figures/volcano_plots.png', dpi=150, bbox_inches='tight')
print("  ✓ Saved: volcano_plots.png")

# Cell type heatmap
print("\nCreating cell-type summary heatmap...")

if len(celltype_results) > 0:
    # Count significant features per cell type
    ct_summary = []
    for ct, de_ct in celltype_results.items():
        sig_genes_up = len(de_ct[(de_ct['feature_type'] == 'gene') & (de_ct['pvals_adj'] < 0.05) & (de_ct['logfoldchanges'] > 0)])
        sig_genes_dn = len(de_ct[(de_ct['feature_type'] == 'gene') & (de_ct['pvals_adj'] < 0.05) & (de_ct['logfoldchanges'] < 0)])
        sig_tes_up = len(de_ct[(de_ct['feature_type'] == 'TE') & (de_ct['pvals_adj'] < 0.05) & (de_ct['logfoldchanges'] > 0)])
        sig_tes_dn = len(de_ct[(de_ct['feature_type'] == 'TE') & (de_ct['pvals_adj'] < 0.05) & (de_ct['logfoldchanges'] < 0)])
        
        ct_summary.append({
            'Cell Type': ct,
            'Genes Up': sig_genes_up,
            'Genes Down': sig_genes_dn,
            'TEs Up': sig_tes_up,
            'TEs Down': sig_tes_dn
        })
    
    ct_df = pd.DataFrame(ct_summary)
    ct_df.to_csv('differential_results/celltype_summary.csv', index=False)
    
    print("\nCell-type summary:")
    print(ct_df.to_string(index=False))

print(f"\n{'='*70}")
print("✓ DIFFERENTIAL ANALYSIS COMPLETE!")
print(f"{'='*70}")
print("\nResults saved to: differential_results/")
print("  - genes_AD_vs_Control_all_cells.csv")
print("  - TEs_AD_vs_Control_all_cells.csv")
print("  - <celltype>_AD_vs_Control.csv (for each cell type)")
print("  - celltype_summary.csv")
print("  - figures/volcano_plots.png")
