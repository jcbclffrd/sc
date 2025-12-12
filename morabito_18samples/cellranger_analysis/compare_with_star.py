#!/usr/bin/env python3
"""
Compare CellRanger gene quantification with STAR+scTE gene quantification.

This script loads:
1. merged_cellranger_18samples.h5ad (CellRanger genes)
2. ../merged_18samples_genes.h5ad (STAR+scTE genes only)

And generates comparison CSVs:
- gene_detection_comparison.csv (which genes detected by each)
- gene_expression_comparison.csv (expression correlation)
- per_cell_comparison.csv (UMI counts, genes detected per cell)
"""

import os
import pandas as pd
import scanpy as sc
import numpy as np
from scipy.stats import pearsonr, spearmanr

print("=" * 70)
print("COMPARING CELLRANGER vs STAR+scTE GENE QUANTIFICATION")
print("=" * 70)

# Load CellRanger data
print("\nLoading CellRanger data...")
if not os.path.exists('merged_cellranger_18samples.h5ad'):
    print("✗ CellRanger data not found: merged_cellranger_18samples.h5ad")
    print("  Run: python3 merge_cellranger_outputs.py")
    exit(1)

adata_cr = sc.read_h5ad('merged_cellranger_18samples.h5ad')
print(f"  CellRanger: {adata_cr.n_obs:,} cells × {adata_cr.n_vars:,} genes")

# Load STAR data
print("\nLoading STAR+scTE data...")
if not os.path.exists('../merged_18samples_genes.h5ad'):
    print("✗ STAR data not found: ../merged_18samples_genes.h5ad")
    exit(1)

adata_star = sc.read_h5ad('../merged_18samples_genes.h5ad')
print(f"  STAR+scTE: {adata_star.n_obs:,} cells × {adata_star.n_vars:,} genes")

# Find common cells and genes
print(f"\n{'='*70}")
print("FINDING COMMON CELLS AND GENES:")
print(f"{'='*70}")

common_cells = list(set(adata_cr.obs_names) & set(adata_star.obs_names))
common_genes = list(set(adata_cr.var_names) & set(adata_star.var_names))

print(f"\nCommon cells: {len(common_cells):,}")
print(f"Common genes: {len(common_genes):,}")

if len(common_cells) == 0:
    print("\n✗ No common cells found! Check barcode formatting.")
    exit(1)

# Subset to common cells and genes
adata_cr_sub = adata_cr[common_cells, common_genes].copy()
adata_star_sub = adata_star[common_cells, common_genes].copy()

print(f"\nSubset shape: {len(common_cells):,} cells × {len(common_genes):,} genes")

# 1. Gene detection comparison
print(f"\n{'='*70}")
print("1. GENE DETECTION COMPARISON")
print(f"{'='*70}")

# Calculate mean expression per gene
cr_mean = np.array(adata_cr_sub.X.mean(axis=0)).flatten()
star_mean = np.array(adata_star_sub.X.mean(axis=0)).flatten()

# Detection threshold: mean > 0
cr_detected = cr_mean > 0
star_detected = star_mean > 0

gene_comparison = pd.DataFrame({
    'gene_name': common_genes,
    'cellranger_mean': cr_mean,
    'star_mean': star_mean,
    'detected_cellranger': cr_detected,
    'detected_star': star_detected,
    'fold_difference': np.where(star_mean > 0, cr_mean / star_mean, np.nan)
})

# Categorize detection
def categorize_detection(row):
    if row['detected_cellranger'] and row['detected_star']:
        return 'Both'
    elif row['detected_cellranger']:
        return 'CellRanger_only'
    elif row['detected_star']:
        return 'STAR_only'
    else:
        return 'Neither'

gene_comparison['detection_category'] = gene_comparison.apply(categorize_detection, axis=1)

print("\nGene detection categories:")
print(gene_comparison['detection_category'].value_counts())

# Calculate correlation for commonly detected genes
both_detected = gene_comparison['detection_category'] == 'Both'
if both_detected.sum() > 0:
    cr_vals = gene_comparison.loc[both_detected, 'cellranger_mean'].values
    star_vals = gene_comparison.loc[both_detected, 'star_mean'].values
    
    pearson_r, pearson_p = pearsonr(cr_vals, star_vals)
    spearman_r, spearman_p = spearmanr(cr_vals, star_vals)
    
    print(f"\nCorrelation (genes detected by both):")
    print(f"  Pearson r: {pearson_r:.4f} (p={pearson_p:.2e})")
    print(f"  Spearman r: {spearman_r:.4f} (p={spearman_p:.2e})")

gene_comparison.to_csv('gene_expression_comparison.csv', index=False)
print(f"\n✓ Saved: gene_expression_comparison.csv")

# 2. Per-cell comparison
print(f"\n{'='*70}")
print("2. PER-CELL COMPARISON")
print(f"{'='*70}")

cell_comparison = pd.DataFrame({
    'cell_barcode': common_cells,
    'cellranger_umi_counts': np.array(adata_cr_sub.X.sum(axis=1)).flatten(),
    'star_umi_counts': np.array(adata_star_sub.X.sum(axis=1)).flatten(),
    'cellranger_genes_detected': np.array((adata_cr_sub.X > 0).sum(axis=1)).flatten(),
    'star_genes_detected': np.array((adata_star_sub.X > 0).sum(axis=1)).flatten(),
})

# Add metadata
for col in ['cell_type', 'diagnosis', 'sample_id']:
    if col in adata_cr_sub.obs.columns:
        cell_comparison[col] = adata_cr_sub.obs[col].values

print(f"\nPer-cell statistics:")
print(f"  CellRanger mean UMIs: {cell_comparison['cellranger_umi_counts'].mean():.0f}")
print(f"  STAR mean UMIs: {cell_comparison['star_umi_counts'].mean():.0f}")
print(f"  CellRanger mean genes: {cell_comparison['cellranger_genes_detected'].mean():.0f}")
print(f"  STAR mean genes: {cell_comparison['star_genes_detected'].mean():.0f}")

# Calculate correlation
pearson_umi, _ = pearsonr(cell_comparison['cellranger_umi_counts'], 
                          cell_comparison['star_umi_counts'])
pearson_genes, _ = pearsonr(cell_comparison['cellranger_genes_detected'], 
                            cell_comparison['star_genes_detected'])

print(f"\nPer-cell correlation:")
print(f"  UMI counts: r={pearson_umi:.4f}")
print(f"  Genes detected: r={pearson_genes:.4f}")

cell_comparison.to_csv('per_cell_comparison.csv', index=False)
print(f"\n✓ Saved: per_cell_comparison.csv")

# 3. Summary statistics
print(f"\n{'='*70}")
print("3. SUMMARY STATISTICS")
print(f"{'='*70}")

summary = {
    'metric': [
        'Total cells',
        'Total genes',
        'Genes detected by both',
        'Genes only in CellRanger',
        'Genes only in STAR',
        'Mean correlation (gene expression)',
        'Mean correlation (UMI counts per cell)',
        'CellRanger mean UMIs/cell',
        'STAR mean UMIs/cell',
    ],
    'value': [
        len(common_cells),
        len(common_genes),
        (gene_comparison['detection_category'] == 'Both').sum(),
        (gene_comparison['detection_category'] == 'CellRanger_only').sum(),
        (gene_comparison['detection_category'] == 'STAR_only').sum(),
        f"{pearson_r:.4f}" if both_detected.sum() > 0 else "N/A",
        f"{pearson_umi:.4f}",
        f"{cell_comparison['cellranger_umi_counts'].mean():.0f}",
        f"{cell_comparison['star_umi_counts'].mean():.0f}",
    ]
}

summary_df = pd.DataFrame(summary)
summary_df.to_csv('summary_cellranger_vs_star.csv', index=False)

print("\nSummary:")
for i, row in summary_df.iterrows():
    print(f"  {row['metric']}: {row['value']}")

print(f"\n✓ Saved: summary_cellranger_vs_star.csv")

print(f"\n{'='*70}")
print("✓ COMPARISON COMPLETE")
print(f"{'='*70}")
print("\nOutput files:")
print("  - gene_expression_comparison.csv")
print("  - per_cell_comparison.csv")
print("  - summary_cellranger_vs_star.csv")
print("")
