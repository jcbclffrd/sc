#!/usr/bin/env python3
"""
Differential expression analysis for SoloTE data: AD vs Control

Performs the same analysis as ../05_differential_analysis.py but for SoloTE TEs.
"""

import os
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

print("=" * 70)
print("SOLOTE DIFFERENTIAL EXPRESSION: AD vs CONTROL")
print("=" * 70)

# Load merged SoloTE data
print("\nLoading merged SoloTE data...")
adata = sc.read_h5ad('merged_soloTE_18samples.h5ad')
print(f"  Shape: {adata.n_obs:,} cells × {adata.n_vars:,} TEs")

# Show diagnosis and cell type counts
print("\nDiagnosis counts:")
print(adata.obs['diagnosis'].value_counts())

print("\nCell type counts:")
print(adata.obs['cell_type'].value_counts())

# Filter cells
print("\nFiltering cells...")
cells_with_data = (adata.obs['diagnosis'].notna()) & (adata.obs['cell_type'].notna())
adata = adata[cells_with_data].copy()
print(f"Filtered to: {adata.n_obs:,} cells with diagnosis and cell type")

# Normalize
print("\nNormalizing data...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Create output directory
os.makedirs('differential_results_soloTE', exist_ok=True)

print(f"\n{'='*70}")
print("1. OVERALL DIFFERENTIAL ANALYSIS (ALL CELL TYPES)")
print(f"{'='*70}")

# Run differential expression
adata.obs['diagnosis_key'] = adata.obs['diagnosis']
sc.tl.rank_genes_groups(
    adata, 
    groupby='diagnosis_key', 
    groups=['AD'], 
    reference='Control',
    method='wilcoxon',
    key_added='de_diagnosis'
)

# Extract results
result = adata.uns['de_diagnosis']
groups = result['names'].dtype.names

de_results = []
for group in groups:
    group_df = pd.DataFrame({
        'names': result['names'][group],
        'scores': result['scores'][group],
        'logfoldchanges': result['logfoldchanges'][group],
        'pvals': result['pvals'][group],
        'pvals_adj': result['pvals_adj'][group]
    })
    
    # Add TE annotations
    group_df = group_df.merge(
        adata.var[['TE_subfamily', 'TE_family', 'TE_class']], 
        left_on='names', 
        right_index=True, 
        how='left'
    )
    
    de_results.append(group_df)
    
    # Save
    output_file = f'differential_results_soloTE/TEs_AD_vs_Control_all_cells.csv'
    group_df.to_csv(output_file, index=False)
    print(f"\nSaved: {output_file}")

# Show top results
print(f"\nTop 20 upregulated TEs in AD:")
top_up = de_results[0].sort_values('pvals_adj').head(20)
for idx, row in top_up.iterrows():
    print(f"  {row['TE_subfamily']:20s}  {row['TE_family']:10s}  logFC={row['logfoldchanges']:6.3f}  padj={row['pvals_adj']:.2e}")

print(f"\n{'='*70}")
print("2. CELL-TYPE SPECIFIC DIFFERENTIAL ANALYSIS")
print(f"{'='*70}")

# Run DE for each cell type
cell_types = adata.obs['cell_type'].unique()
cell_types = [ct for ct in cell_types if pd.notna(ct)]

print(f"\nAnalyzing {len(cell_types)} cell types: {', '.join(cell_types)}")

celltype_summary = []

for ct in cell_types:
    print(f"\n{ct}:")
    
    # Subset to cell type
    adata_ct = adata[adata.obs['cell_type'] == ct].copy()
    
    # Check we have both conditions
    n_ad = (adata_ct.obs['diagnosis'] == 'AD').sum()
    n_control = (adata_ct.obs['diagnosis'] == 'Control').sum()
    
    print(f"  Cells: AD={n_ad}, Control={n_control}")
    
    if n_ad < 10 or n_control < 10:
        print(f"  Skipping (too few cells)")
        continue
    
    # Run DE
    sc.tl.rank_genes_groups(
        adata_ct,
        groupby='diagnosis_key',
        groups=['AD'],
        reference='Control',
        method='wilcoxon',
        key_added='de_diagnosis'
    )
    
    # Extract results
    result = adata_ct.uns['de_diagnosis']
    ct_df = pd.DataFrame({
        'names': result['names']['AD'],
        'scores': result['scores']['AD'],
        'logfoldchanges': result['logfoldchanges']['AD'],
        'pvals': result['pvals']['AD'],
        'pvals_adj': result['pvals_adj']['AD']
    })
    
    # Add TE annotations
    ct_df = ct_df.merge(
        adata_ct.var[['TE_subfamily', 'TE_family', 'TE_class']], 
        left_on='names',
        right_index=True,
        how='left'
    )
    
    # Save
    output_file = f'differential_results_soloTE/{ct}_AD_vs_Control.csv'
    ct_df.to_csv(output_file, index=False)
    
    # Count significant
    sig = ct_df[ct_df['pvals_adj'] < 0.05]
    n_up = (sig['logfoldchanges'] > 0).sum()
    n_down = (sig['logfoldchanges'] < 0).sum()
    
    print(f"  Significant (padj < 0.05): {len(sig)} TEs ({n_up} up, {n_down} down)")
    
    celltype_summary.append({
        'Cell Type': ct,
        'TEs Up': n_up,
        'TEs Down': n_down,
        'Total Sig': len(sig)
    })

# Save summary
summary_df = pd.DataFrame(celltype_summary)
summary_df.to_csv('differential_results_soloTE/celltype_summary.csv', index=False)
print(f"\n✓ Saved: differential_results_soloTE/celltype_summary.csv")

print(f"\n{'='*70}")
print("✓ SOLOTE DIFFERENTIAL ANALYSIS COMPLETE!")
print(f"{'='*70}")
print(f"\nResults saved to: differential_results_soloTE/")
