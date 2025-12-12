#!/usr/bin/env python3
"""
Compare scTE vs SoloTE: TE quantification and differential expression

Creates comprehensive comparison CSV files showing:
1. TE detection comparison (which TEs found by each method)
2. Expression level correlation (counts per TE)
3. Differential expression comparison (AD vs Control)
4. Method concordance analysis

Output CSVs with joined data using TE names as keys.
"""

import os
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

print("=" * 70)
print("COMPARING scTE vs SoloTE")
print("=" * 70)

# Create output directory
os.makedirs('comparison', exist_ok=True)

# =========================================================================
# PART 1: LOAD DATA
# =========================================================================

print("\n" + "=" * 70)
print("LOADING DATA")
print("=" * 70)

# Load scTE data
print("\nLoading scTE data...")
scte_data = sc.read_h5ad('../merged_18samples.h5ad')
print(f"  scTE: {scte_data.n_obs:,} cells × {scte_data.n_vars:,} features")

# Filter to TEs only
scte_tes = scte_data[:, scte_data.var['feature_type'] == 'TE'].copy()
print(f"  scTE TEs: {scte_tes.n_vars:,} TEs")

# Load SoloTE data
print("\nLoading SoloTE data...")
solote_data = sc.read_h5ad('merged_soloTE_18samples.h5ad')
print(f"  SoloTE: {solote_data.n_obs:,} cells × {solote_data.n_vars:,} TEs")

# =========================================================================
# PART 2: TE DETECTION COMPARISON
# =========================================================================

print("\n" + "=" * 70)
print("TE DETECTION COMPARISON")
print("=" * 70)

# Get TE names from each method
scte_te_names = set(scte_tes.var_names)
# For soloTE, the var_names ARE the TE subfamily names
solote_te_names = set(solote_data.var_names)

print(f"\nscTE TEs detected: {len(scte_te_names):,}")
print(f"SoloTE TEs detected: {len(solote_te_names):,}")

# Find overlap
both = scte_te_names & solote_te_names
only_scte = scte_te_names - solote_te_names
only_solote = solote_te_names - scte_te_names

print(f"\nDetected by both: {len(both):,} TEs")
print(f"Only scTE: {len(only_scte):,} TEs")
print(f"Only SoloTE: {len(only_solote):,} TEs")

# Create detection comparison table
detection_df = []

# TEs in both
for te in both:
    detection_df.append({
        'TE_name': te,
        'detected_scTE': True,
        'detected_soloTE': True,
        'detection_method': 'Both'
    })

# Only scTE
for te in only_scte:
    detection_df.append({
        'TE_name': te,
        'detected_scTE': True,
        'detected_soloTE': False,
        'detection_method': 'scTE_only'
    })

# Only SoloTE
for te in only_solote:
    detection_df.append({
        'TE_name': te,
        'detected_scTE': False,
        'detected_soloTE': True,
        'detection_method': 'SoloTE_only'
    })

detection_df = pd.DataFrame(detection_df)
detection_df.to_csv('comparison/01_TE_detection_comparison.csv', index=False)
print(f"\n✓ Saved: comparison/01_TE_detection_comparison.csv")

# =========================================================================
# PART 3: EXPRESSION LEVEL COMPARISON (for TEs in both)
# =========================================================================

print("\n" + "=" * 70)
print("EXPRESSION LEVEL COMPARISON")
print("=" * 70)

print(f"\nComparing expression for {len(both):,} TEs detected by both methods...")

# Calculate mean expression per TE
scte_means = {}
for te in both:
    if te in scte_tes.var_names:
        te_data = scte_tes[:, te].X
        if hasattr(te_data, 'toarray'):
            te_data = te_data.toarray().flatten()
        else:
            te_data = np.array(te_data).flatten()
        scte_means[te] = te_data.mean()

solote_means = {}
for te in both:
    # For soloTE, var_names are the TE names
    if te in solote_data.var_names:
        te_data = solote_data[:, te].X
        if hasattr(te_data, 'toarray'):
            te_data = te_data.toarray().flatten()
        else:
            te_data = np.array(te_data).flatten()
        solote_means[te] = te_data.mean()

# Create comparison dataframe
expression_comparison = []
for te in both:
    if te in scte_means and te in solote_means:
        expression_comparison.append({
            'TE_name': te,
            'scTE_mean_expression': scte_means[te],
            'soloTE_mean_expression': solote_means[te],
            'log2_scTE': np.log2(scte_means[te] + 1),
            'log2_soloTE': np.log2(solote_means[te] + 1),
            'fold_difference': solote_means[te] / (scte_means[te] + 1e-10)
        })

expression_df = pd.DataFrame(expression_comparison)

# Calculate correlation
if len(expression_df) > 0:
    pearson_r, pearson_p = pearsonr(expression_df['log2_scTE'], expression_df['log2_soloTE'])
    spearman_r, spearman_p = spearmanr(expression_df['log2_scTE'], expression_df['log2_soloTE'])
    
    print(f"\nCorrelation (log2 expression):")
    print(f"  Pearson r = {pearson_r:.3f} (p = {pearson_p:.2e})")
    print(f"  Spearman r = {spearman_r:.3f} (p = {spearman_p:.2e})")
    
    expression_df['pearson_r'] = pearson_r
    expression_df['spearman_r'] = spearman_r

expression_df.to_csv('comparison/02_expression_comparison.csv', index=False)
print(f"\n✓ Saved: comparison/02_expression_comparison.csv")

# =========================================================================
# PART 4: DIFFERENTIAL EXPRESSION COMPARISON
# =========================================================================

print("\n" + "=" * 70)
print("DIFFERENTIAL EXPRESSION COMPARISON")
print("=" * 70)

# Load DE results
print("\nLoading differential expression results...")

scte_de = pd.read_csv('../differential_results/TEs_AD_vs_Control_all_cells.csv')
print(f"  scTE DE: {len(scte_de):,} TEs")

solote_de = pd.read_csv('differential_results_soloTE/TEs_AD_vs_Control_all_cells.csv')
print(f"  SoloTE DE: {len(solote_de):,} TEs")

# Merge on TE name
scte_de_subset = scte_de[['names', 'logfoldchanges', 'pvals', 'pvals_adj']].copy()
scte_de_subset.columns = ['TE_name', 'scTE_logFC', 'scTE_pval', 'scTE_padj']

solote_de_subset = solote_de[['names', 'logfoldchanges', 'pvals', 'pvals_adj', 'TE_family', 'TE_class']].copy()
solote_de_subset.columns = ['TE_name', 'soloTE_logFC', 'soloTE_pval', 'soloTE_padj', 'TE_family', 'TE_class']

# Merge
de_comparison = scte_de_subset.merge(solote_de_subset, on='TE_name', how='outer', indicator=True)

# Add significance flags
de_comparison['scTE_significant'] = de_comparison['scTE_padj'] < 0.05
de_comparison['soloTE_significant'] = de_comparison['soloTE_padj'] < 0.05

# Add direction
de_comparison['scTE_direction'] = de_comparison['scTE_logFC'].apply(
    lambda x: 'Up' if x > 0 else 'Down' if x < 0 else 'None'
)
de_comparison['soloTE_direction'] = de_comparison['soloTE_logFC'].apply(
    lambda x: 'Up' if x > 0 else 'Down' if x < 0 else 'None'
)

# Add agreement category
def get_agreement(row):
    if pd.isna(row['scTE_logFC']) or pd.isna(row['soloTE_logFC']):
        return 'Missing_data'
    
    scte_sig = row['scTE_significant']
    solote_sig = row['soloTE_significant']
    
    if scte_sig and solote_sig:
        if row['scTE_direction'] == row['soloTE_direction']:
            return 'Both_significant_same_direction'
        else:
            return 'Both_significant_opposite_direction'
    elif scte_sig and not solote_sig:
        return 'Only_scTE_significant'
    elif not scte_sig and solote_sig:
        return 'Only_soloTE_significant'
    else:
        return 'Neither_significant'

de_comparison['agreement'] = de_comparison.apply(get_agreement, axis=1)

# Calculate correlation for TEs tested by both
both_tested = de_comparison[
    de_comparison['_merge'] == 'both'
].copy()

if len(both_tested) > 10:
    valid = both_tested[
        both_tested['scTE_logFC'].notna() & 
        both_tested['soloTE_logFC'].notna()
    ]
    
    if len(valid) > 10:
        pearson_r, pearson_p = pearsonr(valid['scTE_logFC'], valid['soloTE_logFC'])
        spearman_r, spearman_p = spearmanr(valid['scTE_logFC'], valid['soloTE_logFC'])
        
        print(f"\nCorrelation (log fold changes):")
        print(f"  Pearson r = {pearson_r:.3f} (p = {pearson_p:.2e})")
        print(f"  Spearman r = {spearman_r:.3f} (p = {spearman_p:.2e})")
        
        de_comparison['logFC_pearson_r'] = pearson_r
        de_comparison['logFC_spearman_r'] = spearman_r

# Summary statistics
print(f"\nAgreement summary:")
for agreement, count in de_comparison['agreement'].value_counts().items():
    print(f"  {agreement}: {count}")

de_comparison.to_csv('comparison/03_differential_expression_comparison.csv', index=False)
print(f"\n✓ Saved: comparison/03_differential_expression_comparison.csv")

# =========================================================================
# PART 5: CREATE SUMMARY REPORT
# =========================================================================

print("\n" + "=" * 70)
print("CREATING SUMMARY REPORT")
print("=" * 70)

summary_stats = {
    'Metric': [],
    'scTE': [],
    'SoloTE': [],
    'Overlap': []
}

# Detection
summary_stats['Metric'].append('TEs detected')
summary_stats['scTE'].append(len(scte_te_names))
summary_stats['SoloTE'].append(len(solote_te_names))
summary_stats['Overlap'].append(len(both))

# Significant in DE
scte_sig = scte_de[scte_de['pvals_adj'] < 0.05]
solote_sig = solote_de[solote_de['pvals_adj'] < 0.05]

summary_stats['Metric'].append('Significant TEs (AD vs Control)')
summary_stats['scTE'].append(len(scte_sig))
summary_stats['SoloTE'].append(len(solote_sig))
summary_stats['Overlap'].append(len(de_comparison[
    (de_comparison['scTE_significant']) & (de_comparison['soloTE_significant'])
]))

# Upregulated
scte_up = scte_sig[scte_sig['logfoldchanges'] > 0]
solote_up = solote_sig[solote_sig['logfoldchanges'] > 0]

summary_stats['Metric'].append('Upregulated in AD')
summary_stats['scTE'].append(len(scte_up))
summary_stats['SoloTE'].append(len(solote_up))
summary_stats['Overlap'].append(len(de_comparison[
    (de_comparison['scTE_significant']) & 
    (de_comparison['soloTE_significant']) &
    (de_comparison['scTE_direction'] == 'Up') &
    (de_comparison['soloTE_direction'] == 'Up')
]))

# Downregulated
scte_down = scte_sig[scte_sig['logfoldchanges'] < 0]
solote_down = solote_sig[solote_sig['logfoldchanges'] < 0]

summary_stats['Metric'].append('Downregulated in AD')
summary_stats['scTE'].append(len(scte_down))
summary_stats['SoloTE'].append(len(solote_down))
summary_stats['Overlap'].append(len(de_comparison[
    (de_comparison['scTE_significant']) & 
    (de_comparison['soloTE_significant']) &
    (de_comparison['scTE_direction'] == 'Down') &
    (de_comparison['soloTE_direction'] == 'Down')
]))

summary_df = pd.DataFrame(summary_stats)
summary_df.to_csv('comparison/00_summary_statistics.csv', index=False)
print(f"\n✓ Saved: comparison/00_summary_statistics.csv")

print("\nSummary:")
print(summary_df.to_string(index=False))

print("\n" + "=" * 70)
print("✓ COMPARISON COMPLETE!")
print("=" * 70)
print("\nGenerated comparison files:")
print("  - comparison/00_summary_statistics.csv")
print("  - comparison/01_TE_detection_comparison.csv")
print("  - comparison/02_expression_comparison.csv")
print("  - comparison/03_differential_expression_comparison.csv")
print("\nThese files contain joined data using TE names as keys,")
print("showing counts, DE statistics, and method concordance.")
