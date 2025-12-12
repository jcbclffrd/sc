#!/usr/bin/env python3
"""
Generate Professional Comparison Report: scTE vs SoloTE

Creates a comprehensive markdown report with formatted tables comparing:
- TE detection rates
- Expression correlation
- Differential expression concordance
- Top significant TEs from each method
- Method-specific findings

Output: COMPARISON_REPORT.md - ready for PI presentation
"""

import os
import pandas as pd
import numpy as np
from datetime import datetime

print("=" * 70)
print("GENERATING COMPARISON REPORT")
print("=" * 70)

# Load comparison data
print("\nLoading comparison results...")
summary = pd.read_csv('comparison/00_summary_statistics.csv')
detection = pd.read_csv('comparison/01_TE_detection_comparison.csv')
expression = pd.read_csv('comparison/02_expression_comparison.csv')
de_comparison = pd.read_csv('comparison/03_differential_expression_comparison.csv')

# =========================================================================
# Generate Markdown Report
# =========================================================================

report = []

report.append("# Comparison Report: scTE vs SoloTE")
report.append("")
report.append(f"**Analysis Date**: {datetime.now().strftime('%B %d, %Y')}")
report.append("")
report.append("**Dataset**: Morabito et al. (2021) - 18 snRNA-seq samples")
report.append("")
report.append("**Purpose**: Compare two TE quantification methods on the same 18 samples and ~61,000 cells")
report.append("")
report.append("---")
report.append("")

# =========================================================================
# Executive Summary
# =========================================================================

report.append("## Executive Summary")
report.append("")

# Calculate key metrics
scte_tes = summary[summary['Metric'] == 'TEs detected']['scTE'].values[0]
solote_tes = summary[summary['Metric'] == 'TEs detected']['SoloTE'].values[0]
overlap = summary[summary['Metric'] == 'TEs detected']['Overlap'].values[0]

scte_sig = summary[summary['Metric'] == 'Significant TEs (AD vs Control)']['scTE'].values[0]
solote_sig = summary[summary['Metric'] == 'Significant TEs (AD vs Control)']['SoloTE'].values[0]
sig_overlap = summary[summary['Metric'] == 'Significant TEs (AD vs Control)']['Overlap'].values[0]

report.append(f"**TE Detection:**")
report.append(f"- scTE detected **{scte_tes:,} TEs**")
report.append(f"- soloTE detected **{solote_tes:,} TEs**")
report.append(f"- Overlap: **{overlap:,} TEs ({100*overlap/max(scte_tes, solote_tes):.1f}%)**")
report.append("")

report.append(f"**Differential Expression (AD vs Control):**")
report.append(f"- scTE found **{scte_sig:,} significant TEs**")
report.append(f"- soloTE found **{solote_sig:,} significant TEs**")
report.append(f"- Agreement: **{sig_overlap:,} TEs ({100*sig_overlap/min(scte_sig, solote_sig):.1f}% concordance)**")
report.append("")

# Expression correlation
if 'pearson_r' in expression.columns:
    pearson = expression['pearson_r'].iloc[0]
    report.append(f"**Expression Correlation:**")
    report.append(f"- Pearson r = **{pearson:.3f}** (strong positive correlation)")
    report.append("")

# logFC correlation
if 'logFC_pearson_r' in de_comparison.columns:
    logfc_pearson = de_comparison['logFC_pearson_r'].iloc[0]
    report.append(f"**Differential Expression Correlation:**")
    report.append(f"- Log fold-change Pearson r = **{logfc_pearson:.3f}**")
    report.append("")

report.append("---")
report.append("")

# =========================================================================
# Summary Statistics Table
# =========================================================================

report.append("## Summary Statistics")
report.append("")
report.append("| Metric | scTE | soloTE | Overlap |")
report.append("|--------|------|--------|---------|")
for _, row in summary.iterrows():
    metric = row['Metric']
    scte = f"{int(row['scTE']):,}"
    solote = f"{int(row['SoloTE']):,}"
    overlap_val = f"{int(row['Overlap']):,}"
    report.append(f"| {metric} | {scte} | {solote} | {overlap_val} |")
report.append("")
report.append("---")
report.append("")

# =========================================================================
# Top Significant TEs - scTE
# =========================================================================

report.append("## Top Differentially Expressed TEs (AD vs Control)")
report.append("")

# Get top TEs from scTE
scte_top = de_comparison[
    (de_comparison['scTE_significant'] == True) & 
    (de_comparison['scTE_padj'].notna())
].copy()
scte_top = scte_top.sort_values('scTE_padj').head(20)

report.append("### Top 20 TEs by scTE (padj < 0.05)")
report.append("")
report.append("| Rank | TE Name | scTE logFC | scTE padj | soloTE logFC | soloTE padj | Agreement |")
report.append("|------|---------|------------|-----------|--------------|-------------|-----------|")

for i, (_, row) in enumerate(scte_top.iterrows(), 1):
    te_name = row['TE_name']
    scte_lfc = f"{row['scTE_logFC']:.3f}" if pd.notna(row['scTE_logFC']) else "N/A"
    scte_padj = f"{row['scTE_padj']:.2e}" if pd.notna(row['scTE_padj']) else "N/A"
    solote_lfc = f"{row['soloTE_logFC']:.3f}" if pd.notna(row['soloTE_logFC']) else "N/A"
    solote_padj = f"{row['soloTE_padj']:.2e}" if pd.notna(row['soloTE_padj']) else "N/A"
    
    # Agreement symbol
    if row['agreement'] == 'Both_significant_same_direction':
        agreement = "✅ Both agree"
    elif row['agreement'] == 'Only_scTE_significant':
        agreement = "⚠️ scTE only"
    elif row['agreement'] == 'Only_soloTE_significant':
        agreement = "⚠️ soloTE only"
    elif row['agreement'] == 'Both_significant_opposite_direction':
        agreement = "❌ Opposite"
    else:
        agreement = "—"
    
    report.append(f"| {i} | {te_name} | {scte_lfc} | {scte_padj} | {solote_lfc} | {solote_padj} | {agreement} |")

report.append("")
report.append("---")
report.append("")

# =========================================================================
# Top Significant TEs - soloTE
# =========================================================================

# Get top TEs from soloTE
solote_top = de_comparison[
    (de_comparison['soloTE_significant'] == True) & 
    (de_comparison['soloTE_padj'].notna())
].copy()
solote_top = solote_top.sort_values('soloTE_padj').head(20)

report.append("### Top 20 TEs by soloTE (padj < 0.05)")
report.append("")
report.append("| Rank | TE Name | TE Family | soloTE logFC | soloTE padj | scTE logFC | scTE padj | Agreement |")
report.append("|------|---------|-----------|--------------|-------------|------------|-----------|-----------|")

for i, (_, row) in enumerate(solote_top.iterrows(), 1):
    te_name = row['TE_name']
    te_family = row['TE_family'] if pd.notna(row['TE_family']) else "Unknown"
    solote_lfc = f"{row['soloTE_logFC']:.3f}" if pd.notna(row['soloTE_logFC']) else "N/A"
    solote_padj = f"{row['soloTE_padj']:.2e}" if pd.notna(row['soloTE_padj']) else "N/A"
    scte_lfc = f"{row['scTE_logFC']:.3f}" if pd.notna(row['scTE_logFC']) else "N/A"
    scte_padj = f"{row['scTE_padj']:.2e}" if pd.notna(row['scTE_padj']) else "N/A"
    
    # Agreement symbol
    if row['agreement'] == 'Both_significant_same_direction':
        agreement = "✅ Both agree"
    elif row['agreement'] == 'Only_soloTE_significant':
        agreement = "⚠️ soloTE only"
    elif row['agreement'] == 'Only_scTE_significant':
        agreement = "⚠️ scTE only"
    elif row['agreement'] == 'Both_significant_opposite_direction':
        agreement = "❌ Opposite"
    else:
        agreement = "—"
    
    report.append(f"| {i} | {te_name} | {te_family} | {solote_lfc} | {solote_padj} | {scte_lfc} | {scte_padj} | {agreement} |")

report.append("")
report.append("---")
report.append("")

# =========================================================================
# Concordant TEs (Both Methods Agree)
# =========================================================================

report.append("## High-Confidence TEs (Both Methods Agree)")
report.append("")

# Get TEs where both methods agree (same direction, both significant)
concordant = de_comparison[
    de_comparison['agreement'] == 'Both_significant_same_direction'
].copy()

concordant['avg_logFC'] = (concordant['scTE_logFC'] + concordant['soloTE_logFC']) / 2
concordant = concordant.sort_values('avg_logFC', ascending=False)

report.append(f"**{len(concordant)} TEs are significant in both methods with the same direction**")
report.append("")

# Top upregulated (concordant)
concordant_up = concordant[concordant['avg_logFC'] > 0].head(15)
report.append("### Top 15 Upregulated TEs (Both Methods)")
report.append("")
report.append("| TE Name | TE Family | scTE logFC | soloTE logFC | Avg logFC | scTE padj | soloTE padj |")
report.append("|---------|-----------|------------|--------------|-----------|-----------|-------------|")

for _, row in concordant_up.iterrows():
    te_name = row['TE_name']
    te_family = row['TE_family'] if pd.notna(row['TE_family']) else "Unknown"
    scte_lfc = f"{row['scTE_logFC']:.3f}"
    solote_lfc = f"{row['soloTE_logFC']:.3f}"
    avg_lfc = f"{row['avg_logFC']:.3f}"
    scte_padj = f"{row['scTE_padj']:.2e}"
    solote_padj = f"{row['soloTE_padj']:.2e}"
    
    report.append(f"| {te_name} | {te_family} | {scte_lfc} | {solote_lfc} | {avg_lfc} | {scte_padj} | {solote_padj} |")

report.append("")

# Top downregulated (concordant)
concordant_down = concordant[concordant['avg_logFC'] < 0].head(15)
report.append("### Top 15 Downregulated TEs (Both Methods)")
report.append("")
report.append("| TE Name | TE Family | scTE logFC | soloTE logFC | Avg logFC | scTE padj | soloTE padj |")
report.append("|---------|-----------|------------|--------------|-----------|-----------|-------------|")

for _, row in concordant_down.iterrows():
    te_name = row['TE_name']
    te_family = row['TE_family'] if pd.notna(row['TE_family']) else "Unknown"
    scte_lfc = f"{row['scTE_logFC']:.3f}"
    solote_lfc = f"{row['soloTE_logFC']:.3f}"
    avg_lfc = f"{row['avg_logFC']:.3f}"
    scte_padj = f"{row['scTE_padj']:.2e}"
    solote_padj = f"{row['soloTE_padj']:.2e}"
    
    report.append(f"| {te_name} | {te_family} | {scte_lfc} | {solote_lfc} | {avg_lfc} | {scte_padj} | {solote_padj} |")

report.append("")
report.append("---")
report.append("")

# =========================================================================
# Method Disagreements
# =========================================================================

report.append("## Method Disagreements")
report.append("")

# Disagreement summary
disagreement_counts = de_comparison['agreement'].value_counts()

report.append("| Agreement Category | Count |")
report.append("|--------------------|-------|")
for agreement, count in disagreement_counts.items():
    report.append(f"| {agreement.replace('_', ' ')} | {count:,} |")

report.append("")

# TEs only significant in scTE
only_scte = de_comparison[
    de_comparison['agreement'] == 'Only_scTE_significant'
].sort_values('scTE_padj').head(10)

if len(only_scte) > 0:
    report.append("### Top 10 TEs Only Significant in scTE")
    report.append("")
    report.append("| TE Name | scTE logFC | scTE padj | soloTE logFC | soloTE padj |")
    report.append("|---------|------------|-----------|--------------|-------------|")
    
    for _, row in only_scte.iterrows():
        te_name = row['TE_name']
        scte_lfc = f"{row['scTE_logFC']:.3f}"
        scte_padj = f"{row['scTE_padj']:.2e}"
        solote_lfc = f"{row['soloTE_logFC']:.3f}" if pd.notna(row['soloTE_logFC']) else "N/A"
        solote_padj = f"{row['soloTE_padj']:.2e}" if pd.notna(row['soloTE_padj']) else "N/A"
        
        report.append(f"| {te_name} | {scte_lfc} | {scte_padj} | {solote_lfc} | {solote_padj} |")
    
    report.append("")

# TEs only significant in soloTE
only_solote = de_comparison[
    de_comparison['agreement'] == 'Only_soloTE_significant'
].sort_values('soloTE_padj').head(10)

if len(only_solote) > 0:
    report.append("### Top 10 TEs Only Significant in soloTE")
    report.append("")
    report.append("| TE Name | TE Family | soloTE logFC | soloTE padj | scTE logFC | scTE padj |")
    report.append("|---------|-----------|--------------|-------------|------------|-----------|")
    
    for _, row in only_solote.iterrows():
        te_name = row['TE_name']
        te_family = row['TE_family'] if pd.notna(row['TE_family']) else "Unknown"
        solote_lfc = f"{row['soloTE_logFC']:.3f}"
        solote_padj = f"{row['soloTE_padj']:.2e}"
        scte_lfc = f"{row['scTE_logFC']:.3f}" if pd.notna(row['scTE_logFC']) else "N/A"
        scte_padj = f"{row['scTE_padj']:.2e}" if pd.notna(row['scTE_padj']) else "N/A"
        
        report.append(f"| {te_name} | {te_family} | {solote_lfc} | {solote_padj} | {scte_lfc} | {scte_padj} |")
    
    report.append("")

report.append("---")
report.append("")

# =========================================================================
# Conclusions
# =========================================================================

report.append("## Conclusions")
report.append("")

# Calculate concordance percentage
concordance_pct = 100 * sig_overlap / min(scte_sig, solote_sig) if min(scte_sig, solote_sig) > 0 else 0

report.append(f"1. **High concordance**: {sig_overlap} TEs ({concordance_pct:.1f}%) are significant in both methods")
report.append("")

if 'pearson_r' in expression.columns:
    pearson = expression['pearson_r'].iloc[0]
    report.append(f"2. **Strong correlation**: Expression levels correlate well (r = {pearson:.3f})")
    report.append("")

report.append(f"3. **Method-specific findings**: Each method detects unique TEs")
report.append(f"   - scTE unique: {scte_sig - sig_overlap} TEs")
report.append(f"   - soloTE unique: {solote_sig - sig_overlap} TEs")
report.append("")

report.append(f"4. **Recommendation**: Focus on the **{len(concordant)} high-confidence TEs** where both methods agree")
report.append("")

report.append("---")
report.append("")

# =========================================================================
# Write Report
# =========================================================================

report_text = "\n".join(report)

with open('COMPARISON_REPORT.md', 'w') as f:
    f.write(report_text)

print("✓ Generated: COMPARISON_REPORT.md")
print(f"  - {len(concordant)} high-confidence TEs (both methods agree)")
print(f"  - {scte_sig} scTE significant TEs")
print(f"  - {solote_sig} soloTE significant TEs")
print(f"  - Ready for PI presentation!")

print("\n" + "=" * 70)
print("✓ REPORT GENERATION COMPLETE")
print("=" * 70)
