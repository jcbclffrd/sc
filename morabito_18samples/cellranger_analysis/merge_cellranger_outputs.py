#!/usr/bin/env python3
"""
Merge CellRanger outputs for 18 samples into a single AnnData object.

This script:
1. Loads all 18 CellRanger filtered_feature_bc_matrix.h5 files
2. Merges them into a single AnnData
3. Adds Morabito's cell metadata
4. Filters to the same 60,328 cells used in STAR+scTE analysis
5. Saves merged_cellranger_18samples.h5ad

Output can be directly compared with:
- merged_18samples.h5ad (STAR+scTE genes+TEs)
- merged_18samples_genes.h5ad (STAR+scTE genes only)
- merged_soloTE_18samples.h5ad (soloTE TEs only)
"""

import os
import pandas as pd
import scanpy as sc
import numpy as np

print("=" * 70)
print("MERGING CELLRANGER OUTPUTS: 18 Samples")
print("=" * 70)

# Load sample mapping
mapping = pd.read_csv('../sample_mapping.csv')
print(f"\nLoaded {len(mapping)} sample mappings")

# Load Morabito's metadata
print("\nLoading Morabito's cell metadata...")
meta = pd.read_csv('../../analysis/morabito_reference/GSE174367_snRNA-seq_cell_meta.csv.gz')
print(f"  {len(meta):,} cells")

# Create barcode + sample to metadata mapping
meta['barcode_sample'] = meta['barcode'] + '_' + meta['Sample']
meta_dict = meta.set_index('barcode_sample').to_dict('index')
print(f"  Created lookup for {len(meta_dict):,} cell barcodes")

# Load and merge all samples
adatas = []
output_dir = 'cellranger_output'

print(f"\n{'='*70}")
print("LOADING CELLRANGER OUTPUTS:")
print(f"{'='*70}")

for idx, row in mapping.iterrows():
    sample_id = row['sample_id']
    srr = row['srr']
    
    h5_file = f"{output_dir}/{srr}/outs/filtered_feature_bc_matrix.h5"
    
    if not os.path.exists(h5_file):
        print(f"\n[{idx+1}/{len(mapping)}] {sample_id} ({srr})")
        print(f"  ✗ File not found: {h5_file}")
        continue
    
    print(f"\n[{idx+1}/{len(mapping)}] {sample_id} ({srr})")
    
    # Load CellRanger output
    adata = sc.read_10x_h5(h5_file)
    
    # Add sample ID to barcodes (CellRanger uses -1 suffix)
    adata.obs_names = [bc.replace('-1', f'_{sample_id}') for bc in adata.obs_names]
    
    # Add sample metadata
    adata.obs['sample_id'] = sample_id
    adata.obs['srr'] = srr
    
    print(f"  Cells: {adata.n_obs:>6,}")
    print(f"  Genes: {adata.n_vars:>6,}")
    
    adatas.append(adata)

if len(adatas) == 0:
    print("\n✗ No CellRanger outputs found!")
    exit(1)

print(f"\n{'='*70}")
print("CONCATENATING SAMPLES:")
print(f"{'='*70}")

# Merge all samples
adata = sc.concat(adatas, join='outer', fill_value=0)
print(f"\nMerged shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Filter to Morabito's cells
print(f"\n{'='*70}")
print("FILTERING TO MORABITO'S CELLS:")
print(f"{'='*70}")

adata.obs['in_morabito'] = adata.obs_names.isin(meta_dict.keys())
n_matched = adata.obs['in_morabito'].sum()
print(f"\nMatched cells: {n_matched:,} / {adata.n_obs:,} ({n_matched/adata.n_obs*100:.1f}%)")

adata = adata[adata.obs['in_morabito']].copy()
print(f"After filtering: {adata.n_obs:,} cells")

# Add Morabito's metadata
print(f"\n{'='*70}")
print("ADDING METADATA:")
print(f"{'='*70}")

for col in ['cell_type', 'diagnosis', 'Age', 'Sex', 'PMI', 'Cognitive_Status']:
    if col in meta.columns:
        adata.obs[col] = adata.obs_names.map(lambda x: meta_dict.get(x, {}).get(col, None))
        print(f"  Added: {col}")

# Calculate QC metrics
print(f"\n{'='*70}")
print("CALCULATING QC METRICS:")
print(f"{'='*70}")

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print(f"\nCell type distribution:")
print(adata.obs['cell_type'].value_counts())

print(f"\nDiagnosis distribution:")
print(adata.obs['diagnosis'].value_counts())

print(f"\nQC statistics:")
print(f"  Mean UMIs/cell: {adata.obs['total_counts'].mean():.0f}")
print(f"  Median UMIs/cell: {adata.obs['total_counts'].median():.0f}")
print(f"  Mean genes/cell: {adata.obs['n_genes_by_counts'].mean():.0f}")
print(f"  Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
print(f"  Mean % MT: {adata.obs['pct_counts_mt'].mean():.2f}%")

# Save merged data
output_file = 'merged_cellranger_18samples.h5ad'
print(f"\n{'='*70}")
print(f"SAVING: {output_file}")
print(f"{'='*70}")

adata.write(output_file)
file_size_mb = os.path.getsize(output_file) / (1024**2)
print(f"✓ Saved: {file_size_mb:.1f} MB")

print(f"\n{'='*70}")
print("✓ MERGE COMPLETE")
print(f"{'='*70}")
print(f"\nOutput: {output_file}")
print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"\nNext step: Compare with STAR results")
print(f"  python3 compare_with_star.py")
print("")
