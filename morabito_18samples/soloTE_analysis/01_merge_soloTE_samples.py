#!/usr/bin/env python3
"""
Merge SoloTE output matrices for 18 samples into a single AnnData object.

This script:
1. Loads all 18 SoloTE output matrices
2. Merges them into a single AnnData
3. Adds Morabito's cell metadata
4. Filters to the same 60,328 cells used in scTE analysis
5. Saves merged_soloTE_18samples.h5ad

Similar to ../04_merge_and_qc.py but for SoloTE output format.
"""

import os
import pandas as pd
import scanpy as sc
import numpy as np
import scipy.sparse as sp
from scipy.io import mmread

print("=" * 70)
print("MERGING SOLOTE OUTPUTS: 18 Samples")
print("=" * 70)

# Load sample mapping
mapping = pd.read_csv('../sample_mapping.csv')
print(f"\nLoaded {len(mapping)} sample mappings")

# Load Morabito's metadata
print("\nLoading Morabito's cell metadata...")
meta = pd.read_csv('../../analysis/morabito_reference/GSE174367_snRNA-seq_cell_meta.csv.gz')
print(f"  {len(meta):,} cells")

# Create barcode + sample to metadata mapping
meta['barcode_clean'] = meta['Barcode'].str.split('-').str[0]
meta['barcode_sample_key'] = meta['barcode_clean'] + '_' + meta['SampleID']

print(f"\n{'='*70}")
print("LOADING SOLOTE SAMPLES:")
print(f"{'='*70}")

adatas = []

for idx, row in mapping.iterrows():
    srr = row['srr']
    sample_id = row['sample_id']
    
    # SoloTE creates: soloTE_output/SRR/SRR_SoloTE_output/SRR_subfamilytes_MATRIX/
    output_dir = f"soloTE_output/{srr}/{srr}_SoloTE_output/{srr}_subfamilytes_MATRIX"
    
    # SoloTE outputs (standard 10x format)
    mtx_file = f"{output_dir}/matrix.mtx"
    barcode_file = f"{output_dir}/barcodes.tsv"
    feature_file = f"{output_dir}/features.tsv"
    
    if not os.path.exists(mtx_file):
        print(f"✗ {sample_id:12s} ({srr}): NOT FOUND at {mtx_file}")
        continue
    
    # Load matrix (market matrix format)
    counts = mmread(mtx_file).T.tocsr()  # Transpose to cells x features
    
    # Load barcodes
    barcodes = pd.read_csv(barcode_file, header=None, sep='\t')[0].values
    
    # Load features (TE subfamily names)
    features = pd.read_csv(feature_file, header=None, sep='\t')
    if features.shape[1] >= 2:
        feature_names = features[0].values
        feature_info = features[1].values if features.shape[1] > 1 else None
    else:
        feature_names = features[0].values
        feature_info = None
    
    # Create AnnData with ALL features first
    adata = sc.AnnData(
        X=counts,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=feature_names)
    )
    
    # Filter to TEs only (soloTE prefixes TEs with "SoloTE|")
    te_mask = [name.startswith('SoloTE|') for name in adata.var_names]
    adata = adata[:, te_mask].copy()
    
    # Strip the "SoloTE|" prefix to get actual TE names
    adata.var_names = [name.replace('SoloTE|', '') for name in adata.var_names]
    adata.var_names_make_unique()
    
    # Note: soloTE features are just TE names (e.g., "L1HS", "AluY")
    # We don't have detailed family/class annotations like scTE
    # Mark all as TEs for consistency
    adata.var['feature_type'] = 'TE'
    
    # Add sample info
    adata.obs['sample_id'] = sample_id
    adata.obs['srr'] = srr
    
    print(f"✓ {sample_id:12s} ({srr}): {adata.n_obs:>7,} cells × {adata.n_vars:>5,} TEs")
    adatas.append(adata)

print(f"\nLoaded {len(adatas)}/{len(mapping)} samples")

if len(adatas) == 0:
    print("\n✗ No samples loaded! Run SoloTE first.")
    exit(1)

# Concatenate all samples
print(f"\n{'='*70}")
print("MERGING SAMPLES:")
print(f"{'='*70}")

adata = sc.concat(adatas, join='outer', index_unique=None, merge='same')
print(f"Merged shape: {adata.n_obs:,} cells × {adata.n_vars:,} TEs")

# Add Morabito's metadata
print("\nAdding Morabito's cell metadata...")

# Clean barcodes for matching
adata.obs['barcode_clean'] = [bc.split('-')[0] if '-' in bc else bc for bc in adata.obs_names]
adata.obs['barcode_sample_key'] = adata.obs['barcode_clean'] + '_' + adata.obs['sample_id']

# Prepare metadata for merging
meta_subset = meta[['barcode_sample_key', 'Diagnosis', 'Cell.Type', 'Batch', 'Age', 'Sex', 'PMI', 'Tangle.Stage', 'Plaque.Stage', 'RIN']].copy()
meta_subset.columns = ['barcode_sample_key', 'diagnosis', 'cell_type', 'batch', 'age', 'sex', 'pmi', 'tangle_stage', 'plaque_stage', 'rin']

# Save columns to keep
keep_cols = ['sample_id', 'srr', 'barcode_clean', 'barcode_sample_key']

# Merge metadata
adata.obs = adata.obs[keep_cols].merge(meta_subset, on='barcode_sample_key', how='left')

n_matched = adata.obs['diagnosis'].notna().sum()
match_pct = 100 * n_matched / adata.n_obs
print(f"  Matched {n_matched:,}/{adata.n_obs:,} cells ({match_pct:.1f}%)")

# Filter to cells with metadata (same as scTE analysis)
print(f"\n{'='*70}")
print("FILTERING TO MORABITO'S CELLS:")
print(f"{'='*70}")

cells_with_meta = adata.obs['diagnosis'].notna()
adata = adata[cells_with_meta].copy()
print(f"Filtered to {adata.n_obs:,} cells with metadata")

# Show metadata summary
print(f"\n{'='*70}")
print("METADATA SUMMARY:")
print(f"{'='*70}")

print("\nDiagnosis distribution:")
print(adata.obs['diagnosis'].value_counts())

print("\nCell type distribution:")
print(adata.obs['cell_type'].value_counts())

# Calculate QC metrics
print(f"\n{'='*70}")
print("CALCULATING QC METRICS:")
print(f"{'='*70}")

# UMIs per cell
adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()

# TEs detected per cell
adata.obs['n_tes'] = np.array((adata.X > 0).sum(axis=1)).flatten()

# Note: soloTE subfamilytes output has TE subfamily names only, not family annotations
# We could parse family from subfamily names if needed, but skipping for now
print(f"\nTotal TE subfamilies: {adata.n_vars:,}")

print(f"\nQC statistics:")
print(f"  Mean UMIs/cell: {adata.obs['n_counts'].mean():.0f}")
print(f"  Median UMIs/cell: {adata.obs['n_counts'].median():.0f}")
print(f"  Mean TEs/cell: {adata.obs['n_tes'].mean():.0f}")
print(f"  Median TEs/cell: {adata.obs['n_tes'].median():.0f}")

# Save merged data
output_file = 'merged_soloTE_18samples.h5ad'
print(f"\n{'='*70}")
print(f"SAVING MERGED DATA:")
print(f"{'='*70}")
adata.write_h5ad(output_file)
print(f"✓ Saved: {output_file}")
size_mb = os.path.getsize(output_file) / (1024 * 1024)
print(f"  Size: {size_mb:.1f} MB")

print(f"{'='*70}")
print("✓ MERGING COMPLETE!")
print(f"{'='*70}")
