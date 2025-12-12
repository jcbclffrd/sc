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
    
    output_dir = f"soloTE_output/{srr}"
    
    # SoloTE outputs
    mtx_file = f"{output_dir}/{srr}_TE_counts.mtx"
    barcode_file = f"{output_dir}/{srr}_barcodes.tsv"
    feature_file = f"{output_dir}/{srr}_features.tsv"
    
    if not os.path.exists(mtx_file):
        print(f"✗ {sample_id:12s} ({srr}): NOT FOUND")
        continue
    
    # Load matrix
    counts = mmread(mtx_file).T.tocsr()  # Transpose to cells x features
    
    # Load barcodes
    barcodes = pd.read_csv(barcode_file, header=None, sep='\t')[0].values
    
    # Load features
    features = pd.read_csv(feature_file, header=None, sep='\t')
    if features.shape[1] >= 2:
        feature_names = features[0].values
        feature_info = features[1].values if features.shape[1] > 1 else None
    else:
        feature_names = features[0].values
        feature_info = None
    
    # Create AnnData
    adata = sc.AnnData(
        X=counts,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=feature_names)
    )
    
    # Add feature info if available
    if feature_info is not None:
        adata.var['te_info'] = feature_info
    
    # Extract TE subfamily, family, class from feature names
    # Format: chr1|11505|11675|L1MC5a:L1:LINE|25.1|-
    te_annotations = []
    for feat in feature_names:
        parts = feat.split('|')
        if len(parts) >= 4:
            te_part = parts[3]  # L1MC5a:L1:LINE
            te_fields = te_part.split(':')
            if len(te_fields) >= 3:
                subfamily = te_fields[0]
                family = te_fields[1]
                te_class = te_fields[2]
            else:
                subfamily = te_part
                family = 'Unknown'
                te_class = 'Unknown'
        else:
            subfamily = feat
            family = 'Unknown'
            te_class = 'Unknown'
        
        te_annotations.append({
            'TE_subfamily': subfamily,
            'TE_family': family,
            'TE_class': te_class,
            'feature_type': 'TE'  # All SoloTE features are TEs
        })
    
    adata.var = pd.DataFrame(te_annotations, index=feature_names)
    
    # Add sample info
    adata.obs['sample_id'] = sample_id
    adata.obs['srr'] = srr
    
    print(f"✓ {sample_id:12s} ({srr}): {adata.n_obs:>5,} cells × {adata.n_vars:>6,} TEs")
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

# TE family counts
te_families = adata.var['TE_family'].value_counts()
print(f"\nTE families detected: {len(te_families)}")
print("Top 10 TE families:")
for family, count in te_families.head(10).items():
    print(f"  {family:15s}: {count:>5,} subfamilies")

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
