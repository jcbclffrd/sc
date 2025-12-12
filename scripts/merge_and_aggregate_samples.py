#!/usr/bin/env python3
"""
Merge all scTE output samples and create aggregated count table

This script:
1. Loads all 150 .h5ad files from scTE output
2. Merges them into one AnnData object
3. Saves the merged single-cell data
4. Collapses across cells to get total counts per feature
5. Saves aggregated counts as CSV
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.sparse import issparse
import sys

# Paths
SCTE_OUTPUT = Path("/home/jacobc/sc/scTE_output")
ANALYSIS_DIR = Path("/home/jacobc/sc/analysis")
ANALYSIS_DIR.mkdir(exist_ok=True)

MERGED_FILE = ANALYSIS_DIR / "merged_samples.h5ad"
AGGREGATED_CSV = ANALYSIS_DIR / "combinedCells.csv"

print("=" * 60)
print("scTE Sample Merger and Aggregator")
print("=" * 60)
print()

# Find all h5ad files
print("Finding samples...")
h5ad_files = sorted(SCTE_OUTPUT.glob("*/SRR*.h5ad"))
print(f"Found {len(h5ad_files)} samples")
print()

if len(h5ad_files) == 0:
    print("ERROR: No .h5ad files found!")
    sys.exit(1)

# Load and merge all samples
print("Loading and merging samples...")
adatas = []
failed = []

for i, h5ad_file in enumerate(h5ad_files, 1):
    sample_id = h5ad_file.stem
    try:
        adata = sc.read_h5ad(h5ad_file)
        
        # Add sample ID to cell barcodes to make them unique
        adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]
        
        # Add sample metadata
        adata.obs['sample_id'] = sample_id
        adata.obs['batch'] = sample_id
        
        adatas.append(adata)
        
        if i % 10 == 0:
            print(f"  Loaded {i}/{len(h5ad_files)} samples...")
            
    except Exception as e:
        print(f"  WARNING: Failed to load {sample_id}: {e}")
        failed.append(sample_id)

print(f"Successfully loaded {len(adatas)}/{len(h5ad_files)} samples")
if failed:
    print(f"Failed samples: {', '.join(failed)}")
print()

# Concatenate all samples
print("Concatenating samples...")
merged = sc.concat(adatas, join='outer', merge='same')
print(f"Merged shape: {merged.shape} (cells x features)")
print(f"  Total cells: {merged.n_obs:,}")
print(f"  Total features: {merged.n_vars:,}")
print()

# Save merged data
print(f"Saving merged data to: {MERGED_FILE}")
merged.write_h5ad(MERGED_FILE)
print(f"  File size: {MERGED_FILE.stat().st_size / 1e9:.2f} GB")
print()

# Create aggregated count table (sum across all cells)
print("Aggregating counts across all cells...")

# Get the count matrix
X = merged.X

# Sum across cells (axis=0)
if issparse(X):
    total_counts = np.array(X.sum(axis=0)).flatten()
else:
    total_counts = X.sum(axis=0)

# Create dataframe
agg_df = pd.DataFrame({
    'feature': merged.var_names,
    'total_count': total_counts
})

# Sort by count (highest first)
agg_df = agg_df.sort_values('total_count', ascending=False)

# Add some statistics
agg_df['mean_count_per_cell'] = agg_df['total_count'] / merged.n_obs
agg_df['cells_expressing'] = np.array((X > 0).sum(axis=0)).flatten()
agg_df['pct_cells_expressing'] = (agg_df['cells_expressing'] / merged.n_obs) * 100

print(f"Total non-zero counts: {(total_counts > 0).sum():,} / {len(total_counts):,}")
print()

# Show top expressed features
print("Top 10 most expressed features:")
print(agg_df.head(10).to_string(index=False))
print()

# Save aggregated counts
print(f"Saving aggregated counts to: {AGGREGATED_CSV}")
agg_df.to_csv(AGGREGATED_CSV, index=False)
print(f"  Saved {len(agg_df):,} features")
print()

print("=" * 60)
print("COMPLETE!")
print("=" * 60)
print()
print("Output files:")
print(f"  1. Merged single-cell data: {MERGED_FILE}")
print(f"  2. Aggregated counts: {AGGREGATED_CSV}")
print()
print("Summary:")
print(f"  Total cells: {merged.n_obs:,}")
print(f"  Total features: {merged.n_vars:,}")
print(f"  Total samples: {len(adatas)}")
print(f"  Total UMI counts: {total_counts.sum():,.0f}")
