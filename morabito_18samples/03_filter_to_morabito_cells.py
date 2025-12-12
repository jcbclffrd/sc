#!/usr/bin/env python3
"""
Filter our scTE data to Morabito's exact cell barcodes.

This script loads Morabito's filtered cell barcodes from GEO and filters
our scTE data to match exactly the same cells they used in their analysis.

Requires: 
  - scte_18samples/ (from 02_extract_18samples.py)
  - GSE174367_snRNA-seq_cell_meta.csv.gz (from GEO)
  
Output: scte_18samples_filtered/ with filtered .h5ad files
"""

import os
import pandas as pd
import scanpy as sc
import numpy as np
from collections import defaultdict

print("=" * 70)
print("FILTERING scTE DATA TO MORABITO'S EXACT CELL BARCODES")
print("=" * 70)

# Load Morabito's cell metadata
print("\nLoading Morabito's cell metadata...")
meta = pd.read_csv('../analysis/morabito_reference/GSE174367_snRNA-seq_cell_meta.csv.gz')
print(f"  {len(meta):,} cells in metadata")
print(f"  Columns: {list(meta.columns)}")

# Load sample mapping
mapping = pd.read_csv('sample_mapping.csv')
print(f"\nLoaded {len(mapping)} sample mappings")

# Create Sample-X to SRR mapping
sample_to_srr = dict(zip(mapping['sample_id'], mapping['srr']))

# Group barcodes by sample
print("\nGrouping barcodes by sample...")
sample_barcodes = defaultdict(set)

for idx, row in meta.iterrows():
    barcode = row['Barcode']
    sample_id = row['SampleID']
    
    if sample_id in sample_to_srr:
        srr = sample_to_srr[sample_id]
        # Store just the barcode (without -1 suffix if present)
        bc = barcode.split('-')[0] if '-' in barcode else barcode
        sample_barcodes[srr].add(bc)

print(f"\nBarcodes per sample:")
for sample_id in sorted(mapping['sample_id'].unique()):
    if sample_id in sample_to_srr:
        srr = sample_to_srr[sample_id]
        n_barcodes = len(sample_barcodes[srr])
        print(f"  {sample_id:12s} ({srr}): {n_barcodes:>5,} cells")

# Create output directory
output_dir = 'scte_18samples_filtered'
os.makedirs(output_dir, exist_ok=True)

# Process each sample
print(f"\n{'='*70}")
print("FILTERING SAMPLES:")
print(f"{'='*70}")

total_cells_before = 0
total_cells_after = 0
processed = 0

for idx, row in mapping.iterrows():
    srr = row['srr']
    sample_id = row['sample_id']
    
    # scTE names files as {SRR}.h5ad
    input_file = f"scte_18samples/{srr}/{srr}.h5ad"
    output_subdir = f"{output_dir}/{srr}"
    output_file = f"{output_subdir}/{srr}.h5ad"
    
    if not os.path.exists(input_file):
        print(f"✗ {sample_id} ({srr}): INPUT FILE NOT FOUND")
        continue
    
    # Load our scTE data
    adata = sc.read_h5ad(input_file)
    n_before = adata.n_obs
    total_cells_before += n_before
    
    # Get barcodes for this sample (remove -1 suffix from our data)
    our_barcodes = [bc.split('-')[0] if '-' in bc else bc for bc in adata.obs_names]
    adata.obs['barcode_clean'] = our_barcodes
    
    # Get Morabito's barcodes for this sample
    morabito_bcs = sample_barcodes[srr]
    
    if len(morabito_bcs) == 0:
        print(f"⚠️  {sample_id} ({srr}): NO BARCODES IN MORABITO'S DATA")
        continue
    
    # Filter to matching barcodes
    mask = adata.obs['barcode_clean'].isin(morabito_bcs)
    adata_filtered = adata[mask].copy()
    n_after = adata_filtered.n_obs
    total_cells_after += n_after
    
    # Calculate overlap percentage
    pct = 100 * n_after / len(morabito_bcs)
    
    print(f"{sample_id:12s} ({srr}): {n_before:>5,} → {n_after:>5,} cells ({pct:5.1f}% of Morabito's {len(morabito_bcs):,})")
    
    # Save filtered data
    os.makedirs(output_subdir, exist_ok=True)
    adata_filtered.write_h5ad(output_file)
    processed += 1

print(f"\n{'='*70}")
print("SUMMARY:")
print(f"{'='*70}")
print(f"Processed: {processed}/{len(mapping)} samples")
print(f"Total cells before filtering: {total_cells_before:,}")
print(f"Total cells after filtering: {total_cells_after:,}")
print(f"Retention rate: {100*total_cells_after/total_cells_before:.1f}%")
print(f"Morabito's total cells: {len(meta):,}")

if total_cells_after > 0:
    match_rate = 100 * total_cells_after / len(meta)
    print(f"Match with Morabito: {match_rate:.1f}%")
    
    if match_rate < 90:
        print("\n⚠️  WARNING: Low match rate!")
        print("Possible reasons:")
        print("  - Barcode format mismatch (check -1 suffix)")
        print("  - Different cell calling between CellRanger and STARsolo")
        print("  - Sample mismatch")
    elif match_rate > 110:
        print("\n⚠️  WARNING: More cells than expected!")
        print("Possible reasons:")
        print("  - scTE kept cells that Morabito filtered out")
        print("  - Need to apply additional QC filters")
    else:
        print("\n✓ Good match with Morabito's cell counts!")

print(f"{'='*70}")
