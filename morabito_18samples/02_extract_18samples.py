#!/usr/bin/env python3
"""
Extract scTE output for the 18 samples used by Morabito et al.

This script copies the scTE .h5ad files for the 18 samples into a
separate directory for focused analysis.

Requires: sample_mapping.csv (from 01_map_samples.py)
Output: scte_18samples/ directory with 18 .h5ad files
"""

import os
import shutil
import pandas as pd

print("=" * 70)
print("EXTRACTING scTE OUTPUT FOR MORABITO'S 18 SAMPLES")
print("=" * 70)

# Load sample mapping
mapping = pd.read_csv('sample_mapping.csv')
print(f"\nLoaded mapping: {len(mapping)} samples")

# Create output directory
output_dir = 'scte_18samples'
os.makedirs(output_dir, exist_ok=True)
print(f"Output directory: {output_dir}/")

# Source directory
source_dir = '../scTE_output'

# Extract files
print(f"\n{'='*70}")
print("COPYING scTE FILES:")
print(f"{'='*70}")

copied = 0
missing = 0
total_size = 0

for idx, row in mapping.iterrows():
    srr = row['srr']
    sample_id = row['sample_id']
    
    # scTE names files as {SRR}.h5ad, not scTE_counts.h5ad
    source_file = f"{source_dir}/{srr}/{srr}.h5ad"
    dest_dir = f"{output_dir}/{srr}"
    dest_file = f"{dest_dir}/{srr}.h5ad"
    
    if os.path.exists(source_file):
        # Create destination directory
        os.makedirs(dest_dir, exist_ok=True)
        
        # Copy file
        shutil.copy2(source_file, dest_file)
        
        # Get file size
        size_mb = os.path.getsize(dest_file) / (1024 * 1024)
        total_size += size_mb
        
        print(f"✓ {srr} ({sample_id}): {size_mb:.1f} MB")
        copied += 1
    else:
        print(f"✗ {srr} ({sample_id}): NOT FOUND")
        missing += 1

print(f"\n{'='*70}")
print("SUMMARY:")
print(f"{'='*70}")
print(f"Copied: {copied}/{len(mapping)} samples")
print(f"Missing: {missing}/{len(mapping)} samples")
print(f"Total size: {total_size:.1f} MB ({total_size/1024:.2f} GB)")

if copied == len(mapping):
    print("\n✓ All 18 samples successfully extracted!")
else:
    print(f"\n⚠️  {missing} samples are missing from scTE_output/")
    print("Run the scTE quantification for these samples first.")

print(f"{'='*70}")
