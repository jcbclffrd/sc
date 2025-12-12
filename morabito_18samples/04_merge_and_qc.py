#!/usr/bin/env python3
"""
Merge the 18 filtered samples and add Morabito's cell metadata.

This script combines all 18 filtered samples into a single AnnData object
and adds cell type, diagnosis, and other metadata from Morabito's annotations.

Requires: 
  - scte_18samples_filtered/ (from 03_filter_to_morabito_cells.py)
  - GSE174367_snRNA-seq_cell_meta.csv.gz
  
Output: 
  - merged_18samples.h5ad (genes + TEs)
  - merged_18samples_genes.h5ad (genes only, for comparison)
"""

import os
import re
import pandas as pd
import scanpy as sc
import numpy as np

print("=" * 70)
print("MERGING 18 FILTERED SAMPLES AND ADDING METADATA")
print("=" * 70)

# Load sample mapping
mapping = pd.read_csv('sample_mapping.csv')
print(f"\nLoaded {len(mapping)} sample mappings")

# Load Morabito's metadata
print("\nLoading Morabito's cell metadata...")
meta = pd.read_csv('../analysis/morabito_reference/GSE174367_snRNA-seq_cell_meta.csv.gz')
print(f"  {len(meta):,} cells")
print(f"  Metadata columns: {list(meta.columns)}")

# Create barcode + sample to metadata mapping (barcodes are reused across samples)
meta['barcode_clean'] = meta['Barcode'].str.split('-').str[0]
meta['barcode_sample_key'] = meta['barcode_clean'] + '_' + meta['SampleID']
meta_dict = meta.set_index('barcode_sample_key').to_dict('index')

# Load all samples
print(f"\n{'='*70}")
print("LOADING SAMPLES:")
print(f"{'='*70}")

adatas = []
sample_to_srr = dict(zip(mapping['sample_id'], mapping['srr']))

for idx, row in mapping.iterrows():
    srr = row['srr']
    sample_id = row['sample_id']
    
    # scTE names files as {SRR}.h5ad
    input_file = f"scte_18samples_filtered/{srr}/{srr}.h5ad"
    
    if not os.path.exists(input_file):
        print(f"✗ {sample_id} ({srr}): NOT FOUND")
        continue
    
    adata = sc.read_h5ad(input_file)
    
    # Add sample info
    adata.obs['sample_id'] = sample_id
    adata.obs['srr'] = srr
    
    print(f"✓ {sample_id:12s} ({srr}): {adata.n_obs:>5,} cells × {adata.n_vars:>6,} features")
    adatas.append(adata)

print(f"\nLoaded {len(adatas)}/{len(mapping)} samples")

if len(adatas) == 0:
    print("\n✗ No samples loaded! Run 03_filter_to_morabito_cells.py first.")
    exit(1)

# Add feature_type annotations before merging
print(f"\n{'='*70}")
print("ADDING FEATURE TYPE ANNOTATIONS:")
print(f"{'='*70}")

# Load TE names from RepeatMasker
print("Loading TE annotations from RepeatMasker...")
te_names_set = set()
rmsk_file = '/home/jacobc/sc/annotations/hg38_rmsk.gtf'
if os.path.exists(rmsk_file):
    import re
    with open(rmsk_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 9:
                attrs = fields[8]
                match = re.search(r'gene_id "([^"]+)"', attrs)
                if match:
                    # TE names have format: L1M4c#LINE#L1
                    # scTE uses just the first part: L1M4c
                    te_full_name = match.group(1)
                    te_name = te_full_name.split('#')[0]
                    te_names_set.add(te_name)
    print(f"  Loaded {len(te_names_set):,} TE names")
else:
    print("  WARNING: RepeatMasker GTF not found, cannot annotate TEs!")
    print(f"  Looked for: {rmsk_file}")

# Add feature_type to each sample
for idx, adata in enumerate(adatas):
    # Annotate features as 'gene' or 'TE'
    feature_types = []
    for feature_name in adata.var_names:
        if feature_name in te_names_set:
            feature_types.append('TE')
        else:
            feature_types.append('gene')
    
    adata.var['feature_type'] = feature_types
    
    n_genes = sum(1 for ft in feature_types if ft == 'gene')
    n_tes = sum(1 for ft in feature_types if ft == 'TE')
    
    # Print once for first sample to verify
    if idx == 0:
        print(f"  Example annotation (first sample):")
        print(f"    Genes: {n_genes:,}")
        print(f"    TEs: {n_tes:,}")

# Concatenate all samples
print(f"\n{'='*70}")
print("MERGING SAMPLES:")
print(f"{'='*70}")

adata = sc.concat(adatas, join='outer', index_unique=None, merge='same')
print(f"Merged shape: {adata.n_obs:,} cells × {adata.n_vars:,} features")

# Verify feature_type is preserved
if 'feature_type' in adata.var.columns:
    n_genes = (adata.var['feature_type'] == 'gene').sum()
    n_tes = (adata.var['feature_type'] == 'TE').sum()
    print(f"  ✓ Feature types preserved: {n_genes:,} genes, {n_tes:,} TEs")
else:
    print("  ✗ WARNING: feature_type column was lost during merge!")

# Add Morabito's metadata
print("\nAdding Morabito's cell metadata...")

# Clean barcodes for matching
adata.obs['barcode_clean'] = [bc.split('-')[0] if '-' in bc else bc for bc in adata.obs_names]

# Create unique key for barcode + sample
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

# Show metadata summary
print(f"\n{'='*70}")
print("METADATA SUMMARY:")
print(f"{'='*70}")

print("\nDiagnosis distribution:")
print(adata.obs['diagnosis'].value_counts())

print("\nCell type distribution:")
print(adata.obs['cell_type'].value_counts())

print("\nSample distribution:")
print(adata.obs['sample_id'].value_counts().sort_index())

# Calculate QC metrics
print(f"\n{'='*70}")
print("CALCULATING QC METRICS:")
print(f"{'='*70}")

# UMIs per cell
adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()

# Genes detected per cell
adata.obs['n_genes'] = np.array((adata.X > 0).sum(axis=1)).flatten()

# Separate genes and TEs
if 'feature_type' in adata.var.columns:
    is_te = adata.var['feature_type'] == 'TE'
    is_gene = adata.var['feature_type'] == 'gene'
    
    adata.obs['n_genes_detected'] = np.array((adata[:, is_gene].X > 0).sum(axis=1)).flatten()
    adata.obs['n_tes_detected'] = np.array((adata[:, is_te].X > 0).sum(axis=1)).flatten()
    adata.obs['counts_genes'] = np.array(adata[:, is_gene].X.sum(axis=1)).flatten()
    adata.obs['counts_tes'] = np.array(adata[:, is_te].X.sum(axis=1)).flatten()
    adata.obs['pct_counts_tes'] = 100 * adata.obs['counts_tes'] / adata.obs['n_counts']
    
    print(f"Gene features: {is_gene.sum():,}")
    print(f"TE features: {is_te.sum():,}")

print(f"\nQC statistics:")
print(f"  Mean UMIs/cell: {adata.obs['n_counts'].mean():.0f}")
print(f"  Median UMIs/cell: {adata.obs['n_counts'].median():.0f}")
print(f"  Mean genes/cell: {adata.obs['n_genes'].mean():.0f}")
print(f"  Median genes/cell: {adata.obs['n_genes'].median():.0f}")

if 'pct_counts_tes' in adata.obs.columns:
    print(f"  Mean % TE counts: {adata.obs['pct_counts_tes'].mean():.1f}%")
    print(f"  Median % TE counts: {adata.obs['pct_counts_tes'].median():.1f}%")

# Save merged data (genes + TEs)
output_file = 'merged_18samples.h5ad'
print(f"\n{'='*70}")
print(f"SAVING MERGED DATA:")
print(f"{'='*70}")
adata.write_h5ad(output_file)
print(f"✓ Saved: {output_file}")
size_mb = os.path.getsize(output_file) / (1024 * 1024)
print(f"  Size: {size_mb:.1f} MB")

# Save genes-only version (for comparison with Morabito)
if 'feature_type' in adata.var.columns:
    adata_genes = adata[:, adata.var['feature_type'] == 'gene'].copy()
    genes_file = 'merged_18samples_genes.h5ad'
    adata_genes.write_h5ad(genes_file)
    print(f"✓ Saved genes-only: {genes_file}")
    size_mb = os.path.getsize(genes_file) / (1024 * 1024)
    print(f"  Size: {size_mb:.1f} MB")
    print(f"  Shape: {adata_genes.n_obs:,} cells × {adata_genes.n_vars:,} genes")

print(f"{'='*70}")
print("✓ MERGING COMPLETE!")
print(f"{'='*70}")
