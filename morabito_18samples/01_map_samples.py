#!/usr/bin/env python3
"""
Map Morabito's Sample-X names to SRR accessions.

This script parses the GEO series matrix to create a mapping between
Morabito's sample names (Sample-17, Sample-19, etc.) and the corresponding
SRR accessions in our dataset.

Output: sample_mapping.csv
"""

import re
import pandas as pd
import os

print("=" * 70)
print("MAPPING MORABITO'S SAMPLE NAMES TO SRR ACCESSIONS")
print("=" * 70)

# Parse the GEO series matrix
geo_file = '../analysis/morabito_reference/GSE174367_series_matrix.txt'
print(f"\nReading GEO series matrix: {geo_file}")

with open(geo_file, 'r') as f:
    lines = [l.strip() for l in f.readlines()]

# Extract relevant fields
data = {}
for line in lines:
    if line.startswith('!Sample_title'):
        data['title'] = [t.strip('"') for t in line.split('\t')[1:]]
    elif line.startswith('!Sample_geo_accession'):
        data['gsm'] = [g.strip('"') for g in line.split('\t')[1:]]
    elif line.startswith('!Sample_relation') and 'SRA' in line:
        sra_ids = []
        for item in line.split('\t')[1:]:
            item = item.strip('"')
            if 'SRX' in item:
                srx = re.search(r'SRX\d+', item)
                if srx:
                    sra_ids.append(srx.group())
            else:
                sra_ids.append('')
        data['srx'] = sra_ids

# Create DataFrame
df = pd.DataFrame(data)

# Filter to snRNA-seq only
snrna_df = df[df['title'].str.contains('snRNA-seq', na=False)].copy()

# Clean up the Sample-X name
snrna_df['sample_id'] = snrna_df['title'].str.replace(r' \(snRNA-seq\)', '', regex=True)

print(f"\nFound {len(snrna_df)} snRNA-seq samples in GEO")

# Map SRX to SRR using our SRA metadata
sra = pd.read_csv('../data/sra_metadata.csv')
print(f"Loaded SRA metadata: {len(sra)} records")

# Create SRX to SRR mapping
srx_to_srr = dict(zip(sra['Experiment'], sra['Run']))
snrna_df['srr'] = snrna_df['srx'].map(srx_to_srr)

# Morabito's 18 samples from the paper
morabito_18_samples = [
    'Sample-100', 'Sample-17', 'Sample-19', 'Sample-22', 'Sample-27', 'Sample-33',
    'Sample-37', 'Sample-43', 'Sample-45', 'Sample-46', 'Sample-47', 'Sample-50',
    'Sample-52', 'Sample-58', 'Sample-66', 'Sample-82', 'Sample-90', 'Sample-96'
]

# Filter to the 18 samples
mapping_18 = snrna_df[snrna_df['sample_id'].isin(morabito_18_samples)].copy()
mapping_18 = mapping_18[['sample_id', 'gsm', 'srx', 'srr']].sort_values('sample_id')

print(f"\n{'='*70}")
print("MORABITO'S 18 SAMPLES MAPPING:")
print(f"{'='*70}")
print(mapping_18.to_string(index=False))

# Check which samples we have aligned
aligned_dir = '../starsolo_aligned'
if os.path.exists(aligned_dir):
    aligned_samples = set([d for d in os.listdir(aligned_dir) if d.startswith('SRR')])
    morabito_srrs = set(mapping_18['srr'].dropna())
    overlap = morabito_srrs & aligned_samples
    
    print(f"\n{'='*70}")
    print("AVAILABILITY CHECK:")
    print(f"{'='*70}")
    print(f"Morabito's 18 samples: {len(morabito_srrs)} SRR IDs")
    print(f"Our aligned samples: {len(aligned_samples)}")
    print(f"Overlap (available): {len(overlap)} samples")
    
    if len(overlap) < len(morabito_srrs):
        missing = morabito_srrs - overlap
        print(f"\n⚠️  Missing {len(missing)} samples:")
        for srr in sorted(missing):
            sample_id = mapping_18[mapping_18['srr'] == srr]['sample_id'].iloc[0]
            print(f"  {srr} ({sample_id})")
    else:
        print(f"\n✓ All 18 samples are available!")

# Check scTE output
scte_dir = '../scTE_output'
if os.path.exists(scte_dir):
    scte_samples = set([d for d in os.listdir(scte_dir) if d.startswith('SRR')])
    overlap_scte = morabito_srrs & scte_samples
    
    print(f"\nscTE output availability: {len(overlap_scte)}/{len(morabito_srrs)} samples")
    
    if len(overlap_scte) < len(morabito_srrs):
        missing_scte = morabito_srrs - overlap_scte
        print(f"⚠️  Missing scTE output for {len(missing_scte)} samples:")
        for srr in sorted(missing_scte):
            sample_id = mapping_18[mapping_18['srr'] == srr]['sample_id'].iloc[0]
            print(f"  {srr} ({sample_id})")

# Save mapping
output_file = 'sample_mapping.csv'
mapping_18.to_csv(output_file, index=False)
print(f"\n{'='*70}")
print(f"✓ Saved mapping to: {output_file}")
print(f"{'='*70}")

# Also save the full snRNA-seq mapping for reference
full_output = 'all_snrnaseq_samples.csv'
snrna_df[['sample_id', 'gsm', 'srx', 'srr']].to_csv(full_output, index=False)
print(f"✓ Saved full snRNA-seq mapping to: {full_output}")
