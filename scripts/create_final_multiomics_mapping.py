#!/usr/bin/env python3
"""
Build complete ATAC-seq to snRNA-seq patient mapping by querying GEO
"""

import pandas as pd
import requests
import re
import time

print("Loading base data...")
sra = pd.read_csv('data/sra_metadata.csv')
snrna_18 = pd.read_csv('morabito_18samples/sample_mapping.csv')
snrna_cells = pd.read_csv('data/GSE174367_snRNA-seq_cell_meta.csv.gz')

# Get clinical metadata
clinical = snrna_cells.groupby('SampleID').agg({
    'Diagnosis': 'first',
    'Age': 'first',
    'Sex': 'first',
    'PMI': 'first',
    'Tangle.Stage': 'first',
    'Plaque.Stage': 'first',
    'RIN': 'first',
    'Batch': 'first'
}).reset_index()

# Get all ATAC-seq GSM IDs
atac_gsms = sra[sra['LibraryStrategy'] == 'ATAC-seq']['SampleName'].unique()
print(f"Found {len(atac_gsms)} ATAC-seq GSM IDs")

# Query GEO for each GSM to get Sample-XX mapping
gsm_to_sample = {}
print("\nQuerying GEO for ATAC-seq sample names...")
for i, gsm in enumerate(sorted(atac_gsms), 1):
    try:
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"
        response = requests.get(url)
        if response.status_code == 200:
            # Extract Sample-XX from title
            match = re.search(r'Sample-(\d+) \(snATAC-seq\)', response.text)
            if match:
                sample_id = f"Sample-{match.group(1)}"
                gsm_to_sample[gsm] = sample_id
                print(f"  {i}/{len(atac_gsms)}: {gsm} → {sample_id}")
            else:
                print(f"  {i}/{len(atac_gsms)}: {gsm} → NOT FOUND")
        time.sleep(0.3)  # Be nice to NCBI servers
    except Exception as e:
        print(f"  Error for {gsm}: {e}")

print(f"\nMapped {len(gsm_to_sample)} ATAC GSMs to Sample IDs")

# Build comprehensive mapping
results = []
for _, row in snrna_18.iterrows():
    sample_id = row['sample_id']
    
    # Get clinical info
    clin = clinical[clinical['SampleID'] == sample_id].iloc[0] if sample_id in clinical['SampleID'].values else None
    
    # Get snRNA info
    snrna_gsm = row['gsm']
    snrna_srr = row['srr']
    
    # Find ATAC info - get GSM(s) that map to this sample
    atac_gsms_for_sample = [gsm for gsm, sid in gsm_to_sample.items() if sid == sample_id]
    
    # Get all SRR runs for these GSMs
    atac_srrs = []
    for gsm in atac_gsms_for_sample:
        srrs = sra[(sra['SampleName'] == gsm) & (sra['LibraryStrategy'] == 'ATAC-seq')]['Run'].tolist()
        atac_srrs.extend(srrs)
    
    result = {
        'patient_sample_id': sample_id,
        'sample_number': sample_id.replace('Sample-', ''),
        'diagnosis': clin['Diagnosis'] if clin is not None else None,
        'age': int(clin['Age']) if clin is not None else None,
        'sex': clin['Sex'] if clin is not None else None,
        'pmi': float(clin['PMI']) if clin is not None else None,
        'tangle_stage': clin['Tangle.Stage'] if clin is not None else None,
        'plaque_stage': clin['Plaque.Stage'] if clin is not None else None,
        'rin': float(clin['RIN']) if clin is not None else None,
        'batch': int(clin['Batch']) if clin is not None else None,
        'snRNA_gsm': snrna_gsm,
        'snRNA_srr': snrna_srr,
        'atac_gsm': ','.join(atac_gsms_for_sample) if atac_gsms_for_sample else None,
        'atac_srr': ','.join(sorted(atac_srrs)) if atac_srrs else None,
        'atac_num_runs': len(atac_srrs),
        'has_both_modalities': len(atac_srrs) > 0,
    }
    results.append(result)

df = pd.DataFrame(results)

# Save
output_file = 'data/patient_multiomics_mapping.csv'
df.to_csv(output_file, index=False)

print("\n" + "="*70)
print("COMPLETE MULTI-OMICS MAPPING")
print("="*70)
print(f"✓ Saved to: {output_file}")
print(f"\nTotal patients: {len(df)}")
print(f"With snRNA-seq: {len(df)}")
print(f"With ATAC-seq: {df['has_both_modalities'].sum()}")
print(f"Total ATAC-seq runs: {df['atac_num_runs'].sum()}")

print("\nBy Diagnosis:")
print(df.groupby('diagnosis').agg({
    'patient_sample_id': 'count',
    'has_both_modalities': 'sum',
    'atac_num_runs': 'sum'
}).rename(columns={
    'patient_sample_id': 'n_patients',
    'has_both_modalities': 'with_ATAC',
    'atac_num_runs': 'total_ATAC_runs'
}))

print("\n" + "="*70)
print("Patients with BOTH modalities:")
print("="*70)
both = df[df['has_both_modalities']]
for _, row in both.iterrows():
    print(f"\n{row['patient_sample_id']}: {row['diagnosis']}, Age {row['age']}, {row['sex']}")
    print(f"  snRNA: {row['snRNA_srr']}")
    print(f"  ATAC:  {row['atac_num_runs']} runs - {row['atac_srr']}")

print("\n" + "="*70)
