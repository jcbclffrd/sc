#!/usr/bin/env python3
"""
Create comprehensive metadata for all bulk RNA-seq samples.

This script:
1. Scans for all .cntTable files
2. Loads existing metadata from bulk_rnaseq_18patients.csv
3. Infers diagnosis from GSE174367 sample naming (if possible)
4. Creates a complete metadata file for differential expression
"""

import pandas as pd
from pathlib import Path
import sys
import argparse

def infer_diagnosis_from_gsm(gsm_id):
    """
    Try to infer diagnosis from GSM ID based on Morabito dataset patterns.
    Returns 'AD', 'Control', or 'Unknown'
    """
    # This is a placeholder - would need the full GEO metadata
    # For now, return Unknown for samples not in our metadata
    return 'Unknown'

def create_metadata(tecount_dir, existing_metadata_file, output_file):
    """Create comprehensive metadata for all samples."""
    
    print("=" * 70)
    print("CREATE BULK RNA-SEQ METADATA FOR ALL SAMPLES")
    print("=" * 70)
    print()
    
    # Find all .cntTable files
    print(f"Scanning directory: {tecount_dir}")
    count_files = list(tecount_dir.glob("*.cntTable"))
    print(f"  Found {len(count_files)} .cntTable files")
    print()
    
    # Load existing metadata
    existing_meta = None
    if existing_metadata_file.exists():
        print(f"Loading existing metadata: {existing_metadata_file}")
        existing_meta = pd.read_csv(existing_metadata_file)
        print(f"  Existing metadata: {len(existing_meta)} samples")
        print()
    
    # Create metadata for all samples
    all_metadata = []
    
    for file_path in sorted(count_files):
        srr_id = file_path.stem
        
        # Check if we have existing metadata
        if existing_meta is not None and 'srr_id' in existing_meta.columns:
            existing_row = existing_meta[existing_meta['srr_id'] == srr_id]
            if len(existing_row) > 0:
                # Use existing metadata
                row_dict = existing_row.iloc[0].to_dict()
                all_metadata.append(row_dict)
                continue
        
        # Create minimal metadata for new samples
        sample_dict = {
            'srr_id': srr_id,
            'diagnosis': 'Unknown',  # Will need to be filled in manually
            'te_count_file': f"{srr_id}.cntTable",
            'included_in_combined_matrix': 'Yes',
            'included_in_de_analysis': 'No',  # Only include known diagnosis
            'notes': 'Auto-detected from .cntTable files'
        }
        all_metadata.append(sample_dict)
    
    # Create DataFrame
    metadata_df = pd.DataFrame(all_metadata)
    
    # Set included_in_de_analysis based on diagnosis
    if 'diagnosis' in metadata_df.columns:
        metadata_df['included_in_de_analysis'] = metadata_df['diagnosis'].apply(
            lambda x: 'Yes' if x in ['AD', 'Control'] else 'No'
        )
    
    # Sort by SRR ID
    metadata_df = metadata_df.sort_values('srr_id')
    
    # Save
    print(f"Saving metadata to: {output_file}")
    metadata_df.to_csv(output_file, index=False)
    print(f"  Saved {len(metadata_df)} samples")
    print()
    
    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    if 'diagnosis' in metadata_df.columns:
        print("\nDiagnosis breakdown:")
        for diagnosis, count in metadata_df['diagnosis'].value_counts().items():
            print(f"  {diagnosis}: {count}")
    
    if 'included_in_de_analysis' in metadata_df.columns:
        n_included = (metadata_df['included_in_de_analysis'] == 'Yes').sum()
        n_excluded = (metadata_df['included_in_de_analysis'] == 'No').sum()
        print(f"\nSamples for DE analysis:")
        print(f"  Included: {n_included}")
        print(f"  Excluded (unknown diagnosis): {n_excluded}")
    
    print()
    print("NOTE: Samples with 'Unknown' diagnosis need to be updated manually")
    print("      with diagnosis information before running DE analysis.")
    print()

def main():
    parser = argparse.ArgumentParser(
        description='Create comprehensive metadata for all bulk RNA-seq samples',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('-i', '--tecount-dir', type=Path, required=True,
                        help='Directory containing .cntTable files')
    parser.add_argument('-e', '--existing-metadata', type=Path, required=False,
                        help='Existing metadata CSV to merge with')
    parser.add_argument('-o', '--output', type=Path, required=True,
                        help='Output metadata CSV file')
    
    args = parser.parse_args()
    
    if not args.tecount_dir.exists():
        print(f"Error: TEcount directory not found: {args.tecount_dir}", file=sys.stderr)
        sys.exit(1)
    
    create_metadata(args.tecount_dir, args.existing_metadata, args.output)

if __name__ == '__main__':
    main()
