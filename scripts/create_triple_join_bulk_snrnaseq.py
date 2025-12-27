#!/usr/bin/env python3
"""
Create a comprehensive join file combining:
1. NEW Bulk RNA-seq TE differential expression (18 samples, AD vs Control)
2. Existing join file with snRNA-seq + previous bulk analysis

Uses TE family names (first substring before colon) as the join key.
The new bulk data will be prefixed with 'new_bulk_' to distinguish from
existing 'bulk_' columns from the previous analysis.
"""

import pandas as pd
import argparse
import sys
from pathlib import Path

def extract_te_family(te_name):
    """
    Extract TE family name (first part before colon).
    Examples:
        L1HS:L1:LINE -> L1HS
        MLT1E3-int:ERVL-MaLR:LTR -> MLT1E3-int
        L1M:L1:LINE -> L1M
    """
    if pd.isna(te_name):
        return None
    te_name = str(te_name)
    return te_name.split(':')[0] if ':' in te_name else te_name

def main():
    parser = argparse.ArgumentParser(
        description='Join NEW bulk RNA-seq TE results with existing snRNA-seq + bulk join file'
    )
    parser.add_argument(
        '--bulk-results',
        default='~/sc/tetranscripts_bulk/results/differential_TE_results_AD_vs_Control.csv',
        help='Path to NEW bulk RNA-seq TE differential expression results CSV (18 samples)'
    )
    parser.add_argument(
        '--snrnaseq-join',
        default='~/sc/tetranscripts_bulk/joined_TE_results_singlecell_bulk.csv',
        help='Path to existing snRNA-seq + previous bulk join file'
    )
    parser.add_argument(
        '--output',
        default='~/sc/tetranscripts_bulk/results/triple_join_bulk_snrnaseq_TEs.csv',
        help='Output path for triple-join results'
    )
    
    args = parser.parse_args()
    
    # Expand paths
    bulk_path = Path(args.bulk_results).expanduser()
    snrnaseq_path = Path(args.snrnaseq_join).expanduser()
    output_path = Path(args.output).expanduser()
    
    print("="*80)
    print("CREATE TRIPLE-JOIN FILE")
    print("Combining: snRNA-seq + previous bulk + NEW bulk (18 samples)")
    print("="*80)
    print()
    
    # Load NEW bulk RNA-seq results (18 samples)
    print(f"Loading NEW bulk RNA-seq results (18 samples):")
    print(f"  {bulk_path}")
    if not bulk_path.exists():
        print(f"ERROR: File not found: {bulk_path}")
        sys.exit(1)
    
    bulk_df = pd.read_csv(bulk_path)
    print(f"  ✓ Loaded {len(bulk_df)} TEs")
    
    # Extract TE family from bulk data
    bulk_df['TE_family'] = bulk_df['TE'].apply(extract_te_family)
    print(f"  ✓ Extracted {bulk_df['TE_family'].nunique()} unique TE families")
    
    # Prefix NEW bulk columns with 'new_bulk_' to distinguish from existing 'bulk_' columns
    new_bulk_rename = {}
    for col in bulk_df.columns:
        if col not in ['TE_family']:
            new_bulk_rename[col] = f'new_bulk_{col}'
    bulk_df = bulk_df.rename(columns=new_bulk_rename)
    
    # Load existing snRNA-seq + previous bulk join file
    print(f"\nLoading existing snRNA-seq + previous bulk join file:")
    print(f"  {snrnaseq_path}")
    if not snrnaseq_path.exists():
        print(f"ERROR: File not found: {snrnaseq_path}")
        sys.exit(1)
    
    snrnaseq_df = pd.read_csv(snrnaseq_path)
    print(f"  ✓ Loaded {len(snrnaseq_df)} entries")
    
    # The snRNA-seq join file uses 'TE_subfamily' as the TE family column
    # We'll use this directly for joining
    if 'TE_subfamily' not in snrnaseq_df.columns:
        print("ERROR: Expected 'TE_subfamily' column not found in snRNA-seq join file")
        sys.exit(1)
    
    # Rename for consistency
    snrnaseq_df = snrnaseq_df.rename(columns={'TE_subfamily': 'TE_family'})
    
    # Check existing columns
    bulk_cols = [c for c in snrnaseq_df.columns if c.startswith('bulk_')]
    print(f"  ✓ Found {len(bulk_cols)} columns from previous bulk analysis")
    
    # Merge on TE_family (outer join to keep all TEs from all three sources)
    print(f"\nPerforming outer join on TE_family...")
    
    merged_df = pd.merge(
        snrnaseq_df,
        bulk_df,
        on='TE_family',
        how='outer',
        indicator=True
    )
    
    print(f"  ✓ Merged dataset contains {len(merged_df)} unique TE families")
    
    # Sort by new bulk significance
    if 'new_bulk_padj' in merged_df.columns:
        merged_df = merged_df.sort_values('new_bulk_padj', na_position='last')
    
    # Create output directory if needed
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save output
    print(f"\nSaving triple-join file to:")
    print(f"  {output_path}")
    merged_df.to_csv(output_path, index=False)
    print(f"  ✓ Saved {len(merged_df)} TE families")
    
    # Summary statistics
    print("\n" + "="*80)
    print("TRIPLE-JOIN SUMMARY")
    print("="*80)
    
    # Count TEs by source
    n_snrnaseq_only = (merged_df['_merge'] == 'left_only').sum()
    n_new_bulk_only = (merged_df['_merge'] == 'right_only').sum()
    n_both = (merged_df['_merge'] == 'both').sum()
    
    print(f"\nTE families by source:")
    print(f"  From snRNA-seq/previous bulk only: {n_snrnaseq_only}")
    print(f"  From NEW bulk only (18 samples):   {n_new_bulk_only}")
    print(f"  In multiple sources:                {n_both}")
    print(f"  Total unique TE families:           {len(merged_df)}")
    
    # Significant TEs
    print(f"\nSignificant TEs (padj < 0.05):")
    
    if 'new_bulk_padj' in merged_df.columns:
        n_sig_new = (merged_df['new_bulk_padj'] < 0.05).sum()
        print(f"  NEW bulk (18 samples): {n_sig_new}")
    
    if 'bulk_padj' in merged_df.columns:
        n_sig_old = (merged_df['bulk_padj'] < 0.05).sum()
        print(f"  Previous bulk:         {n_sig_old}")
    
    if 'pvals_adj' in merged_df.columns:
        n_sig_sc = (merged_df['pvals_adj'] < 0.05).sum()
        print(f"  snRNA-seq:             {n_sig_sc}")
    
    # Concordance (if TE is in both analyses)
    if 'new_bulk_padj' in merged_df.columns and 'bulk_padj' in merged_df.columns:
        both_present = merged_df['new_bulk_padj'].notna() & merged_df['bulk_padj'].notna()
        n_both_present = both_present.sum()
        
        if n_both_present > 0:
            both_sig_new = (merged_df.loc[both_present, 'new_bulk_padj'] < 0.05).sum()
            both_sig_old = (merged_df.loc[both_present, 'bulk_padj'] < 0.05).sum()
            
            print(f"\nConcordance between bulk analyses:")
            print(f"  TEs in both analyses:        {n_both_present}")
            print(f"  Sig in NEW bulk:             {both_sig_new}")
            print(f"  Sig in previous bulk:        {both_sig_old}")
    
    print("\n" + "="*80)
    print("COMPLETE")
    print("="*80)
    print()

if __name__ == '__main__':
    main()
