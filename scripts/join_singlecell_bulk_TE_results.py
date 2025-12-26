#!/usr/bin/env python3
"""
Join single-cell and bulk RNA-seq TE differential expression results.

Single-cell file: de_TEs_only.csv from cellranger-morabito-results
- Extract TE_subfamily from 'names' column (first string before semicolon)

Bulk file: differential_TE_results_AD_vs_Control.csv from sc repo
- Extract TE_subfamily from 'TE' column (first string before first colon)

Output: Left join on TE_subfamily, keeping all single-cell TEs
"""

import pandas as pd
import argparse
import sys

def extract_te_subfamily_singlecell(name):
    """Extract TE subfamily from single-cell names column.
    
    Examples:
        'MLT1E3-int:TE' -> 'MLT1E3-int'
        'L1M:TE' -> 'L1M'
        'MIR4500HG' -> 'MIR4500HG' (no semicolon, return as is)
    """
    # Split by semicolon first (for entries like "name;extra")
    if ';' in name:
        name = name.split(';')[0]
    
    # Then split by colon (for entries like "subfamily:TE")
    if ':' in name:
        return name.split(':')[0]
    
    return name

def extract_te_subfamily_bulk(te_name):
    """Extract TE subfamily from bulk TE column.
    
    Examples:
        'L2a:L2:LINE' -> 'L2a'
        'AluJb:Alu:SINE' -> 'AluJb'
        'Charlie1:hAT-Charlie:DNA' -> 'Charlie1'
    """
    if ':' in te_name:
        return te_name.split(':')[0]
    return te_name

def join_te_results(singlecell_file, bulk_file, output_file):
    """Join single-cell and bulk TE differential expression results."""
    
    print("Loading single-cell TE results...")
    sc_df = pd.read_csv(singlecell_file)
    print(f"  Loaded {len(sc_df)} single-cell TEs")
    
    print("\nLoading bulk TE results...")
    bulk_df = pd.read_csv(bulk_file)
    print(f"  Loaded {len(bulk_df)} bulk TEs")
    
    # Extract TE_subfamily from single-cell data
    print("\nExtracting TE subfamily from single-cell 'names' column...")
    sc_df['TE_subfamily'] = sc_df['names'].apply(extract_te_subfamily_singlecell)
    
    # Extract TE_subfamily from bulk data
    print("Extracting TE subfamily from bulk 'TE' column...")
    bulk_df['TE_subfamily'] = bulk_df['TE'].apply(extract_te_subfamily_bulk)
    
    # Rename bulk columns to avoid conflicts
    print("\nRenaming bulk columns with 'bulk_' prefix...")
    bulk_rename = {}
    for col in bulk_df.columns:
        if col != 'TE_subfamily':
            bulk_rename[col] = f'bulk_{col}'
    bulk_df = bulk_df.rename(columns=bulk_rename)
    
    # Perform left join on TE_subfamily
    print("\nPerforming left join on TE_subfamily...")
    print(f"  Single-cell unique TE subfamilies: {sc_df['TE_subfamily'].nunique()}")
    print(f"  Bulk unique TE subfamilies: {bulk_df['TE_subfamily'].nunique()}")
    
    merged_df = sc_df.merge(bulk_df, on='TE_subfamily', how='left', suffixes=('_sc', '_bulk'))
    
    print(f"\nJoined {len(merged_df)} rows")
    print(f"  Matched with bulk data: {merged_df['bulk_TE'].notna().sum()}")
    print(f"  No bulk match: {merged_df['bulk_TE'].isna().sum()}")
    
    # Reorder columns for better readability
    first_cols = ['TE_subfamily', 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj', 
                  'bulk_TE', 'bulk_baseMean', 'bulk_log2FoldChange', 'bulk_pvalue', 'bulk_padj']
    other_cols = [col for col in merged_df.columns if col not in first_cols]
    merged_df = merged_df[first_cols + other_cols]
    
    # Save output
    print(f"\nSaving joined results to {output_file}...")
    merged_df.to_csv(output_file, index=False)
    
    # Summary statistics
    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)
    
    # Count significant TEs in both datasets
    sc_sig = (sc_df['pvals_adj'] < 0.05).sum()
    bulk_sig = (bulk_df['bulk_padj'] < 0.05).sum()
    both_sig = ((merged_df['pvals_adj'] < 0.05) & (merged_df['bulk_padj'] < 0.05)).sum()
    
    print(f"\nSignificant TEs (FDR < 0.05):")
    print(f"  Single-cell only: {sc_sig}")
    print(f"  Bulk only: {bulk_sig}")
    print(f"  Significant in both: {both_sig}")
    
    # Direction concordance for matched TEs
    matched = merged_df[merged_df['bulk_TE'].notna()].copy()
    if len(matched) > 0:
        concordant = ((matched['logfoldchanges'] > 0) == (matched['bulk_log2FoldChange'] > 0)).sum()
        print(f"\nDirection concordance (matched TEs):")
        print(f"  Same direction: {concordant} / {len(matched)} ({100*concordant/len(matched):.1f}%)")
        
        # Top concordant TEs by single-cell significance
        print(f"\nTop 10 TEs by single-cell significance (matched with bulk):")
        top_matched = matched.sort_values('pvals_adj').head(10)
        for idx, row in top_matched.iterrows():
            sc_dir = "↑" if row['logfoldchanges'] > 0 else "↓"
            bulk_dir = "↑" if row['bulk_log2FoldChange'] > 0 else "↓"
            agree = "✓" if sc_dir == bulk_dir else "✗"
            print(f"  {agree} {row['TE_subfamily']:20s} SC: {sc_dir} {row['logfoldchanges']:6.2f} (p={row['pvals_adj']:.2e})  "
                  f"Bulk: {bulk_dir} {row['bulk_log2FoldChange']:6.2f} (p={row['bulk_padj']:.2e})")
    
    print("\n" + "="*70)
    print(f"Output saved to: {output_file}")
    print("="*70)
    
    return merged_df

def main():
    parser = argparse.ArgumentParser(
        description='Join single-cell and bulk TE differential expression results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Join the two TE results files
  %(prog)s -s de_TEs_only.csv -b differential_TE_results_AD_vs_Control.csv -o joined_TE_results.csv
        """
    )
    
    parser.add_argument('-s', '--singlecell', required=True,
                        help='Single-cell TE results CSV (de_TEs_only.csv)')
    parser.add_argument('-b', '--bulk', required=True,
                        help='Bulk TE results CSV (differential_TE_results_AD_vs_Control.csv)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output CSV file for joined results')
    
    args = parser.parse_args()
    
    try:
        join_te_results(args.singlecell, args.bulk, args.output)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
