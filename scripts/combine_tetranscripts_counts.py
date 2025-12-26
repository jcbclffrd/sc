#!/usr/bin/env python3
"""
Combine TEtranscripts count tables from multiple samples into a single expression matrix.
Output: Combined matrix with genes/TEs as rows and samples as columns.
"""

import os
import pandas as pd
import argparse
from pathlib import Path

def read_count_table(filepath):
    """Read a single TEcount .cntTable file."""
    # Read the file, first column is gene/TE name, second is count
    df = pd.read_csv(filepath, sep='\t', header=0, index_col=0)
    # Get the sample name from the filename
    sample_name = Path(filepath).stem.replace('.cntTable', '')
    # Rename the count column to the sample name
    df.columns = [sample_name]
    return df

def combine_count_tables(input_dir, output_file, sample_mapping=None):
    """Combine all .cntTable files in the directory."""
    
    # Find all .cntTable files
    count_files = sorted(Path(input_dir).glob('*.cntTable'))
    
    if len(count_files) == 0:
        print(f"ERROR: No .cntTable files found in {input_dir}")
        return
    
    print(f"Found {len(count_files)} count tables to combine")
    
    # Read all count tables
    dfs = []
    for f in count_files:
        print(f"  Reading {f.name}...")
        df = read_count_table(f)
        dfs.append(df)
    
    # Combine all dataframes
    print("\nCombining count tables...")
    combined = pd.concat(dfs, axis=1, join='outer')
    
    # Fill NaN with 0 (genes/TEs not present in some samples)
    combined = combined.fillna(0).astype(int)
    
    # Sort by row names (gene/TE names)
    combined = combined.sort_index()
    
    # If sample mapping provided, add patient IDs as additional row
    if sample_mapping:
        print(f"\nAdding patient mapping from {sample_mapping}...")
        mapping_df = pd.read_csv(sample_mapping)
        # Create a dictionary mapping SRR to patient
        srr_to_patient = dict(zip(mapping_df['srr'], mapping_df['patient']))
        # Rename columns to include patient ID
        new_cols = []
        for col in combined.columns:
            if col in srr_to_patient:
                new_cols.append(f"{col}_{srr_to_patient[col]}")
            else:
                new_cols.append(col)
        combined.columns = new_cols
    
    # Separate genes and TEs
    genes = combined[combined.index.str.startswith('ENSG')]
    tes = combined[~combined.index.str.startswith('ENSG')]
    
    # Write outputs
    print(f"\nWriting combined matrix to {output_file}...")
    combined.to_csv(output_file, sep='\t')
    
    # Write separate gene and TE matrices
    gene_file = output_file.replace('.tsv', '_genes_only.tsv')
    te_file = output_file.replace('.tsv', '_TEs_only.tsv')
    
    print(f"Writing gene-only matrix to {gene_file}...")
    genes.to_csv(gene_file, sep='\t')
    
    print(f"Writing TE-only matrix to {te_file}...")
    tes.to_csv(te_file, sep='\t')
    
    # Print summary statistics
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total samples: {combined.shape[1]}")
    print(f"Total features: {combined.shape[0]:,}")
    print(f"  - Genes (ENSG*): {genes.shape[0]:,}")
    print(f"  - TEs: {tes.shape[0]:,}")
    print(f"\nOutput files:")
    print(f"  - Combined: {output_file}")
    print(f"  - Genes only: {gene_file}")
    print(f"  - TEs only: {te_file}")
    print("="*60)
    
    # Print top expressed TEs
    print("\nTop 20 most abundant TEs (total counts across all samples):")
    te_sums = tes.sum(axis=1).sort_values(ascending=False)
    for i, (te_name, count) in enumerate(te_sums.head(20).items(), 1):
        print(f"  {i:2d}. {te_name:30s} {count:>12,}")
    
    return combined, genes, tes

def main():
    parser = argparse.ArgumentParser(description='Combine TEtranscripts count tables')
    parser.add_argument('--input-dir', '-i', required=True,
                        help='Directory containing .cntTable files')
    parser.add_argument('--output', '-o', required=True,
                        help='Output TSV file for combined matrix')
    parser.add_argument('--mapping', '-m', default=None,
                        help='Optional: CSV file mapping SRR to patient IDs (columns: srr, patient)')
    
    args = parser.parse_args()
    
    combine_count_tables(args.input_dir, args.output, args.mapping)

if __name__ == '__main__':
    main()
