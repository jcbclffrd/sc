#!/usr/bin/env python3
"""
Combine TEtranscripts count tables from bulk RNA-seq samples.

This script:
1. Reads all .cntTable files from TEcount output
2. Combines them into unified count matrices
3. Generates three matrices:
   - combined_counts_matrix.tsv (all features)
   - combined_counts_matrix_genes_only.tsv (genes only)
   - combined_counts_matrix_TEs_only.tsv (TEs only)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import argparse

def load_tetranscripts_table(file_path):
    """Load a single TEtranscripts .cntTable file."""
    df = pd.read_csv(file_path, sep='\t', index_col=0)
    # TEcount outputs a single column with the sample name
    return df

def combine_count_tables(input_dir, metadata_file, output_dir):
    """
    Combine all TEcount .cntTable files into unified count matrices.
    
    Parameters
    ----------
    input_dir : Path
        Directory containing .cntTable files
    metadata_file : Path or None
        CSV file with sample metadata (optional - will auto-detect if None)
    output_dir : Path
        Directory to write combined matrices
    """
    
    print("=" * 70)
    print("COMBINE TETRANSCRIPTS COUNT TABLES")
    print("=" * 70)
    print()
    
    # Auto-detect all .cntTable files
    print(f"Scanning directory: {input_dir}")
    count_files = list(input_dir.glob("*.cntTable"))
    print(f"  Found {len(count_files)} .cntTable files")
    
    # Load metadata if provided
    metadata = None
    if metadata_file and metadata_file.exists():
        print(f"\nLoading metadata from: {metadata_file}")
        metadata = pd.read_csv(metadata_file)
        print(f"  Found {len(metadata)} samples in metadata")
        
        # Filter to samples included in analysis (if column exists)
        if 'included_in_combined_matrix' in metadata.columns:
            metadata = metadata[metadata['included_in_combined_matrix'] == 'Yes'].copy()
            print(f"  Using {len(metadata)} samples for combined matrix")
    else:
        print("\nNo metadata file provided - using all detected .cntTable files")
    
    print()
    
    # Load all count tables
    count_dfs = {}
    failed_samples = []
    
    print("Loading count tables...")
    for file_path in sorted(count_files):
        # Extract sample ID from filename (SRR*.cntTable -> SRR*)
        sample_id = file_path.stem
        
        # If metadata provided, optionally filter
        if metadata is not None:
            if 'srr_id' in metadata.columns:
                if sample_id not in metadata['srr_id'].values:
                    continue
        
        try:
            df = load_tetranscripts_table(file_path)
            # Rename the count column to the sample ID
            df.columns = [sample_id]
            count_dfs[sample_id] = df
            print(f"  ✓ {sample_id}: {len(df):,} features")
        except Exception as e:
            print(f"  ✗ {sample_id}: Error loading - {e}")
            failed_samples.append(sample_id)
    
    print()
    print(f"Successfully loaded: {len(count_dfs)}/{len(count_files)} samples")
    
    if failed_samples:
        print(f"Failed samples: {', '.join(failed_samples)}")
    
    if len(count_dfs) == 0:
        print("\n✗ No count tables loaded! Check input directory and metadata.")
        sys.exit(1)
    
    # Combine all count tables
    print()
    print("Combining count tables...")
    combined = pd.concat(count_dfs.values(), axis=1, join='outer')
    combined = combined.fillna(0).astype(int)
    
    print(f"  Combined shape: {combined.shape[0]:,} features × {combined.shape[1]} samples")
    print()
    
    # Identify genes vs TEs
    # TEtranscripts uses the naming convention:
    # - Genes: start with "ENS" (ENSEMBL IDs)
    # - TEs: do not start with "ENS"
    is_gene = combined.index.str.startswith('ENS')
    is_te = ~is_gene
    
    n_genes = is_gene.sum()
    n_tes = is_te.sum()
    
    print(f"Feature breakdown:")
    print(f"  Genes: {n_genes:,}")
    print(f"  TEs: {n_tes:,}")
    print()
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save combined matrix (all features)
    output_all = output_dir / "combined_counts_matrix.tsv"
    print(f"Saving combined matrix (all features)...")
    print(f"  -> {output_all}")
    combined.to_csv(output_all, sep='\t')
    file_size_mb = output_all.stat().st_size / (1024 * 1024)
    print(f"  File size: {file_size_mb:.1f} MB")
    print()
    
    # Save genes only
    output_genes = output_dir / "combined_counts_matrix_genes_only.tsv"
    print(f"Saving genes only matrix...")
    print(f"  -> {output_genes}")
    combined[is_gene].to_csv(output_genes, sep='\t')
    file_size_mb = output_genes.stat().st_size / (1024 * 1024)
    print(f"  File size: {file_size_mb:.1f} MB")
    print()
    
    # Save TEs only
    output_tes = output_dir / "combined_counts_matrix_TEs_only.tsv"
    print(f"Saving TEs only matrix...")
    print(f"  -> {output_tes}")
    combined[is_te].to_csv(output_tes, sep='\t')
    file_size_mb = output_tes.stat().st_size / (1024 * 1024)
    print(f"  File size: {file_size_mb:.1f} MB")
    print()
    
    # Summary statistics
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total samples: {combined.shape[1]}")
    print(f"Total features: {combined.shape[0]:,}")
    print(f"  Genes: {n_genes:,}")
    print(f"  TEs: {n_tes:,}")
    print()
    
    # Show most highly expressed features
    total_counts = combined.sum(axis=1).sort_values(ascending=False)
    
    print("Top 10 most expressed genes:")
    top_genes = total_counts[is_gene].head(10)
    for gene, count in top_genes.items():
        print(f"  {gene}: {count:,}")
    print()
    
    print("Top 10 most expressed TEs:")
    top_tes = total_counts[is_te].head(10)
    for te, count in top_tes.items():
        print(f"  {te}: {count:,}")
    print()
    
    print("=" * 70)
    print("COMPLETE!")
    print("=" * 70)
    print()
    print("Output files:")
    print(f"  1. All features: {output_all}")
    print(f"  2. Genes only: {output_genes}")
    print(f"  3. TEs only: {output_tes}")
    print()

def main():
    parser = argparse.ArgumentParser(
        description='Combine TEtranscripts count tables from bulk RNA-seq',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect all .cntTable files
  python3 %(prog)s -i ~/sc/tetranscripts_bulk -o ~/sc/tetranscripts_bulk/data
  
  # Use metadata to filter samples
  python3 %(prog)s -i ~/sc/tetranscripts_bulk -m metadata.csv -o ~/sc/tetranscripts_bulk/data
"""
    )
    
    parser.add_argument('-i', '--input-dir', type=Path, required=True,
                        help='Directory containing .cntTable files')
    parser.add_argument('-m', '--metadata', type=Path, required=False,
                        help='Metadata CSV file with sample information (optional)')
    parser.add_argument('-o', '--output-dir', type=Path, required=True,
                        help='Output directory for combined matrices')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.input_dir.exists():
        print(f"Error: Input directory not found: {args.input_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Run combination
    combine_count_tables(args.input_dir, args.metadata, args.output_dir)

if __name__ == '__main__':
    main()
