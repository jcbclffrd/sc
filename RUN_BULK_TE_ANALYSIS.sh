#!/bin/bash
# Quick commands to run the bulk TE analysis pipeline on spark-bd86
# Analyzes ALL available bulk RNA-seq samples (~156 samples)

# Navigate to working directory
cd ~/sc

# Option 1: Run the complete pipeline (recommended)
# This will automatically detect all .cntTable files and run the full analysis
bash scripts/run_bulk_te_analysis_pipeline.sh

# Option 2: Run steps individually

# Step 0: Create metadata for all samples (optional - pipeline does this automatically)
python3 scripts/create_bulk_metadata_all_samples.py \
    --tecount-dir ~/sc/tetranscripts_bulk \
    --existing-metadata ~/sc/tetranscripts_bulk_metadata.csv \
    --output ~/sc/tetranscripts_bulk/bulk_metadata_all_samples.csv

# Step 1: Combine count tables (auto-detects all .cntTable files)
python3 scripts/combine_tetranscripts_counts_bulk.py \
    --input-dir ~/sc/tetranscripts_bulk \
    --output-dir ~/sc/tetranscripts_bulk/data

# Step 2: Differential expression analysis
Rscript scripts/differential_TE_analysis_bulk.R \
    ~/sc/tetranscripts_bulk/data/combined_counts_matrix_TEs_only.tsv \
    ~/sc/tetranscripts_bulk/bulk_metadata_all_samples.csv \
    ~/sc/tetranscripts_bulk/results

# Check outputs
echo "Count matrices:"
ls -lh ~/sc/tetranscripts_bulk/data/

echo ""
echo "Results:"
ls -lh ~/sc/tetranscripts_bulk/results/

echo ""
echo "Sample counts:"
wc -l ~/sc/tetranscripts_bulk/bulk_metadata_all_samples.csv
