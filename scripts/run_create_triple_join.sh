#!/bin/bash

# Create triple-join file combining NEW bulk RNA-seq (18 samples) with 
# existing snRNA-seq + previous bulk join file

PYTHON_BIN=~/sc/tetranscripts_env/bin/python3

# Input files - UPDATED to correct paths
BULK_RESULTS=~/sc/tetranscripts_bulk/results/differential_TE_results_AD_vs_Control.csv
SNRNASEQ_JOIN=~/sc/tetranscripts_bulk/joined_TE_results_singlecell_bulk.csv
OUTPUT=~/sc/tetranscripts_bulk/results/triple_join_bulk_snrnaseq_TEs.csv

echo "========================================================================"
echo "CREATE TRIPLE-JOIN FILE"
echo "Combining: snRNA-seq + previous bulk + NEW bulk (18 samples)"
echo "========================================================================"
echo ""
echo "Inputs:"
echo "  NEW bulk results (18 samples): $BULK_RESULTS"
echo "  snRNA-seq + old bulk join:     $SNRNASEQ_JOIN"
echo "  Output:                        $OUTPUT"
echo ""

# Check if input files exist
if [ ! -f "$BULK_RESULTS" ]; then
    echo "Error: NEW bulk results file not found: $BULK_RESULTS"
    echo "Please run the bulk TE analysis pipeline first"
    exit 1
fi

if [ ! -f "$SNRNASEQ_JOIN" ]; then
    echo "Warning: snRNA-seq + bulk join file not found at: $SNRNASEQ_JOIN"
    echo ""
    echo "Please specify the correct path:"
    read -p "Path to snRNA-seq join file: " SNRNASEQ_JOIN
    
    if [ ! -f "$SNRNASEQ_JOIN" ]; then
        echo "Error: File not found: $SNRNASEQ_JOIN"
        exit 1
    fi
fi

# Run the join script
$PYTHON_BIN ~/sc/scripts/create_triple_join_bulk_snrnaseq.py \
    --bulk-results "$BULK_RESULTS" \
    --snrnaseq-join "$SNRNASEQ_JOIN" \
    --output "$OUTPUT"

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Triple-join file created successfully"
    echo "  Output: $OUTPUT"
    echo ""
    echo "This file combines three analyses:"
    echo "  1. snRNA-seq TE differential expression"
    echo "  2. Previous bulk RNA-seq analysis (bulk_* columns)"
    echo "  3. NEW bulk RNA-seq analysis - 18 samples (new_bulk_* columns)"
else
    echo ""
    echo "✗ Error creating triple-join file"
    exit 1
fi
