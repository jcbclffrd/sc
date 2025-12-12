#!/bin/bash
#
# Master script to run the complete Morabito 18-sample analysis pipeline.
#
# This script executes all steps in sequence:
# 1. Map Sample-X names to SRR IDs
# 2. Extract scTE output for 18 samples
# 3. Filter to Morabito's exact cell barcodes
# 4. Merge samples and add metadata
# 5. Differential analysis (AD vs Control)

set -e  # Exit on error

echo "========================================================================"
echo "MORABITO 18-SAMPLE ANALYSIS PIPELINE"
echo "========================================================================"
echo ""
echo "This pipeline will:"
echo "  1. Map Sample-X names to SRR accessions"
echo "  2. Extract scTE data for 18 samples"
echo "  3. Filter to Morabito's cell barcodes"
echo "  4. Merge samples and add metadata"
echo "  5. Perform differential analysis"
echo ""
echo "Working directory: $(pwd)"
echo "========================================================================"
echo ""

# Activate virtual environment if needed
if [ -d "../../hcaTE/.venv" ]; then
    echo "Activating virtual environment..."
    source ../../hcaTE/.venv/bin/activate
fi

# Step 1: Map samples
echo ""
echo "========================================================================"
echo "STEP 1: Mapping Sample-X to SRR IDs"
echo "========================================================================"
python3 01_map_samples.py
if [ $? -ne 0 ]; then
    echo "✗ Step 1 failed!"
    exit 1
fi

# Step 2: Extract 18 samples
echo ""
echo "========================================================================"
echo "STEP 2: Extracting scTE output for 18 samples"
echo "========================================================================"
python3 02_extract_18samples.py
if [ $? -ne 0 ]; then
    echo "✗ Step 2 failed!"
    exit 1
fi

# Step 3: Filter to Morabito's cells
echo ""
echo "========================================================================"
echo "STEP 3: Filtering to Morabito's cell barcodes"
echo "========================================================================"
python3 03_filter_to_morabito_cells.py
if [ $? -ne 0 ]; then
    echo "✗ Step 3 failed!"
    exit 1
fi

# Step 4: Merge and QC
echo ""
echo "========================================================================"
echo "STEP 4: Merging samples and adding metadata"
echo "========================================================================"
python3 04_merge_and_qc.py
if [ $? -ne 0 ]; then
    echo "✗ Step 4 failed!"
    exit 1
fi

# Step 5: Differential analysis
echo ""
echo "========================================================================"
echo "STEP 5: Differential expression analysis"
echo "========================================================================"
python3 05_differential_analysis.py
if [ $? -ne 0 ]; then
    echo "✗ Step 5 failed!"
    exit 1
fi

echo ""
echo "========================================================================"
echo "✓ PIPELINE COMPLETE!"
echo "========================================================================"
echo ""
echo "Output files:"
echo "  - sample_mapping.csv"
echo "  - scte_18samples/ (extracted samples)"
echo "  - scte_18samples_filtered/ (filtered to Morabito's cells)"
echo "  - merged_18samples.h5ad (genes + TEs)"
echo "  - merged_18samples_genes.h5ad (genes only)"
echo "  - differential_results/ (DE analysis results)"
echo ""
echo "Check differential_results/ for:"
echo "  - genes_AD_vs_Control_all_cells.csv"
echo "  - TEs_AD_vs_Control_all_cells.csv"
echo "  - <celltype>_AD_vs_Control.csv"
echo "  - figures/volcano_plots.png"
echo ""
echo "========================================================================"
