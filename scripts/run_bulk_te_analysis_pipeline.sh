#!/bin/bash

# Run complete bulk RNA-seq TE analysis pipeline
# 1. Combine TEcount output tables
# 2. Run differential expression analysis (AD vs Control)

set -e

# =============================================================================
# CONFIGURATION
# =============================================================================

# Directories (adjust these for your server paths)
TECOUNT_DIR=~/sc/tetranscripts_bulk
OUTPUT_DIR=~/sc/tetranscripts_bulk
SCRIPT_DIR=~/sc/scripts

# Python from virtual environment
PYTHON_BIN=~/sc/tetranscripts_env/bin/python3

# Metadata files
METADATA_18=~/sc/tetranscripts_bulk_metadata.csv
METADATA_ALL=${OUTPUT_DIR}/bulk_metadata_all_samples.csv

# Output subdirectories
DATA_DIR=${OUTPUT_DIR}/data
RESULTS_DIR=${OUTPUT_DIR}/results

# Create output directories
mkdir -p ${DATA_DIR}
mkdir -p ${RESULTS_DIR}

echo "========================================================================"
echo "BULK RNA-SEQ TE ANALYSIS PIPELINE - ALL SAMPLES"
echo "========================================================================"
echo ""

# Count available .cntTable files
NFILES=$(find ${TECOUNT_DIR} -name "*.cntTable" | wc -l)

echo "Configuration:"
echo "  TEcount directory: ${TECOUNT_DIR}"
echo "  Available samples: ${NFILES}"
echo "  Output directory: ${OUTPUT_DIR}"
echo "  Start time: $(date)"
echo "========================================================================"
echo ""

# =============================================================================
# STEP 0: CREATE COMPREHENSIVE METADATA (IF NEEDED)
# =============================================================================

if [ ! -f "${METADATA_ALL}" ]; then
    echo "========================================================================"
    echo "STEP 0: CREATING COMPREHENSIVE METADATA"
    echo "========================================================================"
    echo ""
    
    ${PYTHON_BIN} ${SCRIPT_DIR}/create_bulk_metadata_all_samples.py \
        --tecount-dir ${TECOUNT_DIR} \
        --existing-metadata ${METADATA_18} \
        --output ${METADATA_ALL}
    
    if [ $? -eq 0 ]; then
        echo ""
        echo "✓ Metadata created successfully"
        echo ""
    else
        echo ""
        echo "✗ Error creating metadata"
        exit 1
    fi
else
    echo "Using existing metadata: ${METADATA_ALL}"
    echo ""
fi

# =============================================================================
# STEP 1: COMBINE COUNT TABLES
# =============================================================================

echo "========================================================================"
echo "STEP 1: COMBINING COUNT TABLES"
echo "========================================================================"
echo ""

${PYTHON_BIN} ${SCRIPT_DIR}/combine_tetranscripts_counts_bulk.py \
    --input-dir ${TECOUNT_DIR} \
    --output-dir ${DATA_DIR}

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Count tables combined successfully"
    echo ""
else
    echo ""
    echo "✗ Error combining count tables"
    exit 1
fi

# =============================================================================
# STEP 2: DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

echo "========================================================================"
echo "STEP 2: DIFFERENTIAL EXPRESSION ANALYSIS"
echo "========================================================================"
echo ""

# Check if combined matrix exists
COMBINED_MATRIX=${DATA_DIR}/combined_counts_matrix_TEs_only.tsv

if [ ! -f "${COMBINED_MATRIX}" ]; then
    echo "✗ Error: Combined count matrix not found: ${COMBINED_MATRIX}"
    exit 1
fi

echo "Running DESeq2 analysis on TEs..."
echo ""

Rscript ${SCRIPT_DIR}/differential_TE_analysis_bulk.R \
    ${COMBINED_MATRIX} \
    ${METADATA_ALL} \
    ${RESULTS_DIR}

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Differential expression analysis complete"
    echo ""
else
    echo ""
    echo "✗ Error in differential expression analysis"
    exit 1
fi

# =============================================================================
# SUMMARY
# =============================================================================

echo "========================================================================"
echo "PIPELINE COMPLETE!"
echo "========================================================================"
echo ""
echo "Output files:"
echo ""
echo "Count Matrices:"
echo "  ${DATA_DIR}/combined_counts_matrix.tsv"
echo "  ${DATA_DIR}/combined_counts_matrix_genes_only.tsv"
echo "  ${DATA_DIR}/combined_counts_matrix_TEs_only.tsv"
echo ""
echo "Differential Expression Results:"
echo "  ${RESULTS_DIR}/differential_TE_results_AD_vs_Control.csv"
echo "  ${RESULTS_DIR}/volcano_plot_TE_AD_vs_Control.pdf"
echo "  ${RESULTS_DIR}/MA_plot_TE_AD_vs_Control.pdf"
echo "  ${RESULTS_DIR}/PCA_plot_TE_expression.pdf"
echo ""
echo "End time: $(date)"
echo "========================================================================"
