#!/bin/bash
#
# Master script to run full SoloTE analysis pipeline
# Steps:
#   1. Run SoloTE on 18 samples (batch_run via run_soloTE_batch.sh)
#   2. Merge SoloTE outputs
#   3. Differential expression analysis
#   4. Compare scTE vs SoloTE

set -e

# Change to script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "======================================================================"
echo "SOLOTE FULL ANALYSIS PIPELINE"
echo "======================================================================"

# Activate venv
source ~/hcaTE/.venv/bin/activate

# Check if SoloTE batch is still running
if pgrep -f "run_soloTE_batch.sh" > /dev/null; then
    echo ""
    echo "⏳ SoloTE batch processing is currently running..."
    echo "   Check progress: tail -f batch_run.log"
    echo "   This script will wait for it to complete..."
    echo ""
    
    # Wait for SoloTE to finish
    while pgrep -f "run_soloTE_batch.sh" > /dev/null; do
        sleep 60
        echo "   Still running... ($(date '+%H:%M:%S'))"
    done
    
    echo "✓ SoloTE batch processing completed!"
else
    # Check if SoloTE outputs exist (soloTE creates subfamilytes_MATRIX/matrix.mtx)
    N_COMPLETED=$(find soloTE_output -name "matrix.mtx" -path "*subfamilytes_MATRIX*" 2>/dev/null | wc -l)
    
    if [ $N_COMPLETED -eq 0 ]; then
        echo "✗ No SoloTE outputs found"
        echo "Run: bash run_soloTE_parallel.sh first"
        exit 1
    elif [ $N_COMPLETED -lt 18 ]; then
        echo "⚠️  Only $N_COMPLETED/18 samples processed"
        echo "   Continuing anyway..."
    else
        echo "✓ All 18 SoloTE outputs found"
    fi
fi

echo ""
echo "======================================================================"
echo "STEP 2: MERGING SOLOTE SAMPLES"
echo "======================================================================"

python3 01_merge_soloTE_samples.py

if [ ! -f "merged_soloTE_18samples.h5ad" ]; then
    echo "✗ Merge failed"
    exit 1
fi

echo ""
echo "======================================================================"
echo "STEP 3: DIFFERENTIAL EXPRESSION ANALYSIS (SoloTE)"
echo "======================================================================"

python3 02_differential_analysis_soloTE.py

if [ ! -d "differential_results_soloTE" ]; then
    echo "✗ Differential analysis failed"
    exit 1
fi

echo ""
echo "======================================================================"
echo "STEP 4: COMPARING scTE vs SoloTE"
echo "======================================================================"

python3 03_compare_scTE_vs_soloTE.py

if [ ! -d "comparison" ]; then
    echo "✗ Comparison failed"
    exit 1
fi

echo ""
echo "======================================================================"
echo "STEP 5: GENERATING COMPARISON REPORT"
echo "======================================================================"

python3 04_generate_comparison_report.py

if [ ! -f "COMPARISON_REPORT.md" ]; then
    echo "✗ Report generation failed"
    exit 1
fi

echo ""
echo "======================================================================"
echo "✓ FULL PIPELINE COMPLETE!"
echo "======================================================================"
echo ""
echo "Generated files:"
echo "  - merged_soloTE_18samples.h5ad"
echo "  - differential_results_soloTE/"
echo "  - comparison/00_summary_statistics.csv"
echo "  - comparison/01_TE_detection_comparison.csv"
echo "  - comparison/02_expression_comparison.csv"
echo "  - COMPARISON_REPORT.md ← Ready for PI presentation!"
echo "  - comparison/03_differential_expression_comparison.csv"
echo ""
