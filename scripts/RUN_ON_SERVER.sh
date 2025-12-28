#!/bin/bash
# Run this script ON THE SPARK SERVER: ssh spark-bd86, then bash this script

cd ~/sc

# Remove tidyverse requirement from R script (if not already done)
sed -i '/library(tidyverse)/d' scripts/differential_TE_analysis_bulk.R 2>/dev/null

echo "========================================================================"
echo "RUNNING BULK TE ANALYSIS - Starting at $(date)"
echo "========================================================================"

# Run the full pipeline
bash scripts/run_bulk_te_analysis_pipeline.sh 2>&1 | tee bulk_analysis_run.log

echo ""
echo "========================================================================"
echo "COMPLETE - Finished at $(date)"
echo "========================================================================"
echo ""
echo "Check results at:"
echo "  - tetranscripts_bulk/data/combined_counts_matrix_TEs_only.tsv"
echo "  - tetranscripts_bulk/results/differential_TE_results_AD_vs_Control.csv"
echo "  - tetranscripts_bulk/results/*.pdf"
echo ""
