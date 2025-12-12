#!/bin/bash
#
# Quick reference for running the Morabito 18-sample analysis
#

# Navigate to the analysis directory
cd /home/jacobc/sc/morabito_18samples

# Run the complete pipeline (all steps)
./run_pipeline.sh

# OR run steps individually:

# Step 1: Map Sample-X to SRR IDs
python3 01_map_samples.py

# Step 2: Extract 18 samples from scTE_output/
python3 02_extract_18samples.py

# Step 3: Filter to Morabito's exact cell barcodes
python3 03_filter_to_morabito_cells.py

# Step 4: Merge samples and add metadata
python3 04_merge_and_qc.py

# Step 5: Differential analysis (AD vs Control)
python3 05_differential_analysis.py

# View results
ls -lh differential_results/
cat differential_results/celltype_summary.csv
