#!/bin/bash

#######################################################
# Build scTE Index for TE Quantification
# This creates an exclusive index for TEs
#######################################################

set -e

# Activate virtual environment
source /home/jacobc/hcaTE/.venv/bin/activate

# Paths
GTF="/home/jacobc/sc/annotations/gencode.v45.primary_assembly.annotation.gtf"
TE_GTF="/home/jacobc/sc/annotations/hg38_rmsk.gtf"
OUTPUT_PREFIX="/home/jacobc/sc/annotations/hg38_gencode45"

echo "================================================"
echo "Building scTE Index"
echo "Started: $(date)"
echo "================================================"
echo ""
echo "Gene annotation: $GTF"
echo "TE annotation: $TE_GTF"
echo "Output prefix: $OUTPUT_PREFIX"
echo ""

# Build exclusive index (TEs not overlapping genes)
# scTE has built-in hg38 TE annotations, just provide gene GTF
scTE_build \
    -g hg38 \
    -gene "$GTF" \
    -m exclusive \
    -o "$OUTPUT_PREFIX"

echo ""
echo "================================================"
echo "scTE Index Build Complete!"
echo "Finished: $(date)"
echo "================================================"
