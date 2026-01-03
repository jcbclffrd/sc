#!/bin/bash

# Simple test: Count TEs from CellRanger BAM using featureCounts
# This tests if the basic TE counting works before complex multimapper pipeline

set -e

SAMPLE_ID="${1:-SRR14514130}"
CELLRANGER_BAM="$HOME/sc/cellranger_atac_output/$SAMPLE_ID/outs/possorted_bam.bam"
TE_GTF="$HOME/sc/annotations/hg38_rmsk_TE.gtf"
OUTPUT_DIR="$HOME/sc/test_te_counts"
THREADS=8

mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "Testing TE Quantification - Simple Approach"
echo "=========================================="
echo "Sample: $SAMPLE_ID"
echo "Input BAM: $CELLRANGER_BAM"
echo "TE GTF: $TE_GTF"
echo ""

# Check inputs
if [ ! -f "$CELLRANGER_BAM" ]; then
    echo "ERROR: CellRanger BAM not found: $CELLRANGER_BAM"
    exit 1
fi

if [ ! -f "$TE_GTF" ]; then
    echo "ERROR: TE GTF not found: $TE_GTF"
    exit 1
fi

# Test 1: Basic featureCounts on whole BAM
echo "[Test 1] Running featureCounts on whole BAM (pseudobulk)..."
featureCounts \
    -a "$TE_GTF" \
    -o "$OUTPUT_DIR/${SAMPLE_ID}_pseudobulk_counts.txt" \
    -F GTF \
    -t exon \
    -g gene_id \
    -p \
    --countReadPairs \
    -M \
    --fraction \
    -T "$THREADS" \
    "$CELLRANGER_BAM" \
    2>&1 | tee "$OUTPUT_DIR/${SAMPLE_ID}_featurecounts.log"

echo ""
echo "✓ Test 1 complete"
echo ""

# Show summary
echo "Results:"
echo "  - Count file: $OUTPUT_DIR/${SAMPLE_ID}_pseudobulk_counts.txt"
echo "  - TEs quantified: $(tail -n +3 $OUTPUT_DIR/${SAMPLE_ID}_pseudobulk_counts.txt | wc -l)"
echo "  - Total counts: $(tail -n +3 $OUTPUT_DIR/${SAMPLE_ID}_pseudobulk_counts.txt | cut -f7 | awk '{sum+=$1} END {print sum}')"
echo ""
echo "This proves TE counting works on the CellRanger BAM!"
echo "Next step: Implement per-cell counting using CB tags"
