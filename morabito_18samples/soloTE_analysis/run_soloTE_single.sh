#!/bin/bash
#
# Run SoloTE on a single sample (for testing)
# Usage: bash run_soloTE_single.sh SRR14513984

if [ $# -eq 0 ]; then
    echo "Usage: $0 <SRR_ID>"
    echo "Example: $0 SRR14513984"
    exit 1
fi

SRR=$1
BAM_DIR="/home/jacobc/sc/starsolo_aligned"
OUTPUT_DIR="soloTE_output"
TE_BED="SoloTE/hg38_rmsk.bed"
THREADS=8

echo "======================================================================"
echo "SOLOTE SINGLE SAMPLE: $SRR"
echo "======================================================================"

# Check if setup was run
if [ ! -f "$TE_BED" ]; then
    echo "✗ TE annotation not found: $TE_BED"
    echo "Run setup first: bash setup_soloTE.sh"
    exit 1
fi

BAM_FILE="$BAM_DIR/$SRR/Aligned.sortedByCoord.out.bam"
SAMPLE_OUTPUT="$OUTPUT_DIR/$SRR"

# Check BAM exists
if [ ! -f "$BAM_FILE" ]; then
    echo "✗ BAM file not found: $BAM_FILE"
    exit 1
fi

echo "Input BAM: $(du -h $BAM_FILE | cut -f1)"
echo "Output dir: $SAMPLE_OUTPUT"
echo "Threads: $THREADS"
echo ""
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"

START_TIME=$(date +%s)

# Run SoloTE
python3 SoloTE/SoloTE_pipeline.py \
    --threads $THREADS \
    --bam "$BAM_FILE" \
    --teannotation "$TE_BED" \
    --outputprefix "$SRR" \
    --outputdir "$SAMPLE_OUTPUT"

END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo ""
echo "======================================================================"
echo "✓ COMPLETE"
echo "======================================================================"
echo "Time: $(($TOTAL_TIME / 60)) min $(($TOTAL_TIME % 60)) sec"
echo ""

if [ -d "$SAMPLE_OUTPUT" ]; then
    echo "Output files:"
    ls -lh "$SAMPLE_OUTPUT" | tail -n +2 | awk '{print "  " $5 " " $9}'
fi

echo ""
