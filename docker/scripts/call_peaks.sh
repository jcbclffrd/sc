#!/bin/bash
# call_peaks.sh - Call peaks with MACS3
set -e

SAMPLE_ID=${SAMPLE_ID:-$1}
BAM_FILE=${BAM_FILE:-/data/${SAMPLE_ID}/Aligned.sortedByCoord.out.bam}
OUTPUT_DIR=${OUTPUT_DIR:-/output}

if [[ -z "$SAMPLE_ID" ]]; then
    echo "ERROR: SAMPLE_ID not set"
    exit 1
fi

if [[ ! -f "$BAM_FILE" ]]; then
    echo "ERROR: BAM file not found: $BAM_FILE"
    exit 1
fi

echo "=========================================="
echo "Calling peaks: $SAMPLE_ID"
echo "BAM: $BAM_FILE"
echo "=========================================="

PEAK_DIR="$OUTPUT_DIR/$SAMPLE_ID"
mkdir -p "$PEAK_DIR"

macs3 callpeak \
    -t "$BAM_FILE" \
    -f BAMPE \
    -g hs \
    -n "$SAMPLE_ID" \
    --outdir "$PEAK_DIR" \
    -q 0.05 \
    --keep-dup all \
    --call-summits

PEAKS=$(wc -l < "$PEAK_DIR/${SAMPLE_ID}_peaks.narrowPeak")
echo "✓ Complete: $PEAKS peaks called"
