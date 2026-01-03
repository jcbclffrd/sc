#!/bin/bash

# Count TE reads using featureCounts (simpler than TEtranscripts)

set -e

BAM_DIR="/home/jacobc/sc/atacseq_aligned"
OUTPUT_DIR="/home/jacobc/sc/atacseq_te_counts"
GTF_TES="/home/jacobc/sc/annotations/hg38_rmsk_tetranscripts.gtf"

mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "TE Counting with featureCounts"
echo "=========================================="
echo "BAM directory: $BAM_DIR"
echo "TE GTF: $GTF_TES"
echo "Output: $OUTPUT_DIR"
echo "=========================================="
echo ""

# Test on SRR14514140
SAMPLE="SRR14514140"
BAM="$BAM_DIR/$SAMPLE/Aligned.sortedByCoord.out.bam"

if [[ ! -f "$BAM" ]]; then
    echo "Error: BAM file not found: $BAM"
    exit 1
fi

echo "Counting TEs for: $SAMPLE"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"

# Use featureCounts with multi-mapping support
featureCounts \
    -p \
    -M \
    -O \
    -T 8 \
    -t exon \
    -g gene_id \
    -a "$GTF_TES" \
    -o "$OUTPUT_DIR/${SAMPLE}_te_counts.txt" \
    "$BAM" \
    2>&1 | tee "$OUTPUT_DIR/${SAMPLE}_featurecounts.log"

echo ""
echo "✓ TE counting complete at $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "Results:"
ls -lh "$OUTPUT_DIR/${SAMPLE}_te_counts.txt"*
echo ""
echo "Top 10 TEs by count:"
tail -n +3 "$OUTPUT_DIR/${SAMPLE}_te_counts.txt" | sort -k7 -rn | head -10 | cut -f1,7
