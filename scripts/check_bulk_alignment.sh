#!/bin/bash

#####################################################################
# Monitor Bulk RNA-seq STAR Alignment Progress
#####################################################################

OUTPUT_DIR="/home/jacobc/sc/star_bulk_aligned"
SAMPLE_LIST="/home/jacobc/sc/data/bulk_rnaseq_samples.txt"

TOTAL=$(wc -l < "$SAMPLE_LIST")

echo "=== Bulk RNA-seq STAR Alignment Status ==="
echo "Time: $(date)"
echo ""

# Count completed alignments
COMPLETED=$(find "$OUTPUT_DIR" -name "Aligned.sortedByCoord.out.bam" -size +1M 2>/dev/null | wc -l)

echo "Progress: $COMPLETED / $TOTAL samples completed"
echo "Percentage: $(awk "BEGIN {printf \"%.1f\", ($COMPLETED/$TOTAL)*100}")%"
echo ""

# Check if alignment is currently running
if pgrep -f "run_star_bulk_alignment.sh" > /dev/null; then
    echo "Status: RUNNING ✓"
    echo ""
    
    # Show most recent sample being processed
    if [[ -f "${OUTPUT_DIR}/bulk_alignment.log" ]]; then
        echo "Recent log entries:"
        tail -20 "${OUTPUT_DIR}/bulk_alignment.log" | grep -E "Processing sample|SUCCESS|FAILED"
    fi
else
    echo "Status: NOT RUNNING"
    echo ""
    
    if [[ $COMPLETED -eq $TOTAL ]]; then
        echo "🎉 All samples completed!"
    elif [[ $COMPLETED -gt 0 ]]; then
        echo "⚠️  Alignment stopped with $COMPLETED/$TOTAL samples completed"
        echo "Run: ./scripts/run_star_bulk_alignment.sh to resume"
    else
        echo "Ready to start alignment"
        echo "Run: ./scripts/run_star_bulk_alignment.sh"
    fi
fi

echo ""
echo "Recent BAM files created:"
find "$OUTPUT_DIR" -name "Aligned.sortedByCoord.out.bam" -type f -printf '%T@ %p\n' 2>/dev/null | \
    sort -rn | head -5 | while read -r timestamp file; do
    size=$(stat -c%s "$file")
    echo "  $(basename $(dirname "$file")): $(numfmt --to=iec-i --suffix=B $size)"
done

echo ""
echo "Total disk usage: $(du -sh $OUTPUT_DIR 2>/dev/null | cut -f1)"
