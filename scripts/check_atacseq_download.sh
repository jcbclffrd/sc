#!/bin/bash

# Monitor ATAC-seq download progress
# Shows status of all 32 ATAC-seq samples

echo "=========================================="
echo "ATAC-seq Download Status Monitor"
echo "=========================================="
echo "Checked: $(date)"
echo

OUTPUT_DIR="sra_downloads"
LOG_DIR="logs_atacseq"

# All 32 ATAC-seq samples
ALL_SAMPLES=($(seq 129 160 | xargs -I{} echo "SRR14514{}"))

TOTAL=${#ALL_SAMPLES[@]}
DOWNLOADED=0
EMPTY=0
MISSING=0

echo "Sample Status:"
echo "----------------------------------------"

for sample in "${ALL_SAMPLES[@]}"; do
    if [ -d "$OUTPUT_DIR/$sample" ]; then
        # Check for FASTQ files
        FASTQ_COUNT=$(find "$OUTPUT_DIR/$sample" -name "*.fastq.gz" 2>/dev/null | wc -l)
        SIZE=$(du -sh "$OUTPUT_DIR/$sample" 2>/dev/null | cut -f1)
        
        if [ "$FASTQ_COUNT" -ge 2 ]; then
            echo "✓ $sample: $FASTQ_COUNT files ($SIZE)"
            DOWNLOADED=$((DOWNLOADED + 1))
        else
            echo "⚠ $sample: Empty or incomplete ($SIZE)"
            EMPTY=$((EMPTY + 1))
        fi
    else
        echo "✗ $sample: Not downloaded"
        MISSING=$((MISSING + 1))
    fi
done

echo
echo "=========================================="
echo "Summary:"
echo "----------------------------------------"
echo "Total ATAC-seq samples:     $TOTAL"
echo "✓ Downloaded with data:     $DOWNLOADED"
echo "⚠ Empty/incomplete:         $EMPTY"
echo "✗ Not downloaded:           $MISSING"
echo "=========================================="

# Calculate total size
if [ -d "$OUTPUT_DIR" ]; then
    echo
    TOTAL_SIZE=$(find "$OUTPUT_DIR" -name "SRR14514[1][2-6][0-9]" -type d -exec du -sh {} \; 2>/dev/null | awk '{sum+=$1} END {print sum}')
    echo "Estimated total size of ATAC-seq data:"
    du -ch "$OUTPUT_DIR"/SRR14514{129..160} 2>/dev/null | tail -1 | awk '{print "  "$1}'
fi

# Show active downloads
if [ -d "$LOG_DIR" ]; then
    echo
    echo "Recent activity (last 5 log files):"
    echo "----------------------------------------"
    ls -lt "$LOG_DIR"/*.log 2>/dev/null | head -5 | awk '{print "  "$9" ("$6, $7, $8")"}'
    
    # Check if any downloads are currently running
    ACTIVE=$(ps aux | grep -E "fasterq-dump.*SRR14514[1][2-6]" | grep -v grep | wc -l)
    if [ "$ACTIVE" -gt 0 ]; then
        echo
        echo "⟳ Active downloads: $ACTIVE"
        ps aux | grep -E "fasterq-dump.*SRR14514[1][2-6]" | grep -v grep | awk '{print "  "$11" (PID "$2")"}'
    fi
fi

echo
echo "=========================================="
echo "To start/resume download:"
echo "  ./scripts/download_atacseq_samples.sh"
echo
echo "To monitor in real-time:"
echo "  watch -n 30 ./scripts/check_atacseq_download.sh"
echo "=========================================="
