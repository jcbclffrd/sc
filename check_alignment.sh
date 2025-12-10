#!/bin/bash

echo "========================================"
echo "STARsolo Alignment Status"
echo "========================================"
echo "Time: $(date)"
echo ""

# Check if process is running
RUNNING=$(ps aux | grep run_starsolo_alignment.sh | grep -v grep | wc -l)
if [ $RUNNING -gt 0 ]; then
    echo "Status: ✓ RUNNING"
    PID=$(ps aux | grep run_starsolo_alignment.sh | grep -v grep | awk '{print $2}')
    echo "PID: $PID"
else
    echo "Status: ✗ NOT RUNNING"
fi
echo ""

# Count completed alignments
COMPLETED=$(find starsolo_aligned -name "Aligned.sortedByCoord.out.bam" 2>/dev/null | wc -l)
TOTAL=$(wc -l < ready_for_alignment.txt)
PERCENT=$((COMPLETED * 100 / TOTAL))

echo "Progress: $COMPLETED / $TOTAL samples ($PERCENT%)"
echo ""

# Show disk usage
echo "Disk Usage:"
echo "  Alignments: $(du -sh starsolo_aligned 2>/dev/null | cut -f1)"
echo "  Total project: $(du -sh ~/sc 2>/dev/null | cut -f1)"
echo ""

# Show recent activity
echo "Recent log output:"
echo "----------------------------------------"
tail -30 starsolo_alignment.log
echo ""
echo "========================================"
echo "To monitor live: tail -f starsolo_alignment.log"
echo "To check specific sample: ls -lh starsolo_aligned/SRR*/Solo.out/Gene/filtered/"
echo "========================================"
