#!/bin/bash
echo "========================================"
echo "PRJNA729525 Download Status"
echo "========================================"
echo "Time: $(date)"
echo ""

# Check if process is running
RUNNING=$(ps aux | grep download_all_samples.sh | grep -v grep | wc -l)
if [ $RUNNING -gt 0 ]; then
    echo "Status: ✓ RUNNING"
else
    echo "Status: ✗ NOT RUNNING (may be completed or stopped)"
fi
echo ""

# Count completed downloads
COMPLETED=$(find sra_downloads -name "*_1.fastq.gz" 2>/dev/null | wc -l)
TOTAL=375
PERCENT=$((COMPLETED * 100 / TOTAL))

echo "Progress: $COMPLETED / $TOTAL samples ($PERCENT%)"
echo ""

# Show disk usage
echo "Disk Usage: $(du -sh sra_downloads 2>/dev/null | cut -f1)"
echo ""

# Show recent activity
echo "Recent log output:"
echo "----------------------------------------"
tail -20 download_all.log
echo ""
echo "========================================"
echo "To monitor live: tail -f download_all.log"
echo "To check individual sample logs: ls -lh logs/"
echo "========================================"
