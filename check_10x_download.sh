#!/bin/bash

echo "========================================"
echo "10x snRNA-seq Download Status"
echo "========================================"
echo "Time: $(date)"
echo ""

# Check if process is running
RUNNING=$(ps aux | grep download_10x_snRNA_samples.sh | grep -v grep | wc -l)
if [ $RUNNING -gt 0 ]; then
    echo "Status: ✓ RUNNING"
else
    echo "Status: ✗ NOT RUNNING (may be completed or stopped)"
fi
echo ""

# Count completed downloads from the 10x list
COMPLETED=0
while read srr; do
    if ls sra_downloads/$srr/*.fastq.gz 1> /dev/null 2>&1; then
        COMPLETED=$((COMPLETED + 1))
    fi
done < snRNA_seq_10x_samples.txt

TOTAL=152
PERCENT=$((COMPLETED * 100 / TOTAL))

echo "Progress: $COMPLETED / $TOTAL 10x snRNA-seq samples ($PERCENT%)"
echo ""

# Show disk usage
echo "Disk Usage: $(du -sh sra_downloads 2>/dev/null | cut -f1)"
echo ""

# Show recent activity
echo "Recent log output:"
echo "----------------------------------------"
tail -20 download_10x.log
echo ""
echo "========================================"
echo "To monitor live: tail -f download_10x.log"
echo "To check individual sample logs: ls -lh logs_10x/"
echo "========================================"
