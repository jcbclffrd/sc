#!/bin/bash
echo "=== Download Status ==="
echo "Process running: $(ps aux | grep download_samples.sh | grep -v grep | wc -l)"
echo ""
echo "=== Recent log output ==="
tail -30 download.log
echo ""
echo "=== Downloaded samples ==="
ls -lh sra_downloads/*/*.fastq.gz 2>/dev/null | wc -l
echo " FASTQ files found"
echo ""
echo "Use: tail -f download.log  (to monitor in real-time)"
