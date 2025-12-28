#!/bin/bash
#
# Main script to run on spark-bd86 server
# This will process all remaining bulk RNA-seq samples with TEtranscripts
#
# IMPORTANT: This will take ~100-150 hours of compute time!
# - Each sample takes ~30-60 minutes
# - With 4 parallel jobs, expect ~1-2 days total runtime
#
# Usage (on server):
#   ssh spark-bd86
#   cd ~/sc
#   bash RUN_TETRANSCRIPTS_ALL_BULK.sh
#

# Run in a screen session (recommended for long-running jobs)
if [ -z "$STY" ]; then
    echo "WARNING: Not running in a screen session!"
    echo ""
    echo "Recommended: Start in screen to avoid interruption:"
    echo "  screen -S tetranscripts"
    echo "  bash $0"
    echo ""
    read -p "Continue anyway? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

cd ~/sc

echo "========================================================================"
echo "TEtranscripts BULK RNA-SEQ - PROCESS ALL REMAINING SAMPLES"
echo "========================================================================"
echo ""
echo "This will:"
echo "  1. Find all BAM files in star_bulk_aligned/"
echo "  2. Skip samples that already have .cntTable files"
echo "  3. Run TEtranscripts on remaining ~138 samples"
echo "  4. Process 4 samples in parallel"
echo ""
echo "Estimated time: 24-48 hours"
echo "Start time: $(date)"
echo ""
read -p "Press Enter to begin..."

# Run the processing script
bash scripts/run_tetranscripts_remaining_bulk.sh 2>&1 | tee logs/tetranscripts_all_bulk_$(date +%Y%m%d_%H%M%S).log

echo ""
echo "========================================================================"
echo "PROCESSING COMPLETE!"
echo "========================================================================"
echo ""
echo "End time: $(date)"
echo ""
echo "Check results:"
echo "  ls -lh tetranscripts_bulk/*.cntTable | wc -l"
echo ""
echo "Next: Run full analysis with all samples:"
echo "  bash RUN_BULK_TE_ANALYSIS.sh"
echo ""
