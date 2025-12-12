#!/bin/bash
#
# Monitor SoloTE batch progress
# Shows completion status and estimated time remaining

echo "======================================================================"
echo "SOLOTE PROGRESS MONITOR"
echo "======================================================================"

# Count completed samples
N_TOTAL=18
N_COMPLETED=$(find soloTE_output -name "*_TE_counts.mtx" 2>/dev/null | wc -l)

echo ""
echo "Completed: $N_COMPLETED/$N_TOTAL samples"
echo ""

if [ $N_COMPLETED -gt 0 ]; then
    echo "Completed samples:"
    find soloTE_output -name "*_TE_counts.mtx" | while read file; do
        dir=$(dirname "$file")
        srr=$(basename "$dir")
        size=$(du -h "$file" | cut -f1)
        echo "  ✓ $srr ($size)"
    done
    echo ""
fi

# Check if still running
if pgrep -f "run_soloTE_batch.sh" > /dev/null || pgrep -f "SoloTE_pipeline.py" > /dev/null; then
    echo "Status: ⏳ RUNNING"
    echo ""
    
    # Show current sample
    CURRENT_SAMPLE=$(pgrep -f "SoloTE_pipeline.py" -a | grep -oP 'SRR\d+' | head -1)
    if [ ! -z "$CURRENT_SAMPLE" ]; then
        echo "Current sample: $CURRENT_SAMPLE"
        echo ""
    fi
    
    # Show recent log
    echo "Recent log (last 15 lines):"
    echo "---"
    tail -15 batch_run.log
    echo "---"
    echo ""
    echo "Full log: tail -f batch_run.log"
else
    if [ $N_COMPLETED -eq $N_TOTAL ]; then
        echo "Status: ✓ COMPLETE"
        echo ""
        echo "Next step: bash run_full_pipeline.sh"
    else
        echo "Status: ✗ STOPPED (incomplete)"
        echo ""
        echo "Restart: bash run_soloTE_batch.sh"
    fi
fi

echo ""
