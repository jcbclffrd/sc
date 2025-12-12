#!/bin/bash
#
# Monitor parallel SoloTE processing
# Shows real-time progress

OUTPUT_DIR="soloTE_output"

echo "======================================================================"
echo "SOLOTE PARALLEL PROCESSING MONITOR"
echo "======================================================================"
echo ""

# Check completed samples
N_COMPLETED=$(find $OUTPUT_DIR -name "*_TE_counts.mtx" 2>/dev/null | wc -l)
N_TOTAL=18

echo "✓ Completed: $N_COMPLETED/$N_TOTAL samples"
echo ""

# Check running processes
N_RUNNING=$(ps aux | grep "SoloTE_pipeline.py" | grep -v grep | wc -l)

if [ $N_RUNNING -gt 0 ]; then
    echo "⏳ Currently running: $N_RUNNING samples"
    echo ""
    
    # Show which samples are running
    echo "Active samples:"
    ps aux | grep "SoloTE_pipeline.py" | grep -v grep | while read line; do
        PID=$(echo $line | awk '{print $2}')
        ETIME=$(ps -p $PID -o etime= 2>/dev/null | tr -d ' ')
        
        # Find sample from outputdir argument
        OUTDIR=$(ps -p $PID -o args= | grep -oP '\-\-outputdir\s+\K[^ ]+' || echo "unknown")
        SAMPLE=$(basename $OUTDIR)
        
        if [ -n "$SAMPLE" ]; then
            echo "  [$SAMPLE] Running for $ETIME"
        fi
    done
    echo ""
fi

# Show recently completed (last 5 minutes)
echo "Recently completed (last 5 min):"
find $OUTPUT_DIR -name "*_TE_counts.mtx" -mmin -5 2>/dev/null | while read file; do
    SAMPLE=$(basename $(dirname $file))
    SIZE=$(ls -lh $file | awk '{print $5}')
    echo "  ✓ $SAMPLE ($SIZE)"
done
echo ""

# Resource usage
echo "System resources:"
free -h | grep "Mem:" | awk '{print "  RAM: " $3 " used / " $2 " total (" $7 " available)"}'
uptime | awk -F'load average:' '{print "  Load: " $2}'
echo ""

# Estimate time remaining
if [ $N_RUNNING -gt 0 ] && [ $N_COMPLETED -gt 0 ]; then
    # Estimate based on completed samples
    ELAPSED=$(ps -p $(pgrep -f run_soloTE_parallel.sh | head -1) -o etime= 2>/dev/null | tr -d ' ' || echo "0:00")
    ELAPSED_MIN=$(echo $ELAPSED | awk -F: '{if (NF==3) print $1*60+$2; else print $1}')
    
    if [ $ELAPSED_MIN -gt 0 ]; then
        AVG_TIME=$(($ELAPSED_MIN / ($N_COMPLETED + 1)))
        REMAINING=$(($N_TOTAL - $N_COMPLETED))
        
        # Account for parallelization
        PARALLEL=$N_RUNNING
        EST_MIN=$(($REMAINING * $AVG_TIME / $PARALLEL))
        EST_HR=$(($EST_MIN / 60))
        EST_MIN_REM=$(($EST_MIN % 60))
        
        echo "Estimated time remaining: ~${EST_HR}h ${EST_MIN_REM}m"
        echo "(Running $PARALLEL samples in parallel)"
        echo ""
    fi
fi

# Show last few log lines
if [ -f "parallel_batch.log" ]; then
    echo "Recent activity (last 5 lines):"
    tail -5 parallel_batch.log | sed 's/^/  /'
    echo ""
fi

echo "======================================================================"
echo "Run this script again to update: bash monitor_parallel.sh"
echo "Watch continuously: watch -n 10 bash monitor_parallel.sh"
echo "======================================================================"
