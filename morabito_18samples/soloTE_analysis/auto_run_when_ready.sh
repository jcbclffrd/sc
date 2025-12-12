#!/bin/bash
#
# Auto-run pipeline once SoloTE completes
# Checks every hour if SoloTE is done, then runs full analysis

echo "======================================================================"
echo "AUTO-RUN: Waiting for SoloTE to Complete"
echo "======================================================================"
echo ""
echo "This script will:"
echo "  1. Check every hour if SoloTE batch is complete"
echo "  2. Automatically run the full pipeline when ready"
echo "  3. Generate all comparison CSV files"
echo ""
echo "Started: $(date)"
echo "Logging to: auto_pipeline.log"
echo ""
echo "To stop: kill $(echo $$)"
echo ""

# Log file
LOG="auto_pipeline.log"
echo "=== Auto-pipeline started: $(date) ===" >> $LOG

# Check every hour
while true; do
    # Count completed
    N_COMPLETED=$(find soloTE_output -name "*_TE_counts.mtx" 2>/dev/null | wc -l)
    
    echo "[$(date '+%H:%M')] Progress: $N_COMPLETED/18 samples" | tee -a $LOG
    
    # Check if SoloTE is still running
    if pgrep -f "run_soloTE_batch.sh" > /dev/null || pgrep -f "SoloTE_pipeline.py" > /dev/null; then
        echo "  Status: SoloTE still running..." | tee -a $LOG
    elif [ $N_COMPLETED -eq 18 ]; then
        echo "  Status: SoloTE COMPLETE! Starting pipeline..." | tee -a $LOG
        echo "" | tee -a $LOG
        
        # Run full pipeline
        source ~/hcaTE/.venv/bin/activate
        bash run_full_pipeline.sh 2>&1 | tee -a $LOG
        
        echo "" | tee -a $LOG
        echo "=== Pipeline complete: $(date) ===" | tee -a $LOG
        echo "" | tee -a $LOG
        echo "✓ DONE! Check comparison/ directory for results." | tee -a $LOG
        
        exit 0
    else
        echo "  Status: SoloTE stopped but incomplete ($N_COMPLETED/18)" | tee -a $LOG
        echo "  Waiting for restart..." | tee -a $LOG
    fi
    
    # Wait 1 hour
    sleep 3600
done
