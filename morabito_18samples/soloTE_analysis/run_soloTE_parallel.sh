#!/bin/bash
#
# Run SoloTE in PARALLEL on multiple samples
# Uses GNU parallel or runs multiple background jobs
# Much faster with 120GB RAM available!

set -e

SAMPLE_MAP="../sample_mapping.csv"
BAM_DIR="/home/jacobc/sc/starsolo_aligned"
OUTPUT_DIR="soloTE_output"
TE_BED="SoloTE/hg38_rmsk.bed"
THREADS_PER_JOB=4  # 4 threads per sample
MAX_PARALLEL=6     # Run 6 samples simultaneously (4*6=24 threads total)

echo "======================================================================"
echo "SOLOTE PARALLEL PROCESSING: 18 Samples"
echo "======================================================================"

# Check if setup was run
if [ ! -f "$TE_BED" ]; then
    echo "✗ TE annotation not found: $TE_BED"
    echo "Run setup first: bash setup_soloTE.sh"
    exit 1
fi

# Count samples
N_SAMPLES=$(tail -n +2 "$SAMPLE_MAP" | wc -l)
echo "Processing $N_SAMPLES samples"
echo "Parallelization: $MAX_PARALLEL samples at once, $THREADS_PER_JOB threads each"
echo "Total threads used: $((MAX_PARALLEL * THREADS_PER_JOB))"
echo ""

# Function to process one sample
process_sample() {
    local sample_id=$1
    local srr=$2
    
    BAM_FILE="$BAM_DIR/$srr/Aligned.sortedByCoord.out.bam"
    SAMPLE_OUTPUT="$OUTPUT_DIR/$srr"
    
    # Check if already processed
    if [ -f "$SAMPLE_OUTPUT/${srr}_TE_counts.mtx" ]; then
        echo "[$srr] Already processed, skipping"
        return 0
    fi
    
    # Check BAM exists
    if [ ! -f "$BAM_FILE" ]; then
        echo "[$srr] BAM file not found: $BAM_FILE"
        return 1
    fi
    
    echo "[$srr] Started: $(date '+%H:%M:%S')"
    
    SAMPLE_START=$(date +%s)
    
    # Run SoloTE
    python3 SoloTE/SoloTE_pipeline.py \
        --threads $THREADS_PER_JOB \
        --bam "$BAM_FILE" \
        --teannotation "$TE_BED" \
        --outputprefix "$srr" \
        --outputdir "$SAMPLE_OUTPUT" \
        > "${SAMPLE_OUTPUT}_solote.log" 2>&1
    
    SAMPLE_END=$(date +%s)
    SAMPLE_TIME=$((SAMPLE_END - SAMPLE_START))
    
    echo "[$srr] Completed in $(($SAMPLE_TIME / 60)) min $(($SAMPLE_TIME % 60)) sec"
    
    return 0
}

export -f process_sample
export BAM_DIR OUTPUT_DIR TE_BED THREADS_PER_JOB

# Create job list
tail -n +2 "$SAMPLE_MAP" | while IFS=, read -r sample_id gsm srx srr; do
    echo "$sample_id $srr"
done > /tmp/solote_jobs.txt

START_TIME=$(date +%s)

# Check if GNU parallel is available
if command -v parallel &> /dev/null; then
    echo "Using GNU parallel for efficient parallelization"
    echo ""
    
    cat /tmp/solote_jobs.txt | parallel --colsep ' ' -j $MAX_PARALLEL --progress \
        process_sample {1} {2}
else
    echo "GNU parallel not found, using background jobs"
    echo "(Install with: sudo apt install parallel for better progress tracking)"
    echo ""
    
    # Use background jobs
    RUNNING=0
    while IFS=' ' read -r sample_id srr; do
        # Wait if we have too many running
        while [ $RUNNING -ge $MAX_PARALLEL ]; do
            sleep 5
            RUNNING=$(jobs -r | wc -l)
        done
        
        # Start job in background
        process_sample "$sample_id" "$srr" &
        RUNNING=$((RUNNING + 1))
        
        echo "Active jobs: $RUNNING/$MAX_PARALLEL"
        sleep 2
    done < /tmp/solote_jobs.txt
    
    # Wait for all jobs to finish
    echo ""
    echo "Waiting for remaining jobs to complete..."
    wait
fi

END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo ""
echo "======================================================================"
echo "✓ PARALLEL PROCESSING COMPLETE"
echo "======================================================================"
echo "Total time: $(($TOTAL_TIME / 3600)) hours $(($TOTAL_TIME % 3600 / 60)) min"
echo "Average per sample: $(($TOTAL_TIME / $N_SAMPLES / 60)) min"
echo ""

# Show summary
N_COMPLETED=$(find $OUTPUT_DIR -name "*_TE_counts.mtx" 2>/dev/null | wc -l)
echo "Completed: $N_COMPLETED/$N_SAMPLES samples"
echo ""

if [ $N_COMPLETED -eq $N_SAMPLES ]; then
    echo "✓ All samples completed successfully!"
    echo ""
    echo "Next steps:"
    echo "  1. Merge count matrices: python3 01_merge_soloTE_samples.py"
    echo "  2. Or run full pipeline: bash run_full_pipeline.sh"
else
    echo "⚠️  Some samples failed or are incomplete"
    echo "Check logs in: soloTE_output/*_solote.log"
fi
echo ""

rm -f /tmp/solote_jobs.txt
