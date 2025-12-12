#!/bin/bash
#
# Run SoloTE on all 18 Morabito samples
# Processes BAM files and generates TE-annotated counts

set -e

SAMPLE_MAP="../sample_mapping.csv"
BAM_DIR="/home/jacobc/sc/starsolo_aligned"
OUTPUT_DIR="soloTE_output"
TE_BED="SoloTE/hg38_rmsk.bed"
THREADS=8

echo "======================================================================"
echo "SOLOTE BATCH PROCESSING: 18 Samples"
echo "======================================================================"

# Check if setup was run
if [ ! -f "$TE_BED" ]; then
    echo "✗ TE annotation not found: $TE_BED"
    echo "Run setup first: bash setup_soloTE.sh"
    exit 1
fi

# Count samples
N_SAMPLES=$(tail -n +2 "$SAMPLE_MAP" | wc -l)
echo "Processing $N_SAMPLES samples with $THREADS threads each"
echo ""

# Track progress
CURRENT=0
START_TIME=$(date +%s)

# Process each sample
while IFS=, read -r sample_id gsm srx srr; do
    if [ "$sample_id" = "sample_id" ]; then continue; fi  # Skip header
    
    CURRENT=$((CURRENT + 1))
    
    echo "======================================================================"
    echo "[$CURRENT/$N_SAMPLES] Processing: $sample_id ($srr)"
    echo "======================================================================"
    
    BAM_FILE="$BAM_DIR/$srr/Aligned.sortedByCoord.out.bam"
    SAMPLE_OUTPUT="$OUTPUT_DIR/$srr"
    
    # Check if already processed
    if [ -f "$SAMPLE_OUTPUT/${srr}_TE_counts.mtx" ]; then
        echo "✓ Already processed, skipping..."
        continue
    fi
    
    # Check BAM exists
    if [ ! -f "$BAM_FILE" ]; then
        echo "✗ BAM file not found: $BAM_FILE"
        continue
    fi
    
    echo "Input BAM: $(du -h $BAM_FILE | cut -f1)"
    echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
    
    SAMPLE_START=$(date +%s)
    
    # Run SoloTE
    python3 SoloTE/SoloTE_pipeline.py \
        --threads $THREADS \
        --bam "$BAM_FILE" \
        --teannotation "$TE_BED" \
        --outputprefix "$srr" \
        --outputdir "$SAMPLE_OUTPUT"
    
    SAMPLE_END=$(date +%s)
    SAMPLE_TIME=$((SAMPLE_END - SAMPLE_START))
    
    echo "✓ Completed in $(($SAMPLE_TIME / 60)) min $(($SAMPLE_TIME % 60)) sec"
    echo ""
    
    # Show output files
    if [ -d "$SAMPLE_OUTPUT" ]; then
        echo "Output files:"
        ls -lh "$SAMPLE_OUTPUT" | tail -n +2 | awk '{print "  " $5 " " $9}'
    fi
    
    echo ""
    
done < "$SAMPLE_MAP"

END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo "======================================================================"
echo "✓ BATCH PROCESSING COMPLETE"
echo "======================================================================"
echo "Total time: $(($TOTAL_TIME / 3600)) hours $(($TOTAL_TIME % 3600 / 60)) min"
echo "Output directory: $OUTPUT_DIR/"
echo ""
echo "Next steps:"
echo "  1. Merge count matrices: python3 merge_soloTE_samples.py"
echo "  2. Add metadata and run DE analysis"
echo "  3. Compare with scTE results"
echo ""
