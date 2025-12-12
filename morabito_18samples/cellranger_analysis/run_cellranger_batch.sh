#!/bin/bash
#
# Run CellRanger count on 18 Morabito samples
# Processes FASTQ files from /app/fastq_data (mounted from host)

set -e

FASTQ_DIR="/app/fastq_data"
OUTPUT_DIR="/app/cellranger_output"
REFERENCE="${CELLRANGER_REFERENCE:-/app/refdata-gex-GRCh38-2024-A}"
SAMPLE_MAP="/app/sample_mapping/sample_mapping.csv"

# CellRanger parameters
THREADS=8
MEM_GB=64
CHEMISTRY="SC3Pv3"  # 10x Chromium Single Cell 3' v3

echo "======================================================================"
echo "CELLRANGER BATCH PROCESSING: 18 Samples"
echo "======================================================================"

# Check reference exists
if [ ! -d "$REFERENCE" ]; then
    echo "✗ Reference genome not found: $REFERENCE"
    exit 1
fi

# Count samples
N_SAMPLES=$(tail -n +2 "$SAMPLE_MAP" | wc -l)
echo "Processing $N_SAMPLES samples"
echo "Threads: $THREADS"
echo "Memory: ${MEM_GB}GB"
echo ""

# Process each sample
CURRENT=0
START_TIME=$(date +%s)

while IFS=, read -r sample_id gsm srx srr; do
    if [ "$sample_id" = "sample_id" ]; then continue; fi  # Skip header
    
    CURRENT=$((CURRENT + 1))
    
    echo "======================================================================"
    echo "[$CURRENT/$N_SAMPLES] Processing: $sample_id ($srr)"
    echo "======================================================================"
    
    SAMPLE_FASTQ_DIR="$FASTQ_DIR/$srr"
    SAMPLE_OUTPUT_DIR="$OUTPUT_DIR/$srr"
    
    # Check if already processed
    if [ -d "$SAMPLE_OUTPUT_DIR/outs" ]; then
        echo "✓ Already processed, skipping..."
        continue
    fi
    
    # Check FASTQ directory exists
    if [ ! -d "$SAMPLE_FASTQ_DIR" ]; then
        echo "✗ FASTQ directory not found: $SAMPLE_FASTQ_DIR"
        continue
    fi
    
    echo "Input FASTQ: $(du -sh $SAMPLE_FASTQ_DIR | cut -f1)"
    echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
    
    SAMPLE_START=$(date +%s)
    
    # Run CellRanger count
    cd $OUTPUT_DIR
    cellranger count \
        --id=$srr \
        --transcriptome=$REFERENCE \
        --fastqs=$SAMPLE_FASTQ_DIR \
        --sample=$srr \
        --chemistry=$CHEMISTRY \
        --localcores=$THREADS \
        --localmem=$MEM_GB \
        --nosecondary  # Skip t-SNE/clustering (we'll do that ourselves)
    
    SAMPLE_END=$(date +%s)
    SAMPLE_TIME=$((SAMPLE_END - SAMPLE_START))
    
    echo "✓ Completed in $(($SAMPLE_TIME / 60)) min $(($SAMPLE_TIME % 60)) sec"
    echo ""
    
    # Show output summary
    if [ -f "$SAMPLE_OUTPUT_DIR/outs/metrics_summary.csv" ]; then
        echo "Output summary:"
        cat "$SAMPLE_OUTPUT_DIR/outs/metrics_summary.csv"
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
echo "  1. Merge count matrices: python3 merge_cellranger_outputs.py"
echo "  2. Compare with STAR+scTE: python3 compare_with_star.py"
echo ""
