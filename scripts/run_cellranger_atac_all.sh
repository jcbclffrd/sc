#!/bin/bash

# Process all 32 ATAC-seq samples with CellRanger ATAC
# Run on R-instance with sufficient memory

set -e

# Configuration
FASTQ_DIR="$HOME/sc/sra_downloads/ATAC-seq"
OUTPUT_DIR="$HOME/sc/cellranger_atac_output"
REF_DIR="$HOME/sc/cellranger_references/refdata-cellranger-arc-GRCh38-2024-A"
LOG_DIR="$HOME/sc/logs_atacseq"

# Resource settings (adjust based on instance)
CORES=16        # Cores per job
MEMORY=64       # GB per job
MAX_PARALLEL=1  # Number of samples to process in parallel (depends on RAM)

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# All 32 ATAC-seq samples
SAMPLES=(
    SRR14514129 SRR14514130 SRR14514131 SRR14514132
    SRR14514133 SRR14514134 SRR14514135 SRR14514136
    SRR14514137 SRR14514138 SRR14514139 SRR14514140
    SRR14514141 SRR14514142 SRR14514143 SRR14514144
    SRR14514145 SRR14514146 SRR14514147 SRR14514148
    SRR14514149 SRR14514150 SRR14514151 SRR14514152
    SRR14514153 SRR14514154 SRR14514155 SRR14514156
    SRR14514157 SRR14514158 SRR14514159 SRR14514160
)

echo "=========================================="
echo "CellRanger ATAC - Batch Processing"
echo "=========================================="
echo "Total samples: ${#SAMPLES[@]}"
echo "Cores per job: $CORES"
echo "Memory per job: ${MEMORY}GB"
echo "Parallel jobs: $MAX_PARALLEL"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Function to process one sample
process_sample() {
    local SAMPLE_ID=$1
    local STAGING_DIR="$OUTPUT_DIR/${SAMPLE_ID}_fastqs"
    
    echo "[$SAMPLE_ID] Started: $(date '+%H:%M:%S')"
    
    # Check if already processed
    if [ -f "$OUTPUT_DIR/${SAMPLE_ID}/outs/possorted_bam.bam" ]; then
        echo "[$SAMPLE_ID] ✓ Already processed, skipping"
        return 0
    fi
    
    # Check if FASTQs exist
    if [ ! -f "${FASTQ_DIR}/${SAMPLE_ID}_1.fastq.gz" ]; then
        echo "[$SAMPLE_ID] ✗ FASTQs not found, skipping"
        return 1
    fi
    
    # Create staging directory with proper naming
    # CellRanger ATAC expects: R1=50bp genomic, R2=16bp barcode, R3=50bp genomic, I1=8bp index
    mkdir -p "$STAGING_DIR"
    ln -sf "${FASTQ_DIR}/${SAMPLE_ID}_2.fastq.gz" "$STAGING_DIR/${SAMPLE_ID}_S1_L001_R1_001.fastq.gz"
    ln -sf "${FASTQ_DIR}/${SAMPLE_ID}_3.fastq.gz" "$STAGING_DIR/${SAMPLE_ID}_S1_L001_R2_001.fastq.gz"
    ln -sf "${FASTQ_DIR}/${SAMPLE_ID}_4.fastq.gz" "$STAGING_DIR/${SAMPLE_ID}_S1_L001_R3_001.fastq.gz"
    ln -sf "${FASTQ_DIR}/${SAMPLE_ID}_1.fastq.gz" "$STAGING_DIR/${SAMPLE_ID}_S1_L001_I1_001.fastq.gz"
    
    # Run CellRanger
    cd "$OUTPUT_DIR"
    if cellranger-atac count \
        --id="${SAMPLE_ID}" \
        --reference="$REF_DIR" \
        --fastqs="$STAGING_DIR" \
        --sample="$SAMPLE_ID" \
        --localcores=$CORES \
        --localmem=$MEMORY \
        >> "$LOG_DIR/${SAMPLE_ID}_cellranger.log" 2>&1; then
        
        BAM_SIZE=$(du -h "$OUTPUT_DIR/${SAMPLE_ID}/outs/possorted_bam.bam" 2>/dev/null | cut -f1 || echo "N/A")
        echo "[$SAMPLE_ID] ✓ Complete (BAM: $BAM_SIZE) - $(date '+%H:%M:%S')"
    else
        echo "[$SAMPLE_ID] ✗ Failed - see $LOG_DIR/${SAMPLE_ID}_cellranger.log"
        return 1
    fi
}

export -f process_sample
export FASTQ_DIR OUTPUT_DIR REF_DIR LOG_DIR CORES MEMORY

# Process all samples (in parallel if MAX_PARALLEL > 1)
printf '%s\n' "${SAMPLES[@]}" | xargs -P "$MAX_PARALLEL" -I {} bash -c 'process_sample "$@"' _ {}

echo ""
echo "=========================================="
echo "Batch Processing Complete"
echo "=========================================="
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Summary
TOTAL_SAMPLES=${#SAMPLES[@]}
COMPLETED=$(find "$OUTPUT_DIR" -name "possorted_bam.bam" -type f 2>/dev/null | wc -l)

echo "Summary:"
echo "  Total samples: $TOTAL_SAMPLES"
echo "  Completed: $COMPLETED"
echo "  Failed/Pending: $((TOTAL_SAMPLES - COMPLETED))"
echo ""

if [ $COMPLETED -eq $TOTAL_SAMPLES ]; then
    echo "✓ All samples processed successfully!"
else
    echo "⚠ Some samples incomplete. Check logs in $LOG_DIR"
    echo ""
    echo "Missing samples:"
    for sample in "${SAMPLES[@]}"; do
        if [ ! -f "$OUTPUT_DIR/$sample/outs/possorted_bam.bam" ]; then
            echo "  - $sample"
        fi
    done
fi
