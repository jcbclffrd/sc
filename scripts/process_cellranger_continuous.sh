#!/bin/bash

# Continuously process ATAC-seq samples with CellRanger ATAC as they become ready
# This script monitors for completed downloads and processes them automatically

# Note: Removed set -e to prevent silent failures
# set -e

# Add CellRanger to PATH
export PATH="$HOME/software/cellranger-atac-2.2.0:$PATH"

# Configuration
FASTQ_DIR="$HOME/sc/sra_downloads/ATAC-seq"
OUTPUT_DIR="$HOME/sc/cellranger_atac_output"
REF_DIR="$HOME/sc/cellranger_references/refdata-cellranger-arc-GRCh38-2024-A"
LOG_DIR="$HOME/sc/logs_atacseq"
MAIN_LOG="$LOG_DIR/continuous_processing.log"

# Resource settings
CORES=16
MEMORY=64
CHECK_INTERVAL=300  # Check every 5 minutes for new samples

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# All 32 ATAC-seq samples
ALL_SAMPLES=(
    SRR14514129 SRR14514130 SRR14514131 SRR14514132
    SRR14514133 SRR14514134 SRR14514135 SRR14514136
    SRR14514137 SRR14514138 SRR14514139 SRR14514140
    SRR14514141 SRR14514142 SRR14514143 SRR14514144
    SRR14514145 SRR14514146 SRR14514147 SRR14514148
    SRR14514149 SRR14514150 SRR14514151 SRR14514152
    SRR14514153 SRR14514154 SRR14514155 SRR14514156
    SRR14514157 SRR14514158 SRR14514159 SRR14514160
)

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$MAIN_LOG"
}

# Function to check if sample is ready (has all 4 fastq.gz files)
is_sample_ready() {
    local SAMPLE_ID=$1
    local COUNT=$(ls -1 "${FASTQ_DIR}/${SAMPLE_ID}_"*.fastq.gz 2>/dev/null | wc -l)
    [ $COUNT -eq 4 ]
}

# Function to check if sample is already processed
is_sample_processed() {
    local SAMPLE_ID=$1
    [ -f "$OUTPUT_DIR/${SAMPLE_ID}/outs/possorted_bam.bam" ]
}

# Function to check if sample is currently processing
is_sample_processing() {
    local SAMPLE_ID=$1
    pgrep -f "cellranger-atac.*${SAMPLE_ID}" > /dev/null
}

# Function to process one sample
process_sample() {
    local SAMPLE_ID=$1
    local STAGING_DIR="$OUTPUT_DIR/${SAMPLE_ID}_fastqs"
    
    log "[$SAMPLE_ID] Starting CellRanger ATAC processing..."
    
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
        log "[$SAMPLE_ID] ✓ Completed successfully! BAM size: $BAM_SIZE"
        return 0
    else
        log "[$SAMPLE_ID] ✗ Failed! Check log: $LOG_DIR/${SAMPLE_ID}_cellranger.log"
        return 1
    fi
}

# Main monitoring loop
log "=========================================="
log "CellRanger ATAC - Continuous Processing"
log "=========================================="
log "Monitoring ${#ALL_SAMPLES[@]} samples"
log "Check interval: ${CHECK_INTERVAL}s"
log "Cores: $CORES, Memory: ${MEMORY}GB"
log ""

TOTAL_PROCESSED=0
TOTAL_FAILED=0

while true; do
    # Count status
    READY_COUNT=0
    PROCESSED_COUNT=0
    PENDING_COUNT=0
    
    for SAMPLE_ID in "${ALL_SAMPLES[@]}"; do
        if is_sample_processed "$SAMPLE_ID"; then
            ((PROCESSED_COUNT++))
        elif is_sample_ready "$SAMPLE_ID"; then
            ((READY_COUNT++))
        else
            ((PENDING_COUNT++))
        fi
    done
    
    log "Status: $PROCESSED_COUNT processed, $READY_COUNT ready, $PENDING_COUNT pending (downloading or missing)"
    
    # Process ready samples (one at a time to avoid memory issues)
    PROCESSED_THIS_ROUND=0
    for SAMPLE_ID in "${ALL_SAMPLES[@]}"; do
        # Skip if already processed
        if is_sample_processed "$SAMPLE_ID"; then
            continue
        fi
        
        # Skip if currently processing
        if is_sample_processing "$SAMPLE_ID"; then
            log "[$SAMPLE_ID] Already processing, skipping..."
            continue
        fi
        
        # Process if ready
        if is_sample_ready "$SAMPLE_ID"; then
            if process_sample "$SAMPLE_ID"; then
                ((TOTAL_PROCESSED++))
                ((PROCESSED_THIS_ROUND++))
            else
                ((TOTAL_FAILED++))
            fi
            
            # Process one at a time (can adjust to parallel if needed)
            break
        fi
    done
    
    # Check if all samples are done
    if [ $PROCESSED_COUNT -eq ${#ALL_SAMPLES[@]} ]; then
        log "=========================================="
        log "All samples processed!"
        log "Total successful: $TOTAL_PROCESSED"
        log "Total failed: $TOTAL_FAILED"
        log "=========================================="
        break
    fi
    
    # If nothing was processed this round and nothing is currently processing, wait before checking again
    if [ $PROCESSED_THIS_ROUND -eq 0 ]; then
        RUNNING=$(ps aux | grep -c "cellranger-atac count" | grep -v grep || echo 0)
        if [ "$RUNNING" -eq 0 ]; then
            log "No samples ready to process. Waiting ${CHECK_INTERVAL}s before next check..."
            sleep $CHECK_INTERVAL
        else
            log "CellRanger jobs running. Waiting 60s before next check..."
            sleep 60
        fi
    fi
done

log "Script completed at $(date '+%Y-%m-%d %H:%M:%S')"
