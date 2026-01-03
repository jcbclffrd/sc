#!/bin/bash

# Download ATAC-seq samples from SRA using prefetch + fasterq-dump
# This is more reliable than direct fasterq-dump streaming

set -e

# Configuration
MAX_JOBS=6  # Parallel downloads (6 samples at once)
THREADS=6   # Threads per fasterq-dump
SRA_DIR="$HOME/sc/sra_downloads/ATAC-seq"
NCBI_DIR="$HOME/ncbi/public/sra"  # Default SRA cache location
LOG_DIR="$HOME/sc/logs_atacseq"

# Add tools to PATH
export PATH=$HOME/sratoolkit.3.3.0-ubuntu64/bin:$HOME/.local/bin:$PATH

mkdir -p "$SRA_DIR" "$LOG_DIR" "$NCBI_DIR"

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
echo "Downloading ATAC-seq Samples (Prefetch Method)"
echo "=========================================="
echo "Samples to download: ${#SAMPLES[@]}"
echo "Parallel jobs: $MAX_JOBS"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Function to download and process one sample
process_sample() {
    local SAMPLE=$1
    local LOG_FILE="$LOG_DIR/${SAMPLE}_prefetch.log"
    
    echo "[$SAMPLE] Started: $(date '+%H:%M:%S')"
    
    # Check if already complete
    if [ -f "$SRA_DIR/${SAMPLE}_1.fastq.gz" ] && \
       [ -f "$SRA_DIR/${SAMPLE}_2.fastq.gz" ] && \
       [ -f "$SRA_DIR/${SAMPLE}_3.fastq.gz" ] && \
       [ -f "$SRA_DIR/${SAMPLE}_4.fastq.gz" ]; then
        echo "[$SAMPLE] ✓ Already complete, skipping"
        return 0
    fi
    
    # Step 1: Prefetch - download complete SRA file
    echo "[$SAMPLE] Prefetching SRA file..." >> "$LOG_FILE"
    if ! prefetch "$SAMPLE" --max-size 200G --output-directory "$NCBI_DIR" >> "$LOG_FILE" 2>&1; then
        echo "[$SAMPLE] ✗ Prefetch failed"
        return 1
    fi
    
    # Step 2: Extract FASTQs from local SRA file
    echo "[$SAMPLE] Extracting FASTQs..." >> "$LOG_FILE"
    cd "$SRA_DIR"
    if ! fasterq-dump "$NCBI_DIR/${SAMPLE}/${SAMPLE}.sra" \
        --split-files \
        --include-technical \
        --threads "$THREADS" \
        --mem 8GB \
        --temp "$SRA_DIR" \
        --outdir "$SRA_DIR" >> "$LOG_FILE" 2>&1; then
        echo "[$SAMPLE] ✗ fasterq-dump failed"
        return 1
    fi
    
    # Step 3: Compress FASTQs
    echo "[$SAMPLE] Compressing..." >> "$LOG_FILE"
    if pigz -p "$THREADS" "${SAMPLE}_1.fastq" "${SAMPLE}_2.fastq" "${SAMPLE}_3.fastq" "${SAMPLE}_4.fastq" 2>> "$LOG_FILE"; then
        echo "[$SAMPLE] ✓ Downloaded and compressed: $(date '+%H:%M:%S')"
    else
        echo "[$SAMPLE] ✗ Compression failed"
        return 1
    fi
    
    # Step 4: Clean up SRA file to save space
    echo "[$SAMPLE] Cleaning up SRA file..." >> "$LOG_FILE"
    rm -rf "$NCBI_DIR/${SAMPLE}" 2>> "$LOG_FILE" || true
    
    return 0
}

export -f process_sample
export SRA_DIR NCBI_DIR LOG_DIR THREADS PATH

# Process samples in parallel
printf '%s\n' "${SAMPLES[@]}" | xargs -P "$MAX_JOBS" -I {} bash -c 'process_sample "$@"' _ {}

echo ""
echo "=========================================="
echo "Download completed: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================="

# Summary
COMPRESSED=$(ls "$SRA_DIR"/*.fastq.gz 2>/dev/null | wc -l)
COMPLETE=$((COMPRESSED / 4))
echo ""
echo "Complete samples: $COMPLETE / 32"
echo "Compressed files: $COMPRESSED / 128"
echo ""
echo "Check logs in: $LOG_DIR"
