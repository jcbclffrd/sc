#!/bin/bash

# Download missing ATAC-seq samples from PRJNA729525
# Total: 23 missing samples out of 32 total ATAC-seq samples

set -e

# Configuration
OUTPUT_DIR="sra_downloads"
LOG_DIR="logs_atacseq"
THREADS=8

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Missing ATAC-seq sample list (23 samples)
SAMPLES=(
    # 129-139 range (missing 140-143 which are already downloaded)
    SRR14514129
    SRR14514130
    SRR14514131
    SRR14514132
    SRR14514133
    SRR14514134
    SRR14514135
    SRR14514136
    SRR14514137
    SRR14514138
    SRR14514139
    # 144-145 (failed/empty downloads)
    SRR14514144
    SRR14514145
    # 151-160 range
    SRR14514151
    SRR14514152
    SRR14514153
    SRR14514154
    SRR14514155
    SRR14514156
    SRR14514157
    SRR14514158
    SRR14514159
    SRR14514160
)

echo "=========================================="
echo "ATAC-seq Sample Download Script"
echo "=========================================="
echo "Total samples to download: ${#SAMPLES[@]}"
echo "Output directory: $OUTPUT_DIR"
echo "Log directory: $LOG_DIR"
echo "Threads per download: $THREADS"
echo "Started: $(date)"
echo "=========================================="
echo

# Function to download a single sample
download_sample() {
    local sample=$1
    local sample_dir="$OUTPUT_DIR/$sample"
    local log_file="$LOG_DIR/${sample}.log"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting download: $sample"
    
    # Create sample directory
    mkdir -p "$sample_dir"
    
    # Download using fasterq-dump (converts SRA to FASTQ)
    # --split-files: split paired-end reads into separate files
    # --threads: number of threads
    # --outdir: output directory
    # --progress: show progress
    {
        echo "Starting fasterq-dump for $sample at $(date)"
        fasterq-dump "$sample" \
            --split-files \
            --threads "$THREADS" \
            --outdir "$sample_dir" \
            --progress \
            --verbose
        
        # Compress FASTQ files
        echo "Compressing FASTQ files for $sample at $(date)"
        cd "$sample_dir"
        for f in *.fastq; do
            if [ -f "$f" ]; then
                echo "Compressing $f..."
                gzip -v "$f"
            fi
        done
        cd - > /dev/null
        
        echo "Completed $sample at $(date)"
        
        # Show file sizes
        ls -lh "$sample_dir/"
        
    } > "$log_file" 2>&1
    
    if [ $? -eq 0 ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✓ Completed: $sample"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✗ Failed: $sample (check $log_file)"
    fi
    
    echo
}

# Main download loop
TOTAL=${#SAMPLES[@]}
CURRENT=0

for sample in "${SAMPLES[@]}"; do
    CURRENT=$((CURRENT + 1))
    echo "=========================================="
    echo "Progress: $CURRENT / $TOTAL"
    echo "=========================================="
    
    # Check if sample already exists and has data
    if [ -d "$OUTPUT_DIR/$sample" ]; then
        FASTQ_COUNT=$(find "$OUTPUT_DIR/$sample" -name "*.fastq.gz" -size +1M 2>/dev/null | wc -l)
        if [ "$FASTQ_COUNT" -ge 2 ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ⊘ Skipping: $sample (already downloaded with $FASTQ_COUNT FASTQ files)"
            echo
            continue
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] → Re-downloading: $sample (incomplete download)"
        fi
    fi
    
    download_sample "$sample"
    
    # Brief pause between downloads to avoid overwhelming the server
    sleep 2
done

echo "=========================================="
echo "Download Complete!"
echo "=========================================="
echo "Finished: $(date)"
echo

# Summary
echo "Summary of downloaded files:"
for sample in "${SAMPLES[@]}"; do
    if [ -d "$OUTPUT_DIR/$sample" ]; then
        COUNT=$(find "$OUTPUT_DIR/$sample" -name "*.fastq.gz" 2>/dev/null | wc -l)
        SIZE=$(du -sh "$OUTPUT_DIR/$sample" 2>/dev/null | cut -f1)
        if [ "$COUNT" -gt 0 ]; then
            echo "  $sample: $COUNT files ($SIZE)"
        else
            echo "  $sample: FAILED or EMPTY"
        fi
    else
        echo "  $sample: NOT DOWNLOADED"
    fi
done

echo
echo "Check individual logs in: $LOG_DIR/"
echo "=========================================="
