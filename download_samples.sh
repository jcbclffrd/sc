#!/bin/bash

# Download script for PRJNA729525 snRNA-seq samples
# This will download multiple samples using prefetch and fasterq-dump

# List of SRR accessions for 10x snRNA-seq samples
# Starting with first 5 samples for initial analysis
SAMPLES=(
    SRR14514292
    SRR14514293
    SRR14514294
    SRR14514295
    SRR14514296
)

echo "================================================"
echo "Starting SRA download for PRJNA729525"
echo "Started: $(date)"
echo "Number of samples: ${#SAMPLES[@]}"
echo "================================================"

for srr in "${SAMPLES[@]}"; do
    echo ""
    echo "------------------------------------------------"
    echo "Downloading $srr..."
    echo "Started: $(date)"
    echo "------------------------------------------------"
    
    mkdir -p sra_downloads/$srr
    
    # Download using prefetch first (more reliable)
    echo "Step 1: Prefetch $srr..."
    prefetch $srr -O sra_downloads/$srr --max-size 50G
    
    # Convert to FASTQ
    echo "Step 2: Converting to FASTQ..."
    fasterq-dump --split-files --threads 8 \
        --outdir sra_downloads/$srr \
        sra_downloads/$srr/$srr/$srr.sra
    
    # Compress FASTQ files
    echo "Step 3: Compressing FASTQ files..."
    gzip sra_downloads/$srr/*.fastq
    
    # Clean up SRA file
    rm -rf sra_downloads/$srr/$srr
    
    echo "✓ Completed $srr at $(date)"
    ls -lh sra_downloads/$srr/
done

echo ""
echo "================================================"
echo "All downloads completed!"
echo "Finished: $(date)"
echo "================================================"
