#!/bin/bash

# Download ONLY the 152 10x Chromium v3 snRNA-seq samples from PRJNA729525
# These have 208bp reads (28bp barcode/UMI + 91bp cDNA, paired-end)

SRR_LIST="snRNA_seq_10x_samples.txt"
LOG_DIR="logs_10x"
DOWNLOAD_DIR="sra_downloads"

mkdir -p $LOG_DIR
mkdir -p $DOWNLOAD_DIR

# Count total samples
TOTAL=$(wc -l < $SRR_LIST)
CURRENT=0

echo "================================================"
echo "PRJNA729525 10x snRNA-seq Download"
echo "10x Chromium Single Cell 3' v3 platform"
echo "Started: $(date)"
echo "Total samples: $TOTAL"
echo "================================================"

while read srr; do
    CURRENT=$((CURRENT + 1))
    
    echo ""
    echo "================================================"
    echo "[$CURRENT/$TOTAL] Processing $srr"
    echo "Started: $(date)"
    echo "================================================"
    
    mkdir -p $DOWNLOAD_DIR/$srr
    SAMPLE_LOG="$LOG_DIR/${srr}.log"
    
    # Check if already downloaded (look for any FASTQ.gz files)
    if ls $DOWNLOAD_DIR/$srr/*.fastq.gz 1> /dev/null 2>&1; then
        echo "✓ $srr already downloaded, skipping..."
        continue
    fi
    
    {
        echo "Download started: $(date)"
        
        # Step 1: Prefetch (download SRA file)
        echo "Step 1/3: Downloading SRA file..."
        if prefetch $srr -O $DOWNLOAD_DIR/$srr --max-size 50G; then
            echo "✓ Prefetch completed"
        else
            echo "ERROR: Prefetch failed for $srr"
            continue
        fi
        
        # Step 2: Convert to FASTQ
        echo "Step 2/3: Converting to FASTQ..."
        if fasterq-dump --split-files --threads 8 \
            --outdir $DOWNLOAD_DIR/$srr \
            $DOWNLOAD_DIR/$srr/$srr/$srr.sra; then
            echo "✓ FASTQ conversion completed"
        else
            echo "ERROR: FASTQ conversion failed for $srr"
            continue
        fi
        
        # Step 3: Compress
        echo "Step 3/3: Compressing FASTQ files..."
        if gzip $DOWNLOAD_DIR/$srr/*.fastq; then
            echo "✓ Compression completed"
        else
            echo "ERROR: Compression failed for $srr"
            continue
        fi
        
        # Clean up SRA file
        rm -rf $DOWNLOAD_DIR/$srr/$srr
        
        echo "✓ COMPLETED: $srr at $(date)"
        ls -lh $DOWNLOAD_DIR/$srr/
        
    } > $SAMPLE_LOG 2>&1
    
    # Show summary
    if ls $DOWNLOAD_DIR/$srr/*.fastq.gz 1> /dev/null 2>&1; then
        SIZE=$(du -sh $DOWNLOAD_DIR/$srr | cut -f1)
        echo "✓ Success - Size: $SIZE"
    else
        echo "✗ Failed - Check $SAMPLE_LOG"
    fi
    
done < $SRR_LIST

echo ""
echo "================================================"
echo "Download Complete!"
echo "Finished: $(date)"
echo "================================================"
echo ""
echo "Summary:"
SUCCESS=$(find $DOWNLOAD_DIR -name "*.fastq.gz" | wc -l)
echo "  Successfully downloaded FASTQ files: $SUCCESS"
echo "  Total disk usage: $(du -sh $DOWNLOAD_DIR | cut -f1)"
echo ""
echo "Check individual logs in: $LOG_DIR/"
