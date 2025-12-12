#!/bin/bash
#
# Create CellRanger-compatible symlinks for FASTQ files
# Converts: SRR_2.fastq.gz / SRR_3.fastq.gz
# To:       SRR_S1_L001_R1_001.fastq.gz / SRR_S1_L001_R2_001.fastq.gz

set -e

SOURCE_DIR="/home/jacobc/sc/sra_downloads"
TARGET_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/fastq_cellranger_format"
SAMPLE_MAP="../sample_mapping.csv"

echo "======================================================================"
echo "PREPARING FASTQ FILES FOR CELLRANGER"
echo "======================================================================"
echo ""
echo "This script creates symlinks with CellRanger-compatible naming"
echo ""
echo "Source: $SOURCE_DIR"
echo "Target: $TARGET_DIR"
echo ""

# Create target directory
mkdir -p "$TARGET_DIR"

# Process each sample
N_PROCESSED=0
N_SKIPPED=0

while IFS=, read -r sample_id gsm srx srr; do
    if [ "$sample_id" = "sample_id" ]; then continue; fi  # Skip header
    
    SOURCE_SAMPLE_DIR="$SOURCE_DIR/$srr"
    TARGET_SAMPLE_DIR="$TARGET_DIR/$srr"
    
    # Check if source files exist
    if [ ! -f "$SOURCE_SAMPLE_DIR/${srr}_2.fastq.gz" ] || [ ! -f "$SOURCE_SAMPLE_DIR/${srr}_3.fastq.gz" ]; then
        echo "⚠️  $srr: Source files not found, skipping"
        N_SKIPPED=$((N_SKIPPED + 1))
        continue
    fi
    
    # Create sample directory
    mkdir -p "$TARGET_SAMPLE_DIR"
    
    # Create symlinks with CellRanger naming
    # _2 = R1 (barcode + UMI, 28bp)
    # _3 = R2 (cDNA sequence, ~100bp)
    ln -sf "$SOURCE_SAMPLE_DIR/${srr}_2.fastq.gz" "$TARGET_SAMPLE_DIR/${srr}_S1_L001_R1_001.fastq.gz"
    ln -sf "$SOURCE_SAMPLE_DIR/${srr}_3.fastq.gz" "$TARGET_SAMPLE_DIR/${srr}_S1_L001_R2_001.fastq.gz"
    
    echo "✓ $srr: Created symlinks"
    N_PROCESSED=$((N_PROCESSED + 1))
    
done < "$SAMPLE_MAP"

echo ""
echo "======================================================================"
echo "COMPLETE"
echo "======================================================================"
echo ""
echo "Processed: $N_PROCESSED samples"
echo "Skipped: $N_SKIPPED samples"
echo ""
echo "Output directory: $TARGET_DIR"
echo ""

if [ $N_PROCESSED -ge 18 ]; then
    echo "✓ All 18 samples ready for CellRanger"
    echo ""
    echo "Update run_docker.sh to use:"
    echo "  HOST_FASTQ_DIR=\"$TARGET_DIR\""
else
    echo "⚠️  Only $N_PROCESSED/18 samples available"
fi
echo ""
