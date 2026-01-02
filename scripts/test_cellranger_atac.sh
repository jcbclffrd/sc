#!/bin/bash

# Test CellRanger ATAC on a single sample
# This script processes one ATAC-seq sample with CellRanger ATAC

set -e

# Configuration
SAMPLE_ID="${1:-SRR14514129}"  # Default to first sample, or pass as argument
FASTQ_DIR="$HOME/sc/sra_downloads/ATAC-seq"
OUTPUT_DIR="$HOME/sc/cellranger_atac_output"
REF_DIR="$HOME/sc/cellranger_references/refdata-cellranger-arc-GRCh38-2024-A"
LOG_DIR="$HOME/sc/logs_atacseq"

# Resource settings (adjust based on instance)
CORES=16        # Adjust for your R-instance
MEMORY=64       # GB of RAM to use

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# Check if CellRanger ATAC is in PATH
if ! command -v cellranger-atac &> /dev/null; then
    echo "ERROR: cellranger-atac not found in PATH"
    echo "Run: export PATH=\$HOME/software/cellranger-atac-2.1.0:\$PATH"
    exit 1
fi

# Check if reference exists
if [ ! -d "$REF_DIR" ]; then
    echo "ERROR: Reference genome not found at $REF_DIR"
    echo "Run: ~/sc/scripts/setup_cellranger_atac.sh"
    exit 1
fi

# Check if FASTQs exist and are compressed
EXPECTED_FILES=(
    "${FASTQ_DIR}/${SAMPLE_ID}_1.fastq.gz"
    "${FASTQ_DIR}/${SAMPLE_ID}_2.fastq.gz"
    "${FASTQ_DIR}/${SAMPLE_ID}_3.fastq.gz"
    "${FASTQ_DIR}/${SAMPLE_ID}_4.fastq.gz"
)

echo "=========================================="
echo "CellRanger ATAC - Test Run"
echo "=========================================="
echo "Sample: $SAMPLE_ID"
echo "FASTQ directory: $FASTQ_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Reference: $REF_DIR"
echo "Cores: $CORES"
echo "Memory: ${MEMORY}GB"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Check for FASTQ files
echo "Checking FASTQ files..."
MISSING=0
for file in "${EXPECTED_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "  ✗ Missing: $file"
        MISSING=1
    else
        SIZE=$(du -h "$file" | cut -f1)
        echo "  ✓ Found: $(basename $file) ($SIZE)"
    fi
done

if [ $MISSING -eq 1 ]; then
    echo ""
    echo "ERROR: Missing FASTQ files. Make sure downloads are complete."
    exit 1
fi

echo ""
echo "All FASTQ files found!"
echo ""

# Important: CellRanger ATAC expects specific FASTQ naming
# The SRA FASTQs from fasterq-dump are organized as:
# _1.fastq.gz = I7 (sample index, 8bp)
# _2.fastq.gz = Read 1 (genomic, 50bp) 
# _3.fastq.gz = I5 (cell barcode, 16bp)
# _4.fastq.gz = Read 2 (genomic, 50bp)

# CellRanger ATAC expects: 
# R1 = genomic read (50bp)
# R2 = cell barcode (16bp) ← KEY: barcode goes in R2, not R1!
# R3 = genomic read (50bp)
# I1 = sample index (8bp)

# We'll create symlinks in proper format
STAGING_DIR="$OUTPUT_DIR/${SAMPLE_ID}_fastqs"
mkdir -p "$STAGING_DIR"

echo "Creating properly named symlinks for CellRanger..."
ln -sf "${FASTQ_DIR}/${SAMPLE_ID}_2.fastq.gz" "$STAGING_DIR/${SAMPLE_ID}_S1_L001_R1_001.fastq.gz"
ln -sf "${FASTQ_DIR}/${SAMPLE_ID}_3.fastq.gz" "$STAGING_DIR/${SAMPLE_ID}_S1_L001_R2_001.fastq.gz"
ln -sf "${FASTQ_DIR}/${SAMPLE_ID}_4.fastq.gz" "$STAGING_DIR/${SAMPLE_ID}_S1_L001_R3_001.fastq.gz"
ln -sf "${FASTQ_DIR}/${SAMPLE_ID}_1.fastq.gz" "$STAGING_DIR/${SAMPLE_ID}_S1_L001_I1_001.fastq.gz"

echo "  ✓ Symlinks created in $STAGING_DIR"
echo ""

# Change to output directory (CellRanger creates output here)
cd "$OUTPUT_DIR"

echo "=========================================="
echo "Running CellRanger ATAC Count"
echo "=========================================="
echo ""

# Run CellRanger ATAC
cellranger-atac count \
    --id="${SAMPLE_ID}" \
    --reference="$REF_DIR" \
    --fastqs="$STAGING_DIR" \
    --sample="$SAMPLE_ID" \
    --localcores=$CORES \
    --localmem=$MEMORY \
    2>&1 | tee "$LOG_DIR/${SAMPLE_ID}_cellranger.log"

echo ""
echo "=========================================="
echo "CellRanger ATAC Complete"
echo "=========================================="
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "Output directory: $OUTPUT_DIR/${SAMPLE_ID}"
echo ""
echo "Key outputs:"
echo "  - BAM: $OUTPUT_DIR/${SAMPLE_ID}/outs/possorted_bam.bam"
echo "  - Fragments: $OUTPUT_DIR/${SAMPLE_ID}/outs/fragments.tsv.gz"
echo "  - Peaks: $OUTPUT_DIR/${SAMPLE_ID}/outs/peaks.bed"
echo "  - Summary: $OUTPUT_DIR/${SAMPLE_ID}/outs/web_summary.html"
echo ""
echo "Check the web summary for QC metrics!"
echo ""

# Check if successful
if [ -f "$OUTPUT_DIR/${SAMPLE_ID}/outs/possorted_bam.bam" ]; then
    BAM_SIZE=$(du -h "$OUTPUT_DIR/${SAMPLE_ID}/outs/possorted_bam.bam" | cut -f1)
    echo "✓ SUCCESS: BAM file created ($BAM_SIZE)"
else
    echo "✗ FAILED: No BAM file found"
    exit 1
fi
