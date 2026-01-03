#!/bin/bash

# Custom ATAC-seq pipeline with multimapper support for TE quantification
# 
# Strategy:
# 1. Use CellRanger ATAC output to identify valid cell barcodes (QC done!)
# 2. Extract reads from CellRanger's BAM (already demultiplexed & QC'd)
# 3. Re-align with BWA-MEM allowing multimappers (-a flag)
# 4. Preserve cell barcode tags (CB:Z:)
# 5. Use for TE quantification with scTE
#
# This hybrid approach gives us:
# - CellRanger's excellent cell calling and QC
# - Multimapper support for TE analysis

set -e

# Configuration
SAMPLE_ID="${1:-SRR14514130}"
FASTQ_DIR="$HOME/sc/sra_downloads/ATAC-seq"
OUTPUT_DIR="$HOME/sc/atacseq_multimapper_output"     # Multimapper BAMs go here (SEPARATE from CellRanger!)
GENOME_DIR="$HOME/sc/genome"
LOG_DIR="$HOME/sc/logs_atacseq"
CELLRANGER_OUT="$HOME/sc/cellranger_atac_output"     # Read CellRanger QC data from here (input only)

# BWA parameters
THREADS=16
BWA_INDEX="$GENOME_DIR/GRCh38.primary_assembly.genome.fa"

# Create directories
mkdir -p "$OUTPUT_DIR/$SAMPLE_ID" "$LOG_DIR"

echo "=========================================="
echo "ATAC-seq Multimapper Pipeline for TE Analysis"
echo "=========================================="
echo "Sample: $SAMPLE_ID"
echo "Strategy: Use CellRanger QC + BWA multimapper alignment"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Step 1: Verify CellRanger output exists
echo "[Step 1] Checking CellRanger output..."
CELLRANGER_BAM="$CELLRANGER_OUT/$SAMPLE_ID/outs/possorted_bam.bam"
CELLRANGER_CSV="$CELLRANGER_OUT/$SAMPLE_ID/outs/singlecell.csv"

if [ ! -f "$CELLRANGER_BAM" ]; then
    echo "  ERROR: CellRanger BAM not found: $CELLRANGER_BAM"
    echo "  Please run CellRanger ATAC first to get cell QC data"
    exit 1
fi

if [ ! -f "$CELLRANGER_CSV" ]; then
    echo "  ERROR: CellRanger singlecell.csv not found"
    echo "  This file contains cell barcode QC information"
    exit 1
fi

echo "  ✓ CellRanger output found"
echo "  BAM: $CELLRANGER_BAM"
echo "  CSV: $CELLRANGER_CSV"

# Step 2: Extract valid cell barcodes from CellRanger
echo ""
echo "[Step 2] Extracting valid cell barcodes from CellRanger QC..."

# Column 3 in singlecell.csv: is__cell_barcode (1 = cell, 0 = background)
awk -F',' 'NR>1 && $3 == 1 {print $1}' "$CELLRANGER_CSV" \
    > "$OUTPUT_DIR/$SAMPLE_ID/valid_barcodes.txt"

VALID_BC_COUNT=$(wc -l < "$OUTPUT_DIR/$SAMPLE_ID/valid_barcodes.txt")
TOTAL_BC_COUNT=$(awk 'NR>1' "$CELLRANGER_CSV" | wc -l)

echo "  ✓ Extracted barcodes"
echo "  Valid cells: $VALID_BC_COUNT / $TOTAL_BC_COUNT barcodes"
echo "  Output: $OUTPUT_DIR/$SAMPLE_ID/valid_barcodes.txt"

# Step 3: Extract reads from CellRanger BAM
# CellRanger has already:
# - Demultiplexed by barcode
# - Added CB:Z: tags with cell barcodes
# - Filtered to valid cells
# We extract these reads and re-align with multimapper support
echo ""
echo "[Step 3] Extracting reads from CellRanger BAM..."
echo "  CellRanger already did the hard work of demultiplexing!"
echo "  We'll extract and preserve cell barcode tags (CB:Z:)"

# Extract to FASTQ, preserving barcode in read name using CB tag
# Using samtools fastq with special handling for CB tag
samtools view -h "$CELLRANGER_BAM" | \
    awk '
    /^@/ {print; next}  # Pass through header
    {
        # Extract CB tag
        cb = ""
        for (i=12; i<=NF; i++) {
            if ($i ~ /^CB:Z:/) {
                cb = substr($i, 6)
                break
            }
        }
        # Add CB to read name
        if (cb != "") {
            sub(/@/, "@", $1)
            $1 = $1 "_" cb
        }
        print
    }
    ' | \
    samtools fastq \
        -1 "$OUTPUT_DIR/$SAMPLE_ID/reads_R1.fastq.gz" \
        -2 "$OUTPUT_DIR/$SAMPLE_ID/reads_R2.fastq.gz" \
        -0 /dev/null -s /dev/null -n -

R1_COUNT=$(zcat "$OUTPUT_DIR/$SAMPLE_ID/reads_R1.fastq.gz" | wc -l)
R1_READS=$((R1_COUNT / 4))

echo "  ✓ Extracted $R1_READS read pairs"
echo "  R1: $OUTPUT_DIR/$SAMPLE_ID/reads_R1.fastq.gz"
echo "  R2: $OUTPUT_DIR/$SAMPLE_ID/reads_R2.fastq.gz"
echo "  Note: Cell barcodes preserved in read names"

# Step 4: Align with BWA-MEM allowing multimappers
echo ""
echo "[Step 4] Aligning with BWA-MEM (multimapper mode)..."
echo "  Key parameter: -a (report ALL alignments for multimappers)"

# Check BWA index
if [ ! -f "$BWA_INDEX" ]; then
    echo "  ERROR: BWA index not found: $BWA_INDEX"
    echo "  Creating BWA index... (this takes ~1 hour)"
    bwa index "$BWA_INDEX"
fi

# BWA-MEM with multimapper support
# -a : Output all alignments (CRITICAL for multimappers!)
# -T 30 : Minimum score (lower = more multimappers kept)
# -t 16 : Threads
echo "  Running BWA-MEM..."
echo "  Command: bwa mem -a -T 30 -t $THREADS"

bwa mem -a -T 30 -t $THREADS \
    "$BWA_INDEX" \
    "$OUTPUT_DIR/$SAMPLE_ID/reads_R1.fastq.gz" \
    "$OUTPUT_DIR/$SAMPLE_ID/reads_R2.fastq.gz" \
    2> "$LOG_DIR/${SAMPLE_ID}_bwa_multimapper.log" \
    | samtools view -bS - \
    | samtools sort -@ 8 -o "$OUTPUT_DIR/$SAMPLE_ID/aligned_multimapper_unsorted.bam" -

echo "  ✓ Alignment complete"
echo "  Output: $OUTPUT_DIR/$SAMPLE_ID/aligned_multimapper_unsorted.bam"

# Step 5: Re-add cell barcode tags (CB:Z:) to BAM
echo ""
echo "[Step 5] Re-adding cell barcode (CB:Z:) tags to BAM..."
echo "  Parsing barcodes from read names and adding as SAM tags"

# Python script to add CB tags back (inline for simplicity)
python3 << 'EOF' "$OUTPUT_DIR/$SAMPLE_ID/aligned_multimapper_unsorted.bam" \
    "$OUTPUT_DIR/$SAMPLE_ID/aligned_multimapper.bam"
import pysam
import sys

input_bam = sys.argv[1]
output_bam = sys.argv[2]

with pysam.AlignmentFile(input_bam, "rb") as infile:
    with pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:
        for read in infile:
            # Extract barcode from read name (format: READID_BARCODE)
            if '_' in read.query_name:
                parts = read.query_name.rsplit('_', 1)
                if len(parts) == 2:
                    barcode = parts[1]
                    # Add CB tag
                    read.set_tag('CB', barcode, value_type='Z')
            outfile.write(read)
EOF

# Index the final BAM
samtools index "$OUTPUT_DIR/$SAMPLE_ID/aligned_multimapper.bam"

echo "  ✓ Cell barcodes re-added as CB:Z: tags"
echo "  Output: $OUTPUT_DIR/$SAMPLE_ID/aligned_multimapper.bam"

# Clean up intermediate files
rm "$OUTPUT_DIR/$SAMPLE_ID/aligned_multimapper_unsorted.bam"
rm "$OUTPUT_DIR/$SAMPLE_ID/reads_R1.fastq.gz" "$OUTPUT_DIR/$SAMPLE_ID/reads_R2.fastq.gz"

# Step 6: Filter for valid cell barcodes
echo ""
echo "[Step 6] Filtering BAM for valid cell barcodes..."

samtools view -h "$OUTPUT_DIR/$SAMPLE_ID/aligned_multimapper.bam" | \
    awk -v bc_file="$OUTPUT_DIR/$SAMPLE_ID/valid_barcodes.txt" '
    BEGIN {
        while (getline < bc_file) valid[$1] = 1
    }
    /^@/ {print; next}  # Print header lines
    {
        # Extract CB tag
        for (i=12; i<=NF; i++) {
            if ($i ~ /^CB:Z:/) {
                split($i, tag, ":")
                barcode = tag[3]
                if (barcode in valid) {
                    print
                    break
                }
            }
        }
    }' | \
    samtools view -b -o "$OUTPUT_DIR/$SAMPLE_ID/aligned_multimapper_filtered.bam"

samtools index "$OUTPUT_DIR/$SAMPLE_ID/aligned_multimapper_filtered.bam"

echo "  ✓ BAM filtered for valid cell barcodes"
echo "  Output: $OUTPUT_DIR/$SAMPLE_ID/aligned_multimapper_filtered.bam"

# Step 6: Statistics
echo ""
echo "=========================================="
echo "Pipeline Statistics"
echo "=========================================="
echo ""
echo "Multimapper BAM (all valid cells):"
samtools flagstat "$OUTPUT_DIR/$SAMPLE_ID/aligned_multimapper_filtered.bam" \
    | tee "$LOG_DIR/${SAMPLE_ID}_multimapper_stats.txt"

# Compare with CellRanger BAM
if [ -f "$CELLRANGER_BAM" ]; then
    echo ""
    echo "CellRanger BAM (filtered, no multimappers):"
    samtools flagstat "$CELLRANGER_BAM" \
        | tee "$LOG_DIR/${SAMPLE_ID}_cellranger_stats.txt"
fi

echo ""
echo "=========================================="
echo "Multimapper alignment completed: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - Multimapper BAM: $OUTPUT_DIR/$SAMPLE_ID/aligned_multimapper_filtered.bam"
echo "  - Valid barcodes: $OUTPUT_DIR/$SAMPLE_ID/valid_barcodes.txt"
echo "  - Log files: $LOG_DIR/${SAMPLE_ID}_*.log"
echo ""
echo "Next steps for TE quantification:"
echo "  1. Install scTE: conda install -c bioconda scte"
echo "  2. Run scTE on filtered BAM:"
echo "     scTE -i $OUTPUT_DIR/$SAMPLE_ID/aligned_multimapper_filtered.bam \\"
echo "          -o $OUTPUT_DIR/$SAMPLE_ID/scTE_output \\"
echo "          -g hg38 -p $THREADS"
echo "  3. Compare TE accessibility between CellRanger and multimapper results"
echo "  4. Analyze cell-type specific TE patterns"
