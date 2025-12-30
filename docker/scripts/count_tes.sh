#!/bin/bash
# count_tes.sh - Quantify TEs with featureCounts
set -e

TE_GTF=${TE_GTF:-/annotations/hg38_rmsk_TE.gtf}
THREADS=${THREADS:-8}
OUTPUT_DIR=${OUTPUT_DIR:-/output}
DATA_DIR=${DATA_DIR:-/data}

if [[ ! -f "$TE_GTF" ]]; then
    echo "ERROR: TE GTF not found: $TE_GTF"
    exit 1
fi

echo "=========================================="
echo "TE Quantification"
echo "GTF: $TE_GTF"
echo "Threads: $THREADS"
echo "=========================================="

# Find all BAM files
BAM_FILES=($(find "$DATA_DIR" -name "Aligned.sortedByCoord.out.bam"))

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files found in $DATA_DIR"
    exit 1
fi

echo "Found ${#BAM_FILES[@]} BAM files"

mkdir -p "$OUTPUT_DIR"

# Run featureCounts
featureCounts \
    -a "$TE_GTF" \
    -o "$OUTPUT_DIR/atacseq_te_counts.txt" \
    -F GTF \
    -t exon \
    -g gene_id \
    -p \
    --countReadPairs \
    -M \
    --fraction \
    -T "$THREADS" \
    "${BAM_FILES[@]}"

# Create clean matrix
cut -f1,7- "$OUTPUT_DIR/atacseq_te_counts.txt" > "$OUTPUT_DIR/te_counts_matrix.tsv"

TE_COUNT=$(tail -n +3 "$OUTPUT_DIR/te_counts_matrix.tsv" | wc -l)
echo "✓ Complete: $TE_COUNT TEs × ${#BAM_FILES[@]} samples"
