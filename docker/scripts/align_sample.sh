#!/bin/bash
# align_sample.sh - Download from SRA and align with STAR
set -e

SAMPLE_ID=${SAMPLE_ID:-$1}
THREADS=${THREADS:-4}
GENOME_DIR=${GENOME_DIR:-/genome}
OUTPUT_DIR=${OUTPUT_DIR:-/output}
FASTQ_DIR=${FASTQ_DIR:-/data}

if [[ -z "$SAMPLE_ID" ]]; then
    echo "ERROR: SAMPLE_ID not set"
    echo "Usage: docker run -e SAMPLE_ID=SRR14514129 atacseq-alignment"
    exit 1
fi

echo "=========================================="
echo "Aligning: $SAMPLE_ID"
echo "Threads: $THREADS"
echo "Genome: $GENOME_DIR"
echo "=========================================="

cd "$FASTQ_DIR"

# Download from SRA
echo "Downloading from SRA..."
fasterq-dump "$SAMPLE_ID" --split-files --threads 2 --mem 4GB

# Compress
pigz -p 2 "${SAMPLE_ID}_2.fastq" &
pigz -p 2 "${SAMPLE_ID}_4.fastq" &
wait

READ1="${SAMPLE_ID}_2.fastq.gz"
READ2="${SAMPLE_ID}_4.fastq.gz"

# Create output directory
SAMPLE_OUT="$OUTPUT_DIR/$SAMPLE_ID"
mkdir -p "$SAMPLE_OUT"

# Align with STAR
echo "Aligning with STAR..."
STAR \
    --runThreadN "$THREADS" \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$READ1" "$READ2" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$SAMPLE_OUT/" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped None \
    --outSAMattributes NH HI AS nM NM MD \
    --outFilterMultimapNmax 100 \
    --outFilterMismatchNmax 3 \
    --alignEndsType EndToEnd \
    --alignIntronMax 1 \
    --alignEndsProtrude 0 ConcordantPair \
    --winAnchorMultimapNmax 100 \
    --outMultimapperOrder Random \
    --outSAMmultNmax 1

# Index BAM
samtools index "$SAMPLE_OUT/Aligned.sortedByCoord.out.bam"

# Cleanup
rm -rf "$SAMPLE_OUT/_STARtmp"
rm -f "$FASTQ_DIR/${SAMPLE_ID}"*.fastq*

echo "✓ Complete: $SAMPLE_OUT/Aligned.sortedByCoord.out.bam"
