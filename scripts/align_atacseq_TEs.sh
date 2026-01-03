#!/bin/bash

# ATAC-seq alignment pipeline for TE analysis
# Uses STAR with multi-mapping (n=100) to capture reads mapping to repetitive TEs

set -e

# Configuration
STAR_INDEX="/home/jacobc/sc/star_index"
FASTQ_DIR="/home/jacobc/sc/sra_downloads/ATAC-seq"
OUTPUT_DIR="/home/jacobc/sc/atacseq_aligned"
THREADS=8

# STAR parameters for TE analysis
MULTIMAP_MAX=100  # Allow up to 100 mapping locations (critical for TEs)

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"

echo "=========================================="
echo "ATAC-seq TE Alignment Pipeline"
echo "=========================================="
echo "STAR index: $STAR_INDEX"
echo "Input FASTQs: $FASTQ_DIR"
echo "Output: $OUTPUT_DIR"
echo "Multi-mapping max: $MULTIMAP_MAX"
echo "Threads: $THREADS"
echo "=========================================="
echo ""

# Find all sample pairs
cd "$FASTQ_DIR"
for R1 in *_1.fastq.gz *_1.fastq; do
    # Skip if no files found
    [[ ! -f "$R1" ]] && continue
    
    # Get sample name and R2 file
    SAMPLE=$(basename "$R1" | sed 's/_1\.fastq.*//')
    R2="${R1/_1/_2}"
    
    # Skip if R2 doesn't exist
    [[ ! -f "$R2" ]] && continue
    
    # Skip if already aligned
    if [[ -f "$OUTPUT_DIR/${SAMPLE}/Aligned.sortedByCoord.out.bam" ]]; then
        echo "✓ Skipping $SAMPLE (already aligned)"
        continue
    fi
    
    echo "=========================================="
    echo "Aligning: $SAMPLE"
    echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "R1: $R1"
    echo "R2: $R2"
    echo "=========================================="
    
    # Create sample output directory
    mkdir -p "$OUTPUT_DIR/${SAMPLE}"
    
    # Decompress if needed
    if [[ "$R1" == *.gz ]]; then
        READ1="<(gunzip -c $FASTQ_DIR/$R1)"
        READ2="<(gunzip -c $FASTQ_DIR/$R2)"
    else
        READ1="$FASTQ_DIR/$R1"
        READ2="$FASTQ_DIR/$R2"
    fi
    
    # Run STAR alignment
    STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$FASTQ_DIR/$R1" "$FASTQ_DIR/$R2" \
        --readFilesCommand $(if [[ "$R1" == *.gz ]]; then echo "zcat"; else echo "cat"; fi) \
        --outFileNamePrefix "$OUTPUT_DIR/${SAMPLE}/" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM NM \
        --outFilterMultimapNmax "$MULTIMAP_MAX" \
        --winAnchorMultimapNmax "$MULTIMAP_MAX" \
        --outMultimapperOrder Random \
        --outSAMmultNmax 1 \
        --alignEndsType EndToEnd \
        --alignIntronMax 1 \
        --alignMatesGapMax 2000 \
        2>&1 | tee "$OUTPUT_DIR/logs/${SAMPLE}_alignment.log"
    
    # Index BAM
    echo "Indexing BAM..."
    samtools index "$OUTPUT_DIR/${SAMPLE}/Aligned.sortedByCoord.out.bam"
    
    # Get alignment stats
    echo "Alignment stats:"
    samtools flagstat "$OUTPUT_DIR/${SAMPLE}/Aligned.sortedByCoord.out.bam" | tee "$OUTPUT_DIR/logs/${SAMPLE}_flagstat.txt"
    
    echo "✓ $SAMPLE complete at $(date '+%Y-%m-%d %H:%M:%S')"
    echo ""
done

echo "=========================================="
echo "All alignments complete!"
echo "=========================================="
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Next steps:"
echo "1. Count TE reads with featureCounts or TEtranscripts"
echo "2. Differential accessibility analysis with DESeq2"
echo ""
echo "Summary:"
ls -d "$OUTPUT_DIR"/SRR* 2>/dev/null | wc -l | xargs echo "Samples aligned:"
du -sh "$OUTPUT_DIR"
