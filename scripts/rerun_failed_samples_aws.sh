#!/bin/bash

# Rerun failed ATAC-seq samples on AWS with lower parallelization (3 jobs)
# This avoids OOM issues that occurred with 9 parallel jobs

set -e

# Configuration
MAX_JOBS=3  # Reduced from 9 to avoid OOM
THREADS=4   # 3 jobs × 4 threads = 12 cores used
GENOME_DIR="$HOME/sc/star_index"
TE_GTF="$HOME/sc/annotations/hg38_rmsk_TE.gtf"
SRA_DIR="$HOME/sc/sra_downloads/ATAC-seq"
ALIGNED_DIR="$HOME/sc/atacseq_aligned"
LOG_DIR="$HOME/sc/logs_atacseq"

# Add tools to PATH
export PATH=$HOME/sratoolkit.3.3.0-ubuntu64/bin:$HOME/.local/bin:$PATH

mkdir -p "$SRA_DIR" "$ALIGNED_DIR" "$LOG_DIR"

# Failed samples to rerun (17 samples)
FAILED_SAMPLES=(
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
    SRR14514140
    SRR14514141
    SRR14514142
    SRR14514143
    SRR14514144
    SRR14514149
)

echo "=========================================="
echo "Rerunning Failed ATAC-seq Samples"
echo "=========================================="
echo "Instance: c5.9xlarge (36 cores, 72GB RAM)"
echo "Failed samples: ${#FAILED_SAMPLES[@]}"
echo "Parallel jobs: $MAX_JOBS"
echo "Threads per job: $THREADS"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Function to align one sample
align_sample() {
    local SAMPLE=$1
    local SAMPLE_DIR="$SRA_DIR"
    local SAMPLE_OUT="$ALIGNED_DIR/$SAMPLE"
    local LOG_FILE="$LOG_DIR/${SAMPLE}_alignment.log"
    
    echo "[$SAMPLE] Started: $(date '+%H:%M:%S')"
    
    mkdir -p "$SAMPLE_OUT"
    cd "$SAMPLE_DIR"
    
    # Download from SRA (include technical reads to get all 4 files: _1=read1, _2=barcode, _3=read2, _4=UMI)
    if ! fasterq-dump "$SAMPLE" --split-files --include-technical --threads 2 --mem 4GB >> "$LOG_FILE" 2>&1; then
        echo "[$SAMPLE] ✗ Download failed"
        return 1
    fi
    
    # Compress all 4 FASTQs immediately to save space
    pigz -p 4 "${SAMPLE}_1.fastq" "${SAMPLE}_2.fastq" "${SAMPLE}_3.fastq" "${SAMPLE}_4.fastq" 2>> "$LOG_FILE" || true
    
    # Technical reads (_1=barcodes length=8, _3=UMIs length=16)
    BARCODE="${SAMPLE}_1.fastq.gz"
    UMI="${SAMPLE}_3.fastq.gz"
    
    # Biological reads for alignment (_2 and _4, both length=50)
    READ2="${SAMPLE}_2.fastq.gz"
    READ4="${SAMPLE}_4.fastq.gz"
    
    if [ ! -f "$READ2" ] || [ ! -f "$READ4" ]; then
        echo "[$SAMPLE] ✗ FASTQ files not found"
        return 1
    fi
    
    # STAR alignment with multi-mapping (n=100) for TE analysis
    # Use biological reads (_2 and _4, both length=50); technical reads (_1=barcodes, _3=UMIs) kept for demux/dedup
    READFILES_CMD="gunzip -c"
    if ! STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$READ2" "$READ4" \
        --readFilesCommand $READFILES_CMD \
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
        --outSAMmultNmax 1 \
        --limitBAMsortRAM 8000000000 \
        >> "$LOG_FILE" 2>&1; then
        echo "[$SAMPLE] ✗ Alignment failed"
        return 1
    fi
    
    # Index BAM
    if ! samtools index "$SAMPLE_OUT/Aligned.sortedByCoord.out.bam" >> "$LOG_FILE" 2>&1; then
        echo "[$SAMPLE] ✗ Indexing failed"
        return 1
    fi
    
    # Clean up temp files and FASTQs to save space
    rm -rf "$SAMPLE_OUT/_STARtmp"
    rm -f "$READ1" "$READ2"
    
    BAM_SIZE=$(du -sh "$SAMPLE_OUT/Aligned.sortedByCoord.out.bam" | cut -f1)
    echo "[$SAMPLE] ✓ Aligned (BAM: $BAM_SIZE)"
}

export -f align_sample
export GENOME_DIR TE_GTF SRA_DIR ALIGNED_DIR LOG_DIR THREADS PATH

# Run alignment in parallel (3 jobs at a time)
echo "=========================================="
echo "Running alignments (3 parallel jobs)..."
echo "=========================================="

printf '%s\n' "${FAILED_SAMPLES[@]}" | xargs -P "$MAX_JOBS" -I {} bash -c 'align_sample "$@"' _ {}

echo ""
echo "=========================================="
echo "Alignment Complete"
echo "=========================================="
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "Checking results..."
TOTAL_BAMS=$(ls "$ALIGNED_DIR"/SRR*/Aligned.sortedByCoord.out.bam 2>/dev/null | wc -l)
echo "Total BAM files: $TOTAL_BAMS / 32"
echo ""

if [ "$TOTAL_BAMS" -eq 32 ]; then
    echo "✓ All samples aligned successfully!"
else
    echo "⚠ Some samples failed. Check logs in $LOG_DIR"
fi
