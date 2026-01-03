#!/bin/bash

# Align remaining ATAC-seq samples from existing FASTQs
# Lower parallelization (2 jobs) to avoid OOM

set -e

# Configuration
MAX_JOBS=2  # Conservative to avoid OOM
THREADS=6   # 2 jobs × 6 threads = 12 cores used
GENOME_DIR="$HOME/sc/star_index"
SRA_DIR="$HOME/sc/sra_downloads/ATAC-seq"
ALIGNED_DIR="$HOME/sc/atacseq_aligned"
LOG_DIR="$HOME/sc/logs_atacseq"

# Add tools to PATH
export PATH=$HOME/sratoolkit.3.3.0-ubuntu64/bin:$HOME/.local/bin:$PATH

mkdir -p "$ALIGNED_DIR" "$LOG_DIR"

# Samples to align (have FASTQs but no valid BAMs)
SAMPLES=(
    SRR14514129 SRR14514130 SRR14514131 SRR14514132
    SRR14514133 SRR14514134 SRR14514135 SRR14514136
    SRR14514137 SRR14514138 SRR14514139 SRR14514140
    SRR14514142 SRR14514143 SRR14514149 SRR14514150
    SRR14514151 SRR14514153 SRR14514154 SRR14514156
    SRR14514157 SRR14514158 SRR14514159
)

echo "=========================================="
echo "Aligning Remaining ATAC-seq Samples"
echo "=========================================="
echo "Samples to align: ${#SAMPLES[@]}"
echo "Parallel jobs: $MAX_JOBS"
echo "Threads per job: $THREADS"
echo "Free space: $(df -h / | tail -1 | awk '{print $4}')"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Function to align one sample
align_sample() {
    local SAMPLE=$1
    local SAMPLE_OUT="$ALIGNED_DIR/$SAMPLE"
    local LOG_FILE="$LOG_DIR/${SAMPLE}_alignment.log"
    
    echo "[$SAMPLE] Started: $(date '+%H:%M:%S')"
    
    mkdir -p "$SAMPLE_OUT"
    
    # Technical reads (_1=barcodes length=8, _3=UMIs length=16)
    BARCODE="$SRA_DIR/${SAMPLE}_1.fastq.gz"
    UMI="$SRA_DIR/${SAMPLE}_3.fastq.gz"
    
    # Biological reads for alignment (_2 and _4, both length=50)
    READ2="$SRA_DIR/${SAMPLE}_2.fastq.gz"
    READ4="$SRA_DIR/${SAMPLE}_4.fastq.gz"
    
    # Check if FASTQs exist, if not download them
    if [ ! -f "$READ2" ] || [ ! -f "$READ4" ]; then
        echo "[$SAMPLE] Downloading from SRA..."
        cd "$SRA_DIR"
        
        # Remove old files with wrong format
        rm -f "${SAMPLE}_2.fastq.gz" "${SAMPLE}_4.fastq" "${SAMPLE}_4.fastq.gz" 2>/dev/null || true
        
        # Download with --include-technical to get all 4 files (_1=read1, _2=barcode, _3=read2, _4=UMI)
        if ! fasterq-dump "$SAMPLE" --split-files --include-technical --threads 2 --mem 4GB >> "$LOG_FILE" 2>&1; then
            echo "[$SAMPLE] ✗ Download failed"
            return 1
        fi
        
        # Compress all 4 FASTQs
        pigz -p 4 "${SAMPLE}_1.fastq" "${SAMPLE}_2.fastq" "${SAMPLE}_3.fastq" "${SAMPLE}_4.fastq" 2>> "$LOG_FILE" || true
        
        if [ ! -f "$READ2" ] || [ ! -f "$READ4" ]; then
            echo "[$SAMPLE] ✗ FASTQ files not found after download"
            return 1
        fi
    fi
    
    # STAR alignment with multi-mapping (n=100) for TE analysis
    # Use biological reads (_2 and _4, both length=50); technical reads (_1=barcodes, _3=UMIs) kept for demux/dedup
    if ! STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$READ2" "$READ4" \
        --readFilesCommand gunzip -c \
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
        --limitBAMsortRAM 10000000000 \
        >> "$LOG_FILE" 2>&1; then
        echo "[$SAMPLE] ✗ Alignment failed"
        return 1
    fi
    
    # Index BAM
    if ! samtools index "$SAMPLE_OUT/Aligned.sortedByCoord.out.bam" >> "$LOG_FILE" 2>&1; then
        echo "[$SAMPLE] ✗ Indexing failed"
        return 1
    fi
    
    # Clean up temp files to save space
    rm -rf "$SAMPLE_OUT/_STARtmp"
    
    BAM_SIZE=$(du -sh "$SAMPLE_OUT/Aligned.sortedByCoord.out.bam" | cut -f1)
    echo "[$SAMPLE] ✓ Aligned (BAM: $BAM_SIZE)"
}

export -f align_sample
export GENOME_DIR SRA_DIR ALIGNED_DIR LOG_DIR THREADS PATH

# Run alignment in parallel (2 jobs at a time)
echo "Running alignments..."
echo ""

printf '%s\n' "${SAMPLES[@]}" | xargs -P "$MAX_JOBS" -I {} bash -c 'align_sample "$@"' _ {}

echo ""
echo "=========================================="
echo "Alignment Complete"
echo "=========================================="
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "Checking results..."
TOTAL_BAMS=$(find "$ALIGNED_DIR" -name "Aligned.sortedByCoord.out.bam" -size +100M 2>/dev/null | wc -l)
TOTAL_INDEXED=$(find "$ALIGNED_DIR" -name "Aligned.sortedByCoord.out.bam.bai" -size +100k 2>/dev/null | wc -l)
echo "Valid BAM files: $TOTAL_BAMS / 32"
echo "Indexed BAMs: $TOTAL_INDEXED / 32"
echo "Free space remaining: $(df -h / | tail -1 | awk '{print $4}')"
echo ""

if [ "$TOTAL_INDEXED" -eq 32 ]; then
    echo "✓ All samples aligned and indexed successfully!"
else
    echo "⚠ Some samples incomplete. Check logs in $LOG_DIR"
    echo ""
    echo "Missing BAMs:"
    for i in $(seq 14514129 14514160); do
        s="SRR$i"
        if [ ! -s "$ALIGNED_DIR/$s/Aligned.sortedByCoord.out.bam" ]; then
            echo "  - $s"
        fi
    done
fi
