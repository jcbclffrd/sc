#!/bin/bash

# AWS ATAC-seq STAR Alignment Script
# Optimized for c5.9xlarge (36 cores, 72GB RAM)
# Downloads from SRA and aligns with STAR in parallel

set -e

# Configuration
FASTQ_DIR="$HOME/sc/sra_downloads/ATAC-seq"
OUTPUT_BASE="$HOME/sc/atacseq_aligned"
GENOME_DIR="$HOME/sc/star_index"
LOG_DIR="$HOME/sc/logs_atacseq"
THREADS=4          # Threads per sample
MAX_JOBS=9         # Run 9 samples in parallel (9 × 4 = 36 cores)

# Create directories
mkdir -p "$OUTPUT_BASE"
mkdir -p "$LOG_DIR"
mkdir -p "$FASTQ_DIR"

# All 32 ATAC-seq samples
SAMPLES=(
    SRR14514129 SRR14514130 SRR14514131 SRR14514132 SRR14514133
    SRR14514134 SRR14514135 SRR14514136 SRR14514137 SRR14514138
    SRR14514139 SRR14514140 SRR14514141 SRR14514142 SRR14514143
    SRR14514144 SRR14514145 SRR14514146 SRR14514147 SRR14514148
    SRR14514149 SRR14514150 SRR14514151 SRR14514152 SRR14514153
    SRR14514154 SRR14514155 SRR14514156 SRR14514157 SRR14514158
    SRR14514159 SRR14514160
)

echo "=========================================="
echo "AWS ATAC-seq STAR Alignment"
echo "=========================================="
echo "Instance: c5.9xlarge (36 cores, 72GB RAM)"
echo "Total samples: ${#SAMPLES[@]}"
echo "Parallel jobs: $MAX_JOBS (4 threads each)"
echo "Multi-mapping: 100 (for TE analysis)"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================="
echo ""

# Function to download and align a single sample
align_sample_aws() {
    local sample=$1
    
    echo "----------------------------------------"
    echo "Sample: $sample"
    echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
    
    # Check if already aligned
    if [[ -f "$OUTPUT_BASE/${sample}/Aligned.sortedByCoord.out.bam" ]]; then
        BAM_SIZE=$(stat -c%s "$OUTPUT_BASE/${sample}/Aligned.sortedByCoord.out.bam" 2>/dev/null || echo 0)
        if [[ $BAM_SIZE -gt 1000000000 ]]; then
            echo "✓ Already aligned (BAM: $(numfmt --to=iec-i --suffix=B $BAM_SIZE))"
            return 0
        fi
    fi
    
    SAMPLE_OUT="$OUTPUT_BASE/${sample}"
    mkdir -p "$SAMPLE_OUT"
    LOG_FILE="$LOG_DIR/${sample}_aws_alignment.log"
    
    # Download with fasterq-dump (streams from SRA in us-east-1)
    echo "Downloading from SRA..."
    cd "$FASTQ_DIR"
    
    if fasterq-dump "$sample" \
        --split-files \
        --threads 2 \
        --progress \
        --mem 4GB \
        --bufsize 1GB \
        > "$LOG_FILE" 2>&1; then
        echo "✓ Download complete"
    else
        echo "✗ Download failed"
        return 1
    fi
    
    # Compress immediately with pigz to save disk space
    echo "Compressing FASTQ files..."
    if [[ -f "${sample}_2.fastq" ]]; then
        pigz -p 2 "${sample}_2.fastq" &
    fi
    if [[ -f "${sample}_4.fastq" ]]; then
        pigz -p 2 "${sample}_4.fastq" &
    fi
    wait
    
    # Find FASTQ files
    if [[ -f "${sample}_2.fastq.gz" ]]; then
        READ1="${sample}_2.fastq.gz"
    elif [[ -f "${sample}_2.fastq" ]]; then
        READ1="${sample}_2.fastq"
    else
        echo "✗ ERROR: ${sample}_2 not found"
        return 1
    fi
    
    if [[ -f "${sample}_4.fastq.gz" ]]; then
        READ2="${sample}_4.fastq.gz"
    elif [[ -f "${sample}_4.fastq" ]]; then
        READ2="${sample}_4.fastq"
    else
        echo "✗ ERROR: ${sample}_4 not found"
        return 1
    fi
    
    # Determine read command
    if [[ "$READ1" == *.gz && "$READ2" == *.gz ]]; then
        READFILES_CMD="zcat"
    else
        READFILES_CMD="cat"
    fi
    
    echo "Aligning with STAR..."
    
    if STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$READ1" "$READ2" \
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
        >> "$LOG_FILE" 2>&1; then
        
        # Get statistics
        MAPPED=$(grep "Uniquely mapped reads number" "$SAMPLE_OUT/Log.final.out" | awk '{print $NF}')
        TOTAL=$(grep "Number of input reads" "$SAMPLE_OUT/Log.final.out" | awk '{print $NF}')
        BAM_SIZE=$(ls -lh "$SAMPLE_OUT/Aligned.sortedByCoord.out.bam" | awk '{print $5}')
        
        echo "✓ Alignment complete"
        echo "  Total reads: $TOTAL"
        echo "  Mapped reads: $MAPPED"
        echo "  BAM size: $BAM_SIZE"
        
        # Cleanup
        rm -rf "$SAMPLE_OUT/_STARtmp"
        
        # Delete FASTQ files to save space
        rm -f "$FASTQ_DIR/${sample}_2.fastq" "$FASTQ_DIR/${sample}_2.fastq.gz"
        rm -f "$FASTQ_DIR/${sample}_4.fastq" "$FASTQ_DIR/${sample}_4.fastq.gz"
        
    else
        echo "✗ Alignment failed"
        return 1
    fi
    
    echo "Complete: $(date '+%Y-%m-%d %H:%M:%S')"
    echo ""
}

export -f align_sample_aws
export FASTQ_DIR OUTPUT_BASE GENOME_DIR LOG_DIR THREADS

# Run alignments in parallel
echo "Starting parallel alignment (9 samples at a time)..."
printf '%s\n' "${SAMPLES[@]}" | xargs -P "$MAX_JOBS" -I {} bash -c 'align_sample_aws "$@"' _ {}

echo "=========================================="
echo "All Alignments Complete"
echo "=========================================="
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "Results summary:"
COMPLETED=$(ls -d "$OUTPUT_BASE"/SRR*/Aligned.sortedByCoord.out.bam 2>/dev/null | wc -l)
echo "Completed: $COMPLETED / ${#SAMPLES[@]} samples"
echo ""
echo "Total disk usage:"
du -sh "$OUTPUT_BASE"
echo ""
echo "Download BAMs to local machine:"
echo "  rsync -avz -e 'ssh -i ~/.ssh/awsWebsite.pem' ubuntu@\$(curl -s http://169.254.169.254/latest/meta-data/public-ipv4):~/sc/atacseq_aligned/ ~/sc/atacseq_aligned/"
echo "=========================================="
