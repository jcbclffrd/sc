#!/bin/bash

#####################################################################
# Bulk RNA-seq STAR Alignment Script
# For PRJNA729525 bulk RNA-seq samples (200bp paired-end)
# Standard STAR alignment (not STARsolo)
#####################################################################

# Configuration
SAMPLE_LIST="/home/jacobc/sc/data/bulk_rnaseq_samples.txt"
FASTQ_DIR="/home/jacobc/sc/sra_downloads"
OUTPUT_DIR="/home/jacobc/sc/star_bulk_aligned"
STAR_INDEX="/home/jacobc/sc/star_index"
THREADS=20

# TE-specific multimapping parameter
MULTIMAP_MAX=100  # Allow up to 100 alignments per read for TE detection

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="${OUTPUT_DIR}/bulk_alignment.log"
echo "=== Bulk RNA-seq STAR Alignment Started: $(date) ===" | tee "$LOG_FILE"
echo "Samples to process: $(wc -l < $SAMPLE_LIST)" | tee -a "$LOG_FILE"
echo "Threads: $THREADS" | tee -a "$LOG_FILE"
echo "Multimapping max: $MULTIMAP_MAX" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Counter for tracking progress
TOTAL=$(wc -l < "$SAMPLE_LIST")
COUNT=0
SUCCESS=0
FAILED=0

# Process each sample
while IFS= read -r SAMPLE; do
    COUNT=$((COUNT + 1))
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing sample $COUNT/$TOTAL: $SAMPLE" | tee -a "$LOG_FILE"
    
    # Define input files
    R1="${FASTQ_DIR}/${SAMPLE}/${SAMPLE}_1.fastq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}/${SAMPLE}_2.fastq.gz"
    
    # Check if FASTQ files exist
    if [[ ! -f "$R1" ]] || [[ ! -f "$R2" ]]; then
        echo "  ⚠️  FASTQ files not found for $SAMPLE - SKIPPING" | tee -a "$LOG_FILE"
        FAILED=$((FAILED + 1))
        continue
    fi
    
    # Define output directory for this sample
    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
    mkdir -p "$SAMPLE_OUT"
    
    # Check if alignment already exists
    if [[ -f "${SAMPLE_OUT}/Aligned.sortedByCoord.out.bam" ]]; then
        BAM_SIZE=$(stat -c%s "${SAMPLE_OUT}/Aligned.sortedByCoord.out.bam")
        if [[ $BAM_SIZE -gt 1000000 ]]; then
            echo "  ✓ Alignment already exists for $SAMPLE ($(numfmt --to=iec-i --suffix=B $BAM_SIZE)) - SKIPPING" | tee -a "$LOG_FILE"
            SUCCESS=$((SUCCESS + 1))
            continue
        else
            echo "  ⚠️  Existing BAM file is too small ($BAM_SIZE bytes), re-aligning..." | tee -a "$LOG_FILE"
            rm -rf "$SAMPLE_OUT"
            mkdir -p "$SAMPLE_OUT"
        fi
    fi
    
    # Run STAR alignment
    echo "  Running STAR alignment..." | tee -a "$LOG_FILE"
    
    STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${SAMPLE_OUT}/" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes NH HI AS nM NM MD \
        --outFilterMultimapNmax "$MULTIMAP_MAX" \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1 \
        --limitBAMsortRAM 30000000000 \
        2>&1 | tee -a "$LOG_FILE"
    
    # Check if alignment was successful
    if [[ -f "${SAMPLE_OUT}/Aligned.sortedByCoord.out.bam" ]]; then
        BAM_SIZE=$(stat -c%s "${SAMPLE_OUT}/Aligned.sortedByCoord.out.bam")
        if [[ $BAM_SIZE -gt 1000000 ]]; then
            echo "  ✓ SUCCESS: $SAMPLE aligned ($(numfmt --to=iec-i --suffix=B $BAM_SIZE))" | tee -a "$LOG_FILE"
            SUCCESS=$((SUCCESS + 1))
            
            # Index the BAM file
            echo "  Indexing BAM file..." | tee -a "$LOG_FILE"
            samtools index "${SAMPLE_OUT}/Aligned.sortedByCoord.out.bam"
        else
            echo "  ✗ FAILED: $SAMPLE - BAM file too small ($BAM_SIZE bytes)" | tee -a "$LOG_FILE"
            FAILED=$((FAILED + 1))
        fi
    else
        echo "  ✗ FAILED: $SAMPLE - BAM file not created" | tee -a "$LOG_FILE"
        FAILED=$((FAILED + 1))
    fi
    
    echo "" | tee -a "$LOG_FILE"
    
done < "$SAMPLE_LIST"

# Final summary
echo "=== Bulk RNA-seq STAR Alignment Completed: $(date) ===" | tee -a "$LOG_FILE"
echo "Total samples: $TOTAL" | tee -a "$LOG_FILE"
echo "Successful: $SUCCESS" | tee -a "$LOG_FILE"
echo "Failed: $FAILED" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

if [[ $FAILED -eq 0 ]]; then
    echo "🎉 All samples aligned successfully!" | tee -a "$LOG_FILE"
else
    echo "⚠️  Some samples failed. Check the log for details." | tee -a "$LOG_FILE"
fi
