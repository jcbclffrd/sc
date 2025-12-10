#!/bin/bash

# STARsolo alignment script for 10x Chromium v3 snRNA-seq data
# PRJNA729525 - Alzheimer's disease brain tissue
# 
# Read structure: _2 file = R1 (100bp with 16bp barcode + 12bp UMI embedded)
#                 _3 file = R2 (100bp cDNA)

# Configuration
STAR_INDEX="$HOME/sc/star_index"
WHITELIST="$HOME/sc/annotations/3M-february-2018.txt"
THREADS=20
SAMPLE_LIST="ready_for_alignment.txt"
OUTPUT_DIR="starsolo_aligned"

# 10x Chromium v3 chemistry parameters
CB_START=1
CB_LEN=16
UMI_START=17
UMI_LEN=12

echo "================================================"
echo "STARsolo Alignment Pipeline"
echo "10x Chromium Single Cell 3' v3"
echo "Started: $(date)"
echo "================================================"
echo ""
echo "Configuration:"
echo "  STAR index: $STAR_INDEX"
echo "  Whitelist: $WHITELIST"
echo "  Threads: $THREADS"
echo "  Cell barcode: ${CB_LEN}bp (start: $CB_START)"
echo "  UMI: ${UMI_LEN}bp (start: $UMI_START)"
echo ""

# Check if STAR is available
if ! command -v STAR &> /dev/null; then
    echo "ERROR: STAR not found in PATH"
    exit 1
fi

# Check STAR version
echo "STAR version:"
STAR --version
echo ""

# Count samples
TOTAL=$(wc -l < $SAMPLE_LIST)
CURRENT=0

echo "Processing $TOTAL samples"
echo "================================================"
echo ""

while read srr; do
    CURRENT=$((CURRENT + 1))
    
    echo "[$CURRENT/$TOTAL] Processing $srr"
    echo "Started: $(date)"
    echo "---"
    
    # Define file paths
    # For 10x data in SRA: _2 = barcode/UMI (R1), _3 = cDNA (R2)
    R1_FILE="sra_downloads/${srr}/${srr}_2.fastq.gz"  # Barcode + UMI
    R2_FILE="sra_downloads/${srr}/${srr}_3.fastq.gz"  # cDNA
    
    if [[ ! -f "$R1_FILE" ]] || [[ ! -f "$R2_FILE" ]]; then
        echo "ERROR: FASTQ files not found for $srr"
        echo "  R1: $R1_FILE"
        echo "  R2: $R2_FILE"
        continue
    fi
    
    SAMPLE_OUT="$OUTPUT_DIR/$srr"
    mkdir -p $SAMPLE_OUT
    
    # Skip if already processed
    if [[ -f "$SAMPLE_OUT/Aligned.sortedByCoord.out.bam" ]]; then
        echo "Already processed, skipping..."
        echo ""
        continue
    fi
    
    # Run STARsolo
    STAR --runThreadN $THREADS \
         --genomeDir $STAR_INDEX \
         --readFilesCommand zcat \
         --readFilesIn $R2_FILE $R1_FILE \
         --soloType CB_UMI_Simple \
         --soloCBwhitelist $WHITELIST \
         --soloCBstart $CB_START \
         --soloCBlen $CB_LEN \
         --soloUMIstart $UMI_START \
         --soloUMIlen $UMI_LEN \
         --soloBarcodeReadLength 100 \
         --soloFeatures Gene GeneFull \
         --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
         --soloUMIfiltering MultiGeneUMI_CR \
         --soloUMIdedup 1MM_CR \
         --outSAMattributes NH HI AS nM CR UR CB UB GX GN \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix $SAMPLE_OUT/ \
         --limitBAMsortRAM 60000000000 \
         --outFilterMultimapNmax 100 \
         --winAnchorMultimapNmax 100    # Check if alignment succeeded
    if [[ -f "$SAMPLE_OUT/Aligned.sortedByCoord.out.bam" ]]; then
        # Index BAM file
        samtools index $SAMPLE_OUT/Aligned.sortedByCoord.out.bam
        
        # Get cell count
        if [[ -f "$SAMPLE_OUT/Solo.out/Gene/filtered/barcodes.tsv" ]]; then
            CELL_COUNT=$(wc -l < $SAMPLE_OUT/Solo.out/Gene/filtered/barcodes.tsv)
            echo "✓ Cells detected: $CELL_COUNT"
        fi
        
        # Show summary stats
        if [[ -f "$SAMPLE_OUT/Log.final.out" ]]; then
            echo ""
            echo "Alignment summary:"
            grep "Uniquely mapped reads %" $SAMPLE_OUT/Log.final.out
            grep "Number of reads mapped to multiple loci" $SAMPLE_OUT/Log.final.out | head -1
        fi
        
        echo "✓ Completed $srr"
    else
        echo "✗ FAILED: $srr"
    fi
    
    echo ""
    
done < $SAMPLE_LIST

echo "================================================"
echo "STARsolo Alignment Complete!"
echo "Finished: $(date)"
echo "================================================"
echo ""
echo "Summary:"
SUCCESS=$(find $OUTPUT_DIR -name "Aligned.sortedByCoord.out.bam" 2>/dev/null | wc -l)
echo "  Successfully aligned: $SUCCESS / $TOTAL samples"
echo ""
echo "Next step: Run scTE for TE quantification"
