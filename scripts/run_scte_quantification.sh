#!/bin/bash

#######################################################
# Run scTE for TE Quantification on All Samples
# Single-cell TE expression quantification
#######################################################

set -e

# Activate virtual environment
source /home/jacobc/hcaTE/.venv/bin/activate

# Paths
SCTE_INDEX="/home/jacobc/sc/annotations/hg38_gencode45"
ALIGNED_DIR="/home/jacobc/sc/starsolo_aligned"
OUTPUT_DIR="/home/jacobc/sc/scTE_output"
THREADS=10  # Reduced from 20 to avoid silent crashes

# Check if index exists
if [ ! -f "${SCTE_INDEX}.exclusive.idx" ]; then
    echo "ERROR: scTE index not found at ${SCTE_INDEX}.exclusive.idx"
    echo "Please run build_scte_index.sh first"
    exit 1
fi

echo "================================================"
echo "scTE TE Quantification Pipeline"
echo "Started: $(date)"
echo "================================================"
echo ""
echo "Configuration:"
echo "  scTE index: ${SCTE_INDEX}"
echo "  Input directory: ${ALIGNED_DIR}"
echo "  Output directory: ${OUTPUT_DIR}"
echo "  Threads: ${THREADS}"
echo ""

# Get list of samples
SAMPLES=($(ls -d ${ALIGNED_DIR}/SRR* | xargs -n 1 basename))
TOTAL=${#SAMPLES[@]}

echo "Processing ${TOTAL} samples"
echo "================================================"
echo ""

# Process each sample
COUNT=0
SUCCESS=0
FAILED=0

for SAMPLE in "${SAMPLES[@]}"; do
    COUNT=$((COUNT + 1))
    echo "[${COUNT}/${TOTAL}] Processing ${SAMPLE}"
    echo "Started: $(date)"
    echo "---"
    
    BAM_FILE="${ALIGNED_DIR}/${SAMPLE}/Aligned.sortedByCoord.out.bam"
    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
    
    # Check if BAM file exists
    if [ ! -f "$BAM_FILE" ]; then
        echo "ERROR: BAM file not found: $BAM_FILE"
        FAILED=$((FAILED + 1))
        continue
    fi
    
    # Create output directory
    mkdir -p "$SAMPLE_OUT"
    
    # Run scTE (need to cd into output dir due to scTE path handling bug)
    # Note: scTE expects cell barcodes in CB tag (which STARsolo provides)
    (cd "$SAMPLE_OUT" && \
     scTE \
        -i "$(realpath $BAM_FILE)" \
        -o "${SAMPLE}" \
        -x "$(realpath ${SCTE_INDEX}.exclusive.idx)" \
        -CB CR \
        -UMI UR \
        -p $THREADS \
        --hdf5 True)
    
    if [ $? -eq 0 ]; then
        echo "✓ Completed ${SAMPLE}"
        SUCCESS=$((SUCCESS + 1))
        
        # Show summary if available
        if [ -f "${SAMPLE_OUT}/scTE.stats.txt" ]; then
            echo ""
            echo "Summary:"
            head -20 "${SAMPLE_OUT}/scTE.stats.txt"
        fi
    else
        echo "✗ Failed ${SAMPLE}"
        FAILED=$((FAILED + 1))
    fi
    
    echo ""
done

echo "================================================"
echo "scTE TE Quantification Complete!"
echo "Finished: $(date)"
echo "================================================"
echo ""
echo "Summary:"
echo "  Successfully processed: ${SUCCESS} / ${TOTAL} samples"
echo "  Failed: ${FAILED} samples"
echo ""
echo "Next step: Differential TE expression analysis"
