#!/bin/bash

# Run TEtranscripts on ALL bulk RNA-seq samples (156 aligned samples)
# Expands analysis beyond 18 multi-omics patients to full cohort

set -e

# Paths
BULK_DIR=~/sc/star_bulk_aligned
OUTPUT_DIR=~/sc/tetranscripts_bulk_all
GENE_GTF=~/sc/annotations/gencode.v45.primary_assembly.annotation.gtf
TE_GTF=~/sc/annotations/GRCh38_GENCODE_rmsk_TE.gtf
VENV=~/sc/tetranscripts_env

# Number of parallel jobs (DGX has 80 cores)
NJOBS=20

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Activate virtual environment
source ${VENV}/bin/activate

echo "========================================================================"
echo "TEtranscripts Analysis: ALL Bulk RNA-seq Samples"
echo "========================================================================"
echo "Total samples to process: 156"
echo "Parallel jobs: ${NJOBS}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Start time: $(date)"
echo "========================================================================"
echo ""

# Function to process one sample
process_sample() {
    local sample=$1
    local bam_file="${BULK_DIR}/${sample}/Aligned.sortedByCoord.out.bam"
    local output_prefix="${OUTPUT_DIR}/${sample}"
    
    # Skip if already processed
    if [ -f "${output_prefix}.cntTable" ]; then
        echo "✓ ${sample} already processed, skipping"
        return 0
    fi
    
    if [ -f "$bam_file" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing ${sample}..."
        
        TEcount \
            --BAM $bam_file \
            --GTF $GENE_GTF \
            --TE $TE_GTF \
            --mode multi \
            --stranded reverse \
            --sortByPos \
            --format BAM \
            --project ${output_prefix} \
            2>&1 | tee ${output_prefix}.log
            
        if [ $? -eq 0 ]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✓ Completed ${sample}"
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✗ Failed ${sample}"
            return 1
        fi
    else
        echo "⚠ Warning: BAM file not found for ${sample}"
        return 1
    fi
}

export -f process_sample
export BULK_DIR OUTPUT_DIR GENE_GTF TE_GTF VENV

# Get list of all aligned samples
SAMPLE_LIST="${OUTPUT_DIR}/sample_list.txt"
ls ${BULK_DIR} | grep '^SRR' > ${SAMPLE_LIST}

TOTAL_SAMPLES=$(wc -l < ${SAMPLE_LIST})
echo "Found ${TOTAL_SAMPLES} aligned samples"
echo ""

# Run in parallel using GNU parallel or xargs
echo "Starting parallel processing..."
echo ""

if command -v parallel &> /dev/null; then
    # GNU parallel with progress bar
    cat ${SAMPLE_LIST} | parallel --progress --joblog ${OUTPUT_DIR}/parallel.log -j ${NJOBS} process_sample {}
else
    # Fallback to xargs
    cat ${SAMPLE_LIST} | xargs -P ${NJOBS} -I {} bash -c 'process_sample "$@"' _ {}
fi

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Processing Complete"
echo "========================================================================"
echo "End time: $(date)"
echo ""

# Count completed samples
COMPLETED=$(ls ${OUTPUT_DIR}/*.cntTable 2>/dev/null | wc -l)
echo "Completed samples: ${COMPLETED} / ${TOTAL_SAMPLES}"

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✓ All samples processed successfully"
else
    echo "⚠ Some samples failed - check logs in ${OUTPUT_DIR}/*.log"
fi

echo ""
echo "Results directory: ${OUTPUT_DIR}"
echo "========================================================================"

exit ${EXIT_CODE}
