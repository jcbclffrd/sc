#!/bin/bash
#
# Run TEtranscripts using bash job control (no GNU parallel needed)
#

set -e

# Configuration
BAM_DIR=~/sc/star_bulk_aligned
OUTPUT_DIR=~/sc/tetranscripts_bulk
ANNOTATION_DIR=~/sc/annotations
GTF=${ANNOTATION_DIR}/gencode.v45.primary_assembly.annotation.gtf
TE_GTF=${ANNOTATION_DIR}/GRCh38_GENCODE_rmsk_TE.gtf
PYTHON_ENV=~/sc/tetranscripts_env/bin/activate
MAX_PARALLEL_JOBS=16
LOG_DIR=~/sc/logs/tetranscripts_bulk

mkdir -p ${LOG_DIR}

echo "========================================================================"
echo "TEtranscripts BULK RNA-SEQ - Using bash job control"
echo "========================================================================"
echo ""

# Activate environment
source ${PYTHON_ENV}

# Find samples to process
ALL_BAMS=($(find ${BAM_DIR} -name "Aligned.sortedByCoord.out.bam" | sort))
SAMPLES_TO_PROCESS=()

for bam_file in "${ALL_BAMS[@]}"; do
    sample_id=$(basename $(dirname ${bam_file}))
    if [ ! -f "${OUTPUT_DIR}/${sample_id}.cntTable" ]; then
        SAMPLES_TO_PROCESS+=("${sample_id}")
    fi
done

echo "Samples to process: ${#SAMPLES_TO_PROCESS[@]}"
echo "Max parallel jobs: ${MAX_PARALLEL_JOBS}"
echo "Start: $(date)"
echo ""

# Process samples with job control
RUNNING_JOBS=0

for sample_id in "${SAMPLES_TO_PROCESS[@]}"; do
    # Wait if we have too many jobs running
    while [ ${RUNNING_JOBS} -ge ${MAX_PARALLEL_JOBS} ]; do
        sleep 10
        # Count running jobs
        RUNNING_JOBS=$(jobs -r | wc -l)
    done
    
    # Start new job in background
    (
        echo "[$(date +%H:%M:%S)] Starting: ${sample_id}"
        TEcount \
            --BAM ${BAM_DIR}/${sample_id}/Aligned.sortedByCoord.out.bam \
            --GTF ${GTF} \
            --TE ${TE_GTF} \
            --mode multi \
            --stranded reverse \
            --project ${OUTPUT_DIR}/${sample_id} \
            --sortByPos \
            > ${LOG_DIR}/${sample_id}_tetranscripts.log 2>&1
        
        if [ $? -eq 0 ]; then
            echo "[$(date +%H:%M:%S)] ✓ Completed: ${sample_id}"
        else
            echo "[$(date +%H:%M:%S)] ✗ FAILED: ${sample_id}"
        fi
    ) &
    
    RUNNING_JOBS=$((RUNNING_JOBS + 1))
    echo "Started ${sample_id} (${RUNNING_JOBS} running)"
done

# Wait for all jobs to complete
echo ""
echo "Waiting for all jobs to complete..."
wait

echo ""
echo "========================================================================"
echo "COMPLETE - $(date)"
echo "========================================================================"
echo ""
find ${OUTPUT_DIR} -name "*.cntTable" | wc -l
echo "samples quantified"
