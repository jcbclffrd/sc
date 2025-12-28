#!/bin/bash
#
# Run TEtranscripts on all remaining bulk RNA-seq BAM files
# This processes the ~138 samples that haven't been quantified yet
#
# Usage: bash run_tetranscripts_remaining_bulk_FIXED.sh
#
# Run on: spark-bd86 server

set -e

# =============================================================================
# CONFIGURATION
# =============================================================================

# Directories
BAM_DIR=~/sc/star_bulk_aligned
OUTPUT_DIR=~/sc/tetranscripts_bulk
ANNOTATION_DIR=~/sc/annotations

# Annotation files - CORRECTED PATHS
GTF=${ANNOTATION_DIR}/gencode.v45.primary_assembly.annotation.gtf
TE_GTF=${ANNOTATION_DIR}/GRCh38_GENCODE_rmsk_TE.gtf

# Python environment with TEtranscripts
PYTHON_ENV=~/sc/tetranscripts_env/bin/activate

# Parallel processing (adjusted for 20-core server)
N_THREADS=8
MAX_PARALLEL_JOBS=2  # Run 2 samples in parallel (2×8=16 cores, leaving 4 for system)

# Log directory
LOG_DIR=~/sc/logs/tetranscripts_bulk
mkdir -p ${LOG_DIR}

# =============================================================================
# VALIDATION
# =============================================================================

echo "========================================================================"
echo "TEtranscripts BULK RNA-SEQ QUANTIFICATION - REMAINING SAMPLES"
echo "========================================================================"
echo ""
echo "Configuration:"
echo "  BAM directory: ${BAM_DIR}"
echo "  Output directory: ${OUTPUT_DIR}"
echo "  GTF: ${GTF}"
echo "  TE GTF: ${TE_GTF}"
echo "  Parallel jobs: ${MAX_PARALLEL_JOBS}"
echo "  Threads per job: ${N_THREADS}"
echo ""

# Check directories exist
if [ ! -d "${BAM_DIR}" ]; then
    echo "ERROR: BAM directory not found: ${BAM_DIR}"
    exit 1
fi

if [ ! -f "${GTF}" ]; then
    echo "ERROR: GTF file not found: ${GTF}"
    exit 1
fi

if [ ! -f "${TE_GTF}" ]; then
    echo "ERROR: TE GTF file not found: ${TE_GTF}"
    exit 1
fi

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Activate TEtranscripts environment
echo "Activating TEtranscripts environment..."
source ${PYTHON_ENV}

# =============================================================================
# FIND SAMPLES TO PROCESS
# =============================================================================

echo ""
echo "Scanning for BAM files and checking completion status..."
echo ""

# Find all BAM files
ALL_BAMS=($(find ${BAM_DIR} -name "Aligned.sortedByCoord.out.bam" | sort))
TOTAL_BAMS=${#ALL_BAMS[@]}

echo "Found ${TOTAL_BAMS} total BAM files"

# Check which ones need processing
SAMPLES_TO_PROCESS=()
SAMPLES_COMPLETED=()

for bam_file in "${ALL_BAMS[@]}"; do
    # Extract sample ID from path: .../SRR14514XXX/Aligned.sortedByCoord.out.bam -> SRR14514XXX
    sample_id=$(basename $(dirname ${bam_file}))
    
    # Skip if .cntTable already exists
    if [ -f "${OUTPUT_DIR}/${sample_id}.cntTable" ]; then
        SAMPLES_COMPLETED+=("${sample_id}")
    else
        SAMPLES_TO_PROCESS+=("${sample_id}")
    fi
done

N_COMPLETED=${#SAMPLES_COMPLETED[@]}
N_TO_PROCESS=${#SAMPLES_TO_PROCESS[@]}

echo ""
echo "Sample status:"
echo "  Already completed: ${N_COMPLETED}"
echo "  Remaining to process: ${N_TO_PROCESS}"
echo ""

if [ ${N_TO_PROCESS} -eq 0 ]; then
    echo "✓ All samples already processed!"
    exit 0
fi

# =============================================================================
# PROCESS SAMPLES
# =============================================================================

echo "========================================================================"
echo "PROCESSING ${N_TO_PROCESS} SAMPLES"
echo "========================================================================"
echo ""
echo "Start time: $(date)"
echo ""

# Function to process a single sample
process_sample() {
    local sample_id=$1
    local bam_file="${BAM_DIR}/${sample_id}/Aligned.sortedByCoord.out.bam"
    local output_prefix="${OUTPUT_DIR}/${sample_id}"
    local log_file="${LOG_DIR}/${sample_id}_tetranscripts.log"
    
    echo "[$(date +%H:%M:%S)] Starting: ${sample_id}"
    
    # Run TEtranscripts
    TEcount \
        --BAM ${bam_file} \
        --GTF ${GTF} \
        --TE ${TE_GTF} \
        --mode multi \
        --stranded reverse \
        --project ${output_prefix} \
        --sortByPos \
        > ${log_file} 2>&1
    
    if [ $? -eq 0 ]; then
        echo "[$(date +%H:%M:%S)] ✓ Completed: ${sample_id}"
        return 0
    else
        echo "[$(date +%H:%M:%S)] ✗ FAILED: ${sample_id} (see ${log_file})"
        return 1
    fi
}

export -f process_sample
export BAM_DIR OUTPUT_DIR GTF TE_GTF LOG_DIR

# Process samples in parallel using GNU parallel or xargs
if command -v parallel &> /dev/null; then
    # Use GNU parallel if available (better)
    echo "Using GNU parallel for processing..."
    printf '%s\n' "${SAMPLES_TO_PROCESS[@]}" | \
        parallel -j ${MAX_PARALLEL_JOBS} --eta process_sample {}
else
    # Fallback to sequential processing
    echo "GNU parallel not found, processing sequentially..."
    for sample_id in "${SAMPLES_TO_PROCESS[@]}"; do
        process_sample ${sample_id}
    done
fi

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "========================================================================"
echo "PROCESSING COMPLETE"
echo "========================================================================"
echo ""
echo "End time: $(date)"
echo ""

# Count completed samples
N_SUCCESS=$(find ${OUTPUT_DIR} -name "*.cntTable" | wc -l)
echo "Successfully quantified: ${N_SUCCESS} / ${TOTAL_BAMS} total samples"
echo ""

# List any that failed
FAILED_SAMPLES=()
for sample_id in "${SAMPLES_TO_PROCESS[@]}"; do
    if [ ! -f "${OUTPUT_DIR}/${sample_id}.cntTable" ]; then
        FAILED_SAMPLES+=("${sample_id}")
    fi
done

if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo "Failed samples (${#FAILED_SAMPLES[@]}):"
    printf '  %s\n' "${FAILED_SAMPLES[@]}"
    echo ""
    echo "Check logs in: ${LOG_DIR}"
else
    echo "✓ All samples processed successfully!"
fi

echo ""
echo "Output files:"
echo "  Count tables: ${OUTPUT_DIR}/*.cntTable"
echo "  Logs: ${LOG_DIR}/"
echo ""

# =============================================================================
# NEXT STEPS
# =============================================================================

echo "Next steps:"
echo "  1. Create metadata for all ${N_SUCCESS} samples:"
echo "     python3 scripts/create_bulk_metadata_all_samples.py \\"
echo "         --tecount-dir ${OUTPUT_DIR} \\"
echo "         --output ${OUTPUT_DIR}/bulk_metadata_all_samples.csv"
echo ""
echo "  2. Run complete analysis pipeline:"
echo "     bash scripts/run_bulk_te_analysis_pipeline.sh"
echo ""
