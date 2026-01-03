#!/bin/bash
#
# Run TEtranscripts on a single bulk RNA-seq sample (for testing/debugging)
#
# Usage: bash run_tetranscripts_single_sample.sh SRR14514XXX
#

if [ $# -ne 1 ]; then
    echo "Usage: $0 <sample_id>"
    echo "Example: $0 SRR14514192"
    exit 1
fi

SAMPLE_ID=$1

# Directories
BAM_DIR=~/sc/star_bulk_aligned
OUTPUT_DIR=~/sc/tetranscripts_bulk
ANNOTATION_DIR=~/sc/annotations

# Files
BAM_FILE=${BAM_DIR}/${SAMPLE_ID}/Aligned.sortedByCoord.out.bam
GTF=${ANNOTATION_DIR}/gencode.v38.annotation.gtf
TE_GTF=${ANNOTATION_DIR}/GRCh38_GENCODE_rmsk_TE.gtf

# Check inputs
if [ ! -f "${BAM_FILE}" ]; then
    echo "ERROR: BAM file not found: ${BAM_FILE}"
    exit 1
fi

if [ ! -f "${GTF}" ] || [ ! -f "${TE_GTF}" ]; then
    echo "ERROR: Annotation files not found"
    exit 1
fi

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Activate environment
source ~/sc/tetranscripts_env/bin/activate

echo "Processing sample: ${SAMPLE_ID}"
echo "BAM: ${BAM_FILE}"
echo ""

# Run TEcount
TEcount \
    --BAM ${BAM_FILE} \
    --GTF ${GTF} \
    --TE ${TE_GTF} \
    --mode multi \
    --stranded reverse \
    --project ${OUTPUT_DIR}/${SAMPLE_ID} \
    --sortByPos

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Success! Output: ${OUTPUT_DIR}/${SAMPLE_ID}.cntTable"
else
    echo ""
    echo "✗ Failed!"
    exit 1
fi
