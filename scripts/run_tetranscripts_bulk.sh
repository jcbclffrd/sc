#!/bin/bash

# Run TEtranscripts on bulk RNA-seq data from Morabito et al.
# This script quantifies both genes and TEs from STAR alignments in parallel

set -e

# Paths
BULK_DIR=~/sc/star_bulk_aligned
OUTPUT_DIR=~/sc/tetranscripts_bulk
GENE_GTF=~/sc/annotations/gencode.v45.primary_assembly.annotation.gtf
TE_GTF=~/sc/annotations/GRCh38_GENCODE_rmsk_TE.gtf
VENV=~/sc/tetranscripts_env
SAMPLE_LIST=~/sc/data/bulk_rnaseq_18patients_srr.txt

# Number of parallel jobs
NJOBS=6

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Activate virtual environment
source ${VENV}/bin/activate

echo "Starting TEtranscripts analysis on 18 patients with $NJOBS parallel jobs..."

# Function to process one sample
process_sample() {
    local sample=$1
    local bam_file="${BULK_DIR}/${sample}/Aligned.sortedByCoord.out.bam"
    local output_prefix="${OUTPUT_DIR}/${sample}"
    
    if [ -f "$bam_file" ]; then
        echo "Processing $sample..."
        
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
            
        echo "Completed $sample"
    else
        echo "Warning: BAM file not found for $sample"
    fi
}

export -f process_sample
export BULK_DIR OUTPUT_DIR GENE_GTF TE_GTF

# Run in parallel using GNU parallel or xargs
if command -v parallel &> /dev/null; then
    cat ${SAMPLE_LIST} | parallel -j ${NJOBS} process_sample {}
else
    cat ${SAMPLE_LIST} | xargs -P ${NJOBS} -I {} bash -c 'process_sample "$@"' _ {}
fi

echo "All samples completed!"

echo "All samples processed!"
echo "Results in: ${OUTPUT_DIR}"
