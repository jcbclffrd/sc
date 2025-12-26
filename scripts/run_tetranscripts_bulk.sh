#!/bin/bash

# Run TEtranscripts on bulk RNA-seq data from Morabito et al.
# This script quantifies both genes and TEs from STAR alignments

set -e

# Paths
BULK_DIR=~/sc/star_bulk_aligned
OUTPUT_DIR=~/sc/tetranscripts_bulk
GENE_GTF=~/sc/annotations/gencode.v45.primary_assembly.annotation.gtf
TE_GTF=~/sc/annotations/hg38_rmsk.gtf
VENV=~/sc/tetranscripts_env

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Activate virtual environment
source ${VENV}/bin/activate

# Get list of 18 patient samples only
SAMPLE_LIST="/home/jacobc/sc/data/bulk_rnaseq_18patients_srr.txt"

cd ${BULK_DIR}

echo "Starting TEtranscripts analysis on 18 patients..."

# Run TEtranscripts only on the 18 patient samples
while read sample; do
    bam_file="${BULK_DIR}/${sample}/Aligned.sortedByCoord.out.bam"
    output_prefix="${OUTPUT_DIR}/${sample}"
    
    if [ -f "$bam_file" ]; then
        echo "Processing $sample..."
        
        # Run TEtranscripts
        TEcount \
            --BAM $bam_file \
            --GTF $GENE_GTF \
            --TE $TE_GTF \
            --mode multi \
            --stranded reverse \
            --format BAM \
            --project ${output_prefix} \
            2>&1 | tee ${output_prefix}.log
            
        echo "Completed $sample"
    else
        echo "Warning: BAM file not found for $sample"
    fi
done < ${SAMPLE_LIST}

echo "All samples processed!"
echo "Results in: ${OUTPUT_DIR}"
