#!/bin/bash

# Complete missing bulk RNA-seq samples for TEtranscripts analysis
# SRR14514235 (Sample-17) and SRR14514163 (Sample-100)

set -e  # Exit on error

# Directories
SRA_DIR=~/sc/sra_downloads
ALIGN_DIR=~/sc/star_bulk_aligned
GENOME_DIR=~/sc/star_index
TE_DIR=~/sc/tetranscripts_bulk
GENE_GTF=~/sc/annotations/gencode.v45.primary_assembly.annotation.gtf
TE_GTF=~/sc/annotations/GRCh38_GENCODE_rmsk_TE.gtf

# Samples
SAMPLES=("SRR14514235" "SRR14514163")

# Activate tetranscripts environment
source ~/sc/tetranscripts_env/bin/activate

echo "Starting download and alignment of missing samples..."
echo "Samples: ${SAMPLES[@]}"
echo "Start time: $(date)"

for SRR in "${SAMPLES[@]}"; do
    echo ""
    echo "=========================================="
    echo "Processing $SRR"
    echo "=========================================="
    
    # Create directories
    mkdir -p ${SRA_DIR}/${SRR}
    mkdir -p ${ALIGN_DIR}/${SRR}
    
    # Download SRA file
    echo "[$(date)] Downloading $SRR..."
    cd ${SRA_DIR}/${SRR}
    if [ ! -f ${SRR}/${SRR}.sra ]; then
        prefetch ${SRR} -O ${SRA_DIR}/${SRR}
    else
        echo "SRA file already exists, skipping download"
    fi
    
    # Convert to FASTQ
    echo "[$(date)] Converting to FASTQ..."
    if [ ! -f ${SRR}_1.fastq ] || [ ! -f ${SRR}_2.fastq ]; then
        fasterq-dump ${SRR}/${SRR}.sra -O ${SRA_DIR}/${SRR} --split-files --threads 8
    else
        echo "FASTQ files already exist, skipping conversion"
    fi
    
    # Run STAR alignment
    echo "[$(date)] Running STAR alignment..."
    cd ${ALIGN_DIR}/${SRR}
    STAR --runThreadN 8 \
         --genomeDir ${GENOME_DIR} \
         --readFilesIn ${SRA_DIR}/${SRR}/${SRR}_1.fastq ${SRA_DIR}/${SRR}/${SRR}_2.fastq \
         --outFileNamePrefix ${ALIGN_DIR}/${SRR}/ \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes Standard
    
    # Check if BAM file was created
    if [ ! -f ${ALIGN_DIR}/${SRR}/Aligned.sortedByCoord.out.bam ]; then
        echo "ERROR: BAM file not created for $SRR"
        exit 1
    fi
    
    echo "[$(date)] Alignment complete for $SRR"
    
    # Delete FASTQ files to save space
    echo "[$(date)] Cleaning up FASTQ files..."
    rm -f ${SRA_DIR}/${SRR}/${SRR}_*.fastq
    
    # Run TEcount
    echo "[$(date)] Running TEcount..."
    TEcount --BAM ${ALIGN_DIR}/${SRR}/Aligned.sortedByCoord.out.bam \
            --GTF ${GENE_GTF} \
            --TE ${TE_GTF} \
            --mode multi \
            --stranded reverse \
            --sortByPos \
            --format BAM \
            --project ${TE_DIR}/${SRR} \
            2>&1 | tee ${TE_DIR}/${SRR}.log
    
    echo "[$(date)] Completed $SRR"
    echo "Output: ${TE_DIR}/${SRR}.cntTable"
done

echo ""
echo "=========================================="
echo "All missing samples completed!"
echo "End time: $(date)"
echo "=========================================="
echo ""
echo "Results:"
ls -lh ${TE_DIR}/*.cntTable | wc -l
echo "count tables generated"
