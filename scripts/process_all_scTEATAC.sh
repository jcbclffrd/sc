#!/bin/bash

# Complete pipeline: CellRanger BAM → STAR multi-mapper → CB transfer → scTEATAC TE counts
# Processes all samples except SRR14514130 (test sample)

set -e

export PATH=~/.local/bin:$PATH

# Configuration
SAMPLES=(
    SRR14514131 SRR14514132 SRR14514133 SRR14514134
    SRR14514135 SRR14514136 SRR14514137 SRR14514138
    SRR14514139 SRR14514140 SRR14514141 SRR14514142
    SRR14514143 SRR14514144 SRR14514145 SRR14514146
    SRR14514147 SRR14514148 SRR14514149 SRR14514150
    SRR14514151 SRR14514152 SRR14514153 SRR14514154
    SRR14514155 SRR14514156 SRR14514157 SRR14514158
    SRR14514159 SRR14514160 SRR14514129
)

# Paths
GENOME=~/sc/genome/GRCh38.primary_assembly.genome.fa
STAR_INDEX=~/sc/star_index
TE_GTF=~/sc/annotations/GRCh38_GENCODE_rmsk_TE.gtf
FASTQ_DIR=~/sc/sra_downloads/ATAC-seq
CELLRANGER_DIR=~/sc/cellranger_atac_output
STAR_DIR=~/sc/atacseq_aligned
OUTPUT_DIR=~/sc/atacseq_te_counts
LOG_DIR=~/sc/logs_atacseq
THREADS=8

mkdir -p ${LOG_DIR}

# Build scTEATAC index once (if not exists)
INDEX_DIR=${OUTPUT_DIR}/scTEATAC_index
if [ ! -f "${INDEX_DIR}/hg38_te.idx" ]; then
    echo "Building scTEATAC index..."
    mkdir -p ${INDEX_DIR}
    cd ${INDEX_DIR}
    
    # Convert GTF to BED
    awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $4-1, $5, $10, ".", $7}' ${TE_GTF} | sed 's/[\";\]//g' > te_annotation.bed
    
    # Build index
    scTEATAC_build -g te_annotation.bed -o hg38_te
    
    echo "✓ Index built: ${INDEX_DIR}/hg38_te.idx"
fi

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_LOG="${LOG_DIR}/${SAMPLE}_scTEATAC_pipeline.log"
    
    echo "=========================================="
    echo "Processing: ${SAMPLE}"
    echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "=========================================="
    
    {
        # Check prerequisites
        CELLRANGER_BAM="${CELLRANGER_DIR}/${SAMPLE}/outs/possorted_bam.bam"
        FASTQ_R1="${FASTQ_DIR}/${SAMPLE}_1.fastq.gz"
        FASTQ_R2="${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"
        
        if [ ! -f "${CELLRANGER_BAM}" ]; then
            echo "⏸️  Skipping ${SAMPLE}: CellRanger BAM not found"
            continue
        fi
        
        if [ ! -f "${FASTQ_R1}" ] || [ ! -f "${FASTQ_R2}" ]; then
            echo "⏸️  Skipping ${SAMPLE}: FASTQ files not found"
            continue
        fi
        
        SAMPLE_STAR_DIR="${STAR_DIR}/${SAMPLE}"
        SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE}_scTEATAC"
        
        mkdir -p ${SAMPLE_STAR_DIR} ${SAMPLE_OUTPUT_DIR}
        
        # Step 1: STAR alignment with multi-mapper support
        STAR_BAM="${SAMPLE_STAR_DIR}/Aligned.sortedByCoord.out.bam"
        if [ ! -f "${STAR_BAM}" ]; then
            echo "[${SAMPLE}] Step 1/4: Running STAR multi-mapper alignment..."
            cd ${SAMPLE_STAR_DIR}
            
            STAR --runThreadN ${THREADS} \
                --genomeDir ${STAR_INDEX} \
                --readFilesIn ${FASTQ_R1} ${FASTQ_R2} \
                --readFilesCommand zcat \
                --outFilterMultimapNmax 100 \
                --winAnchorMultimapNmax 100 \
                --outSAMmultNmax 1 \
                --outFilterMismatchNmax 3 \
                --outSAMtype BAM SortedByCoordinate \
                --outFileNamePrefix ${SAMPLE_STAR_DIR}/
            
            echo "✓ STAR alignment complete"
        else
            echo "[${SAMPLE}] Step 1/4: STAR BAM exists, skipping alignment"
        fi
        
        # Step 2: Transfer CB tags from CellRanger to STAR BAM
        FINAL_BAM="${SAMPLE_STAR_DIR}/star_multimapper_final.bam"
        if [ ! -f "${FINAL_BAM}" ]; then
            echo "[${SAMPLE}] Step 2/4: Transferring CB tags..."
            
            # Extract ReadID → CB mapping
            echo "  Extracting CB tags from CellRanger BAM..."
            samtools view ${CELLRANGER_BAM} | \
                awk '{for(i=12;i<=NF;i++){if($i~/^CB:Z:/){split($i,a,":"); print $1"\t"a[3]; break}}}' | \
                gzip > ${SAMPLE_STAR_DIR}/readid_to_cb.txt.gz
            
            # Add CB tags to STAR BAM
            echo "  Adding CB tags to STAR multi-mapper BAM..."
            python3 ~/sc/scripts/add_cb_tags.py \
                ${SAMPLE_STAR_DIR}/readid_to_cb.txt.gz \
                ${STAR_BAM} \
                ${SAMPLE_STAR_DIR}/star_with_cb.bam
            
            # Filter to only reads with CB tags
            echo "  Filtering for reads with CB tags..."
            samtools view -h -d CB ${SAMPLE_STAR_DIR}/star_with_cb.bam | \
                samtools sort -@ ${THREADS} -o ${FINAL_BAM}
            
            samtools index ${FINAL_BAM}
            
            # Cleanup intermediate files
            rm -f ${SAMPLE_STAR_DIR}/star_with_cb.bam
            rm -f ${SAMPLE_STAR_DIR}/readid_to_cb.txt.gz
            
            echo "✓ CB transfer complete"
        else
            echo "[${SAMPLE}] Step 2/4: Final BAM with CB tags exists, skipping"
        fi
        
        # Step 3: Run scTEATAC
        OUTPUT_H5AD="${SAMPLE_OUTPUT_DIR}/${SAMPLE}_scTEATAC.h5ad"
        if [ ! -f "${OUTPUT_H5AD}" ]; then
            echo "[${SAMPLE}] Step 3/4: Running scTEATAC TE quantification..."
            cd ${SAMPLE_OUTPUT_DIR}
            
            scTEATAC -i ${FINAL_BAM} \
                -o ${SAMPLE}_scTEATAC \
                -x ${INDEX_DIR}/hg38_te.idx \
                -p ${THREADS} \
                -CB True \
                --hdf5 True
            
            echo "✓ scTEATAC quantification complete"
        else
            echo "[${SAMPLE}] Step 3/4: scTEATAC output exists, skipping"
        fi
        
        # Step 4: Verify output
        if [ -f "${OUTPUT_H5AD}" ]; then
            FILE_SIZE=$(ls -lh ${OUTPUT_H5AD} | awk '{print $5}')
            echo "[${SAMPLE}] Step 4/4: ✓ Complete! Output: ${OUTPUT_H5AD} (${FILE_SIZE})"
        else
            echo "[${SAMPLE}] Step 4/4: ⚠️  Warning: Expected output file not found"
        fi
        
        echo "Finished ${SAMPLE}: $(date '+%Y-%m-%d %H:%M:%S')"
        echo ""
        
    } 2>&1 | tee -a ${SAMPLE_LOG}
    
done

echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "Summary:"
find ${OUTPUT_DIR} -name "*.h5ad" -type f | while read f; do
    SIZE=$(ls -lh "$f" | awk '{print $5}')
    echo "  ✓ $f ($SIZE)"
done
