#!/bin/bash
set -e

SAMPLE=SRR14514130
#BAM=~/sc/atacseq_aligned/${SAMPLE}/star_multimapper_final.bam
BAM=~/sc/cellranger_atac_output/${SAMPLE}/outs/possorted_bam.bam
TE_GTF=~/sc/annotations/hg38_rmsk_TE.gtf
OUTPUT_DIR=~/sc/atacseq_te_counts/${SAMPLE}_scTEATAC_cellranger
THREADS=16

export PATH=~/.local/bin:$PATH

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

echo "[1] Converting GTF to BED format..."
awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $4-1, $5, $10, ".", $7}' ${TE_GTF} | sed 's/[\";\]//g' > te_annotation.bed

echo "[2] Building scTEATAC index..."
scTEATAC_build -g te_annotation.bed -o hg38_te

echo "[3] Running scTEATAC on BAM with CB tags..."
# Redirect stderr to /dev/null to suppress bamToBed paired-read warnings (they fill up logs)
scTEATAC -i ${BAM} -o ${SAMPLE}_scTEATAC -x hg38_te.idx -p ${THREADS} -CB True --hdf5 True 2>/dev/null

echo "Done! Output: ${OUTPUT_DIR}/${SAMPLE}_scTEATAC.h5ad"
ls -lh ${SAMPLE}_scTEATAC*
