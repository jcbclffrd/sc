#!/bin/bash
# Build and run the complete ATAC-seq Docker pipeline

set -e

echo "=========================================="
echo "ATAC-seq Docker Pipeline"
echo "=========================================="

# Configuration
SAMPLES=(
    SRR14514129 SRR14514130 SRR14514131 SRR14514132 SRR14514133
    SRR14514134 SRR14514135 SRR14514136 SRR14514137 SRR14514138
    SRR14514139 SRR14514140 SRR14514141 SRR14514142 SRR14514143
    SRR14514144 SRR14514145 SRR14514146 SRR14514147 SRR14514148
    SRR14514149 SRR14514150 SRR14514151 SRR14514152 SRR14514153
    SRR14514154 SRR14514155 SRR14514156 SRR14514157 SRR14514158
    SRR14514159 SRR14514160
)

MAX_PARALLEL=9  # Adjust based on CPU cores available

#============================================
# Step 0: Build Docker images
#============================================

echo "Building Docker images..."
docker build -t atacseq-base:latest -f docker/Dockerfile.base .
docker build -t atacseq-alignment:latest -f docker/Dockerfile.alignment .
docker build -t atacseq-peaks:latest -f docker/Dockerfile.peaks .
docker build -t atacseq-quantify:latest -f docker/Dockerfile.quantify .
docker build -t atacseq-deseq2:latest -f docker/Dockerfile.deseq2 .

echo "✓ Docker images built"
echo ""

#============================================
# Step 1: Alignment (parallel)
#============================================

echo "=========================================="
echo "Step 1: STAR Alignment"
echo "=========================================="

align_sample() {
    local sample=$1
    echo "[$sample] Starting alignment..."
    
    docker run --rm \
        -v "$(pwd)/star_index:/genome:ro" \
        -v "$(pwd)/atacseq_aligned:/output" \
        -v "$(pwd)/sra_downloads/ATAC-seq:/data" \
        -e SAMPLE_ID="$sample" \
        -e THREADS=4 \
        --cpus=4 \
        --memory=16g \
        atacseq-alignment:latest
    
    echo "[$sample] ✓ Complete"
}

export -f align_sample

# Run alignments in parallel
printf '%s\n' "${SAMPLES[@]}" | xargs -P "$MAX_PARALLEL" -I {} bash -c 'align_sample "$@"' _ {}

echo "✓ All alignments complete"
echo ""

#============================================
# Step 2: Peak calling (parallel)
#============================================

echo "=========================================="
echo "Step 2: Peak Calling"
echo "=========================================="

call_peaks() {
    local sample=$1
    echo "[$sample] Calling peaks..."
    
    docker run --rm \
        -v "$(pwd)/atacseq_aligned:/data:ro" \
        -v "$(pwd)/atacseq_peaks:/output" \
        -e SAMPLE_ID="$sample" \
        --cpus=2 \
        --memory=4g \
        atacseq-peaks:latest
    
    echo "[$sample] ✓ Complete"
}

export -f call_peaks

# Run peak calling in parallel (lighter, can run more)
printf '%s\n' "${SAMPLES[@]}" | xargs -P 16 -I {} bash -c 'call_peaks "$@"' _ {}

echo "✓ All peaks called"
echo ""

#============================================
# Step 3: TE quantification
#============================================

echo "=========================================="
echo "Step 3: TE Quantification"
echo "=========================================="

docker run --rm \
    -v "$(pwd)/atacseq_aligned:/data:ro" \
    -v "$(pwd)/annotations:/annotations:ro" \
    -v "$(pwd)/atacseq_te_counts:/output" \
    -e TE_GTF=/annotations/hg38_rmsk_TE.gtf \
    -e THREADS=8 \
    --cpus=8 \
    --memory=16g \
    atacseq-quantify:latest

echo "✓ TE quantification complete"
echo ""

#============================================
# Step 4: DESeq2 analysis
#============================================

echo "=========================================="
echo "Step 4: DESeq2 Analysis"
echo "=========================================="

docker run --rm \
    -v "$(pwd)/atacseq_te_counts:/data:ro" \
    -v "$(pwd)/metadata:/metadata:ro" \
    -v "$(pwd)/deseq2_results:/output" \
    -e COUNT_MATRIX=/data/te_counts_matrix.tsv \
    -e METADATA=/metadata/atacseq_metadata.csv \
    -e CONDITION_COL=condition \
    -e CONTROL=Control \
    -e TREATMENT=AD \
    --cpus=4 \
    --memory=16g \
    atacseq-deseq2:latest

echo "✓ DESeq2 analysis complete"
echo ""

#============================================
# Summary
#============================================

echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
echo "Results:"
echo "  Alignments: $(ls atacseq_aligned/SRR*/Aligned.sortedByCoord.out.bam 2>/dev/null | wc -l) samples"
echo "  Peaks: $(ls atacseq_peaks/SRR*/*_peaks.narrowPeak 2>/dev/null | wc -l) samples"
echo "  TE counts: atacseq_te_counts/te_counts_matrix.tsv"
echo "  DESeq2 results: deseq2_results/deseq2_results.csv"
echo "=========================================="
