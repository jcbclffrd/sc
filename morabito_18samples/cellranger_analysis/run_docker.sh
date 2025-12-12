#!/bin/bash
#
# Run CellRanger Docker container with proper volume mounts
# Maps host directories to container paths

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
IMAGE_NAME="cellranger-morabito:latest"

# Host paths (your actual data locations)
HOST_FASTQ_DIR="$SCRIPT_DIR/fastq_cellranger_format"
HOST_OUTPUT_DIR="$SCRIPT_DIR/cellranger_output"
HOST_SAMPLE_MAP_DIR="$SCRIPT_DIR/.."  # Points to morabito_18samples/

# Create output directory if it doesn't exist
mkdir -p "$HOST_OUTPUT_DIR"

echo "======================================================================"
echo "RUNNING CELLRANGER DOCKER CONTAINER"
echo "======================================================================"
echo ""
echo "Image: $IMAGE_NAME"
echo "Platform: linux/amd64 (emulated on ARM)"
echo ""
echo "Volume mounts:"
echo "  FASTQ data: $HOST_FASTQ_DIR"
echo "              → /app/fastq_data"
echo "  Output:     $HOST_OUTPUT_DIR"
echo "              → /app/cellranger_output"
echo "  Mapping:    $HOST_SAMPLE_MAP_DIR/sample_mapping.csv"
echo "              → /app/sample_mapping/sample_mapping.csv"
echo ""

# Check that FASTQ directory exists
if [ ! -d "$HOST_FASTQ_DIR" ]; then
    echo "✗ ERROR: FASTQ directory not found: $HOST_FASTQ_DIR"
    echo ""
    echo "Make sure you have downloaded the FASTQ files to /home/jacobc/sc/data/"
    exit 1
fi

# Check that sample mapping exists
if [ ! -f "$HOST_SAMPLE_MAP_DIR/sample_mapping.csv" ]; then
    echo "✗ ERROR: Sample mapping not found: $HOST_SAMPLE_MAP_DIR/sample_mapping.csv"
    exit 1
fi

# Count available samples
N_SAMPLES=$(find "$HOST_FASTQ_DIR" -maxdepth 1 -type d -name "SRR*" | wc -l)
echo "Found $N_SAMPLES FASTQ sample directories"
echo ""

if [ $N_SAMPLES -lt 18 ]; then
    echo "⚠️  WARNING: Expected 18 samples, found $N_SAMPLES"
    echo "Some samples may be missing from $HOST_FASTQ_DIR"
    echo ""
fi

echo "Starting container..."
echo ""
echo "⚠️  NOTE: Running x86_64 container on ARM may be SLOW (2-3x slower)"
echo "         Consider running on x86_64 machine for faster processing"
echo ""

# Run container interactively
docker run -it --rm \
    --platform linux/amd64 \
    -v "$HOST_FASTQ_DIR:/app/fastq_data:ro" \
    -v "$HOST_OUTPUT_DIR:/app/cellranger_output" \
    -v "$HOST_SAMPLE_MAP_DIR:/app/sample_mapping:ro" \
    -w /app \
    "$IMAGE_NAME" \
    /bin/bash

echo ""
echo "Container exited."
echo ""
