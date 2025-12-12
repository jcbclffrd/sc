#!/bin/bash
#
# Build CellRanger Docker container for x86_64 (amd64) platform
# Required because CellRanger does not support ARM architecture

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
IMAGE_NAME="cellranger-morabito"
IMAGE_TAG="8.0.1"

echo "======================================================================"
echo "BUILDING CELLRANGER DOCKER IMAGE"
echo "======================================================================"
echo ""
echo "Image: ${IMAGE_NAME}:${IMAGE_TAG}"
echo "Platform: linux/amd64 (x86_64 only)"
echo "Build context: $SCRIPT_DIR"
echo ""
echo "⚠️  WARNING: This is a LARGE build (~15GB)"
echo "   - CellRanger 8.0.1: ~2GB"
echo "   - GRCh38 reference: ~11GB"
echo "   - Base system: ~2GB"
echo ""
echo "Estimated build time: 20-30 minutes"
echo ""

read -p "Continue with build? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Build cancelled."
    exit 1
fi

echo ""
echo "Starting build..."
echo ""

# Check if multiplatform builder exists, create if not
if ! docker buildx ls | grep -q "multiplatform"; then
    echo "Creating multiplatform builder..."
    docker buildx create --name multiplatform --driver docker-container --bootstrap --use || true
fi

# Use multiplatform builder
docker buildx use multiplatform || docker buildx use default

# Build with amd64 platform using buildx
docker buildx build \
    --platform linux/amd64 \
    --load \
    -t ${IMAGE_NAME}:${IMAGE_TAG} \
    -t ${IMAGE_NAME}:latest \
    "$SCRIPT_DIR"

BUILD_EXIT=$?

if [ $BUILD_EXIT -eq 0 ]; then
    echo ""
    echo "======================================================================"
    echo "✓ BUILD SUCCESSFUL"
    echo "======================================================================"
    echo ""
    echo "Image: ${IMAGE_NAME}:${IMAGE_TAG}"
    
    # Show image size
    IMAGE_SIZE=$(docker images ${IMAGE_NAME}:${IMAGE_TAG} --format "{{.Size}}")
    echo "Size: $IMAGE_SIZE"
    echo ""
    
    echo "Next steps:"
    echo "  1. Run container: bash run_docker.sh"
    echo "  2. Inside container: bash run_cellranger_batch.sh"
    echo ""
else
    echo ""
    echo "======================================================================"
    echo "✗ BUILD FAILED"
    echo "======================================================================"
    echo ""
    echo "Check errors above. Common issues:"
    echo "  - Docker not running"
    echo "  - No internet connection (needs to download CellRanger + reference)"
    echo "  - Insufficient disk space (needs ~20GB)"
    echo ""
    exit 1
fi
