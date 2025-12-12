#!/bin/bash
#
# Check if CellRanger download is complete and ready to build

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CELLRANGER_FILE="$SCRIPT_DIR/cellranger-10.0.0.tar.gz"
EXPECTED_SIZE_GB=2

echo "======================================================================"
echo "CELLRANGER DOWNLOAD CHECK"
echo "======================================================================"
echo ""

if [ ! -f "$CELLRANGER_FILE" ]; then
    echo "✗ CellRanger file not found: $CELLRANGER_FILE"
    echo ""
    echo "Expected: cellranger-10.0.0.tar.gz (~2GB)"
    exit 1
fi

CURRENT_SIZE=$(du -h "$CELLRANGER_FILE" | cut -f1)
CURRENT_SIZE_BYTES=$(stat -f%z "$CELLRANGER_FILE" 2>/dev/null || stat -c%s "$CELLRANGER_FILE")
CURRENT_SIZE_GB=$(echo "scale=2; $CURRENT_SIZE_BYTES / 1024 / 1024 / 1024" | bc)

echo "File: $(basename $CELLRANGER_FILE)"
echo "Size: $CURRENT_SIZE ($CURRENT_SIZE_GB GB)"
echo ""

# Check if file is still being written (size changing)
INITIAL_SIZE=$CURRENT_SIZE_BYTES
sleep 2
NEW_SIZE=$(stat -f%z "$CELLRANGER_FILE" 2>/dev/null || stat -c%s "$CELLRANGER_FILE")

if [ "$INITIAL_SIZE" != "$NEW_SIZE" ]; then
    echo "⏳ Download still in progress (file size changing)"
    PERCENT=$(echo "scale=0; $CURRENT_SIZE_GB * 100 / $EXPECTED_SIZE_GB" | bc)
    echo "   Progress: ~${PERCENT}%"
    echo ""
    echo "Wait for download to complete, then run:"
    echo "  bash build_docker.sh"
    exit 1
fi

# Check if file is complete (>1.8GB)
if (( $(echo "$CURRENT_SIZE_GB < 1.8" | bc -l) )); then
    echo "⚠️  File seems incomplete (expected ~2GB)"
    echo "   Current: $CURRENT_SIZE_GB GB"
    echo "   Expected: ~2.0 GB"
    echo ""
    echo "If download is still running, wait for it to complete."
    exit 1
fi

# Verify it's a valid tar.gz
echo "Verifying file integrity..."
if tar -tzf "$CELLRANGER_FILE" >/dev/null 2>&1; then
    echo "✓ Valid tar.gz archive"
else
    echo "✗ File is corrupted (not a valid tar.gz)"
    echo ""
    echo "Try downloading again."
    exit 1
fi

# Check what's inside
CELLRANGER_DIR=$(tar -tzf "$CELLRANGER_FILE" | head -1 | cut -d/ -f1)
echo "✓ Contains: $CELLRANGER_DIR"

echo ""
echo "======================================================================"
echo "✓ READY TO BUILD"
echo "======================================================================"
echo ""
echo "CellRanger 10.0.0 downloaded successfully!"
echo ""
echo "Next step: Build Docker image"
echo "  cd $SCRIPT_DIR"
echo "  bash build_docker.sh"
echo ""
echo "Estimated build time: 20-30 minutes"
echo "Image size: ~15GB"
echo ""
