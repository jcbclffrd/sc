#!/bin/bash
#
# Download CellRanger 8.0.1 from 10x Genomics
# This requires accepting their End User License Agreement

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CELLRANGER_VERSION="8.0.1"
CELLRANGER_FILE="cellranger-${CELLRANGER_VERSION}.tar.gz"
DOWNLOAD_URL="https://cf.10xgenomics.com/releases/cell-exp/${CELLRANGER_FILE}"

echo "======================================================================"
echo "DOWNLOADING CELLRANGER ${CELLRANGER_VERSION}"
echo "======================================================================"
echo ""

# Check if already downloaded
if [ -f "$SCRIPT_DIR/$CELLRANGER_FILE" ]; then
    FILE_SIZE=$(du -h "$SCRIPT_DIR/$CELLRANGER_FILE" | cut -f1)
    echo "✓ CellRanger already downloaded: $CELLRANGER_FILE ($FILE_SIZE)"
    echo ""
    exit 0
fi

echo "CellRanger is proprietary software from 10x Genomics."
echo "By downloading, you agree to their End User License Agreement:"
echo "https://www.10xgenomics.com/legal/end-user-software-license-agreement"
echo ""
echo "Download URL: $DOWNLOAD_URL"
echo "Expected size: ~2.0 GB"
echo ""

read -p "Continue with download? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Download cancelled."
    exit 1
fi

echo ""
echo "Downloading CellRanger ${CELLRANGER_VERSION}..."
echo "This may take 5-10 minutes depending on your connection."
echo ""

# Try wget first
if command -v wget &> /dev/null; then
    wget -O "$SCRIPT_DIR/$CELLRANGER_FILE" "$DOWNLOAD_URL"
    DOWNLOAD_EXIT=$?
# Fallback to curl
elif command -v curl &> /dev/null; then
    curl -L -o "$SCRIPT_DIR/$CELLRANGER_FILE" "$DOWNLOAD_URL"
    DOWNLOAD_EXIT=$?
else
    echo "✗ ERROR: Neither wget nor curl found"
    echo "Please install wget or curl to download CellRanger"
    exit 1
fi

if [ $DOWNLOAD_EXIT -eq 0 ]; then
    FILE_SIZE=$(du -h "$SCRIPT_DIR/$CELLRANGER_FILE" | cut -f1)
    echo ""
    echo "======================================================================"
    echo "✓ DOWNLOAD SUCCESSFUL"
    echo "======================================================================"
    echo ""
    echo "File: $CELLRANGER_FILE"
    echo "Size: $FILE_SIZE"
    echo "Location: $SCRIPT_DIR/"
    echo ""
    
    # Verify it's a valid tar.gz file
    if tar -tzf "$SCRIPT_DIR/$CELLRANGER_FILE" >/dev/null 2>&1; then
        echo "✓ File integrity verified (valid tar.gz)"
    else
        echo "⚠️  WARNING: File may be corrupted (not a valid tar.gz)"
        echo "This could be an HTML error page instead of the actual file."
        echo ""
        echo "If download failed with 403 Forbidden, you need to:"
        echo "  1. Visit: https://www.10xgenomics.com/support/software/cell-ranger/downloads"
        echo "  2. Sign in or create account"
        echo "  3. Accept EULA"
        echo "  4. Download cellranger-8.0.1.tar.gz manually"
        echo "  5. Move it to: $SCRIPT_DIR/"
        rm -f "$SCRIPT_DIR/$CELLRANGER_FILE"
        exit 1
    fi
    
    echo ""
    echo "Next step: Build Docker image"
    echo "  bash build_docker.sh"
    echo ""
else
    echo ""
    echo "======================================================================"
    echo "✗ DOWNLOAD FAILED"
    echo "======================================================================"
    echo ""
    echo "The download requires authentication with 10x Genomics."
    echo ""
    echo "Manual download instructions:"
    echo "  1. Visit: https://www.10xgenomics.com/support/software/cell-ranger/downloads"
    echo "  2. Sign in (or create a free account)"
    echo "  3. Accept the End User License Agreement"
    echo "  4. Download: Cell Ranger 8.0.1 for Linux (tar.gz)"
    echo "  5. Move the file to: $SCRIPT_DIR/"
    echo ""
    echo "Or use wget with authentication:"
    echo "  wget --user=YOUR_EMAIL --ask-password \\"
    echo "    -O $CELLRANGER_FILE \\"
    echo "    $DOWNLOAD_URL"
    echo ""
    exit 1
fi
