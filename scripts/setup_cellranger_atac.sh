#!/bin/bash

# Download and install CellRanger ATAC
# Run this on the R-instance or wherever you'll do alignment

set -e

INSTALL_DIR="$HOME/software"
CELLRANGER_VERSION="2.2.0"

mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

echo "=========================================="
echo "Installing CellRanger ATAC"
echo "=========================================="
echo "Version: $CELLRANGER_VERSION"
echo "Install directory: $INSTALL_DIR"
echo ""

# Download CellRanger ATAC
if [ ! -f "cellranger-atac-${CELLRANGER_VERSION}.tar.gz" ]; then
    echo "Downloading CellRanger ATAC ${CELLRANGER_VERSION}..."
    curl -o cellranger-atac-${CELLRANGER_VERSION}.tar.gz "https://cf.10xgenomics.com/releases/cell-atac/cellranger-atac-2.2.0.tar.gz?Expires=1767289327&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=nE4eMEQZuTh2OwH3NfxfQLHBzFrUaYdWHSodq227~to21VHqVTkzOTbVUpgKL3W6mkOBaMcrLBaL~G0LqKbbe9~DAsUAbqgm8MQKYQUkHaiOEd-vkMWQb2If0JpiIqnXDqTJxMABMIjRPdJrfCd9TRuiz77WiCxbxpDRu4JJNnYDIf04HMYGPA9PHJNs9nKMO3sm1P5RLq9HcJ4qe5sIsZTqMsMxiFByg6HtTJ8OKf7lspQ5DspVg3MqTNtPWG-gjAET~irr9Ibij6ssyhZtY21Mp9af2J4PBj9htlPlrhfDsd0H~rxx547rgDzNWozTWacRGYr~AW0ZTyAqiMlYcw__"
else
    echo "CellRanger ATAC tarball already exists, skipping download"
fi

# Extract
if [ ! -d "cellranger-atac-${CELLRANGER_VERSION}" ]; then
    echo "Extracting..."
    tar -xzvf cellranger-atac-${CELLRANGER_VERSION}.tar.gz
else
    echo "CellRanger ATAC already extracted"
fi

# Add to PATH
echo ""
echo "=========================================="
echo "Installation Complete"
echo "=========================================="
echo ""
echo "Add this to your ~/.bashrc:"
echo "export PATH=$INSTALL_DIR/cellranger-atac-${CELLRANGER_VERSION}:\$PATH"
echo ""
echo "Or run this now:"
echo "export PATH=$INSTALL_DIR/cellranger-atac-${CELLRANGER_VERSION}:\$PATH"
echo ""
echo "Verify installation:"
echo "cellranger-atac --version"
echo ""

# Download reference genome (hg38)
echo "=========================================="
echo "Downloading Reference Genome (hg38)"
echo "=========================================="
echo ""

REF_DIR="$HOME/sc/cellranger_references"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

REF_VERSION="refdata-cellranger-arc-GRCh38-2024-A"

if [ ! -f "${REF_VERSION}.tar.gz" ]; then
    echo "Downloading hg38 reference (this will take a while, ~14GB)..."
    curl -O "https://cf.10xgenomics.com/supp/cell-arc/${REF_VERSION}.tar.gz"
else
    echo "Reference genome tarball already exists"
fi

if [ ! -d "$REF_VERSION" ]; then
    echo "Extracting reference genome..."
    tar -xzvf ${REF_VERSION}.tar.gz
else
    echo "Reference genome already extracted"
fi

echo ""
echo "=========================================="
echo "Setup Complete!"
echo "=========================================="
echo ""
echo "Reference genome: $REF_DIR/$REF_VERSION"
echo ""
