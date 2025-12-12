#!/bin/bash
#
# Setup script for SoloTE analysis
# Prepares TE annotation and checks dependencies

set -e

echo "======================================================================"
echo "SOLOTE SETUP: Morabito 18 Samples"
echo "======================================================================"

# Check dependencies
echo -e "\n=== Checking Dependencies ==="

check_command() {
    if command -v $1 &> /dev/null; then
        echo "✓ $1 found: $(which $1)"
        if [ "$1" = "samtools" ]; then
            samtools --version | head -1
        elif [ "$1" = "bedtools" ]; then
            bedtools --version
        elif [ "$1" = "R" ]; then
            R --version | head -1
        elif [ "$1" = "python3" ]; then
            python3 --version
        fi
    else
        echo "✗ $1 NOT FOUND"
        return 1
    fi
}

# Check all dependencies
MISSING=""
check_command samtools || MISSING="$MISSING samtools"
check_command bedtools || MISSING="$MISSING bedtools"
check_command R || MISSING="$MISSING R"
check_command python3 || MISSING="$MISSING python3"

if [ ! -z "$MISSING" ]; then
    echo -e "\n✗ Missing dependencies:$MISSING"
    echo -e "\nInstall with:"
    echo "  sudo apt install bedtools  # for bedtools"
    echo "  # samtools and R should already be installed"
    exit 1
fi

# Check Python packages
echo -e "\n=== Checking Python Packages ==="
python3 -c "import pysam; print('✓ pysam installed')" 2>/dev/null || echo "✗ pysam not installed (pip install pysam)"
python3 -c "import pandas; print('✓ pandas installed:', pandas.__version__)" 2>/dev/null || echo "✗ pandas not installed (pip install pandas)"

# Download and prepare TE annotation
echo -e "\n=== Preparing TE Annotation ==="

cd SoloTE

if [ -f "hg38_rmsk.bed" ]; then
    echo "✓ hg38_rmsk.bed already exists"
    echo "  Size: $(du -h hg38_rmsk.bed | cut -f1)"
    echo "  Lines: $(wc -l < hg38_rmsk.bed)"
else
    echo "Downloading RepeatMasker annotation for hg38..."
    python3 SoloTE_RepeatMasker_to_BED.py -g hg38
    
    if [ -f "hg38_rmsk.bed" ]; then
        echo "✓ Downloaded hg38_rmsk.bed"
        echo "  Size: $(du -h hg38_rmsk.bed | cut -f1)"
        echo "  Lines: $(wc -l < hg38_rmsk.bed)"
    else
        echo "✗ Failed to create hg38_rmsk.bed"
        exit 1
    fi
fi

cd ..

# Create output directory
echo -e "\n=== Creating Output Directory ==="
mkdir -p soloTE_output
echo "✓ Created: soloTE_output/"

# Check sample mapping
echo -e "\n=== Checking Sample Files ==="

SAMPLE_MAP="../sample_mapping.csv"
if [ ! -f "$SAMPLE_MAP" ]; then
    echo "✗ Sample mapping not found: $SAMPLE_MAP"
    exit 1
fi

N_SAMPLES=$(tail -n +2 "$SAMPLE_MAP" | wc -l)
echo "✓ Sample mapping found: $N_SAMPLES samples"

# Check BAM files exist
echo -e "\n=== Checking BAM Files ==="

BAM_DIR="/home/jacobc/sc/starsolo_aligned"
MISSING_BAMS=0

while IFS=, read -r sample_id gsm srx srr; do
    if [ "$sample_id" = "sample_id" ]; then continue; fi  # Skip header
    
    BAM_FILE="$BAM_DIR/$srr/Aligned.sortedByCoord.out.bam"
    if [ -f "$BAM_FILE" ]; then
        SIZE=$(du -h "$BAM_FILE" | cut -f1)
        echo "  ✓ $srr: $SIZE"
    else
        echo "  ✗ $srr: NOT FOUND"
        MISSING_BAMS=$((MISSING_BAMS + 1))
    fi
done < "$SAMPLE_MAP"

if [ $MISSING_BAMS -gt 0 ]; then
    echo -e "\n✗ Missing $MISSING_BAMS BAM files"
    echo "Make sure all samples have been aligned with STAR"
    exit 1
fi

echo -e "\n======================================================================"
echo "✓ SETUP COMPLETE!"
echo "======================================================================"
echo ""
echo "Next steps:"
echo "  1. Run batch processing: bash run_soloTE_batch.sh"
echo "  2. Or run single sample: bash run_soloTE_single.sh SRR14513984"
echo ""
echo "Estimated time: 1-2 hours per sample (~18-36 hours total for 18 samples)"
echo ""
