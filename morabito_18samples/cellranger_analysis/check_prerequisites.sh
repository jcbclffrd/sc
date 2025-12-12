#!/bin/bash
#
# Check prerequisites for CellRanger Docker analysis

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
HOST_FASTQ_DIR="/home/jacobc/sc/sra_downloads"
SAMPLE_MAP="../sample_mapping.csv"

echo "======================================================================"
echo "CELLRANGER PREREQUISITES CHECK"
echo "======================================================================"
echo ""

# Check 1: Docker
echo "1. Checking Docker..."
if command -v docker &> /dev/null; then
    DOCKER_VERSION=$(docker --version)
    echo "   ✓ Docker installed: $DOCKER_VERSION"
    
    # Check if Docker is running
    if docker info &> /dev/null; then
        echo "   ✓ Docker daemon is running"
    else
        echo "   ✗ Docker daemon is not running"
        echo "     Start with: sudo systemctl start docker"
        exit 1
    fi
else
    echo "   ✗ Docker not found"
    echo "     Install: https://docs.docker.com/engine/install/"
    exit 1
fi
echo ""

# Check 2: Platform support
echo "2. Checking platform support..."
ARCH=$(uname -m)
echo "   Host architecture: $ARCH"
if [ "$ARCH" = "aarch64" ] || [ "$ARCH" = "arm64" ]; then
    echo "   ⚠️  ARM architecture detected"
    echo "     CellRanger requires x86_64 - will use Docker emulation (slower)"
else
    echo "   ✓ x86_64 architecture - native performance"
fi
echo ""

# Check 3: FASTQ files
echo "3. Checking FASTQ files..."
if [ ! -d "$HOST_FASTQ_DIR" ]; then
    echo "   ✗ FASTQ directory not found: $HOST_FASTQ_DIR"
    exit 1
fi

N_SAMPLES=$(find "$HOST_FASTQ_DIR" -maxdepth 1 -type d -name "SRR*" | wc -l)
echo "   Found: $N_SAMPLES sample directories in $HOST_FASTQ_DIR"

if [ $N_SAMPLES -ge 18 ]; then
    echo "   ✓ All 18 samples available"
else
    echo "   ⚠️  Expected 18 samples, found $N_SAMPLES"
fi

# Check sample structure
echo "   Checking sample structure..."
SAMPLE_CHECK=$(find "$HOST_FASTQ_DIR/SRR14513984" -name "*_R1_*.fastq.gz" -o -name "*_R2_*.fastq.gz" 2>/dev/null | wc -l)
if [ $SAMPLE_CHECK -ge 2 ]; then
    echo "   ✓ FASTQ files have correct naming (*_R1_*.fastq.gz, *_R2_*.fastq.gz)"
else
    echo "   ⚠️  FASTQ naming may need adjustment"
    echo "     CellRanger expects: *_S1_L001_R1_001.fastq.gz format"
fi
echo ""

# Check 4: Sample mapping
echo "4. Checking sample mapping..."
if [ -f "$SAMPLE_MAP" ]; then
    N_MAPPED=$(tail -n +2 "$SAMPLE_MAP" | wc -l)
    echo "   ✓ Sample mapping found: $N_MAPPED samples"
else
    echo "   ✗ Sample mapping not found: $SAMPLE_MAP"
    exit 1
fi
echo ""

# Check 5: Disk space
echo "5. Checking disk space..."
AVAIL_GB=$(df -BG "$SCRIPT_DIR" | tail -1 | awk '{print $4}' | sed 's/G//')
echo "   Available space: ${AVAIL_GB}GB"

if [ $AVAIL_GB -lt 20 ]; then
    echo "   ⚠️  Less than 20GB available"
    echo "     Recommended: 20GB for Docker image + 10GB for outputs"
else
    echo "   ✓ Sufficient disk space"
fi
echo ""

# Check 6: Memory
echo "6. Checking system memory..."
TOTAL_MEM_GB=$(free -g | grep "Mem:" | awk '{print $2}')
echo "   Total RAM: ${TOTAL_MEM_GB}GB"

if [ $TOTAL_MEM_GB -lt 64 ]; then
    echo "   ⚠️  Less than 64GB RAM"
    echo "     CellRanger may run slower or need to process samples sequentially"
else
    echo "   ✓ Sufficient memory for parallel processing"
fi
echo ""

# Check 7: Docker image
echo "7. Checking Docker image..."
if docker images | grep -q "cellranger-morabito"; then
    IMAGE_SIZE=$(docker images cellranger-morabito:latest --format "{{.Size}}")
    echo "   ✓ Docker image already built: $IMAGE_SIZE"
else
    echo "   ✗ Docker image not built yet"
    echo "     Run: bash build_docker.sh"
fi
echo ""

echo "======================================================================"
echo "SUMMARY"
echo "======================================================================"
echo ""

if [ $N_SAMPLES -ge 18 ] && [ $AVAIL_GB -ge 20 ]; then
    echo "✓ Ready to proceed!"
    echo ""
    echo "Next steps:"
    if ! docker images | grep -q "cellranger-morabito"; then
        echo "  1. Build Docker image: bash build_docker.sh"
        echo "  2. Run container: bash run_docker.sh"
    else
        echo "  1. Run container: bash run_docker.sh"
        echo "  2. Inside container: bash run_cellranger_batch.sh"
    fi
else
    echo "⚠️  Some issues detected - review messages above"
fi
echo ""
