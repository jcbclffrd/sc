#!/bin/bash

# AWS EC2 Setup Script for ATAC-seq STAR Alignment
# Run this on a fresh c5.9xlarge Ubuntu 22.04 instance

set -e

echo "=========================================="
echo "AWS ATAC-seq Alignment Setup"
echo "Instance: c5.9xlarge (36 cores, 72GB RAM)"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================="

# Update system
echo "Updating system packages..."
sudo apt-get update -y
sudo apt-get upgrade -y

# Install required packages
echo "Installing dependencies..."
sudo apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    pigz \
    samtools \
    htop \
    iftop \
    unzip \
    python3 \
    python3-pip

# Install SRA Toolkit
echo "Installing SRA Toolkit..."
cd ~
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
rm sratoolkit.current-ubuntu64.tar.gz
echo 'export PATH=$HOME/sratoolkit.3.1.1-ubuntu64/bin:$PATH' >> ~/.bashrc
export PATH=$HOME/sratoolkit.3.1.1-ubuntu64/bin:$PATH

# Configure SRA toolkit to not prefetch (stream directly)
vdb-config --prefetch-to-user-repo off

# Install STAR aligner
echo "Installing STAR aligner..."
cd ~
wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz
tar -xzf 2.7.11b.tar.gz
cd STAR-2.7.11b/source
make STAR
sudo cp STAR /usr/local/bin/
cd ~
rm -rf STAR-2.7.11b 2.7.11b.tar.gz

# Verify installations
echo "Verifying installations..."
which fasterq-dump
which STAR
which samtools
which pigz

fasterq-dump --version
STAR --version
samtools --version
pigz --version

# Clone the sc repository
echo "Cloning sc repository..."
cd ~
if [[ -d "sc" ]]; then
    echo "sc directory already exists, pulling latest..."
    cd sc
    git pull
else
    git clone https://github.com/jcbclffrd/sc.git
    cd sc
fi

# Create necessary directories
echo "Creating directories..."
mkdir -p ~/sc/sra_downloads/ATAC-seq
mkdir -p ~/sc/atacseq_aligned
mkdir -p ~/sc/logs_atacseq
mkdir -p ~/sc/star_index

# Configure git
git config --global user.email "jcbclffrd@users.noreply.github.com"
git config --global user.name "Jacob Clifford"

# Check available resources
echo ""
echo "=========================================="
echo "System Resources:"
echo "=========================================="
echo "CPUs: $(nproc)"
echo "RAM: $(free -h | grep Mem | awk '{print $2}')"
echo "Disk: $(df -h / | tail -1 | awk '{print $4}') available"
echo ""

echo "=========================================="
echo "Setup Complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Upload genome index: rsync -avz -e 'ssh -i ~/.ssh/awsWebsite.pem' ~/sc/star_index/ ubuntu@\$AWS_IP:~/sc/star_index/"
echo "2. Or rebuild genome index: cd ~/sc && ./scripts/build_star_index_aws.sh"
echo "3. Run alignment: cd ~/sc && nohup bash scripts/align_atacseq_aws.sh > logs_atacseq/aws_alignment.log 2>&1 &"
echo "4. Monitor: tail -f logs_atacseq/aws_alignment.log"
echo ""
