#!/bin/bash

# Build STAR genome index on AWS
# This rebuilds the index if you don't want to transfer 27GB

set -e

echo "=========================================="
echo "Building STAR Genome Index"
echo "=========================================="
echo "This will take ~30 minutes"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Configuration
GENOME_DIR="$HOME/sc/star_index"
FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"

mkdir -p "$GENOME_DIR"
cd "$GENOME_DIR"

# Download genome FASTA
echo "Downloading genome FASTA..."
wget -q --show-progress -O GRCh38.fa.gz "$FASTA_URL"
gunzip GRCh38.fa.gz

# Download gene annotations
echo "Downloading GTF annotations..."
wget -q --show-progress -O gencode.v44.gtf.gz "$GTF_URL"
gunzip gencode.v44.gtf.gz

# Build STAR index
echo "Building STAR index (using all 36 cores)..."
STAR \
    --runMode genomeGenerate \
    --runThreadN 36 \
    --genomeDir "$GENOME_DIR" \
    --genomeFastaFiles GRCh38.fa \
    --sjdbGTFfile gencode.v44.gtf \
    --sjdbOverhang 99 \
    --genomeSAindexNbases 14

# Cleanup
echo "Cleaning up..."
rm -f GRCh38.fa gencode.v44.gtf

echo ""
echo "=========================================="
echo "Index Build Complete"
echo "=========================================="
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Index location: $GENOME_DIR"
echo "Index size: $(du -sh $GENOME_DIR | cut -f1)"
echo "=========================================="
