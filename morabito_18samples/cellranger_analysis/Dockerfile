FROM ubuntu:22.04

# CellRanger for Morabito 18-sample analysis
# This container runs CellRanger on the same 18 samples for comparison with STAR+scTE
# MUST use amd64 platform since CellRanger doesn't support ARM
# Platform is specified in build command: docker buildx build --platform linux/amd64

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget \
    curl \
    build-essential \
    python3 \
    python3-pip \
    python3-venv \
    pigz \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Create virtual environment
RUN python3 -m venv .venv
ENV PATH="/app/.venv/bin:$PATH"

# Install Python dependencies
RUN pip install --upgrade pip && \
    pip install pandas numpy scanpy h5py anndata scipy

# Copy and install CellRanger (x86_64 only)
# CellRanger must be downloaded manually from:
# https://www.10xgenomics.com/support/software/cell-ranger/downloads
# Accepts any version: cellranger-8.0.1.tar.gz, cellranger-10.0.0.tar.gz, etc.
COPY cellranger-*.tar.gz /tmp/cellranger.tar.gz
RUN tar -xzf /tmp/cellranger.tar.gz -C /app && \
    rm /tmp/cellranger.tar.gz && \
    ln -s /app/cellranger-* /app/cellranger

ENV PATH="/app/cellranger:$PATH"

# Create directory structure
RUN mkdir -p cellranger_output

# Download human reference (GRCh38-2024-A, ~11GB)
# This matches the reference we used for STAR alignment
RUN wget -O refdata-gex-GRCh38-2024-A.tar.gz \
    "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz" && \
    tar -xzf refdata-gex-GRCh38-2024-A.tar.gz && \
    rm refdata-gex-GRCh38-2024-A.tar.gz

ENV CELLRANGER_REFERENCE="/app/refdata-gex-GRCh38-2024-A"

# Copy processing scripts
COPY run_cellranger_batch.sh .
COPY merge_cellranger_outputs.py .
COPY compare_with_star.py .
COPY README.md .

# Make scripts executable
RUN chmod +x run_cellranger_batch.sh

# Volume mounts for data persistence
# Mount your FASTQ data from /home/jacobc/sc/data/
# Mount output to /home/jacobc/sc/morabito_18samples/cellranger_analysis/cellranger_output/
VOLUME ["/app/fastq_data", "/app/cellranger_output", "/app/sample_mapping"]

CMD ["/bin/bash"]
