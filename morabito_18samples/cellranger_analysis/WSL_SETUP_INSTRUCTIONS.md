# CellRanger Analysis Setup for WSL (Windows x86_64)

## Overview
This guide is for running CellRanger on a **Windows WSL Ubuntu machine with LIMITED resources** (16GB RAM, 50GB storage). CellRanger requires x86_64 architecture and won't work on ARM. This setup uses **network-mounted FASTQ files** to save disk space.

## System Requirements (WSL Machine)
- **Architecture**: x86_64 (confirmed by `uname -m`)
- **RAM**: 16GB (minimum, will process 1-2 samples at a time)
- **Storage**: ~50GB available
  - Docker image: ~21GB
  - CellRanger download: ~1GB
  - Output (all 18 samples): ~4GB
  - Workspace: ~2GB
  - **Total needed**: ~28GB (FASTQ files via network mount)
- **OS**: WSL2 with Ubuntu
- **Docker**: Docker Desktop for Windows with WSL2 backend

## Prerequisites

✅ **Already have**:
- WSL2 with Ubuntu
- VS Code with Remote-WSL extension
- Git

🔧 **Need to install**:
- Docker Desktop for Windows (WSL2 backend)
- sshfs (for mounting remote FASTQ files)

---

## IMPORTANT: Resource-Constrained Setup Strategy

Given **16GB RAM and 50GB storage**, we will:

1. ✅ **Mount FASTQ files via SSH** (saves ~60GB disk space)
2. ✅ **Process samples sequentially** (1 at a time, max 12GB RAM per job)
3. ✅ **Use `--nosecondary`** (skip clustering, saves RAM and disk)
4. ⚠️ **Cannot run parallel jobs** (would need 48GB+ RAM for 3 jobs)
5. ⚠️ **Total runtime**: 18-36 hours for all 18 samples (vs 6-12 hours parallel)

---

## Step-by-Step Setup Instructions

---

## Step-by-Step Setup Instructions

### PHASE 1: Initial Environment Setup (on WSL machine)

#### 1.1 Install Docker Desktop for Windows

**Download and Install**:
1. Download Docker Desktop from: https://www.docker.com/products/docker-desktop/
2. Install on Windows (NOT inside WSL)
3. During installation, ensure "Use WSL 2 instead of Hyper-V" is selected
4. After installation, open Docker Desktop
5. Go to Settings → Resources → WSL Integration
6. Enable integration with your Ubuntu distribution
7. Click "Apply & Restart"

**Verify in WSL**:
```bash
# Open WSL terminal
docker --version
# Should show: Docker version 24.x.x or higher

docker run hello-world
# Should download and run successfully
```

#### 1.2 Install SSH and SSHFS (for network-mounted FASTQ files)

```bash
# Update package list
sudo apt-get update

# Install sshfs for mounting remote directories
sudo apt-get install -y sshfs

# Verify installation
sshfs --version
```

#### 1.3 Check System Resources

```bash
# Check architecture (must be x86_64)
uname -m
# Expected: x86_64

# Check available RAM
free -h
# Should show ~16GB total

# Check available disk space
df -h ~
# Should show ~50GB available

# Check CPU cores
nproc
# Note this number for later configuration
```

---

### PHASE 2: Clone Repository and Setup Workspace

#### 2.1 Clone the Repository

```bash
# Navigate to home directory
cd ~

# Clone repository (replace with your actual repo URL)
git clone https://github.com/YOUR_USERNAME/YOUR_REPO.git sc

# Or if already cloned, pull latest changes
cd ~/sc
git pull origin main

# Navigate to cellranger analysis directory
cd ~/sc/morabito_18samples/cellranger_analysis
```

#### 2.2 Verify All Files Present

```bash
# Check for required scripts
ls -lh *.sh
# Should see: build_docker.sh, run_docker.sh, run_cellranger_batch.sh, etc.

# Check for Dockerfile
ls -lh Dockerfile

# Check for Python scripts
ls -lh *.py
# Should see: merge_cellranger_outputs.py, compare_with_star.py
```

---

### PHASE 3: Mount FASTQ Files from ARM Machine

#### 3.1 Create Mount Point

```bash
# Create directory for mounted FASTQ files
mkdir -p ~/remote_fastq

# Create directory for SRA downloads mount
mkdir -p ~/sc/sra_downloads
```

#### 3.2 Setup SSH Key Authentication (RECOMMENDED)

To avoid entering password repeatedly:

```bash
# Generate SSH key if you don't have one
ssh-keygen -t ed25519 -C "wsl-machine"
# Press Enter to accept defaults

# Copy public key to ARM machine
# Replace USER and ARM_MACHINE_IP with actual values
ssh-copy-id USER@ARM_MACHINE_IP

# Test passwordless login
ssh USER@ARM_MACHINE_IP "echo 'SSH working'"
```

#### 3.3 Mount the FASTQ Files

**Replace placeholders with actual values:**
- `USER`: your username on ARM machine (probably 'jacobc')
- `ARM_MACHINE_IP`: IP address or hostname of ARM machine

```bash
# Mount the sra_downloads directory from ARM machine
sshfs USER@ARM_MACHINE_IP:/home/jacobc/sc/sra_downloads ~/sc/sra_downloads

# Verify mount worked
ls ~/sc/sra_downloads
# Should show: SRR14513984_2.fastq.gz, SRR14513984_3.fastq.gz, etc.

# Count FASTQ files (should be 36: 18 samples × 2 files)
ls ~/sc/sra_downloads/*.fastq.gz | wc -l
# Expected: 36
```

**IMPORTANT: Keep this terminal/connection alive** during CellRanger processing, or add to `/etc/fstab` for persistent mount.

#### 3.4 Create CellRanger-Compatible Symlinks

```bash
cd ~/sc/morabito_18samples/cellranger_analysis

# Run the symlink preparation script
bash prepare_fastq_symlinks.sh

# Verify symlinks created
ls -lh fastq_cellranger_format/SRR14513984/
# Should show: SRR14513984_S1_L001_R1_001.fastq.gz, SRR14513984_S1_L001_R2_001.fastq.gz
```

---

### PHASE 4: Download CellRanger

#### 4.1 Get CellRanger Download URL

1. Open browser and go to: https://www.10xgenomics.com/support/software/cell-ranger/downloads
2. Sign in or create account (free)
3. Find "Cell Ranger 10.0.0"
4. Right-click "Download" → Copy link address
5. The URL will look like: `https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz?Expires=...`

#### 4.2 Download CellRanger

```bash
cd ~/sc/morabito_18samples/cellranger_analysis

# Replace YOUR_SIGNED_URL with the URL you copied
curl -o cellranger-10.0.0.tar.gz "YOUR_SIGNED_URL"

# This will take ~5-10 minutes (~828MB download)
```

#### 4.3 Verify Download

```bash
# Check file size
ls -lh cellranger-10.0.0.tar.gz
# Should be ~828MB

# Verify it's a valid tar.gz
tar -tzf cellranger-10.0.0.tar.gz | head -5
# Should show cellranger-10.0.0/ directory structure
```

---

### PHASE 5: Build Docker Image

⚠️ **This is the LONGEST step** (~20-30 minutes)

#### 5.1 Check Prerequisites

```bash
cd ~/sc/morabito_18samples/cellranger_analysis

# Run prerequisite check (this will fail on some checks, that's OK)
bash check_prerequisites.sh
```

#### 5.2 Build Docker Image

```bash
# Start the build
bash build_docker.sh

# When prompted "Continue with build? (y/n)", type: y

# Expected output:
# - Step 1-5: Install system packages (~2 min)
# - Step 6-7: Extract CellRanger (~1 min)
# - Step 8-9: Download reference genome (~15-20 min, ~11GB download)
# - Step 10-14: Final setup (~2 min)
```

**IMPORTANT**: 
- This downloads ~11GB reference genome
- Total image size: ~21GB
- Monitor disk space: `df -h ~`
- If download fails, re-run `bash build_docker.sh` (it will resume)

#### 5.3 Verify Image Built

```bash
# Check Docker images
docker images | grep cellranger-morabito

# Expected output:
# cellranger-morabito  8.0.1   <IMAGE_ID>   X minutes ago   21.4GB
# cellranger-morabito  latest  <IMAGE_ID>   X minutes ago   21.4GB
```

---

### PHASE 6: Configure for Low-Memory System

#### 6.1 Modify Batch Script for Sequential Processing

The default `run_cellranger_batch.sh` tries to use 64GB RAM per job. We need to reduce this:

```bash
cd ~/sc/morabito_18samples/cellranger_analysis

# Edit the batch script
nano run_cellranger_batch.sh
```

**Find this line**:
```bash
    --localmem=64 \
```

**Change to**:
```bash
    --localmem=12 \
```

Also find this line:
```bash
    --localcores=8 \
```

**Change to** (use output from `nproc` earlier, or use 4 if unsure):
```bash
    --localcores=4 \
```

**Save and exit** (Ctrl+O, Enter, Ctrl+X)

---

### PHASE 7: Test with Single Sample

Before running all 18 samples (18-36 hours), test with ONE sample:

#### 7.1 Launch Docker Container

```bash
cd ~/sc/morabito_18samples/cellranger_analysis

# Start container
bash run_docker.sh

# You should now be INSIDE the container
# Prompt will change to: root@<container-id>:/app#
```

#### 7.2 Verify CellRanger Inside Container

```bash
# Check CellRanger version
cellranger --version
# Expected: cellranger cellranger-10.0.0

# Check reference genome
ls -lh $CELLRANGER_REFERENCE
# Should show refdata-gex-GRCh38-2024-A directory

# Check mounted FASTQ files
ls /app/fastq_data
# Should show 18 sample directories: SRR14513984, SRR14513985, etc.

# Check one sample's FASTQ files
ls -lh /app/fastq_data/SRR14513984/
# Should show: SRR14513984_S1_L001_R1_001.fastq.gz (symlink)
#              SRR14513984_S1_L001_R2_001.fastq.gz (symlink)
```

#### 7.3 Run Single Test Sample

Pick a small sample for testing. Let's use SRR14513984:

```bash
# INSIDE CONTAINER - Run CellRanger on ONE sample
cellranger count \
  --id=SRR14513984_test \
  --transcriptome=$CELLRANGER_REFERENCE \
  --fastqs=/app/fastq_data/SRR14513984 \
  --sample=SRR14513984 \
  --chemistry=SC3Pv3 \
  --localcores=4 \
  --localmem=12 \
  --nosecondary

# This will take 1-2 hours
# Monitor progress - you'll see stages:
# - Aligning reads
# - Counting UMIs
# - Generating summary
```

**Expected Output Location**:
```bash
# Check output directory
ls -lh /app/cellranger_output/SRR14513984_test/outs/

# Key files:
# - filtered_feature_bc_matrix.h5 (~50-200MB)
# - metrics_summary.csv
# - web_summary.html
```

#### 7.4 Verify Test Results

```bash
# INSIDE CONTAINER - Check if output exists
ls -lh /app/cellranger_output/SRR14513984_test/outs/filtered_feature_bc_matrix.h5

# Check metrics
cat /app/cellranger_output/SRR14513984_test/outs/metrics_summary.csv

# Exit container (keep it running in background)
# Press Ctrl+P then Ctrl+Q (detaches without stopping)
```

From WSL (outside container):
```bash
# Check output on host
ls -lh ~/sc/morabito_18samples/cellranger_analysis/cellranger_output/SRR14513984_test/

# If successful, you're ready for full batch!
```

---

### PHASE 8: Process All 18 Samples (LONG-RUNNING)

⚠️ **This will take 18-36 hours** (1-2 hours per sample × 18 samples)

#### 8.1 Start Batch Processing

```bash
# If you exited container, restart it
cd ~/sc/morabito_18samples/cellranger_analysis
bash run_docker.sh

# INSIDE CONTAINER - Start batch processing
bash /app/scripts/run_cellranger_batch.sh

# Or run in background with nohup:
nohup bash /app/scripts/run_cellranger_batch.sh > /app/cellranger_output/batch.log 2>&1 &

# Detach from container: Ctrl+P then Ctrl+Q
```

#### 8.2 Monitor Progress

From WSL (outside container):

```bash
# Check which samples have completed
ls ~/sc/morabito_18samples/cellranger_analysis/cellranger_output/

# Check log file (if using nohup)
tail -f ~/sc/morabito_18samples/cellranger_analysis/cellranger_output/batch.log

# Check disk space periodically
df -h ~

# Count completed samples (should eventually be 18)
ls -d ~/sc/morabito_18samples/cellranger_analysis/cellranger_output/SRR*/outs/ 2>/dev/null | wc -l
```

#### 8.3 Resume if Interrupted

If processing stops (computer sleeps, network drops, etc.):

```bash
# Re-mount FASTQ files if needed
sshfs USER@ARM_MACHINE_IP:/home/jacobc/sc/sra_downloads ~/sc/sra_downloads

# Restart container
cd ~/sc/morabito_18samples/cellranger_analysis
bash run_docker.sh

# INSIDE CONTAINER - The script will skip completed samples
bash /app/scripts/run_cellranger_batch.sh
```

---

### PHASE 9: Merge Outputs and Compare with STAR

Once all 18 samples complete:

#### 9.1 Verify All Samples Completed

```bash
# INSIDE CONTAINER
cd /app/cellranger_output

# Check each sample has output
for sample in SRR14513984 SRR14513985 SRR14513986 SRR14513987 SRR14513988 SRR14513989 SRR14513990 SRR14513991 SRR14513992 SRR14513993 SRR14513994 SRR14513995 SRR14514096 SRR14514124 SRR14514125 SRR14514126 SRR14514127 SRR14514128; do
  if [ -f "$sample/outs/filtered_feature_bc_matrix.h5" ]; then
    echo "✓ $sample completed"
  else
    echo "✗ $sample MISSING"
  fi
done
```

#### 9.2 Merge CellRanger Outputs

This combines all 18 samples into one AnnData file:

```bash
# INSIDE CONTAINER
cd /app

# Run merge script
python merge_cellranger_outputs.py

# Expected output:
# - Loads 18 × filtered_feature_bc_matrix.h5
# - Filters to 60,328 cells (Morabito subset)
# - Adds metadata (cell_type, diagnosis, etc.)
# - Saves: /app/cellranger_output/merged_cellranger_18samples.h5ad

# Verify output
ls -lh /app/cellranger_output/merged_cellranger_18samples.h5ad
# Should be ~500MB-1GB
```

#### 9.3 Compare CellRanger vs STAR

This generates comparison CSV files:

```bash
# INSIDE CONTAINER
cd /app

# Run comparison script
python compare_with_star.py

# Expected output files:
# - gene_expression_comparison.csv (per-gene correlations)
# - per_cell_comparison.csv (per-cell UMI counts)
# - summary_cellranger_vs_star.csv (overall statistics)

# Verify outputs
ls -lh /app/cellranger_output/*.csv
```

#### 9.4 Exit Container

```bash
# INSIDE CONTAINER
exit

# Back in WSL
```

---

### PHASE 10: Copy Results Back to ARM Machine

#### 10.1 Copy Comparison Results

```bash
# From WSL, copy comparison CSVs to ARM machine
scp ~/sc/morabito_18samples/cellranger_analysis/cellranger_output/*.csv \
    USER@ARM_MACHINE_IP:/home/jacobc/sc/morabito_18samples/cellranger_analysis/cellranger_output/

# Copy merged AnnData file
scp ~/sc/morabito_18samples/cellranger_analysis/cellranger_output/merged_cellranger_18samples.h5ad \
    USER@ARM_MACHINE_IP:/home/jacobc/sc/morabito_18samples/cellranger_analysis/cellranger_output/
```

#### 10.2 Verify on ARM Machine

On ARM machine:
```bash
cd /home/jacobc/sc/morabito_18samples/cellranger_analysis/cellranger_output

ls -lh *.csv
# Should show:
# - gene_expression_comparison.csv
# - per_cell_comparison.csv
# - summary_cellranger_vs_star.csv

ls -lh *.h5ad
# Should show:
# - merged_cellranger_18samples.h5ad
```

---

## Troubleshooting

### SSHFS Mount Disconnects

If network mount drops during processing:

```bash
# Unmount
fusermount -u ~/sc/sra_downloads

# Re-mount
sshfs USER@ARM_MACHINE_IP:/home/jacobc/sc/sra_downloads ~/sc/sra_downloads

# Add reconnect option for better stability
sshfs -o reconnect,ServerAliveInterval=15,ServerAliveCountMax=3 \
      USER@ARM_MACHINE_IP:/home/jacobc/sc/sra_downloads ~/sc/sra_downloads
```

### Out of Memory Errors

If CellRanger crashes with OOM:

1. Check memory usage: `free -h`
2. Reduce `--localmem=12` to `--localmem=10` or `--localmem=8`
3. Close other applications
4. Ensure no other Docker containers running: `docker ps`

### Out of Disk Space

If disk fills up:

```bash
# Check disk usage
df -h ~

# Remove test output
rm -rf ~/sc/morabito_18samples/cellranger_analysis/cellranger_output/SRR14513984_test

# Clean Docker cache
docker system prune -a

# Remove intermediate files (ONLY if all samples complete)
# This removes BAM files, keeping only summary files
cd ~/sc/morabito_18samples/cellranger_analysis/cellranger_output
for sample in SRR*; do
  rm -rf $sample/outs/possorted_genome_bam.bam
  rm -rf $sample/outs/possorted_genome_bam.bam.bai
done
```

### Symlinks Broken

If CellRanger can't find FASTQ files:

```bash
# Check if mount is active
ls ~/sc/sra_downloads/SRR14513984_2.fastq.gz

# Check symlinks
ls -lh ~/sc/morabito_18samples/cellranger_analysis/fastq_cellranger_format/SRR14513984/

# Recreate symlinks if needed
cd ~/sc/morabito_18samples/cellranger_analysis
bash prepare_fastq_symlinks.sh
```

### Docker Container Won't Start

```bash
# Check Docker is running
docker ps

# If error "Cannot connect to Docker daemon"
# Restart Docker Desktop on Windows

# Check image exists
docker images | grep cellranger-morabito

# If missing, rebuild
bash build_docker.sh
```

---

## Expected Timeline

| Phase | Task | Duration |
|-------|------|----------|
| 1-2 | Setup environment + clone repo | 30 min |
| 3 | Mount FASTQ files | 5 min |
| 4 | Download CellRanger | 10 min |
| 5 | Build Docker image | 25 min |
| 6-7 | Test single sample | 2 hours |
| 8 | **Process all 18 samples** | **18-36 hours** |
| 9 | Merge + compare | 30 min |
| 10 | Copy results back | 20 min |
| **TOTAL** | **~22-40 hours** |

**Recommendation**: Start overnight or over a weekend.

---

## Resource Usage Summary

### Disk Space
- Docker image: 21GB
- CellRanger tar.gz: 1GB
- Output (18 samples): ~4GB
- Workspace: ~2GB
- **Total**: ~28GB (within 50GB limit ✓)

### RAM
- Per sample: 12GB
- Docker overhead: 2GB
- System: 2GB
- **Total**: ~16GB (matches your RAM ✓)

### Network
- Reference download: 11GB (one-time)
- FASTQ streaming: ~60GB (spread over 18-36 hours)
- Results upload: ~1GB

---

## Quick Command Reference

```bash
# Mount FASTQ files
sshfs -o reconnect USER@ARM_IP:/home/jacobc/sc/sra_downloads ~/sc/sra_downloads

# Check mount
ls ~/sc/sra_downloads/*.fastq.gz | wc -l  # Should be 36

# Create symlinks
cd ~/sc/morabito_18samples/cellranger_analysis
bash prepare_fastq_symlinks.sh

# Build Docker
bash build_docker.sh

# Run container
bash run_docker.sh

# Inside container - test one sample
cellranger count --id=TEST --transcriptome=$CELLRANGER_REFERENCE \
  --fastqs=/app/fastq_data/SRR14513984 --sample=SRR14513984 \
  --chemistry=SC3Pv3 --localcores=4 --localmem=12 --nosecondary

# Inside container - run all samples
bash /app/scripts/run_cellranger_batch.sh

# Inside container - merge and compare
python merge_cellranger_outputs.py
python compare_with_star.py

# Unmount FASTQ
fusermount -u ~/sc/sra_downloads

# Copy results back
scp cellranger_output/*.csv USER@ARM_IP:/home/jacobc/sc/morabito_18samples/cellranger_analysis/cellranger_output/
```

---

## Notes for AI Agent

When working through this guide:

1. **Verify each step completes** before moving to next
2. **Check disk space** after Phases 4, 5, and 8: `df -h ~`
3. **Monitor RAM usage** during test (Phase 7): `free -h`
4. **Keep SSHFS mount alive** - test with: `ls ~/sc/sra_downloads`
5. **Adapt core/memory settings** based on `nproc` and `free -h` output
6. **Save logs** from long-running processes (Phase 8)
7. **Be patient** - Phase 8 takes 18-36 hours sequentially

Good luck! This setup is optimized for your 16GB/50GB system.

