## Docker-based ATAC-seq Pipeline

Complete containerized pipeline from SRA download to differential analysis.

## Architecture

**Modular Docker containers:**
1. **Base image** - Common tools (STAR, SRA toolkit, R, etc.)
2. **Alignment** - Download from SRA + STAR alignment
3. **Peaks** - MACS3 peak calling
4. **Quantify** - featureCounts for TE quantification
5. **DESeq2** - Differential analysis

## Quick Start

### Prerequisites
```bash
# Install Docker
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh

# Add user to docker group (logout/login required)
sudo usermod -aG docker $USER
```

### Build Images
```bash
cd ~/sc
docker build -t atacseq-base:latest -f docker/Dockerfile.base .
docker build -t atacseq-alignment:latest -f docker/Dockerfile.alignment .
docker build -t atacseq-peaks:latest -f docker/Dockerfile.peaks .
docker build -t atacseq-quantify:latest -f docker/Dockerfile.quantify .
docker build -t atacseq-deseq2:latest -f docker/Dockerfile.deseq2 .
```

### Run Complete Pipeline
```bash
# Prepare genome index and annotations
# (These should already exist in your ~/sc directory)

# Run full pipeline
bash docker/run_pipeline.sh
```

## Manual Step-by-Step

### 1. Align One Sample
```bash
docker run --rm \
  -v ~/sc/star_index:/genome:ro \
  -v ~/sc/atacseq_aligned:/output \
  -v ~/sc/sra_downloads/ATAC-seq:/data \
  -e SAMPLE_ID=SRR14514129 \
  -e THREADS=4 \
  --cpus=4 \
  --memory=16g \
  atacseq-alignment:latest
```

### 2. Call Peaks
```bash
docker run --rm \
  -v ~/sc/atacseq_aligned:/data:ro \
  -v ~/sc/atacseq_peaks:/output \
  -e SAMPLE_ID=SRR14514129 \
  --cpus=2 \
  --memory=4g \
  atacseq-peaks:latest
```

### 3. Quantify TEs (All Samples)
```bash
docker run --rm \
  -v ~/sc/atacseq_aligned:/data:ro \
  -v ~/sc/annotations:/annotations:ro \
  -v ~/sc/atacseq_te_counts:/output \
  -e TE_GTF=/annotations/hg38_rmsk_TE.gtf \
  -e THREADS=8 \
  --cpus=8 \
  --memory=16g \
  atacseq-quantify:latest
```

### 4. Run DESeq2
```bash
docker run --rm \
  -v ~/sc/atacseq_te_counts:/data:ro \
  -v ~/sc/metadata:/metadata:ro \
  -v ~/sc/deseq2_results:/output \
  -e COUNT_MATRIX=/data/te_counts_matrix.tsv \
  -e METADATA=/metadata/atacseq_metadata.csv \
  -e CONDITION_COL=condition \
  -e CONTROL=Control \
  -e TREATMENT=AD \
  --cpus=4 \
  --memory=16g \
  atacseq-deseq2:latest
```

## Run on AWS EC2

### Launch c5.9xlarge instance
```bash
# See AWS_SETUP.md for instance launch instructions
```

### Install Docker on EC2
```bash
ssh -i ~/.ssh/awsWebsite.pem ubuntu@$AWS_IP

# Install Docker
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
sudo usermod -aG docker ubuntu
exit

# Reconnect
ssh -i ~/.ssh/awsWebsite.pem ubuntu@$AWS_IP
```

### Transfer and run
```bash
# On AWS instance
git clone https://github.com/jcbclffrd/sc.git
cd sc

# Upload genome index
rsync -avz -e "ssh -i ~/.ssh/awsWebsite.pem" \
  ~/sc/star_index/ ubuntu@$AWS_IP:~/sc/star_index/

# Upload TE annotations
rsync -avz -e "ssh -i ~/.ssh/awsWebsite.pem" \
  ~/sc/annotations/hg38_rmsk_TE.gtf ubuntu@$AWS_IP:~/sc/annotations/

# Run pipeline
bash docker/run_pipeline.sh
```

## Advantages of Docker

✅ **Reproducible** - Exact same environment everywhere  
✅ **Portable** - Run on local, AWS, GCP, Azure, HPC  
✅ **Isolated** - No conflicts with system packages  
✅ **Scalable** - Easy to parallelize with Kubernetes/Swarm  
✅ **Shareable** - Push to Docker Hub, others can reproduce  
✅ **Testable** - Test each module independently  

## Test Individual Modules

```bash
# Test alignment on one sample
SAMPLE_ID=SRR14514129 bash docker/run_pipeline.sh

# Test peak calling only
docker run --rm -v ~/sc/atacseq_aligned:/data:ro -v ~/sc/atacseq_peaks:/output \
  -e SAMPLE_ID=SRR14514129 atacseq-peaks:latest

# Test TE quantification on subset
docker run --rm -v ~/sc/atacseq_aligned:/data:ro -v ~/sc/annotations:/annotations:ro \
  -v ~/sc/atacseq_te_counts:/output atacseq-quantify:latest
```

## Resource Requirements

| Module | CPU | RAM | Time/sample |
|--------|-----|-----|-------------|
| Alignment | 4 cores | 16GB | ~30 min |
| Peaks | 2 cores | 4GB | ~5 min |
| Quantify | 8 cores | 16GB | ~20 min (all) |
| DESeq2 | 4 cores | 16GB | ~10 min |

**c5.9xlarge (36 cores, 72GB):**
- Run 9 alignments in parallel
- Total time: ~2 hours for 32 samples

## Troubleshooting

**Out of memory:**
```bash
# Increase memory limit
--memory=32g
```

**Permission errors:**
```bash
# Add user to docker group
sudo usermod -aG docker $USER
# Logout and login again
```

**Image not found:**
```bash
# Rebuild images
docker build -t atacseq-base:latest -f docker/Dockerfile.base .
```

**Container cleanup:**
```bash
# Remove stopped containers
docker container prune -f

# Remove unused images
docker image prune -a -f
```
