# ATAC-seq Download Instructions for DGX Spark

## Location
Scripts are now on **spark-bd86** at: `~/sc/scripts/`

## Current Status
- **9/32 samples downloaded** (68 GB)
- **23 samples missing** (need ~138-160 GB more)

## How to Run

### 1. SSH to Spark Server
```bash
ssh spark-bd86
cd ~/sc
```

### 2. Check Current Status
```bash
./scripts/check_atacseq_download.sh
```

### 3. Start Downloading Missing Samples
```bash
# Run in a screen/tmux session (recommended for long downloads)
screen -S atacseq_download
./scripts/download_atacseq_samples.sh

# Detach: Ctrl+A, then D
# Reattach: screen -r atacseq_download
```

Or with nohup:
```bash
nohup ./scripts/download_atacseq_samples.sh > atacseq_download.log 2>&1 &
```

### 4. Monitor Progress (in another terminal)
```bash
# One-time check
./scripts/check_atacseq_download.sh

# Real-time monitoring (updates every 30 seconds)
watch -n 30 ./scripts/check_atacseq_download.sh

# Check logs of individual samples
tail -f logs_atacseq/SRR14514129.log
```

## What Will Be Downloaded

**23 missing samples:**
- SRR14514129-139 (11 samples, except 140-143 which exist)
- SRR14514144-145 (2 failed samples will be re-downloaded)
- SRR14514151-160 (10 samples)

## Expected Time & Space

- **Size**: ~138-160 GB (23 samples × 6-7 GB)
- **Time**: 6-12 hours depending on network speed
- **Output**: `~/sc/sra_downloads/SRR14514###/`
- **Logs**: `~/sc/logs_atacseq/SRR14514###.log`

## Script Features

- ✓ Automatically skips already-downloaded samples
- ✓ Re-downloads failed/incomplete samples
- ✓ Compresses FASTQ files with gzip
- ✓ Multi-threaded (8 threads per sample)
- ✓ Individual logs for each sample
- ✓ Progress tracking and summaries

## After Download Completes

Check the final status:
```bash
./scripts/check_atacseq_download.sh
```

All 32 samples should show "✓ Downloaded with data"

---
**Last Updated**: December 25, 2025  
**Server**: spark-bd86:~/sc/
