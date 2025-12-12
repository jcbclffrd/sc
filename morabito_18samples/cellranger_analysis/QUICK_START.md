# CellRanger Quick Start Guide

## Summary

You now have a complete CellRanger analysis setup to compare with your STAR+scTE results!

## What's Been Set Up

✅ **Docker configuration** - Platform: linux/amd64 (for ARM compatibility)
✅ **FASTQ symlinks** - 18 samples with CellRanger-compatible naming
✅ **Batch processing script** - Automated processing of all samples
✅ **Comparison scripts** - Compare CellRanger vs STAR gene counts

## File Structure

```
cellranger_analysis/
├── Dockerfile                        # CellRanger 8.0.1 + GRCh38
├── build_docker.sh                   # Build the Docker image
├── run_docker.sh                     # Start container with volume mounts
├── check_prerequisites.sh            # Verify system requirements
├── prepare_fastq_symlinks.sh         # ✓ DONE - Created symlinks
├── fastq_cellranger_format/          # ✓ READY - Symlinked FASTQs (18 samples)
├── run_cellranger_batch.sh           # Process all 18 samples
├── merge_cellranger_outputs.py       # Merge into AnnData
├── compare_with_star.py              # Compare with STAR results
└── cellranger_output/                # Output directory (will be created)
```

## Step-by-Step Instructions

### Step 1: Build Docker Image (~30 minutes)

```bash
cd /home/jacobc/sc/morabito_18samples/cellranger_analysis
bash build_docker.sh
```

This downloads:
- CellRanger 8.0.1 (~2GB)
- GRCh38-2024-A reference (~11GB)
- Total image size: ~15GB

### Step 2: Run Container

```bash
bash run_docker.sh
```

This starts an interactive container with:
- FASTQ data mounted (read-only)
- Output directory mounted (read-write)
- Sample mapping mounted (read-only)

### Step 3: Inside Container, Run CellRanger

#### Option A: Test Single Sample First (Recommended)

```bash
# Test on one sample (~1-2 hours on ARM)
cellranger count \
    --id=SRR14513984_test \
    --transcriptome=$CELLRANGER_REFERENCE \
    --fastqs=/app/fastq_data/SRR14513984 \
    --sample=SRR14513984 \
    --chemistry=SC3Pv3 \
    --localcores=8 \
    --localmem=64

# Check results
ls -lh /app/cellranger_output/SRR14513984_test/outs/
cat /app/cellranger_output/SRR14513984_test/outs/metrics_summary.csv
```

#### Option B: Run All 18 Samples

```bash
# Sequential processing (18-36 hours on ARM)
bash run_cellranger_batch.sh

# Or edit the script to run in parallel (6 at a time)
# Reduces time to 3-6 hours but uses more RAM
```

### Step 4: Exit Container and Merge Results

```bash
# Exit Docker
exit

# Back on host, merge outputs
source ~/hcaTE/.venv/bin/activate
python3 merge_cellranger_outputs.py

# Compare with STAR
python3 compare_with_star.py
```

## Expected Outputs

### Per-Sample Outputs
Each sample produces:
- `filtered_feature_bc_matrix.h5` - UMI count matrix
- `metrics_summary.csv` - QC metrics
- `web_summary.html` - Interactive report

### Merged Outputs
- `merged_cellranger_18samples.h5ad` - All samples, filtered to 60,328 cells

### Comparison CSVs
- `gene_expression_comparison.csv` - Per-gene: CellRanger vs STAR
- `per_cell_comparison.csv` - Per-cell: UMI counts, genes detected
- `summary_cellranger_vs_star.csv` - Overall statistics

## Performance Notes

⚠️ **ARM Emulation Warning**

Your system (DGX Spark ARM) will run x86_64 containers through emulation:
- **Expected slowdown**: 2-3x slower than native x86_64
- **Single sample**: 2-6 hours (vs 1-2 hours native)
- **All 18 samples**: 
  - Sequential: 36-108 hours
  - Parallel (6 jobs): 6-18 hours

💡 **Recommendations**:
1. Test with 1-2 samples first
2. Monitor memory usage
3. Consider running on x86_64 machine if available
4. Run overnight/over weekend for full batch

## Troubleshooting

### "No cells detected"
- Check FASTQ files: `ls /app/fastq_data/SRR14513984/`
- Verify chemistry: Should auto-detect SC3Pv3
- Check barcode format in R1 file

### Memory errors
- Reduce `--localmem` parameter (try 32 instead of 64)
- Run fewer samples in parallel
- Monitor with `docker stats`

### Slow performance
- Normal for ARM emulation (2-3x slower)
- Consider x86_64 machine for faster processing
- Monitor CPU usage: Should be near 100% per core

## What This Comparison Will Show

1. **Gene quantification accuracy**: Does STAR match CellRanger?
2. **UMI counting differences**: Different algorithms (STAR vs CellRanger)
3. **Cell filtering impact**: EmptyDrops vs manual threshold
4. **Quality validation**: Are your STAR results reliable?

Expected correlation: **r > 0.95** for gene expression if everything is working correctly!

## Current Status

✅ **FASTQ symlinks created** - All 18 samples ready
✅ **Scripts configured** - Paths set to your data
✅ **Prerequisites checked** - Docker working, 361GB space, 119GB RAM
⏸️ **Docker image** - Not built yet (run `build_docker.sh`)
⏸️ **Processing** - Waiting to start

## Comparison with SoloTE

While CellRanger is running (or if you skip it), you also have:
- **soloTE** - Currently processing (3/18 samples done)
- **scTE** - Already complete (318 TEs upregulated in AD)

After all three complete, you'll have:
1. **scTE vs soloTE** - TE quantification comparison
2. **CellRanger vs STAR** - Gene quantification validation
3. **Complete validation** - Both TEs and genes verified!

## Next Steps

Choose one:
1. **Wait for soloTE to finish** (~2-3 hours) - Then do CellRanger
2. **Start CellRanger now** - Run in parallel with soloTE
3. **Skip CellRanger** - Focus on TE analysis (scTE vs soloTE)

The choice is yours! 🚀
