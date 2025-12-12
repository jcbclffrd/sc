# CellRanger Analysis - Morabito 18 Samples

This folder contains CellRanger analysis for direct comparison with STAR+scTE gene quantification.

## Overview

Run CellRanger on the same 18 samples to compare:
- **Gene quantification**: CellRanger vs STAR+scTE
- **UMI counting methods**: CellRanger's algorithm vs STARsolo
- **Cell filtering**: CellRanger's EmptyDrops vs manual filtering

## Why This Comparison?

1. **Validation**: Ensure our STAR+scTE gene counts are accurate
2. **Benchmark**: Compare computational efficiency (STAR vs CellRanger)
3. **Quality**: See if CellRanger's filtering improves cell quality
4. **Reference**: Match Morabito's original CellRanger-based analysis

## Directory Structure

```
cellranger_analysis/
├── Dockerfile                          # Docker container with CellRanger
├── run_cellranger_batch.sh            # Process all 18 samples
├── merge_cellranger_outputs.py        # Merge into single AnnData
├── compare_with_star.py               # Compare with STAR results
├── cellranger_output/                 # CellRanger outputs (18 samples)
│   └── SRR*/
│       └── outs/
│           ├── filtered_feature_bc_matrix.h5
│           ├── metrics_summary.csv
│           └── web_summary.html
└── merged_cellranger_18samples.h5ad   # Final merged data
```

## Quick Start (Docker - Required for ARM)

**Important**: CellRanger only supports x86_64 (amd64) architecture. Since you're on ARM (DGX Spark), you **must** use Docker with platform emulation.

### Prerequisites

1. **Docker installed** with BuildKit support
2. **FASTQ files downloaded** to `/home/jacobc/sc/data/`
3. **~20GB free disk space** for Docker image
4. **120GB+ RAM** recommended (or run samples sequentially)

### 1. Build the Docker image

```bash
cd /home/jacobc/sc/morabito_18samples/cellranger_analysis

# Make scripts executable
chmod +x build_docker.sh run_docker.sh run_cellranger_batch.sh

# Build image (takes 20-30 minutes, downloads ~13GB)
bash build_docker.sh
```

This downloads and installs:
- CellRanger 8.0.1 (~2GB)
- GRCh38-2024-A reference (~11GB)
- Python environment with scanpy

### 2. Run the container

```bash
# Start container with volume mounts to your FASTQ data
bash run_docker.sh
```

This mounts:
- `/home/jacobc/sc/data/` → `/app/fastq_data` (read-only)
- `./cellranger_output/` → `/app/cellranger_output` (output)
- `../sample_mapping.csv` → `/app/sample_mapping/` (read-only)

### 3. Inside the container, run CellRanger

```bash
# Test single sample first (optional)
cellranger count \
    --id=SRR14513984_test \
    --transcriptome=$CELLRANGER_REFERENCE \
    --fastqs=/app/fastq_data/SRR14513984 \
    --sample=SRR14513984 \
    --chemistry=SC3Pv3 \
    --localcores=8 \
    --localmem=64

# Run all 18 samples (takes 9-18 hours on ARM emulation)
bash run_cellranger_batch.sh
```

### 4. Exit container and analyze results

```bash
# Exit Docker container
exit

# Back on host, merge outputs
source ~/hcaTE/.venv/bin/activate
python3 merge_cellranger_outputs.py
python3 compare_with_star.py
```

## Important Notes for ARM Users

⚠️ **Performance**: Running x86_64 containers on ARM is **2-3x slower** due to emulation
- Expected time per sample: **1-2 hours** on x86_64
- Expected time per sample: **2-6 hours** on ARM (emulated)
- Total for 18 samples: **36-108 hours** (sequential) or **6-18 hours** (6 parallel)

💡 **Recommendation**: If possible, run on an x86_64 machine for faster processing

## Using Docker

### Build the container

```bash
cd /home/jacobc/sc/morabito_18samples/cellranger_analysis
docker build -t cellranger-morabito .
```

### Run the container

```bash
docker run -it \
    -v $(pwd)/fastq_data:/app/fastq_data \
    -v $(pwd)/cellranger_output:/app/cellranger_output \
    cellranger-morabito
```

### Inside the container

```bash
# Run batch processing
bash run_cellranger_batch.sh

# Merge outputs
python3 merge_cellranger_outputs.py

# Compare with STAR
python3 compare_with_star.py
```

## Expected Outputs

### 1. Per-sample outputs (18 samples)

Each sample produces:
- `filtered_feature_bc_matrix.h5` - Filtered UMI count matrix (~50-200 MB)
- `metrics_summary.csv` - QC metrics (cells detected, median genes, etc.)
- `web_summary.html` - Interactive QC report
- `molecule_info.h5` - Per-molecule information (for aggregation)

### 2. Merged data

- `merged_cellranger_18samples.h5ad` - All 18 samples, filtered to Morabito's 60,328 cells
  - Cells × Genes matrix
  - Cell metadata (diagnosis, cell type, etc.)
  - QC metrics

### 3. Comparison files

- `gene_expression_comparison.csv` - Per-gene comparison (CellRanger vs STAR)
  - Columns: gene_name, cellranger_mean, star_mean, detection_category
  
- `per_cell_comparison.csv` - Per-cell comparison
  - Columns: cell_barcode, cellranger_umi_counts, star_umi_counts, genes_detected
  
- `summary_cellranger_vs_star.csv` - Overall statistics
  - Gene detection overlap
  - Expression correlation
  - Mean UMIs per cell

## Computational Requirements

### Single Sample
- **Time**: 30-60 minutes per sample
- **RAM**: 16-32 GB
- **Disk**: ~500 MB output per sample

### All 18 Samples
- **Time**: 9-18 hours (sequential) or 1.5-3 hours (6 parallel)
- **RAM**: 16 GB per job (96 GB for 6 parallel)
- **Disk**: ~10 GB output total

## CellRanger vs STAR Comparison

| Feature | CellRanger | STAR+scTE |
|---------|------------|-----------|
| **Algorithm** | Proprietary | Open-source |
| **Speed** | Slower (~1 hr/sample) | Faster (~20 min/sample) |
| **Memory** | Higher (~32 GB) | Lower (~16 GB) |
| **Cell calling** | EmptyDrops algorithm | Manual threshold |
| **Multi-mapping** | Default (10) | Custom (100 for TEs) |
| **TE support** | No | Yes (with scTE/soloTE) |

## Expected Results

Based on published comparisons, we expect:

1. **High gene correlation** (r > 0.95) for commonly detected genes
2. **Similar UMI counts** per cell (within 10-20%)
3. **CellRanger detects ~5-10% more genes** (better UMI deduplication)
4. **STAR detects more reads** (more sensitive alignment)

If correlation is low (<0.8), investigate:
- FASTQ file quality
- Reference genome versions (GRCh38 vs hg38)
- Chemistry detection (SC3Pv3 vs auto-detect)

## Troubleshooting

### CellRanger fails with "No cells detected"

Check:
1. FASTQ file structure: `*_R1_*.fastq.gz` and `*_R2_*.fastq.gz`
2. Chemistry version: Use `--chemistry=SC3Pv3` for 10x Chromium v3
3. Barcode whitelist: Should auto-detect for v3 chemistry

### Low gene detection

Check:
1. Alignment rate in `metrics_summary.csv` (should be >70%)
2. Median genes per cell (should be >1,000)
3. Fraction reads in cells (should be >50%)

### Memory errors

Reduce `--localmem` or process fewer samples in parallel

## References

- **CellRanger**: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
- **STAR**: https://github.com/alexdobin/STAR
- **Comparison**: https://doi.org/10.1186/s13059-021-02414-0 (Benchmarking study)

## Notes

- CellRanger uses GRCh38 (GENCODE), same as our STAR analysis
- Both use 10x barcode whitelist (3M-february-2018.txt)
- CellRanger cannot quantify TEs, so this is gene-only comparison
- For complete validation, we still need scTE vs soloTE TE comparison
