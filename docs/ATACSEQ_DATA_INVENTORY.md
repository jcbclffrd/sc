# ATAC-seq Data Inventory - Morabito PRJNA729525

## Location
**Server**: DGX Spark (spark-bd86)  
**Path**: `~/sc/sra_downloads/`

## Found Samples
Found **11 out of 32** ATAC-seq samples on the Spark server:

| Sample ID | Location | Status |
|-----------|----------|--------|
| SRR14514140 | spark-bd86:~/sc/sra_downloads/ | ✅ Present (~14GB) |
| SRR14514141 | spark-bd86:~/sc/sra_downloads/ | ✅ Present |
| SRR14514142 | spark-bd86:~/sc/sra_downloads/ | ✅ Present |
| SRR14514143 | spark-bd86:~/sc/sra_downloads/ | ✅ Present |
| SRR14514144 | spark-bd86:~/sc/sra_downloads/ | ✅ Present |
| SRR14514145 | spark-bd86:~/sc/sra_downloads/ | ✅ Present |
| SRR14514146 | spark-bd86:~/sc/sra_downloads/ | ✅ Present |
| SRR14514147 | spark-bd86:~/sc/sra_downloads/ | ✅ Present |
| SRR14514148 | spark-bd86:~/sc/sra_downloads/ | ✅ Present |
| SRR14514149 | spark-bd86:~/sc/sra_downloads/ | ✅ Present |
| SRR14514150 | spark-bd86:~/sc/sra_downloads/ | ✅ Present |

## Missing Samples (21 samples)
The following ATAC-seq samples are cataloged but NOT found on Spark:
- SRR14514129
- SRR14514130-139 (10 samples)
- SRR14514151-160 (10 samples)

## Data Details
- **Total size on Spark**: ~68 GB (11 samples)
- **File format**: Paired-end FASTQ.gz (e.g., *_2.fastq.gz, *_4.fastq.gz)
- **Average size per sample**: ~6 GB
- **Estimated total if all downloaded**: ~192 GB (32 samples × 6 GB)

## File Structure Example
```
~/sc/sra_downloads/SRR14514140/
├── SRR14514140_2.fastq.gz (6.5 GB)
└── SRR14514140_4.fastq.gz (6.6 GB)
```

## Sample Characteristics (from metadata)
- **Read lengths**: 
  - SRR14514140-144: 124bp paired-end
  - SRR14514145-160: 224bp paired-end
- **Platform**: Illumina NovaSeq 6000
- **Assay**: Single-nucleus ATAC-seq (snATAC-seq)

## Next Steps
1. Check if remaining 21 samples exist elsewhere on Spark
2. Decide if additional ATAC-seq samples should be downloaded
3. Consider ATAC-seq analysis pipeline:
   - Alignment with Bowtie2 or BWA
   - Peak calling with MACS2
   - TE accessibility analysis
   - Integration with snRNA-seq TE expression data

## Commands to Access
```bash
# SSH to Spark
ssh spark-bd86

# List ATAC-seq samples
ls ~/sc/sra_downloads/SRR14514{140..150}

# Check sample size
du -sh ~/sc/sra_downloads/SRR14514140

# List FASTQ files in sample
ls -lh ~/sc/sra_downloads/SRR14514140/
```

---
**Last Updated**: December 25, 2025  
**Project**: PRJNA729525 - Morabito et al. (2021) Alzheimer's Disease TE Analysis
