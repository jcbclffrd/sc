# ATAC-seq Data Inventory - PRJNA729525 Morabito Dataset

## Summary
- **Total ATAC-seq samples in metadata**: 32 (SRR14514129-SRR14514160)
- **Samples with FASTQ data on Spark**: 9 samples
- **Total size**: ~68 GB
- **Location**: spark-bd86:~/sc/sra_downloads/

## Available ATAC-seq Samples on Spark (9/32)

| Sample ID | FASTQ Files | Size | Read Numbers | Status |
|-----------|-------------|------|--------------|--------|
| SRR14514140 | 2 | 14 GB | _2, _4 | ✅ Complete |
| SRR14514141 | 2 | 14 GB | _2, _4 | ✅ Complete |
| SRR14514142 | 2 | 13 GB | _2, _4 | ✅ Complete |
| SRR14514143 | 2 | 13 GB | _2, _4 | ✅ Complete |
| SRR14514146 | 2 | 3.0 GB | _2, _4 | ✅ Complete |
| SRR14514147 | 2 | 2.7 GB | _2, _4 | ✅ Complete |
| SRR14514148 | 2 | 2.3 GB | _2, _4 | ✅ Complete |
| SRR14514149 | 2 | 3.2 GB | _2, _4 | ✅ Complete |
| SRR14514150 | 2 | 6.4 GB | _2, _4 | ✅ Complete |

## Missing ATAC-seq Samples (23/32)
The following samples exist in metadata but NOT downloaded:
- SRR14514129-139 (except 140-143) - 9 samples
- SRR14514144-145 (empty directories)
- SRR14514151-160 - 10 samples

## File Structure & Naming

Each sample has 2 paired-end FASTQ files:
```bash
~/sc/sra_downloads/SRR14514140/
├── SRR14514140_2.fastq.gz  (Read 1 - 6.5 GB)
└── SRR14514140_4.fastq.gz  (Read 2 - 6.6 GB)
```

**Note**: The read numbers are _2 and _4 (not _1 and _2), which is the SRA format for paired-end data.

## Sample Characteristics

From sra_metadata.csv:
- **Platform**: Illumina NovaSeq 6000
- **Library Strategy**: ATAC-seq
- **Library Source**: GENOMIC
- **Layout**: PAIRED
- **Read Lengths**: 
  - Samples 140-144: 124bp paired-end
  - Samples 145-160: 224bp paired-end

## Complete List of All 32 ATAC-seq Samples

```
SRR14514129  ❌ Not downloaded
SRR14514130  ❌ Not downloaded
SRR14514131  ❌ Not downloaded
SRR14514132  ❌ Not downloaded
SRR14514133  ❌ Not downloaded
SRR14514134  ❌ Not downloaded
SRR14514135  ❌ Not downloaded
SRR14514136  ❌ Not downloaded
SRR14514137  ❌ Not downloaded
SRR14514138  ❌ Not downloaded
SRR14514139  ❌ Not downloaded
SRR14514140  ✅ Downloaded (14 GB)
SRR14514141  ✅ Downloaded (14 GB)
SRR14514142  ✅ Downloaded (13 GB)
SRR14514143  ✅ Downloaded (13 GB)
SRR14514144  ❌ Empty folder
SRR14514145  ❌ Empty folder
SRR14514146  ✅ Downloaded (3.0 GB)
SRR14514147  ✅ Downloaded (2.7 GB)
SRR14514148  ✅ Downloaded (2.3 GB)
SRR14514149  ✅ Downloaded (3.2 GB)
SRR14514150  ✅ Downloaded (6.4 GB)
SRR14514151  ❌ Not downloaded
SRR14514152  ❌ Not downloaded
SRR14514153  ❌ Not downloaded
SRR14514154  ❌ Not downloaded
SRR14514155  ❌ Not downloaded
SRR14514156  ❌ Not downloaded
SRR14514157  ❌ Not downloaded
SRR14514158  ❌ Not downloaded
SRR14514159  ❌ Not downloaded
SRR14514160  ❌ Not downloaded
```

## Access Commands

```bash
# SSH to Spark server
ssh spark-bd86

# Navigate to ATAC-seq data
cd ~/sc/sra_downloads

# List all available ATAC-seq samples
ls -d SRR14514{140..150} 2>/dev/null | xargs -I{} sh -c 'if [ $(ls {}/*.fastq.gz 2>/dev/null | wc -l) -gt 0 ]; then echo {}; fi'

# Check a specific sample
ls -lh ~/sc/sra_downloads/SRR14514140/

# Verify FASTQ format
zcat ~/sc/sra_downloads/SRR14514140/SRR14514140_2.fastq.gz | head -4

# Get total size
du -sh ~/sc/sra_downloads/SRR14514{140..150} 2>/dev/null | grep -v "^4.0K"
```

## Data Quality Notes

- Files are gzip-compressed FASTQ format
- Average ~6-7 GB per FASTQ file
- Read lengths vary (124bp vs 224bp between sample groups)
- Some samples (144, 145) have empty directories suggesting failed downloads

## Potential Next Steps

1. **Complete download** of remaining 23 samples (~138 GB estimated)
2. **ATAC-seq alignment pipeline**:
   - Align with Bowtie2 or BWA-MEM
   - Remove duplicates with Picard
   - Peak calling with MACS2
   - TE-specific accessibility analysis
3. **Integration with snRNA-seq data**:
   - Correlate TE chromatin accessibility with expression
   - Identify TEs with open chromatin in AD vs control
   - Cell-type-specific TE accessibility analysis

## Reference
- **Project**: PRJNA729525
- **Paper**: Morabito et al. (2021) Nature Genetics
- **GEO**: GSE175952
- **Metadata**: ~/sc/data/sra_metadata.csv (also at /home/j/scWSL/data/sra_metadata.csv)

---
**Last Updated**: December 25, 2025  
**Generated from**: spark-bd86:~/sc/sra_downloads/
