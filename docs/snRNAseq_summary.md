# snRNA-seq Analysis Summary

**Date**: December 11, 2025

## Dataset Overview

- **Total cells**: 574,298 nuclei
- **Total features**: 62,651 (genes + TEs)
- **Samples**: 150 (10x Chromium v3)
- **Total UMI counts**: 2.3 billion

## Cell Distribution

- **Mean cells per sample**: 3,829
- **Median cells per sample**: 3,462
- **Range**: 2,666 - 8,616 cells per sample

## TE Expression Summary

- **TE features detected**: 4,763
- **Total TE counts**: 1.44 billion (62.7% of all UMIs)
- **Top TE families**: Alu (dominant), L2, MIR, L1

## Top 15 Most Expressed TEs

1. AluJb - 69.6M
2. AluY - 68.8M  
3. AluSx1 - 68.1M
4. AluSx - 65.5M
5. AluSz - 59.5M
6. L2a - 48.1M
7. AluJr - 42.9M
8. AluJo - 38.3M
9. AluSq2 - 35.6M
10. MIRb - 33.3M
11. L2c - 32.0M
12. AluSp - 31.6M
13. AluSz6 - 27.4M
14. MIR - 26.9M
15. AluSg - 21.3M

## Top Non-TE Features

1. MALAT1 (lncRNA) - 73.2M
2. TALAM1 - 65.4M
3. NEAT1 - 3.4M
4. PLP1 (oligodendrocyte marker) - 2.7M
5. MEG3 (imprinted gene) - 2.7M

## Sample Metadata

- **GEO Series**: GSE175952
- **GSM IDs**: Multiple samples per GSM (technical replicates)
- **Need to download**: Sample metadata to identify AD vs Control

## Next Steps

1. Download GEO sample metadata (GSE175952)
2. Map SRR → GSM → Condition (AD/Control)
3. Basic QC (genes per cell, UMI counts, mitochondrial %)
4. Cell type clustering
5. Compare with paper's reported cell types
6. Differential TE expression: AD vs Control

## Paper Reference

Morabito et al. (2021) "Single-nucleus chromatin accessibility and transcriptomic characterization of Alzheimer's disease"
Nature Genetics. https://doi.org/10.1038/s41588-021-00894-z
