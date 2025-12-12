# Comparison Report: scTE vs SoloTE

**Analysis Date**: December 12, 2025

**Dataset**: Morabito et al. (2021) - 18 snRNA-seq samples

**Purpose**: Compare two TE quantification methods on the same 18 samples and ~61,000 cells

---

## Executive Summary

**TE Detection:**
- scTE detected **984 TEs**
- soloTE detected **1,084 TEs**
- Overlap: **953 TEs (87.9%)**

**Differential Expression (AD vs Control):**
- scTE found **374 significant TEs**
- soloTE found **564 significant TEs**
- Agreement: **365 TEs (97.6% concordance)**

**Expression Correlation:**
- Pearson r = **0.983** (strong positive correlation)

**Differential Expression Correlation:**
- Log fold-change Pearson r = **0.994**

---

## Summary Statistics

| Metric | scTE | soloTE | Overlap |
|--------|------|--------|---------|
| TEs detected | 984 | 1,084 | 953 |
| Significant TEs (AD vs Control) | 374 | 564 | 365 |
| Upregulated in AD | 318 | 468 | 313 |
| Downregulated in AD | 56 | 96 | 51 |

---

## Top Differentially Expressed TEs (AD vs Control)

### Top 20 TEs by scTE (padj < 0.05)

| Rank | TE Name | scTE logFC | scTE padj | soloTE logFC | soloTE padj | Agreement |
|------|---------|------------|-----------|--------------|-------------|-----------|
| 1 | L1M4c | 0.397 | 0.00e+00 | 0.110 | 1.97e-35 | ✅ Both agree |
| 2 | L1M5 | 0.129 | 3.92e-296 | 0.145 | 0.00e+00 | ✅ Both agree |
| 3 | L1MEd | 0.178 | 8.54e-261 | 0.174 | 0.00e+00 | ✅ Both agree |
| 4 | L1MB7 | 0.148 | 7.02e-225 | 0.119 | 6.72e-199 | ✅ Both agree |
| 5 | L1ME1 | 0.121 | 1.50e-178 | 0.137 | 8.30e-300 | ✅ Both agree |
| 6 | L1MEg | 0.157 | 2.73e-157 | 0.159 | 5.65e-203 | ✅ Both agree |
| 7 | L1MEc | 0.135 | 5.41e-149 | 0.141 | 1.20e-204 | ✅ Both agree |
| 8 | L1MEf | 0.172 | 3.34e-141 | 0.180 | 9.29e-199 | ✅ Both agree |
| 9 | L1M4 | 0.121 | 2.60e-137 | 0.129 | 6.43e-209 | ✅ Both agree |
| 10 | HAL1 | 0.108 | 5.79e-134 | 0.111 | 8.76e-210 | ✅ Both agree |
| 11 | L1ME2 | 0.163 | 2.04e-129 | 0.171 | 1.56e-181 | ✅ Both agree |
| 12 | L1MC1 | 0.146 | 3.23e-123 | 0.153 | 1.28e-175 | ✅ Both agree |
| 13 | L1MC4 | 0.100 | 2.64e-108 | 0.107 | 2.29e-165 | ✅ Both agree |
| 14 | L1MA3 | 0.125 | 9.71e-90 | 0.129 | 2.22e-130 | ✅ Both agree |
| 15 | L1MB3 | 0.106 | 4.39e-87 | 0.106 | 5.64e-118 | ✅ Both agree |
| 16 | MamRep4096 | 0.537 | 8.34e-82 | 0.081 | 3.14e-03 | ✅ Both agree |
| 17 | L1MD2 | 0.125 | 5.79e-79 | 0.137 | 2.53e-118 | ✅ Both agree |
| 18 | L1M4b | 0.223 | 1.36e-76 | 0.102 | 1.50e-25 | ✅ Both agree |
| 19 | L1ME3Cz | 0.145 | 1.91e-76 | 0.144 | 2.57e-99 | ✅ Both agree |
| 20 | L1MC5a | 0.107 | 9.36e-73 | 0.116 | 3.46e-108 | ✅ Both agree |

---

### Top 20 TEs by soloTE (padj < 0.05)

| Rank | TE Name | TE Family | soloTE logFC | soloTE padj | scTE logFC | scTE padj | Agreement |
|------|---------|-----------|--------------|-------------|------------|-----------|-----------|
| 1 | L1M5 | LINE-1 | 0.145 | 0.00e+00 | 0.129 | 3.92e-296 | ✅ Both agree |
| 2 | L1MEd | LINE-1 | 0.174 | 0.00e+00 | 0.178 | 8.54e-261 | ✅ Both agree |
| 3 | L1ME1 | LINE-1 | 0.137 | 8.30e-300 | 0.121 | 1.50e-178 | ✅ Both agree |
| 4 | AluY | Alu | -0.115 | 2.15e-290 | -0.037 | 1.76e-37 | ✅ Both agree |
| 5 | HAL1 | Unknown | 0.111 | 8.76e-210 | 0.108 | 5.79e-134 | ✅ Both agree |
| 6 | L1M4 | LINE-1 | 0.129 | 6.43e-209 | 0.121 | 2.60e-137 | ✅ Both agree |
| 7 | L1MEc | LINE-1 | 0.141 | 1.20e-204 | 0.135 | 5.41e-149 | ✅ Both agree |
| 8 | L1MEg | LINE-1 | 0.159 | 5.65e-203 | 0.157 | 2.73e-157 | ✅ Both agree |
| 9 | L1MB7 | LINE-1 | 0.119 | 6.72e-199 | 0.148 | 7.02e-225 | ✅ Both agree |
| 10 | L1MEf | LINE-1 | 0.180 | 9.29e-199 | 0.172 | 3.34e-141 | ✅ Both agree |
| 11 | L1ME2 | LINE-1 | 0.171 | 1.56e-181 | 0.163 | 2.04e-129 | ✅ Both agree |
| 12 | L1MC1 | LINE-1 | 0.153 | 1.28e-175 | 0.146 | 3.23e-123 | ✅ Both agree |
| 13 | L1MC4 | LINE-1 | 0.107 | 2.29e-165 | 0.100 | 2.64e-108 | ✅ Both agree |
| 14 | L1MA3 | LINE-1 | 0.129 | 2.22e-130 | 0.125 | 9.71e-90 | ✅ Both agree |
| 15 | L1MD2 | LINE-1 | 0.137 | 2.53e-118 | 0.125 | 5.79e-79 | ✅ Both agree |
| 16 | L1MB3 | LINE-1 | 0.106 | 5.64e-118 | 0.106 | 4.39e-87 | ✅ Both agree |
| 17 | L1ME4a | LINE-1 | 0.098 | 9.32e-117 | 0.085 | 1.86e-60 | ✅ Both agree |
| 18 | L1MC5a | LINE-1 | 0.116 | 3.46e-108 | 0.107 | 9.36e-73 | ✅ Both agree |
| 19 | L2 | LINE-2 | 0.065 | 1.81e-107 | 0.056 | 1.58e-44 | ✅ Both agree |
| 20 | L2c | LINE-2 | 0.047 | 4.87e-103 | 0.037 | 6.03e-31 | ✅ Both agree |

---

## High-Confidence TEs (Both Methods Agree)

**364 TEs are significant in both methods with the same direction**

### Top 15 Upregulated TEs (Both Methods)

| TE Name | TE Family | scTE logFC | soloTE logFC | Avg logFC | scTE padj | soloTE padj |
|---------|-----------|------------|--------------|-----------|-----------|-------------|
| MLT1E3-int | LTR | 0.506 | 0.536 | 0.521 | 3.25e-02 | 2.88e-03 |
| MLT1-int | LTR | 0.366 | 0.399 | 0.382 | 1.44e-09 | 3.62e-11 |
| L1M3b | LINE-1 | 0.386 | 0.377 | 0.381 | 3.61e-11 | 1.18e-11 |
| L1M2b | LINE-1 | 0.360 | 0.371 | 0.365 | 9.68e-03 | 7.33e-04 |
| L1P4d | LINE-1 | 0.365 | 0.360 | 0.362 | 1.96e-07 | 3.55e-08 |
| LTR37-int | LTR | 0.302 | 0.365 | 0.333 | 7.17e-55 | 1.46e-72 |
| MamRep4096 | Unknown | 0.537 | 0.081 | 0.309 | 8.34e-82 | 3.14e-03 |
| L1M3d | LINE-1 | 0.297 | 0.284 | 0.290 | 4.20e-15 | 1.87e-14 |
| LTR80B | LTR | 0.279 | 0.293 | 0.286 | 7.84e-05 | 6.06e-06 |
| L1MEg1 | LINE-1 | 0.259 | 0.294 | 0.277 | 2.23e-18 | 4.99e-21 |
| L1M | LINE-1 | 0.275 | 0.262 | 0.268 | 5.09e-14 | 1.29e-13 |
| LTR81A | LTR | 0.268 | 0.264 | 0.266 | 3.83e-09 | 1.25e-09 |
| LTR24C | LTR | 0.265 | 0.262 | 0.264 | 1.04e-09 | 9.43e-11 |
| MERX | LTR | 0.245 | 0.271 | 0.258 | 6.99e-07 | 1.39e-08 |
| L1M3a | LINE-1 | 0.244 | 0.265 | 0.254 | 1.41e-10 | 1.63e-12 |

### Top 15 Downregulated TEs (Both Methods)

| TE Name | TE Family | scTE logFC | soloTE logFC | Avg logFC | scTE padj | soloTE padj |
|---------|-----------|------------|--------------|-----------|-----------|-------------|
| AluSz | Alu | -0.026 | -0.031 | -0.028 | 6.58e-17 | 8.28e-06 |
| AluSz6 | Alu | -0.023 | -0.039 | -0.031 | 8.64e-16 | 2.99e-13 |
| AluSx3 | Alu | -0.014 | -0.049 | -0.031 | 6.09e-04 | 4.54e-15 |
| AluJo | Alu | -0.033 | -0.036 | -0.035 | 3.43e-29 | 1.01e-06 |
| AluJr | Alu | -0.042 | -0.030 | -0.036 | 5.34e-57 | 5.30e-06 |
| AluJb | Alu | -0.043 | -0.033 | -0.038 | 1.27e-68 | 6.27e-03 |
| AluSg4 | Alu | -0.019 | -0.061 | -0.040 | 1.76e-03 | 4.25e-11 |
| AluSx | Alu | -0.038 | -0.044 | -0.041 | 7.66e-28 | 4.36e-20 |
| AluSx4 | Alu | -0.026 | -0.061 | -0.043 | 5.45e-09 | 7.59e-17 |
| AluSg | Alu | -0.027 | -0.061 | -0.044 | 5.00e-14 | 8.25e-38 |
| AluSq | Alu | -0.022 | -0.066 | -0.044 | 1.40e-06 | 5.13e-26 |
| AluSp | Alu | -0.032 | -0.064 | -0.048 | 1.85e-20 | 2.52e-50 |
| AluSc | Alu | -0.031 | -0.068 | -0.050 | 1.09e-23 | 7.94e-44 |
| AluSq2 | Alu | -0.038 | -0.063 | -0.051 | 4.99e-32 | 1.21e-44 |
| AluSg7 | Alu | -0.037 | -0.065 | -0.051 | 1.23e-11 | 4.56e-16 |

---

## Method Disagreements

| Agreement Category | Count |
|--------------------|-------|
| Neither significant | 428 |
| Both significant same direction | 364 |
| Missing data | 162 |
| Only soloTE significant | 153 |
| Only scTE significant | 7 |
| Both significant opposite direction | 1 |

### Top 10 TEs Only Significant in scTE

| TE Name | scTE logFC | scTE padj | soloTE logFC | soloTE padj |
|---------|------------|-----------|--------------|-------------|
| AluJr4 | -0.034 | 3.72e-15 | -0.026 | 1.68e-01 |
| MLT1J-int | -0.221 | 9.21e-04 | -0.097 | 3.84e-01 |
| Tigger13a | -0.042 | 9.77e-04 | -0.014 | 2.03e-01 |
| AluYh3a3 | 0.081 | 9.31e-03 | -0.017 | 8.98e-01 |
| L1PA11 | -0.044 | 1.99e-02 | -0.044 | 2.82e-01 |
| MER110 | 0.116 | 2.18e-02 | 0.044 | 5.55e-02 |
| Charlie10 | 0.058 | 2.54e-02 | 0.039 | 5.90e-02 |

### Top 10 TEs Only Significant in soloTE

| TE Name | TE Family | soloTE logFC | soloTE padj | scTE logFC | scTE padj |
|---------|-----------|--------------|-------------|------------|-----------|
| L1PA4 | LINE-1 | -0.147 | 2.96e-54 | -0.010 | 1.00e+00 |
| L1PA3 | LINE-1 | -0.131 | 7.30e-46 | -0.014 | 1.00e+00 |
| AluYk4 | Alu | -0.209 | 1.45e-40 | -0.032 | 7.41e-01 |
| L1PA6 | LINE-1 | -0.129 | 3.14e-29 | -0.026 | 1.00e+00 |
| AluYh3 | Alu | -0.154 | 1.25e-22 | -0.003 | 1.00e+00 |
| AluYk3 | Alu | -0.135 | 5.15e-20 | 0.046 | 2.69e-01 |
| L1P2 | LINE-1 | -0.166 | 7.41e-18 | -0.038 | 1.00e+00 |
| AluYf1 | Alu | -0.121 | 9.50e-18 | -0.020 | 5.71e-01 |
| AluYk2 | Alu | -0.145 | 1.20e-17 | 0.014 | 1.00e+00 |
| AluYe5 | Alu | -0.134 | 7.57e-16 | -0.047 | 8.55e-02 |

---

## Conclusions

1. **High concordance**: 365 TEs (97.6%) are significant in both methods

2. **Strong correlation**: Expression levels correlate well (r = 0.983)

3. **Method-specific findings**: Each method detects unique TEs
   - scTE unique: 9 TEs
   - soloTE unique: 199 TEs

4. **Recommendation**: Focus on the **364 high-confidence TEs** where both methods agree

---
