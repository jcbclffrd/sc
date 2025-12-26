# Multi-Omics Patient Mapping for Integrated Analysis

## Overview
**All 18 patients have BOTH snRNA-seq AND snATAC-seq data!**

This enables integrated multi-omics analysis of TE expression and chromatin accessibility in the same individuals.

## File Location
**Primary mapping file**: [data/patient_multiomics_mapping.csv](../data/patient_multiomics_mapping.csv)

## Summary Statistics

| Category | Count |
|----------|-------|
| **Total patients** | 18 |
| **With snRNA-seq** | 18 (100%) |
| **With ATAC-seq** | 18 (100%) |
| **Total ATAC runs** | 30 (some patients have multiple technical replicates) |

### By Diagnosis

| Diagnosis | N Patients | ATAC Runs |
|-----------|-----------|-----------|
| **Alzheimer's Disease (AD)** | 11 | 17 |
| **Control** | 7 | 13 |

## Data Columns

| Column | Description |
|--------|-------------|
| `patient_sample_id` | Patient identifier (Sample-17, Sample-52, etc.) |
| `sample_number` | Numeric ID |
| `diagnosis` | AD or Control |
| `age` | Age at death (years) |
| `sex` | M or F |
| `pmi` | Post-mortem interval (hours) |
| `tangle_stage` | Braak staging for neurofibrillary tangles |
| `plaque_stage` | CERAD staging for neuritic plaques |
| `rin` | RNA integrity number |
| `batch` | Sequencing batch (1, 2, or 3) |
| `snRNA_gsm` | GEO accession for snRNA-seq |
| `snRNA_srr` | SRA run accession for snRNA-seq |
| `atac_gsm` | GEO accession(s) for ATAC-seq |
| `atac_srr` | SRA run accession(s) for ATAC-seq (comma-separated if multiple) |
| `atac_num_runs` | Number of ATAC-seq technical replicates |
| `has_both_modalities` | TRUE for all 18 patients |

## Key Findings

### Multiple ATAC Technical Replicates
Some patients have 4 ATAC-seq runs (technical replicates that should be merged):
- **Sample-100** (Control): 4 runs
- **Sample-43** (AD): 4 runs  
- **Sample-45** (AD): 4 runs
- **Sample-96** (Control): 4 runs

Most patients have 1 ATAC-seq run.

### Clinical Characteristics

**AD Patients (n=11):**
- Age range: 80-90 years
- Advanced pathology (mostly Stage 5-6 tangles, Stage B-C plaques)

**Control Patients (n=7):**
- Age range: 79-90 years
- Low pathology (Stage 1-2 tangles, Stage A-B plaques)

## Example Integration Workflows

### 1. TE Expression vs Chromatin Accessibility
For each patient, compare:
- TE expression from snRNA-seq (already processed in morabito_18samples)
- TE accessibility from ATAC-seq (need to align and call peaks)

### 2. Cell-Type-Specific TE Analysis
- Use snRNA-seq cell type annotations
- Pseudo-bulk ATAC-seq by cell type from snATAC-seq
- Correlate TE expression with accessibility per cell type

### 3. AD vs Control Comparison
- 11 AD vs 7 Control matched for age
- Multi-omic differential analysis:
  - Which TEs are differentially expressed?
  - Which TEs have altered accessibility?
  - Do changes correlate?

## ATAC-seq Data Status

Currently downloading on spark-bd86. Status from latest check:

```
✓ Downloaded: 15/32 total ATAC runs
⏳ In progress: Remaining 15 runs
📊 Size: ~199 GB so far
```

**For our 18 patients:**
- SRR14514130-160 cover all our samples
- Need to track which ones are downloaded

## Next Steps

1. **Complete ATAC-seq download** (ongoing on Spark)
2. **Align ATAC-seq reads** with Bowtie2/BWA
3. **Call peaks** with MACS2
4. **TE-specific peak calling** around TE loci
5. **Integrate with snRNA-seq TE expression**
6. **Multi-omic differential analysis** (AD vs Control)

## Usage Example

```python
import pandas as pd

# Load mapping
mapping = pd.read_csv('data/patient_multiomics_mapping.csv')

# Get AD patients with both modalities
ad_patients = mapping[mapping['diagnosis'] == 'AD']

# For each patient, load their data
for _, patient in ad_patients.iterrows():
    sample_id = patient['patient_sample_id']
    snrna_file = f"morabito_18samples/soloTE_analysis/{sample_id}_TE_counts.tsv"
    atac_runs = patient['atac_srr'].split(',')
    
    # Process multi-omics data...
```

---
**Created**: December 25, 2025  
**Last Updated**: December 25, 2025  
**Source**: GEO GSE174367, SRA PRJNA729525
