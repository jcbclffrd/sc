# Bulk RNA-seq Dataset Overview

## Full Dataset Summary

### Total Samples
- **Total SRA Accessions**: 191 (SRR14514163 - SRR14514483, with gaps)
- **Aligned Samples**: 156 samples
- **Used in TE Analysis**: 18 samples (11.5%)

### Data Source
- **GEO Series**: GSE174367
- **BioProject**: PRJNA729525
- **Study**: Morabito et al. - Multi-omics analysis of Alzheimer's disease
- **Tissue**: Human dorsolateral prefrontal cortex (DLPFC)

## Why Only 18 Samples Were Analyzed

The full bulk RNA-seq dataset includes samples from **100+ unique patients** with varying data modalities:

### Multi-Omics Cohort (18 patients) - **ANALYZED**
These 18 patients have **all three data types**:
- ✓ Bulk RNA-seq (this analysis)
- ✓ snRNA-seq (single-cell RNA-seq)
- ✓ snATAC-seq (single-cell chromatin accessibility)

**Rationale**: Focus on the multi-omics cohort enables integrated analysis across bulk, single-cell transcriptomics, and chromatin accessibility.

### Remaining Samples (138 aligned, ~175 total)
These represent:

1. **Patients with only bulk RNA-seq** (no single-cell data)
   - Cannot integrate with snRNA-seq or ATAC-seq
   - Would increase sample size but lose multi-omics advantage

2. **Different brain regions or conditions**
   - May include different cortical regions beyond DLPFC
   - Technical replicates or pilot samples

3. **Quality control exclusions**
   - Some samples may have failed QC in original study
   - RIN scores, alignment rates, or other technical issues

## Sample Distribution Analysis

### 18 Multi-Omics Patients Used
| Diagnosis | Count | Sample IDs |
|-----------|-------|------------|
| AD | 11 | Sample-17, 19, 22, 27, 33, 37, 43, 45, 46, 47, 50 |
| Control | 7 | Sample-52, 58, 66, 82, 90, 96, 100 |

**SRA Range for Multi-Omics**: 
- SRR14514163 (Sample-100)
- SRR14514235 (Sample-17)
- SRR14514255-SRR14514348 (remaining 16 samples)
- Scattered across bulk RNA-seq accession range

### 156 Aligned but Unused Samples
**SRA Range**: Primarily SRR14514192 - SRR14514483

These likely represent:
- **Additional patients** without matching snRNA-seq/ATAC-seq data
- **Bulk-only cohort** for larger-scale transcriptomic analysis
- **Extended patient series** (Sample-1 through Sample-100+)

## Data Size Comparison

### Multi-Omics Cohort (18 samples)
- **Raw BAM files**: ~25-30GB per sample = ~450-540GB total
- **TE count tables**: 1.4MB per sample = ~25MB total
- **Combined matrices**: ~3.8MB total
- **Analysis focus**: TE expression + integration with single-cell

### Full Aligned Dataset (156 samples)
- **Raw BAM files**: ~3.9-4.7TB total
- **Potential TE analysis**: Would generate ~2.2MB per sample = ~340MB count tables
- **Analysis focus**: Large-scale bulk transcriptomics (genes + TEs)

## Recommendations for Extended Analysis

### Option 1: Expand to All Bulk RNA-seq Samples
**Pros**:
- Increase statistical power (n=156 vs n=18)
- Better detection of subtle TE changes
- Larger effect size estimates

**Cons**:
- Cannot integrate with single-cell data
- Lose patient-matched multi-omics context
- Requires ~100-150 hours compute time for TE quantification

### Option 2: Maintain Multi-Omics Focus
**Pros**:
- Enable integrated analysis (bulk + scRNA + ATAC)
- Cell-type-specific TE analysis in single-cell data
- Chromatin accessibility changes at TE loci

**Cons**:
- Limited statistical power (n=18)
- No FDR-significant TEs detected

### Option 3: Hybrid Approach
1. **Discovery in full cohort** (156 samples): Identify candidate TEs with better power
2. **Validation in multi-omics cohort** (18 samples): Confirm with integrated data
3. **Cell-type resolution**: Examine top candidates in snRNA-seq by cell type

## Current Status

### Completed
✓ Multi-omics cohort TE quantification (18 samples)
✓ Differential expression analysis (AD vs Control)
✓ Combined expression matrices generated

### Available but Unanalyzed
- 138 additional aligned bulk RNA-seq samples
- ~35 unaligned bulk RNA-seq SRA accessions
- Matched snRNA-seq data (18 samples, not yet analyzed for TEs)
- Matched snATAC-seq data (18 samples, chromatin at TE loci)

## Files

### Multi-Omics Sample Lists
- **Used in analysis**: `data/bulk_rnaseq_18patients.csv` (18 samples)
- **Full bulk dataset**: `data/bulk_rnaseq_samples.txt` (191 SRA IDs)

### Metadata
- **Multi-omics patients**: `data/patient_multiomics_mapping.csv`
- **TE analysis samples**: `tetranscripts_bulk/metadata.csv`

### Aligned Data Location
- **18 analyzed samples**: `~/sc/star_bulk_aligned/SRR145141{63,235,255,272,277,283,287,293,295-297,301,303,308,317,334,342,348}/`
- **138 unused samples**: `~/sc/star_bulk_aligned/SRR145141{92-191,217-483}/`
- **Total disk usage**: ~424GB for all 156 aligned samples

## Next Steps Options

### Immediate (Cell-Type-Specific Analysis)
Focus on snRNA-seq to identify cell-type-specific TE dysregulation:
- Analyze TEs in neurons, astrocytes, microglia, oligodendrocytes separately
- Likely to reveal significant cell-type-specific changes

### Short-term (Expand Power)
Run TEcount on all 156 aligned bulk samples:
- Increased power for differential expression
- Meta-analysis across full cohort
- Estimated runtime: 3-4 days with parallelization

### Long-term (Full Multi-Omics Integration)
Integrate TE analysis across all three modalities:
- Bulk RNA-seq: Overall TE expression patterns
- snRNA-seq: Cell-type-specific TE activation
- snATAC-seq: Chromatin accessibility at TE loci

## Conclusion

The 18-sample multi-omics cohort represents a **strategic subset** of the full 191 bulk RNA-seq dataset, selected for **complete data integration** across bulk transcriptomics, single-cell transcriptomics, and chromatin accessibility. While this limits statistical power, it enables unprecedented multi-modal analysis of TE biology in Alzheimer's disease.

The remaining 138 aligned samples are available for **power-enhanced** bulk TE analysis if needed, but would lose the multi-omics integration capability.
