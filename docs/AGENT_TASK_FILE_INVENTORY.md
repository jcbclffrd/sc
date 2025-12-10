# Agent Task: Create FILE_INVENTORY.md for ~/sc Repository

## Objective
Create a comprehensive FILE_INVENTORY.md document that lists and describes every file in the ~/sc repository, organized by directory structure.

## Reference Example
See `/home/jacobc/hcaTE/FILE_INVENTORY.md` for the format and structure to follow.

## Requirements

### 1. Document Structure
```markdown
# File Inventory - [Project Name] Repository

> Complete listing of all files in the repository with descriptions

**Last Updated**: [Current Date]
**Repository**: [Repository name and purpose]

---

## Root Directory
[Table of root-level files with descriptions]

## [Directory Name]/
[Table of files in this directory with descriptions]

### [Subdirectory Name]/
[Table of files in subdirectory]

---

## Key File Relationships
[Workflow diagram showing how files connect]

---

## Statistics Summary
[Total samples, scripts, data size, etc.]

---

## Notes
[Important information about gitignored files, etc.]
```

### 2. What to Include

For each file, provide:
- **File name** (with relative path if in subdirectory)
- **Brief description** (one line, or a few lines for complex files)
- **Purpose** (what does this file do?)
- **Input/Output** (for scripts)
- **File size** (for large data files)

### 3. How to Organize

Group files by:
1. **Root directory** - README, config files, top-level scripts
2. **Data directories** - sra_downloads/, aligned_bams/, results/, etc.
3. **Scripts directory** - All analysis scripts organized by purpose:
   - Data download & preprocessing
   - Genome & annotation setup
   - Alignment pipeline
   - Quantification (scTE, etc.)
   - Analysis scripts
   - Utilities
4. **Results directories** - Output from analyses
5. **Documentation** - All .md files
6. **Logs and temporary files** - Note what can be deleted

### 4. Special Sections to Include

#### File Relationships Section
Show the workflow connections:
```
FASTQ files (sra_downloads/)
    ↓ [script_name.py]
BAM files (aligned_bams/)
    ↓ [script_name.py]
Quantification matrices
    ↓ [analysis_script.py]
Results & Figures
```

#### Statistics Summary
- Total samples downloaded
- Total cells analyzed
- Number of scripts
- Data size breakdown
- Key results numbers

#### Notes Section
- What files are gitignored (if git repo)
- What files can be safely deleted
- What files are critical
- Any important caveats

### 5. Commands to Help You

```bash
# List all files in repository
find ~/sc -type f -not -path '*/\.git/*' -not -name '*.pyc' | sort

# Count files by type
find ~/sc -type f | grep -E '\.(py|R|sh|md)$' | wc -l

# Check disk usage by directory
du -sh ~/sc/*/

# List largest files
find ~/sc -type f -exec ls -lh {} \; | sort -k5 -hr | head -20

# List all markdown files
find ~/sc -name "*.md" | sort

# List all scripts
ls -lh ~/sc/scripts/*.{py,R,sh} 2>/dev/null

# Check for README files
find ~/sc -name "README*" -o -name "*.md"
```

### 6. Information to Gather

Before writing, collect:
- Project name and purpose (check README.md if exists)
- Data source (BioProject, GEO accession, paper citation)
- Sequencing platform (10x Chromium, Smart-seq2, etc.)
- Number of samples
- Sample types (conditions, tissues, etc.)
- Analysis workflow (alignment → quantification → analysis)
- Key results or outputs
- Software/tools used (STAR, scTE, STARsolo, etc.)

### 7. Style Guidelines

**Good descriptions:**
- `download_sra.py` - "Downloads FASTQ files from SRA using prefetch for all 152 samples"
- `alignment.log` - "STARsolo alignment log with mapping statistics and cell detection"
- `SRR14514109_2.fastq.gz` - "R2: cDNA sequences (91bp) for 10x Chromium sample"

**Bad descriptions:**
- `script.py` - "Python script"
- `data.csv` - "Data file"
- `output/` - "Output directory"

**Be specific about:**
- What platform/technology (10x Chromium v3, Smart-seq2, etc.)
- What organism/genome (human hg38, mouse mm10)
- What analysis method (STARsolo, UMI-tools, scTE)
- File formats and structure
- Whether it's input, output, or intermediate data

### 8. Template Sections

Use these headers as appropriate:
- Root Directory
- scripts/ (with subsections by purpose)
- sra_downloads/ or data/
- aligned_bams/ or alignments/
- genome/ or reference/
- annotations/
- results/ or analysis/
- logs/ or log_files/
- figures/ or plots/
- documentation/
- File Relationships
- Statistics Summary
- Notes

## Example Entry Format

### Table Format (for many files):
```markdown
| File | Description |
|------|-------------|
| `script.py` | Brief description of what it does, inputs, outputs |
| `data.csv` | Description of data contents and structure |
```

### List Format (for directories):
```markdown
| Structure | Description |
|-----------|-------------|
| `SRR14514109/` | Sample directory with FASTQ files for one 10x run |
| `SRR*/SRR*_2.fastq.gz` | R2: cDNA reads (91bp) for alignment |
```

## Deliverable

Create `/home/jacobc/sc/FILE_INVENTORY.md` with:
1. Complete file listing organized by directory
2. One-line descriptions for every file
3. File relationships showing workflow
4. Statistics summary
5. Useful notes about the repository

The document should serve as a "file dictionary" - anyone can look up any file and immediately understand what it does and where it fits in the analysis.

## Quality Checklist

Before finishing, verify:
- [ ] Every directory is documented
- [ ] Every script has a description
- [ ] Data files are explained (format, content, size)
- [ ] Workflow connections are clear
- [ ] Statistics are accurate
- [ ] Notes section covers important caveats
- [ ] File is well-organized and easy to navigate
- [ ] Descriptions are specific, not generic

## Tips

1. Start by exploring the repository structure with `ls -R ~/sc/` or `tree ~/sc/`
2. Read the README.md first (if it exists) to understand the project
3. Check script headers/docstrings for purposes
4. Look at log files to understand what processes ran
5. Use `head` to peek at data file formats
6. Group similar files together
7. Be concise but informative
8. Think about what a new person would need to know

---

**Reference**: This task is based on the FILE_INVENTORY.md created for `/home/jacobc/hcaTE/` which you can examine as a template.
