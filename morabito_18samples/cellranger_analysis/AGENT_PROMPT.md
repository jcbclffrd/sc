# AI Agent Instructions for CellRanger WSL Setup

## Context

I need help setting up and running a CellRanger analysis pipeline on a **Windows WSL Ubuntu machine with limited resources** (16GB RAM, 50GB storage). This analysis will process 18 single-cell RNA-seq samples and compare gene expression quantification between CellRanger and STAR alignment methods.

## Your Mission

Please follow the **WSL_SETUP_INSTRUCTIONS.md** file in this directory to:

1. ✅ Set up Docker Desktop for Windows with WSL2 integration
2. ✅ Clone/verify the repository on WSL
3. ✅ Mount FASTQ files from a remote ARM machine via SSHFS (to save disk space)
4. ✅ Download CellRanger 10.0.0 from 10x Genomics
5. ✅ Build the Docker image (~21GB, includes GRCh38-2024-A reference)
6. ✅ Configure for low-memory sequential processing (12GB per job)
7. ✅ Test with a single sample first
8. ✅ Process all 18 samples (18-36 hours runtime)
9. ✅ Merge outputs and compare with STAR results
10. ✅ Copy comparison results back to ARM machine

## Important Constraints

**System Limitations:**
- RAM: 16GB (can only process 1-2 samples at a time)
- Storage: 50GB available (must use network-mounted FASTQ files)
- Architecture: x86_64 (WSL on Windows, NOT ARM)

**Resource Strategy:**
- Use SSHFS to mount FASTQ files from ARM machine (saves ~60GB disk)
- Process samples sequentially, NOT in parallel (avoid OOM errors)
- Use `--localmem=12` and `--localcores=4` for CellRanger
- Expected total runtime: 18-36 hours for all 18 samples

## Key Files in This Directory

- **WSL_SETUP_INSTRUCTIONS.md** ← **YOUR PRIMARY GUIDE** (follow this step-by-step)
- `Dockerfile` - Docker image definition (CellRanger + Python environment)
- `build_docker.sh` - Builds the Docker image
- `run_docker.sh` - Launches the container with volume mounts
- `run_cellranger_batch.sh` - Processes all 18 samples (MUST edit for low-memory)
- `merge_cellranger_outputs.py` - Combines outputs into AnnData
- `compare_with_star.py` - Generates comparison CSV files
- `prepare_fastq_symlinks.sh` - Creates CellRanger-compatible FASTQ symlinks

## Information You'll Need From Me

Before starting, please ask me for:

1. **ARM machine connection details:**
   - IP address or hostname
   - Username (probably 'jacobc')
   - Whether SSH key authentication is set up

2. **Git repository URL:**
   - If this isn't already cloned on WSL

3. **CellRanger download:**
   - I'll need to get a fresh signed URL from https://www.10xgenomics.com/support/software/cell-ranger/downloads
   - (The URL expires, so I'll get it when you're ready)

## Expected Workflow

### Phase 1: Environment Setup (~30 min)
- Verify Docker Desktop is installed and integrated with WSL
- Install sshfs for network mounting
- Check system resources (RAM, disk, CPU cores)

### Phase 2: Repository & FASTQ Files (~15 min)
- Clone or verify git repository
- Set up SSH key authentication to ARM machine (if needed)
- Mount FASTQ files via sshfs: `/home/jacobc/sc/sra_downloads` → `~/sc/sra_downloads`
- Create CellRanger-compatible symlinks

### Phase 3: CellRanger & Docker (~35 min)
- Download CellRanger 10.0.0 (~828MB)
- Build Docker image (~21GB, takes 20-30 min)
- Verify image built successfully

### Phase 4: Test Single Sample (~2 hours)
- Launch container
- Verify CellRanger and FASTQ mounts
- Process ONE sample (SRR14513984) to validate setup
- Check outputs and disk/RAM usage

### Phase 5: Full Batch Processing (18-36 hours)
- **CRITICAL**: Edit `run_cellranger_batch.sh` to use `--localmem=12` and `--localcores=4`
- Start batch processing (can run in background with nohup)
- Monitor progress periodically
- Handle any interruptions (network drops, etc.)

### Phase 6: Analysis & Transfer (~1 hour)
- Merge 18 sample outputs into single AnnData file
- Compare CellRanger vs STAR gene expression
- Generate comparison CSV files
- Copy results back to ARM machine

## How to Help Me

### Verification Checkpoints
After each major phase, please:
- ✅ Confirm the step completed successfully
- 📊 Show disk space: `df -h ~`
- 💾 Show RAM usage: `free -h`
- ⚠️ Alert me to any warnings or errors
- 📝 Summarize what was done and what's next

### Problem Solving
If you encounter issues:
- 🔍 Check the "Troubleshooting" section in WSL_SETUP_INSTRUCTIONS.md
- 🔄 Suggest solutions based on the error
- 💡 Ask for clarification if needed (IP address, credentials, etc.)

### Long-Running Tasks
For the 18-36 hour batch processing:
- 🕒 Set it up to run in background with logging
- 📈 Provide commands to monitor progress
- 🔔 Explain how to resume if interrupted
- 💤 Make it safe to close terminal/disconnect

## Success Criteria

At the end, we should have:

1. ✅ All 18 samples processed by CellRanger
2. ✅ Merged AnnData file: `merged_cellranger_18samples.h5ad`
3. ✅ Three comparison CSV files:
   - `gene_expression_comparison.csv` (per-gene correlations)
   - `per_cell_comparison.csv` (per-cell UMI counts)
   - `summary_cellranger_vs_star.csv` (overall statistics)
4. ✅ Results copied back to ARM machine at:
   - `/home/jacobc/sc/morabito_18samples/cellranger_analysis/cellranger_output/`

## Communication Style

Please:
- 📋 Show me the exact commands before running them
- 🔍 Explain what each command does
- ⏱️ Give time estimates for long operations
- ⚠️ Warn me before operations that consume lots of disk/RAM
- ✅ Confirm each phase completion before moving on

## Ready to Start?

Once you've read this prompt and the WSL_SETUP_INSTRUCTIONS.md file, please:

1. Confirm you understand the task and constraints
2. Check which phase we should start at (environment setup, or is Docker already installed?)
3. Ask me for any information you need (ARM machine IP, git repo URL, etc.)
4. Show me the first commands you'd like to run

Let's get CellRanger running on WSL! 🚀

---

## Quick Reference

**ARM Machine (where data currently lives):**
- Location: `/home/jacobc/sc/sra_downloads/` (FASTQ files)
- Location: `/home/jacobc/sc/morabito_18samples/` (comparison files)

**WSL Machine (where we're working):**
- Target: `~/sc/morabito_18samples/cellranger_analysis/`
- Mount point: `~/sc/sra_downloads/` (via sshfs)

**Key Constraints:**
- RAM: 16GB → use `--localmem=12`
- Disk: 50GB → use network-mounted FASTQ
- Time: 18-36 hours → run in background

**Main Guide:** See WSL_SETUP_INSTRUCTIONS.md for detailed step-by-step instructions.
