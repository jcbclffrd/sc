# Manual CellRanger Download Instructions

## Why Manual Download is Required

10x Genomics requires users to:
1. Create a free account
2. Accept their End User License Agreement (EULA)
3. Download through their authenticated portal

Direct wget/curl downloads fail with "403 Forbidden" without authentication.

## Download Steps

### 1. Visit 10x Genomics Download Page

Go to: **https://www.10xgenomics.com/support/software/cell-ranger/downloads**

### 2. Sign In or Create Account

- Click "Sign In" (top right)
- Or create a free account if you don't have one

### 3. Accept EULA

- Read and accept the End User License Agreement
- This is required to download any 10x software

### 4. Download CellRanger 8.0.1

Find "Cell Ranger 8.0.1" in the list:
- Select: **"Linux (tar.gz)"** 
- Click "Download"
- File: `cellranger-8.0.1.tar.gz` (~2.0 GB)

### 5. Move File to Build Directory

```bash
# Move downloaded file to cellranger_analysis directory
mv ~/Downloads/cellranger-8.0.1.tar.gz /home/jacobc/sc/morabito_18samples/cellranger_analysis/

# Verify file size (should be ~2GB)
ls -lh /home/jacobc/sc/morabito_18samples/cellranger_analysis/cellranger-8.0.1.tar.gz
```

### 6. Verify File Integrity

```bash
cd /home/jacobc/sc/morabito_18samples/cellranger_analysis
tar -tzf cellranger-8.0.1.tar.gz | head -10
```

Should show files like:
```
cellranger-8.0.1/
cellranger-8.0.1/bin/
cellranger-8.0.1/bin/cellranger
...
```

### 7. Build Docker Image

Once the file is in place:

```bash
cd /home/jacobc/sc/morabito_18samples/cellranger_analysis
bash build_docker.sh
```

## Alternative: Use Existing CellRanger Installation

If you already have CellRanger installed on an x86_64 machine, you can:

1. Skip Docker entirely
2. Run `run_cellranger_batch.sh` directly (edit paths)
3. Or copy the tarball from that machine:

```bash
# From machine with CellRanger
scp /path/to/cellranger-8.0.1.tar.gz user@spark-bd86:/home/jacobc/sc/morabito_18samples/cellranger_analysis/
```

## Troubleshooting

### File is HTML instead of tar.gz

If the file opens as text/HTML, the download failed:
- Make sure you're signed in
- Accept EULA before downloading
- Use the "Download" button, not right-click save

### File is too small (<2GB)

Expected size: ~2.0 GB (2,000,000,000 bytes)

If smaller, download was incomplete:
- Check your internet connection
- Try downloading again

### Cannot access 10x website

If you cannot access the website:
- Check if your institution blocks it
- Try from a different network
- Contact 10x support: support@10xgenomics.com

## Need Help?

10x Genomics Support: https://www.10xgenomics.com/support
