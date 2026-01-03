#!/bin/bash

# Transfer CB tags from CellRanger BAM to STAR multi-mapper BAM
# Matches reads by read ID, not position
#
# This gives us:
# - Multi-mapper support from STAR
# - Cell barcode QC from CellRanger

set -e

SAMPLE_ID="${1:-SRR14514130}"
CELLRANGER_BAM="$HOME/sc/cellranger_atac_output/$SAMPLE_ID/outs/possorted_bam.bam"
STAR_BAM="$HOME/sc/atacseq_aligned/$SAMPLE_ID/Aligned.sortedByCoord.out.bam"
OUTPUT_DIR="$HOME/sc/atacseq_aligned/$SAMPLE_ID"
THREADS=16

echo "=========================================="
echo "Transfer CB Tags from CellRanger to STAR BAM"
echo "=========================================="
echo "Sample: $SAMPLE_ID"
echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Step 1: Extract read ID → CB mapping from CellRanger BAM
echo "[Step 1] Extracting CB tags from CellRanger BAM..."
echo "  This creates a lookup: ReadID → CB tag"

samtools view $CELLRANGER_BAM | \
    awk '{
        readid = $1
        cb = ""
        for(i=12; i<=NF; i++) {
            if($i ~ /^CB:Z:/) {
                cb = substr($i, 6)
                break
            }
        }
        if(cb != "") print readid, cb
    }' | \
    gzip > "$OUTPUT_DIR/readid_to_cb.txt.gz"

CB_COUNT=$(zcat "$OUTPUT_DIR/readid_to_cb.txt.gz" | wc -l)
echo "  ✓ Extracted CB tags for $(printf "%'d" $CB_COUNT) reads"
echo ""

# Step 2: Transfer CB tags to STAR BAM
echo "[Step 2] Adding CB tags to STAR BAM..."
echo "  Matching reads by ID..."
echo "  This will take 30-60 minutes..."

python3 << PYTHON_SCRIPT
#!/usr/bin/env python3
import pysam
import gzip
import sys

readid_cb_file = "$OUTPUT_DIR/readid_to_cb.txt.gz"
star_bam = "$STAR_BAM"
output_bam = "$OUTPUT_DIR/star_with_cb.bam"

print("  Loading read ID → CB lookup table...")
readid_to_cb = {}
with gzip.open(readid_cb_file, 'rt') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            readid_to_cb[parts[0]] = parts[1]

print(f"  Loaded {len(readid_to_cb):,} read IDs with CB tags")
print(f"  Processing STAR BAM...")

matched = 0
unmatched = 0

with pysam.AlignmentFile(star_bam, "rb") as infile:
    with pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:
        for i, read in enumerate(infile):
            # Look up CB tag by read ID
            if read.query_name in readid_to_cb:
                cb = readid_to_cb[read.query_name]
                read.set_tag('CB', cb, value_type='Z')
                matched += 1
            else:
                unmatched += 1
            
            outfile.write(read)
            
            if (i + 1) % 5000000 == 0:
                print(f"  Processed {i+1:,} reads... (matched: {matched:,}, unmatched: {unmatched:,})", flush=True)

print(f"  ✓ Total reads processed: {i+1:,}")
print(f"  ✓ Matched with CB tags: {matched:,} ({100*matched/(i+1):.1f}%)")
print(f"  ✓ No CB tag found: {unmatched:,} ({100*unmatched/(i+1):.1f}%)")
PYTHON_SCRIPT

echo ""
echo "[Step 3] Filtering for reads with CB tags..."
samtools view -h -d CB "$OUTPUT_DIR/star_with_cb.bam" | \
    samtools view -b - | \
    samtools sort -@ $THREADS -o "$OUTPUT_DIR/star_multimapper_final.bam" -

samtools index "$OUTPUT_DIR/star_multimapper_final.bam"

echo ""
echo "[Step 4] Final statistics..."
STAR_TOTAL=$(samtools view -c "$STAR_BAM")
FINAL_TOTAL=$(samtools view -c "$OUTPUT_DIR/star_multimapper_final.bam")
FINAL_SIZE=$(stat -c%s "$OUTPUT_DIR/star_multimapper_final.bam" | numfmt --to=iec-i --suffix=B)

echo "  STAR multi-mapper reads: $(printf "%'d" $STAR_TOTAL)"
echo "  With valid CB tags: $(printf "%'d" $FINAL_TOTAL)"
echo "  Retention: $(awk "BEGIN {printf \"%.1f%%\", 100*$FINAL_TOTAL/$STAR_TOTAL}")"
echo "  Output: star_multimapper_final.bam ($FINAL_SIZE)"
echo ""

echo "=========================================="
echo "✓ Complete!"
echo "=========================================="
echo "Output BAM: $OUTPUT_DIR/star_multimapper_final.bam"
echo "This BAM has:"
echo "  ✓ Multi-mapper support (STAR alignment)"
echo "  ✓ Cell barcode tags from CellRanger QC"
echo ""
echo "Ready for scTE/scTEATAC!"
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"
