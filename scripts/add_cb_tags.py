#!/usr/bin/env python3
"""
Add cell barcode (CB:Z:) tags to STAR-aligned BAM file
using barcodes from CellRanger ATAC FASTQ
"""

import pysam
import gzip
import sys

if len(sys.argv) != 4:
    print("Usage: python3 add_cb_tags.py <barcode_fastq.gz> <input.bam> <output.bam>")
    sys.exit(1)

barcode_file = sys.argv[1]
input_bam = sys.argv[2]
output_bam = sys.argv[3]

print(f'Reading barcodes from {barcode_file}...')
barcodes = []
with gzip.open(barcode_file, 'rt') as f:
    while True:
        header = f.readline()
        if not header:
            break
        seq = f.readline().strip()
        f.readline()  # +
        f.readline()  # quality
        barcodes.append(seq)

print(f'Loaded {len(barcodes):,} barcodes')
print(f'Adding CB:Z: tags to BAM...')

with pysam.AlignmentFile(input_bam, 'rb') as infile:
    with pysam.AlignmentFile(output_bam, 'wb', header=infile.header) as outfile:
        for i, read in enumerate(infile):
            if i < len(barcodes):
                read.set_tag('CB', barcodes[i], value_type='Z')
            outfile.write(read)
            
            if (i + 1) % 1000000 == 0:
                print(f'Processed {i+1:,} reads...', flush=True)

print(f'✓ Added CB tags to {i+1:,} reads')
