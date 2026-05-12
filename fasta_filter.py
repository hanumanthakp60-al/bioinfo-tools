#!/usr/bin/env python
"""
fasta_filter.py — Filter FASTA sequences by length and GC content
Author: Hanumantha K P
Usage: python fasta_filter.py -i input.fasta -l 200 -g 40 -o filtered.fasta
"""

import argparse
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import pandas as pd

# ── Argument parser ───────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Filter FASTA sequences by length and GC content"
)
parser.add_argument("-i", "--input",   required=True,  help="Input FASTA file")
parser.add_argument("-o", "--output",  default="filtered.fasta", help="Output FASTA file")
parser.add_argument("-l", "--length",  type=int,   default=0,   help="Minimum sequence length")
parser.add_argument("-g", "--gc_min",  type=float, default=0.0, help="Minimum GC content %%")
parser.add_argument("-G", "--gc_max",  type=float, default=100.0, help="Maximum GC content %%")
parser.add_argument("-r", "--report",  default="report.csv", help="Output report CSV")
args = parser.parse_args()

# ── Filter sequences ──────────────────────────────────────────
print(f"Reading {args.input}...")
records = list(SeqIO.parse(args.input, "fasta"))
print(f"Total sequences loaded: {len(records)}")

passed  = []
failed  = []
report  = []

for record in records:
    seq    = str(record.seq).upper()
    length = len(seq)
    gc     = gc_fraction(record.seq) * 100

    passes_length = length >= args.length
    passes_gc     = args.gc_min <= gc <= args.gc_max

    report.append({
        "id":     record.id,
        "length": length,
        "gc":     round(gc, 2),
        "passed": passes_length and passes_gc
    })

    if passes_length and passes_gc:
        passed.append(record)
    else:
        failed.append(record)

# ── Save results ──────────────────────────────────────────────
SeqIO.write(passed, args.output, "fasta")

report_df = pd.DataFrame(report)
report_df.to_csv(args.report, index=False)

# ── Print summary ─────────────────────────────────────────────
print(f"\n{'='*50}")
print(f"FILTER RESULTS")
print(f"{'='*50}")
print(f"Total input:    {len(records)}")
print(f"Passed filter:  {len(passed)}")
print(f"Failed filter:  {len(failed)}")
print(f"Filters applied:")
print(f"  Min length:   {args.length} bp")
print(f"  GC range:     {args.gc_min}% — {args.gc_max}%")
print(f"\nOutput saved:   {args.output}")
print(f"Report saved:   {args.report}")
print(f"{'='*50}")
