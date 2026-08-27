#!/usr/bin/env python3
"""
sample_pol3_fasta.py

Takes a CSV of Pol III gene annotations + sequences, groups rows by
`derived_pol3_class`, randomly samples up to N sequences per group
(default 1000), excludes unclassified entries, and writes one FASTA
file per class.

Expected CSV columns (header row required):
    source,chr,start,end,name,seq,pol3_gene_type,gene_symbol,
    rna_class,derived_pol3_gene_type,derived_pol3_class,pol3_type_match

Usage:
    python sample_pol3_fasta.py input.csv -o output_dir -n 1000 --seed 42

Output:
    <output_dir>/<derived_pol3_class>.fasta   (one file per class)

FASTA header format:
    >gene_symbol|chr:start-end|derived_pol3_class

    If gene_symbol is missing/NA, falls back to the `name` column.
"""

import argparse
import csv
import random
import sys
from collections import defaultdict
from pathlib import Path

# Values in `derived_pol3_class` that count as "unclassified" and should be excluded.
UNCLASSIFIED_VALUES = {
    "", "na", "n/a", "none", "nan", "unclassified", "unknown", "unclass",
}


def is_unclassified(value: str) -> bool:
    return value is None or value.strip().lower() in UNCLASSIFIED_VALUES


def sanitize_filename(name: str) -> str:
    """Make a string safe to use as a filename."""
    keep = "-_." 
    return "".join(c if c.isalnum() or c in keep else "_" for c in name).strip("_") or "unnamed"


def build_header(row: dict) -> str:
    gene_symbol = (row.get("gene_symbol") or "").strip()
    if not gene_symbol or gene_symbol.lower() in ("na", "n/a", "none", "nan"):
        gene_symbol = (row.get("name") or "unknown").strip()

    chrom = (row.get("chr") or "NA").strip()
    start = (row.get("start") or "NA").strip()
    end = (row.get("end") or "NA").strip()
    gene_class = (row.get("derived_pol3_class") or "NA").strip()

    location = f"{chrom}:{start}-{end}"
    return f">{gene_symbol}|{location}|{gene_class}"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv_path", help="Path to input CSV file")
    parser.add_argument("-o", "--outdir", default="data/pol3_fasta_by_class", help="Output directory (default: fasta_by_class)")
    parser.add_argument("-n", "--num-seqs", type=int, default=1000, help="Number of random sequences per class (default: 1000)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    parser.add_argument("--seq-col", default="seq", help="Column name containing the sequence (default: seq)")
    parser.add_argument("--class-col", default="derived_pol3_class", help="Column name to group by (default: derived_pol3_class)")
    parser.add_argument("--wrap", type=int, default=0, help="Wrap sequence lines at N characters (0 = no wrap, single line per record)")
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    csv_path = Path(args.csv_path)
    if not csv_path.exists():
        sys.exit(f"Error: input file not found: {csv_path}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    groups = defaultdict(list)
    total_rows = 0
    excluded_unclassified = 0
    excluded_no_seq = 0

    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)

        # Sanity check columns
        missing_cols = [c for c in (args.class_col, args.seq_col) if c not in reader.fieldnames]
        if missing_cols:
            sys.exit(f"Error: expected column(s) not found in CSV: {missing_cols}\nFound columns: {reader.fieldnames}")

        for row in reader:
            total_rows += 1
            gene_class = row.get(args.class_col)

            if is_unclassified(gene_class):
                excluded_unclassified += 1
                continue

            seq = (row.get(args.seq_col) or "").strip()
            if not seq:
                excluded_no_seq += 1
                continue

            groups[gene_class.strip()].append(row)

    if not groups:
        sys.exit("Error: no classified rows with sequences found after filtering.")

    print(f"Read {total_rows} rows from {csv_path}")
    print(f"Excluded {excluded_unclassified} unclassified rows, {excluded_no_seq} rows with empty sequence")
    print(f"Found {len(groups)} classes: {', '.join(sorted(groups.keys()))}\n")

    for gene_class, rows in sorted(groups.items()):
        n_available = len(rows)
        n_sample = min(args.num_seqs, n_available)

        if n_available < args.num_seqs:
            print(f"  [WARNING] class '{gene_class}' has only {n_available} sequences "
                  f"(< requested {args.num_seqs}); using all of them.")

        sampled = random.sample(rows, n_sample)

        out_path = outdir / f"{sanitize_filename(gene_class)}.fasta"
        with open(out_path, "w", encoding="utf-8") as out_fh:
            for row in sampled:
                header = build_header(row)
                seq = row[args.seq_col].strip()

                out_fh.write(header + "\n")
                if args.wrap and args.wrap > 0:
                    for i in range(0, len(seq), args.wrap):
                        out_fh.write(seq[i:i + args.wrap] + "\n")
                else:
                    out_fh.write(seq + "\n")

        print(f"  Class '{gene_class}': wrote {n_sample} sequences -> {out_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()