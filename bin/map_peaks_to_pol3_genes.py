#!/usr/bin/env python3
"""
map_peaks_to_pol3_genes.py

Pipeline:
  1. Read a CSV of Pol III genes with HG18 coordinates.
  2. LiftOver each gene's coordinates from HG18 -> HG38 (via pyliftover +
     a UCSC "hg18ToHg38.over.chain(.gz)" chain file).
  3. Read one or more ChIP-seq peak TSV files (HG38, BED-like: chr, start,
     end, name, ...extra columns ignored).
  4. Read a featureCounts-style RNA-seq CSV/TSV (HG38). Each gene row has
     semicolon-joined per-exon Chr/Start/End/Strand lists; the gene's
     genomic footprint is taken as (min start, max end) on its chromosome.
  5. For every lifted-over Pol III gene locus, count how many ChIP peaks
     and how many RNA-seq gene footprints overlap it.
  6. Write:
       - an annotated genes CSV (original fields + CHIP_peak_count +
         RNAseq_peak_count, coordinates updated to HG38)
       - an "unmapped peaks" CSV listing every ChIP peak / RNA-seq entry
         that did not overlap ANY gene in the input list.

Assumptions (adjust if your files differ):
  * All coordinates (gene CSV, ChIP TSV, RNA-seq CSV) are 0-based,
    half-open BED-style intervals ([start, end)). If your inputs are
    1-based, subtract 1 from `start` after loading, or pass --one-based.
  * ChIP TSV files have NO header row and are tab-delimited, with the
    first four columns being: chr, start, end, name. Any additional
    columns are ignored except that they are scanned for an Ensembl
    gene ID (e.g. "ENSG00000237613.3") to use as `expected_gene_id`
    for peaks that end up unmapped.
  * The RNA-seq file is featureCounts output: comment lines starting
    with '#', then a header row with columns including
    Geneid, Chr, Start, End, Strand, gene_name (extra count columns
    are ignored for the purposes of this script).

Dependencies:
    pip install pyliftover intervaltree pandas --break-system-packages

Chain file:
    Download hg18ToHg38.over.chain.gz from UCSC, e.g.:
    https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz

Usage:
    python map_peaks_to_pol3_genes.py \
        --genes pol3_genes_hg18.csv \
        --chip chip_peaks_1.tsv chip_peaks_2.tsv \
        --rnaseq gene_counts.tsv \
        --chain hg18ToHg38.over.chain.gz \
        --out-genes pol3_genes_hg38_annotated.csv \
        --out-unmapped unmapped_peaks.csv
"""

import argparse
import csv
import re
import sys
from collections import defaultdict

import pandas as pd
from intervaltree import IntervalTree
from pyliftover import LiftOver

GENE_FIELDS = [
    "source", "chr", "start", "end", "name", "seq",
    "pol3_gene_type", "gene_symbol", "rna_class",
    "derived_pol3_gene_type", "derived_pol3_class", "pol3_type_match",
]

ENSG_RE = re.compile(r"ENSG\d+(\.\d+)?")


# --------------------------------------------------------------------------
# 1. Load Pol III gene list (HG18)
# --------------------------------------------------------------------------
def load_genes(genes_csv):
    genes = []
    with open(genes_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [f for f in GENE_FIELDS if f not in reader.fieldnames]
        if missing:
            sys.exit(f"ERROR: genes CSV is missing expected columns: {missing}")
        for row in reader:
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            genes.append(row)
    return genes


# --------------------------------------------------------------------------
# 2. LiftOver HG18 -> HG38
# --------------------------------------------------------------------------
def liftover_genes(genes, chain_file):
    lo = LiftOver(chain_file)
    lifted = []
    n_failed = 0
    for gene in genes:
        chrom, start, end = gene["chr"], gene["start"], gene["end"]
        res_start = lo.convert_coordinate(chrom, start)
        res_end = lo.convert_coordinate(chrom, max(start, end - 1))
        if not res_start or not res_end:
            n_failed += 1
            continue
        c1, p1, strand1, _ = res_start[0]
        c2, p2, strand2, _ = res_end[0]
        if c1 != c2 or strand1 != strand2:
            n_failed += 1
            continue
        new_gene = dict(gene)
        new_gene["chr"] = c1
        new_gene["start"] = min(p1, p2)
        new_gene["end"] = max(p1, p2) + 1
        lifted.append(new_gene)
    if n_failed:
        print(f"WARNING: {n_failed} gene(s) failed to lift over and were dropped.",
              file=sys.stderr)
    return lifted


# --------------------------------------------------------------------------
# Build per-chromosome interval trees for the lifted gene list
# --------------------------------------------------------------------------
def build_gene_trees(genes):
    trees = defaultdict(IntervalTree)
    for idx, gene in enumerate(genes):
        start, end = gene["start"], gene["end"]
        if end <= start:
            end = start + 1
        trees[gene["chr"]][start:end] = idx
    return trees


def overlapping_gene_indices(trees, chrom, start, end):
    if end <= start:
        end = start + 1
    tree = trees.get(chrom)
    if tree is None:
        return []
    return [iv.data for iv in tree.overlap(start, end)]


# --------------------------------------------------------------------------
# 3. Load ChIP-seq peak TSV(s): chr, start, end, name, ...
#
# Peak calls (e.g. narrowPeak) occupy the first 10 columns:
#   chr, start, end, peak_name, score, strand, signalValue, pValue, qValue, offset
# Some peak files have an additional "nearest/overlapping gene" annotation
# appended after that (e.g. from `bedtools closest`), occupying 8 more
# columns:
#   chr, start, end, gene_symbol, score(.), strand, gene_id, distance
# If that annotation block is present, its gene_symbol/gene_id are used as
# `expected_gene` / `expected_gene_id` for peaks that end up unmapped.
# Peaks with no annotation block simply get an empty expected_gene/_id.
# --------------------------------------------------------------------------
def load_chip_peaks(chip_files):
    peaks = []
    for path in chip_files:
        with open(path, newline="") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 4:
                    continue
                chrom, start, end = fields[0], fields[1], fields[2]
                # skip an accidental header row
                if start.lower() in ("start", "chromstart"):
                    continue

                gene_symbol = ""
                gene_id = ""
                if len(fields) > 10:
                    annot = fields[10:]
                    if len(annot) >= 4 and annot[3] not in (".", ""):
                        gene_symbol = annot[3]
                    for extra in annot:
                        m = ENSG_RE.search(extra)
                        if m:
                            gene_id = m.group(0)
                            break

                peaks.append({
                    "chr": chrom,
                    "start": int(start),
                    "end": int(end),
                    "expected_gene": gene_symbol,
                    "expected_gene_id": gene_id,
                })
    return peaks


# --------------------------------------------------------------------------
# 4. Load featureCounts-style RNA-seq CSV/TSV
# --------------------------------------------------------------------------
def load_rnaseq_peaks(rnaseq_file):
    df = pd.read_csv(rnaseq_file, sep=None, engine="python", comment="#")
    required = {"Geneid", "Chr", "Start", "End"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: RNA-seq file is missing expected columns: {missing}")
    gene_name_col = "gene_name" if "gene_name" in df.columns else None

    peaks = []
    for row in df.itertuples(index=False):
        row = row._asdict() if hasattr(row, "_asdict") else dict(zip(df.columns, row))
        chrom_field = str(row["Chr"]).split(";")
        chrom = chrom_field[0]
        starts = [int(x) for x in str(row["Start"]).split(";")]
        ends = [int(x) for x in str(row["End"]).split(";")]
        start = min(starts)
        end = max(ends)
        peaks.append({
            "chr": chrom,
            "start": start,
            "end": end,
            "expected_gene": row.get(gene_name_col, "") if gene_name_col else "",
            "expected_gene_id": row["Geneid"],
        })
    return peaks


# --------------------------------------------------------------------------
# 5. Map peaks onto genes
# --------------------------------------------------------------------------
def annotate_genes_with_peaks(genes, chip_peaks, rnaseq_peaks):
    trees = build_gene_trees(genes)
    chip_counts = [0] * len(genes)
    rnaseq_counts = [0] * len(genes)
    unmapped = []

    for peak in chip_peaks:
        hits = overlapping_gene_indices(trees, peak["chr"], peak["start"], peak["end"])
        if hits:
            for idx in hits:
                chip_counts[idx] += 1
        else:
            unmapped.append({
                "peak_category": "chip_seq",
                "chr": peak["chr"],
                "start": peak["start"],
                "end": peak["end"],
                "expected_gene": peak["expected_gene"],
                "expected_gene_id": peak["expected_gene_id"],
            })

    for peak in rnaseq_peaks:
        hits = overlapping_gene_indices(trees, peak["chr"], peak["start"], peak["end"])
        if hits:
            for idx in hits:
                rnaseq_counts[idx] += 1
        else:
            unmapped.append({
                "peak_category": "rna_seq",
                "chr": peak["chr"],
                "start": peak["start"],
                "end": peak["end"],
                "expected_gene": peak["expected_gene"],
                "expected_gene_id": peak["expected_gene_id"],
            })

    for idx, gene in enumerate(genes):
        gene["CHIP_peak_count"] = chip_counts[idx]
        gene["RNAseq_peak_count"] = rnaseq_counts[idx]

    return genes, unmapped


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------
def write_genes_csv(genes, out_path):
    fieldnames = GENE_FIELDS + ["CHIP_peak_count", "RNAseq_peak_count"]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for gene in genes:
            writer.writerow({k: gene.get(k, "") for k in fieldnames})


def write_unmapped_csv(unmapped, out_path):
    fieldnames = ["peak_category", "chr", "start", "end",
                  "expected_gene", "expected_gene_id"]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in unmapped:
            writer.writerow(row)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", required=True, help="Pol III gene CSV (HG18 coords)")
    ap.add_argument("--chip", required=True, nargs="+",
                     help="One or more ChIP-seq peak TSV files (HG38)")
    ap.add_argument("--rnaseq", required=True,
                     help="featureCounts-style RNA-seq CSV/TSV (HG38)")
    ap.add_argument("--chain", required=True,
                     help="hg18ToHg38.over.chain(.gz) LiftOver chain file")
    ap.add_argument("--out-genes", default="pol3_genes_hg38_annotated.csv",
                     help="Output annotated genes CSV")
    ap.add_argument("--out-unmapped", default="unmapped_peaks.csv",
                     help="Output CSV of peaks that did not overlap any gene")
    args = ap.parse_args()

    print("Loading gene list (HG18)...")
    genes_hg18 = load_genes(args.genes)
    print(f"  {len(genes_hg18)} genes loaded")

    print("Lifting over genes HG18 -> HG38...")
    genes_hg38 = liftover_genes(genes_hg18, args.chain)
    print(f"  {len(genes_hg38)} genes successfully lifted over")

    print("Loading ChIP-seq peaks...")
    chip_peaks = load_chip_peaks(args.chip)
    print(f"  {len(chip_peaks)} ChIP-seq peaks loaded")

    print("Loading RNA-seq gene footprints...")
    rnaseq_peaks = load_rnaseq_peaks(args.rnaseq)
    print(f"  {len(rnaseq_peaks)} RNA-seq entries loaded")

    print("Mapping peaks onto gene loci...")
    genes_annotated, unmapped = annotate_genes_with_peaks(
        genes_hg38, chip_peaks, rnaseq_peaks
    )

    write_genes_csv(genes_annotated, args.out_genes)
    write_unmapped_csv(unmapped, args.out_unmapped)

    n_chip_unmapped = sum(1 for u in unmapped if u["peak_category"] == "chip_seq")
    n_rna_unmapped = sum(1 for u in unmapped if u["peak_category"] == "rna_seq")
    print(f"Done.\n  Genes written: {len(genes_annotated)} -> {args.out_genes}\n"
          f"  Unmapped ChIP peaks: {n_chip_unmapped}\n"
          f"  Unmapped RNA-seq entries: {n_rna_unmapped}\n"
          f"  -> {args.out_unmapped}")


if __name__ == "__main__":
    main()