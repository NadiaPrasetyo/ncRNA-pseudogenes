#!/usr/bin/env python3
"""
fetch_transcript_seq.py
========================

Given a CSV of genomic Pol III loci, this script pulls sequence from a
local GENCODE genome FASTA, annotates gene symbols from a local GENCODE
GTF, and assigns each locus a Pol III promoter class (I / II / III) by a
simple, deterministic rule based on the RNA class/name -- NOT by sequence
clustering or reference-track overlap.

Classification rule
--------------------
Based on Oler et al. 2010 (Nat Struct Mol Biol 17:620-628), Table 1, which
groups every Pol III-transcribed RNA class under one of three promoter
types:

    Class I   (Type 1): 5S rRNA (and 5S-related)
    Class II  (Type 2): tRNA (and tRNA-related / tRNA pseudogene, but NOT
                         tRNA-SeC), Alu, MIR, snaR (and snaR-related),
                         HVG or Vault (and HVG-related), 7SL or SRP
                         (and 7SL-related), BC200 (and BC200-related)
    Class III (Type 3): Y (and Y-related, including the "HY#" naming
                         convention, e.g. HY1/HY3/HY5), U6 (and U6-related),
                         7SK (and 7SK-related), RNase P (and RNase
                         P-related), RNase MRP (and RNase MRP-related),
                         U6atac (and U6atac-related), tRNA-SeC

The one explicit exception: tRNA-SeC (the selenocysteine tRNA) is Class
III, not Class II, despite being a "tRNA".

Each locus's `name` is matched, in priority order, against regex patterns
for each RNA class (see CLASS_RULES below) to assign an `rna_class` label
and a `derived_pol3_class` (I/II/III). This is compared against the
pre-existing `pol3_gene_type` column (if present) to flag any mismatch.

Inputs expected
----------------
  - Genome FASTA: "Genome sequence (GRCh38.p14), ALL regions" from GENCODE
    release 50. Contig headers look like ">chr1 1" -- pyfaidx indexes on
    the first whitespace-delimited token, so it keys on "chr1", matching
    the `chr` column in the input CSV and GENCODE's own GTF/GFF3 naming.
  - GTF (--gtf): GENCODE main annotation GTF, for gene_symbol lookup.

What it does
------------
  1. SEQUENCE: random-access lookup into the genome FASTA via pyfaidx.
     Coordinates are treated as 1-based inclusive, matching the input CSV.
     Loci for which no real sequence could be recovered -- fetch failure,
     an empty/NaN result, or a sequence made up entirely of N's (e.g. an
     assembly gap) -- are dropped from the output.
  2. GENE SYMBOL: IntervalTree built once per chromosome from the GTF's
     `gene` features; each locus is overlap-queried against it.
  3. CLASSIFICATION: `name` is matched against CLASS_RULES to derive
     `rna_class` and `derived_pol3_class` (I/II/III), compared against
     any existing `pol3_gene_type` column.
  4. OUTPUT: all input columns minus `classification_basis`, plus: seq,
     gene_symbol, rna_class, derived_pol3_class, derived_pol3_gene_type
     (numeric 1/2/3), pol3_type_match (bool, only if the input had a
     pol3_gene_type column).

Usage
-----
    python pol3_class_annotate.py input.csv output.csv \\
        --genome-fasta data/downloads/GRCh38.p14.genome.fa \\
        --gtf data/downloads/gencode.v50.annotation.gtf.gz

Requirements
------------
    pip install pandas numpy pyfaidx intervaltree
"""

import argparse
import gzip
import re
import sys

import pandas as pd

# ---------------------------------------------------------------------
# Sequence retrieval (local GENCODE genome FASTA)
# ---------------------------------------------------------------------

def open_genome(fasta_path: str):
    """pyfaidx keys on the first whitespace-delimited token of each header,
    e.g. '>chr1 1' -> 'chr1', matching GENCODE's own chromosome naming."""
    from pyfaidx import Fasta
    return Fasta(fasta_path, sequence_always_upper=True)


def fetch_sequence_local(fa, chrom: str, start: int, end: int) -> str:
    """1-based inclusive coordinates, matching the input CSV convention.
    Returns an empty string (rather than raising) if the chromosome can't
    be found, so the caller can drop the row instead of crashing the run."""
    key = chrom if chrom in fa.keys() else chrom.replace("chr", "")
    if key not in fa.keys():
        print(
            f"  WARNING: chromosome '{chrom}' not found in genome FASTA "
            f"-- no sequence recovered for this locus.",
            file=sys.stderr,
        )
        return ""
    try:
        return str(fa[key][start - 1:end].seq).upper()
    except Exception as e:
        print(
            f"  WARNING: failed to fetch sequence for {chrom}:{start}-{end} ({e}) "
            f"-- no sequence recovered for this locus.",
            file=sys.stderr,
        )
        return ""


# ---------------------------------------------------------------------
# Gene symbol lookup (local GTF, interval-tree indexed for speed)
# ---------------------------------------------------------------------

ATTR_KEYS_GTF = re.compile(r'(\w+) "([^"]+)"')
ATTR_KEYS_GFF3 = re.compile(r'(\w+)=([^;]+)')

def parse_attributes(attr_field: str) -> dict:
    d = dict(ATTR_KEYS_GTF.findall(attr_field))
    if not d:
        d = dict(ATTR_KEYS_GFF3.findall(attr_field))
    return d


def build_gene_interval_index(gtf_path: str):
    from intervaltree import IntervalTree

    trees = {}
    opener = gzip.open if gtf_path.endswith(".gz") else open
    n = 0
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom = fields[0]
            try:
                f_start, f_end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            attrs = parse_attributes(fields[8])
            gene_name = attrs.get("gene_name")
            if gene_name is None:
                continue
            trees.setdefault(chrom, IntervalTree())
            trees[chrom].addi(f_start, f_end + 1, gene_name)
            n += 1

    print(f"  Indexed {n} gene features across {len(trees)} chromosomes from {gtf_path}",
          file=sys.stderr)
    return trees


def lookup_gene_symbol(trees: dict, chrom: str, start: int, end: int) -> str:
    tree = trees.get(chrom)
    if tree is None:
        return "NA"
    hits = tree.overlap(start, end + 1)
    if not hits:
        return "NA"
    names = sorted({iv.data for iv in hits})
    return ";".join(names)


# ---------------------------------------------------------------------
# Pol III class rules, derived from Oler et al. 2010 Table 1
# ---------------------------------------------------------------------
# Ordered list of (regex, rna_class_label, pol3_class_numeric). Matched
# against the `name` column, case-insensitively, in order -- first match
# wins. tRNA-SeC MUST be checked before the generic tRNA rule, since it is
# the one explicit exception (Class III, not Class II).
#
# Patterns cover:
#   - the informal TRNA<isotype><copy#> naming used in this dataset
#     (e.g. TRNAK20), including names where "tRNA" is embedded rather
#     than at the start of the string (e.g. "chr1_15_tRNAGln_CTG+")
#   - common official GENCODE/HGNC symbols (e.g. RNA5S1, RNU6-1, RN7SL1,
#     BCYRN1)
#   - "5S_rRNA-related"-style names where an underscore/hyphen follows
#     "5S" rather than a plain word boundary
#   - the "HY#" naming convention for Y RNA (e.g. HY1, HY3, HY5), in
#     addition to RNY# and Y-RNA/Y_RNA style names
# so the rule set generalizes beyond tRNA-only inputs.

CLASS_RULES = [
    # --- Class III exception: tRNA-SeC must be checked before generic tRNA ---
    (re.compile(r'TRNA.{0,3}SEC', re.IGNORECASE), "tRNA-SeC", 3),

    # --- Class I: 5S rRNA (and 5S-related / pseudogene) ---
    (re.compile(r'^RNA5S', re.IGNORECASE), "5S rRNA", 1),
    # Matches "5S", "5S_rRNA-related", "5S-rRNA", "5S rRNA", etc. -- any
    # name starting with "5S" is treated as 5S rRNA / 5S-related.
    (re.compile(r'^5S', re.IGNORECASE), "5S rRNA", 1),

    # --- Class II: HVG or Vault --- (checked BEFORE the generic tRNA rule
    # below, since "VTRNA..." contains "TRNA" as a substring and would
    # otherwise be misclassified as plain tRNA)
    (re.compile(r'VTRNA|HVG', re.IGNORECASE), "HVG or Vault", 2),

    # --- Class II: tRNA (all other isotypes, consensus/related/pseudogene) ---
    # Matches "TRNA..." at the start (e.g. TRNAK20) as well as "tRNA"
    # embedded elsewhere in the name (e.g. "chr1_15_tRNAGln_CTG+").
    (re.compile(r'TRNA', re.IGNORECASE), "tRNA", 2),

    # --- Class II: Alu --- (repeat names like AluY, AluSx, AluYa5 have a
    # letter immediately after "Alu", so only a left boundary is checked)
    (re.compile(r'(^|[^A-Za-z])ALU', re.IGNORECASE), "Alu", 2),

    # --- Class II: MIR ---
    (re.compile(r'^MIR', re.IGNORECASE), "MIR", 2),

    # --- Class II: snaR ---
    (re.compile(r'SNAR', re.IGNORECASE), "snaR", 2),

    # --- Class II: 7SL or SRP ---
    (re.compile(r'7SL', re.IGNORECASE), "7SL or SRP", 2),

    # --- Class II: BC200 ---
    (re.compile(r'^BCYRN1$|BC200', re.IGNORECASE), "BC200", 2),

    # --- Class III: Y RNA --- (RNY#, "Y-RNA"/"Y_RNA" style names, and the
    # "HY#" naming convention, e.g. HY1, HY3, HY5)
    (re.compile(r'^RNY\d|Y[-_]?RNA|^HY\d', re.IGNORECASE), "Y", 3),

    # --- Class III: U6atac (check before generic U6) ---
    (re.compile(r'U6ATAC', re.IGNORECASE), "U6atac", 3),

    # --- Class III: U6 ---
    (re.compile(r'(^|[^A-Z])U6([^A-Z]|$)|^RNU6', re.IGNORECASE), "U6", 3),

    # --- Class III: 7SK ---
    (re.compile(r'7SK', re.IGNORECASE), "7SK", 3),

    # --- Class III: RNase P ---
    (re.compile(r'^RPPH1$|RNASE.?P\b', re.IGNORECASE), "RNase P", 3),

    # --- Class III: RNase MRP ---
    (re.compile(r'^RMRP$|RNASE.?MRP', re.IGNORECASE), "RNase MRP", 3),
]

POL3_CLASS_ROMAN = {1: "I", 2: "II", 3: "III"}


def classify_rna(name: str):
    """Return (rna_class_label, pol3_class_numeric) or ('Unclassified', None)
    if no rule matches."""
    name = str(name).strip()
    for pattern, label, class_num in CLASS_RULES:
        if pattern.search(name):
            return label, class_num
    return "Unclassified", None


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input_csv", default="data/literature_sourced/dedup_literature_polIII_transcripts_classified.csv",)
    ap.add_argument("-o", "--output_csv", default="data/literature_sourced/literature_polIII_transcript_seq.csv",)
    ap.add_argument("--genome-fasta", required=True,
                     help="Local GENCODE genome FASTA (GRCh38.p14.genome.fa)")
    ap.add_argument("--gtf", default=None,
                     help="Local GENCODE GTF, plain or .gz -- enables gene_symbol lookup")
    args = ap.parse_args()

    df = pd.read_csv(args.input_csv)
    required = {"source", "chr", "start", "end", "name", "seq"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"Input CSV is missing required columns: {missing}")

    print(f"Loading genome FASTA: {args.genome_fasta} ...", file=sys.stderr)
    fa = open_genome(args.genome_fasta)

    gene_trees = None
    if args.gtf:
        print(f"Indexing gene GTF: {args.gtf} ...", file=sys.stderr)
        gene_trees = build_gene_interval_index(args.gtf)

    # --- fetch sequences ---
    print(f"Fetching sequence for {len(df)} loci ...", file=sys.stderr)
    seqs = []
    for i, row in enumerate(df.itertuples(index=False)):
        existing = getattr(row, "seq", None)
        if isinstance(existing, str) and existing.strip():
            seqs.append(existing.strip().upper())
            continue
        seq = fetch_sequence_local(fa, row.chr, int(row.start), int(row.end))
        seqs.append(seq)
        if (i + 1) % 1000 == 0:
            print(f"  ... {i + 1}/{len(df)} sequences fetched", file=sys.stderr)
    df["seq"] = seqs

    # --- drop rows with no sequence recovered (empty/NaN, or all-N) ---
    n_before = len(df)
    seq_str = df["seq"].astype(str).str.strip()
    empty_mask = df["seq"].isna() | (seq_str == "")
    # all-N sequences (any run length, any case) indicate an assembly gap
    # or other unrecovered region, not real sequence
    all_n_mask = seq_str.str.fullmatch(r'[Nn]+', na=False)
    no_seq_mask = empty_mask | all_n_mask
    n_dropped = int(no_seq_mask.sum())
    if n_dropped:
        n_empty = int(empty_mask.sum())
        n_all_n = int((all_n_mask & ~empty_mask).sum())
        print(f"Dropping {n_dropped}/{n_before} loci with no sequence recovered "
              f"({n_empty} empty/NaN, {n_all_n} all-N).", file=sys.stderr)
        df = df.loc[~no_seq_mask].reset_index(drop=True)

    # --- gene symbols ---
    if gene_trees is not None:
        print("Looking up overlapping gene symbols ...", file=sys.stderr)
        df["gene_symbol"] = [
            lookup_gene_symbol(gene_trees, row.chr, int(row.start), int(row.end))
            for row in df.itertuples(index=False)
        ]
    else:
        df["gene_symbol"] = "NA"
        print("No --gtf given: gene_symbol will be 'NA' for all rows", file=sys.stderr)

    # --- rule-based Pol III class assignment ---
    print("Classifying loci by RNA name/class ...", file=sys.stderr)
    classes = [classify_rna(n) for n in df["name"]]
    df["rna_class"] = [c[0] for c in classes]
    df["derived_pol3_gene_type"] = [c[1] for c in classes]
    df["derived_pol3_class"] = [
        POL3_CLASS_ROMAN.get(c[1], "Unclassified") for c in classes
    ]

    if "pol3_gene_type" in df.columns:
        df["pol3_type_match"] = [
            (pd.notna(orig) and derived is not None and int(orig) == derived)
            for orig, derived in zip(df["pol3_gene_type"], df["derived_pol3_gene_type"])
        ]

    df = df.drop(columns=["classification_basis"], errors="ignore")

    df.to_csv(args.output_csv, index=False)
    print(f"\nWrote {len(df)} rows to {args.output_csv}", file=sys.stderr)

    # --- summary ---
    print("\n--- Summary ---", file=sys.stderr)
    print(df["rna_class"].value_counts().to_string(), file=sys.stderr)
    n_unclassified = (df["rna_class"] == "Unclassified").sum()
    if n_unclassified:
        print(f"\n{n_unclassified} loci did not match any rule (rna_class = 'Unclassified') "
              f"-- inspect their `name` values and extend CLASS_RULES if needed.",
              file=sys.stderr)
    if "pol3_type_match" in df.columns:
        n_mismatch = (~df["pol3_type_match"]).sum()
        print(f"\n{n_mismatch}/{len(df)} loci have a derived class that disagrees with "
              f"(or couldn't be compared to) the input pol3_gene_type column.", file=sys.stderr)


if __name__ == "__main__":
    main()