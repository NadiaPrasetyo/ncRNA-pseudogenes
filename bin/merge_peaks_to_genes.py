#!/usr/bin/env python3
"""
merge_peaks_to_genes.py
========================

1) Merges ChIP-seq peaks that are within ~N bp of each other (default
   100 bp) into larger "peak regions".
2) Assigns each merged region to the gene(s) it overlaps -- or, if it
   doesn't overlap any gene, the nearest gene -- using a GENCODE
   "comprehensive" GTF annotation.
3) Optionally attaches RNA-seq expression (featureCounts output) for the
   assigned gene(s).

--------------------------------------------------------------------------
INPUT FILES
--------------------------------------------------------------------------

1) ChIP-seq peak file(s) (TSV, no header). The script auto-detects two
   layouts:

   a) Plain narrowPeak-style, 10 columns:
      chr  start  end  name  score  strand  signalValue  pValue  qValue  summit

   b) Peak-joined-to-gene, 18 columns (e.g. bedtools closest output) --
      only the first 10 columns (the peak itself) are used; the embedded
      gene columns are ignored, since gene assignment is now done
      properly from the GTF instead.

2) GENCODE GTF file (e.g. "gencode.v50.annotation.gtf", comprehensive /
   ALL). Only "gene" feature rows are used. gene_id and gene_name are
   pulled from the attributes column.

3) (optional) RNA-seq featureCounts output TSV:
   Geneid  Chr  Start  End  Strand  Length  gene_name  <sample1>  <sample2> ...

--------------------------------------------------------------------------
OUTPUT
--------------------------------------------------------------------------

TSV, one row per merged peak region:

   region_chr, region_start, region_end, n_peaks, peak_names,
   max_signalValue, gene_name, gene_id, gene_relation, gene_distance,
   RNAseq_mapped_count

   gene_relation is "overlap" if the region falls inside/over a gene body,
   or "nearest" if no gene overlapped and the closest one(s) within
   --max-gene-distance were used instead. gene_distance is 0 for overlaps,
   otherwise the bp distance to the nearest gene boundary. gene_name/
   gene_id can hold multiple ';'-separated values if a region overlaps
   more than one gene, or is equidistant from more than one.

   RNAseq_mapped_count sums the RNA-seq counts (across all sample/BAM
   columns) for whatever gene(s) got assigned to the region. Pass
   --keep-per-sample-counts to get one column per sample instead.

   Chromosome names are normalized to GENCODE/UCSC "chr..." style where
   possible (e.g. "1" -> "chr1", "MT" -> "chrM"). Accession-style
   scaffold/patch names (e.g. "GL000008.2") can't be resolved by a simple
   rule -- supply --chrom-alias-map for those, or use
   --drop-unmapped-chroms to discard them.

--------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------

    python merge_peaks_to_genes.py \
        --peaks peaks.tsv \
        --gtf gencode.v50.annotation.gtf \
        --rnaseq gene_counts.tsv \
        --distance 100 \
        --max-gene-distance 100000 \
        --out merged_regions.tsv

Multiple peak files can be supplied (concatenated first):

    python merge_peaks_to_genes.py --peaks chip1.tsv chip2.tsv --gtf gencode.gtf --out out.tsv
"""

import argparse
import re
import sys

import numpy as np
import pandas as pd

PEAK_COLUMNS_10 = [
    "peak_chr", "peak_start", "peak_end", "peak_name", "score", "strand",
    "signalValue", "pValue", "qValue", "summit",
]

# Standard set of "main" chromosomes/contigs to keep by default. GENCODE
# "comprehensive/ALL" GTFs also include scaffolds, patches, and alternate
# haplotype loci (e.g. chr1_KI270706v1_random, chr19_GL000209v2_alt) which
# you usually don't want competing for gene assignment unless your peaks
# were actually called against those contigs too.
MAIN_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM", "chrMT"}

# Matches Ensembl/NCBI-style "bare" chromosome names for the 24 primary
# chromosomes + mitochondrion, e.g. "1", "17", "X", "MT" -- as opposed to
# accession-style scaffold/patch names like "GL000008.2" or "KI270706.1",
# which are NOT primary chromosomes and can't be turned into "chrN" by a
# simple rule (see normalize_chrom / --chrom-alias-map below).
PRIMARY_CHROM_RE = re.compile(r"^(\d{1,2}|[XY]|MT)$", re.IGNORECASE)


def load_chrom_alias_map(path):
    """
    Load a simple 2-column TSV mapping arbitrary chromosome/contig names to
    their UCSC/GENCODE "chr..." equivalents, e.g.:

        GL000008.2    chr4_GL000008v2_random
        KI270706.1    chr1_KI270706v1_random

    Such a table can be built from an NCBI assembly report (columns
    "GenBank-Accn" and "UCSC-style-name") or from UCSC's own
    hg38.chromAlias.txt. Lines starting with '#' are ignored.
    """
    alias_map = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            alias_map[parts[0]] = parts[1]
    return alias_map


def normalize_chrom(chrom: str, alias_map: dict | None = None) -> str:
    """
    Turn a chromosome/contig name into GENCODE/UCSC "chr..." style where
    possible:
      - already "chr..." -> unchanged
      - found in a supplied alias map -> mapped value
      - bare Ensembl/NCBI primary chromosome name ("1", "17", "X", "MT")
        -> "chr1", "chr17", "chrX", "chrM"
      - anything else (unrecognized scaffold/patch accession) -> unchanged,
        and the caller is responsible for tracking/handling these
    """
    if chrom.startswith("chr"):
        return chrom
    if alias_map and chrom in alias_map:
        return alias_map[chrom]
    m = PRIMARY_CHROM_RE.match(chrom)
    if m:
        val = m.group(1).upper()
        return "chrM" if val == "MT" else f"chr{val}"
    return chrom


# --------------------------------------------------------------------------
# Peak loading + merging
# --------------------------------------------------------------------------

def load_peaks(paths, alias_map=None, drop_unmapped_chroms=False):
    """Load and concatenate one or more peak TSVs (10- or 18-column layout)."""
    frames = []
    for path in paths:
        if not path.lower().endswith((".tsv", ".bed", ".txt", ".narrowpeak")):
            print(f"[warn] {path}: unexpected file extension for a peak file "
                  f"(expected .tsv/.bed/.txt) -- double check this is the right file.",
                  file=sys.stderr)
        try:
            df = pd.read_csv(path, sep=r"\s+", header=None, comment="#")
        except Exception as e:
            sys.exit(f"[error] failed to parse {path} as a whitespace/tab-delimited peak file: {e}")
        if df.shape[1] not in (10, 18):
            sys.exit(
                f"[error] {path}: expected 10 or 18 columns, found {df.shape[1]}. "
                f"Edit PEAK_COLUMNS_10 in the script if your layout differs, or "
                f"double-check this is actually a peak file and not passed in by mistake."
            )
        df = df.iloc[:, :10]
        df.columns = PEAK_COLUMNS_10
        frames.append(df)
    peaks = pd.concat(frames, ignore_index=True)

    for col in ("peak_start", "peak_end", "signalValue"):
        peaks[col] = pd.to_numeric(peaks[col], errors="coerce")
    peaks = peaks.dropna(subset=["peak_start", "peak_end"])
    peaks["peak_start"] = peaks["peak_start"].astype(int)
    peaks["peak_end"] = peaks["peak_end"].astype(int)

    # Normalize chromosome/contig names to "chr..." style so they match the
    # GTF (e.g. "1" -> "chr1", "MT" -> "chrM"). Anything that still isn't
    # "chr..."-prefixed afterwards is an accession-style scaffold/patch
    # name (e.g. "GL000008.2") that a simple rule can't resolve -- these
    # need an explicit --chrom-alias-map entry, or get flagged/dropped.
    peaks["peak_chr"] = peaks["peak_chr"].apply(lambda c: normalize_chrom(c, alias_map))
    unmapped = sorted(peaks.loc[~peaks["peak_chr"].str.startswith("chr"), "peak_chr"].unique())
    if unmapped:
        n_unmapped_peaks = (~peaks["peak_chr"].str.startswith("chr")).sum()
        preview = ", ".join(unmapped[:10]) + (", ..." if len(unmapped) > 10 else "")
        print(f"[warn] {n_unmapped_peaks} peaks are on {len(unmapped)} contig(s) that could not be "
              f"normalized to 'chr...' style (likely unplaced/unlocalized scaffolds or patches, "
              f"not primary chromosomes): {preview}\n"
              f"       These will not match genes in the GTF unless you provide --chrom-alias-map, "
              f"or pass --include-alt-scaffolds together with a correct alias so the scaffold's own "
              f"GTF entry can be used.",
              file=sys.stderr)
        if drop_unmapped_chroms:
            peaks = peaks[peaks["peak_chr"].str.startswith("chr")].reset_index(drop=True)
            print(f"[info] --drop-unmapped-chroms: removed {n_unmapped_peaks} peaks on unmapped contigs",
                  file=sys.stderr)

    return peaks


def merge_peaks(peaks: pd.DataFrame, distance: int) -> pd.DataFrame:
    """Merge peaks within `distance` bp of each other (per chromosome)."""
    merged_rows = []

    for chrom, chrom_df in peaks.groupby("peak_chr"):
        chrom_df = chrom_df.sort_values("peak_start").reset_index(drop=True)

        cur_start = cur_end = None
        cur_peaks, cur_max_signal = [], None

        def flush():
            if cur_start is None:
                return
            merged_rows.append({
                "region_chr": chrom,
                "region_start": cur_start,
                "region_end": cur_end,
                "n_peaks": len(cur_peaks),
                "peak_names": ";".join(cur_peaks),
                "max_signalValue": cur_max_signal,
            })

        for _, row in chrom_df.iterrows():
            start, end, signal = row["peak_start"], row["peak_end"], row["signalValue"]
            if cur_start is None:
                cur_start, cur_end = start, end
                cur_peaks = [str(row["peak_name"])]
                cur_max_signal = signal
            elif start - cur_end <= distance:
                cur_end = max(cur_end, end)
                cur_peaks.append(str(row["peak_name"]))
                if pd.notna(signal) and (pd.isna(cur_max_signal) or signal > cur_max_signal):
                    cur_max_signal = signal
            else:
                flush()
                cur_start, cur_end = start, end
                cur_peaks = [str(row["peak_name"])]
                cur_max_signal = signal
        flush()

    result = pd.DataFrame(merged_rows)
    if not result.empty:
        result = result.sort_values(["region_chr", "region_start"]).reset_index(drop=True)
    return result


# --------------------------------------------------------------------------
# GTF gene annotation loading + assignment
# --------------------------------------------------------------------------

ATTR_RE = re.compile(r'(\w+) "([^"]+)"')


def load_gtf_genes(path, primary_chroms_only=True):
    """
    Parse a GTF and return a DataFrame of gene features:
    gene_chr, gene_start, gene_end, gene_strand, gene_id, gene_name
    (GTF coordinates are 1-based inclusive; kept as-is here.)
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom, _, _, start, end, _, strand = fields[0], *fields[1:7]
            attrs = dict(ATTR_RE.findall(fields[8]))
            rows.append({
                "gene_chr": chrom,
                "gene_start": int(start),
                "gene_end": int(end),
                "gene_strand": strand,
                "gene_id": attrs.get("gene_id", ""),
                "gene_name": attrs.get("gene_name", attrs.get("gene_id", "")),
            })

    genes = pd.DataFrame(rows)
    if primary_chroms_only:
        genes = genes[genes["gene_chr"].isin(MAIN_CHROMS)].reset_index(drop=True)
    return genes


def build_gene_index(genes: pd.DataFrame):
    """Per-chromosome numpy arrays of gene start/end for fast vectorized lookup."""
    index = {}
    for chrom, g in genes.groupby("gene_chr"):
        g = g.sort_values("gene_start").reset_index(drop=True)
        index[chrom] = {
            "start": g["gene_start"].to_numpy(),
            "end": g["gene_end"].to_numpy(),
            "gene_id": g["gene_id"].to_numpy(),
            "gene_name": g["gene_name"].to_numpy(),
        }
    return index


def assign_genes(regions: pd.DataFrame, gene_index: dict, max_gene_distance: int):
    """
    For each merged region, find overlapping genes; if none, the nearest
    gene(s) within max_gene_distance bp. Adds gene_names, gene_ids,
    gene_relation, gene_distance columns.
    """
    names_col, ids_col, relation_col, dist_col = [], [], [], []

    for _, row in regions.iterrows():
        chrom, start, end = row["region_chr"], row["region_start"], row["region_end"]
        idx = gene_index.get(chrom)

        if idx is None or len(idx["start"]) == 0:
            names_col.append("")
            ids_col.append("")
            relation_col.append("none")
            dist_col.append(np.nan)
            continue

        g_start, g_end = idx["start"], idx["end"]

        # overlap: gene_start <= region_end AND gene_end >= region_start
        overlap_mask = (g_start <= end) & (g_end >= start)

        if overlap_mask.any():
            names_col.append(";".join(sorted(set(idx["gene_name"][overlap_mask]))))
            ids_col.append(";".join(sorted(set(idx["gene_id"][overlap_mask]))))
            relation_col.append("overlap")
            dist_col.append(0)
            continue

        # no overlap -> nearest gene boundary distance
        downstream_dist = np.where(g_start > end, g_start - end, np.inf)
        upstream_dist = np.where(g_end < start, start - g_end, np.inf)
        dist = np.minimum(downstream_dist, upstream_dist)
        min_dist = dist.min()

        if min_dist <= max_gene_distance:
            nearest_mask = dist == min_dist
            names_col.append(";".join(sorted(set(idx["gene_name"][nearest_mask]))))
            ids_col.append(";".join(sorted(set(idx["gene_id"][nearest_mask]))))
            relation_col.append("nearest")
            dist_col.append(int(min_dist))
        else:
            names_col.append("")
            ids_col.append("")
            relation_col.append("none")
            dist_col.append(int(min_dist))

    regions = regions.copy()
    regions["gene_name"] = names_col
    regions["gene_id"] = ids_col
    regions["gene_relation"] = relation_col
    regions["gene_distance"] = dist_col
    return regions


# --------------------------------------------------------------------------
# RNA-seq expression
# --------------------------------------------------------------------------

def load_rnaseq(path):
    with open(path) as fh:
        first_line = fh.readline()
    skiprows = 1 if first_line.startswith("#") else 0

    df = pd.read_csv(path, sep="\t", skiprows=skiprows)
    required = {"Geneid", "gene_name"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"[error] {path}: missing expected column(s) {missing}")

    non_sample_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length", "gene_name"]
    sample_cols = [c for c in df.columns if c not in non_sample_cols]
    return df[["Geneid", "gene_name"] + sample_cols], sample_cols


def attach_expression(regions: pd.DataFrame, rnaseq: pd.DataFrame, sample_cols,
                       combine_samples=True, column_name="RNAseq_mapped_count"):
    """
    Sum RNA-seq counts across the assigned gene(s) for each region. By
    default all sample/BAM columns are also summed together into a single
    `column_name` column; pass combine_samples=False to keep one output
    column per original sample instead.
    """
    if regions.empty:
        if combine_samples:
            regions[column_name] = []
        else:
            for col in sample_cols:
                regions[col] = []
        return regions

    # gene_name is not guaranteed unique (paralogs, readthrough genes,
    # duplicate symbols across loci, etc.) -- collapse duplicates by
    # summing their counts so the lookup index below is unique.
    n_before = rnaseq["gene_name"].nunique()
    n_rows = len(rnaseq)
    if n_rows != n_before:
        print(f"[warn] {n_rows - n_before} duplicate gene_name entries in RNA-seq "
              f"file collapsed by summing counts (e.g. paralogs sharing a symbol).",
              file=sys.stderr)
    name_to_counts = rnaseq.groupby("gene_name")[sample_cols].sum()

    def sum_for_region(gene_names_str):
        names = [g for g in gene_names_str.split(";") if g]
        if not names:
            return pd.Series({c: 0 for c in sample_cols})
        matched = name_to_counts.reindex(names).fillna(0)
        return matched.sum(axis=0)

    expr = regions["gene_name"].apply(sum_for_region)
    regions = regions.copy()
    if combine_samples:
        regions[column_name] = expr.sum(axis=1)
        return regions
    return pd.concat([regions, expr], axis=1)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--peaks", nargs="+", required=True, help="One or more ChIP-seq peak TSV files")
    parser.add_argument("--gtf", required=True, help="GENCODE comprehensive GTF annotation file")
    parser.add_argument("--rnaseq", help="Optional featureCounts gene_counts TSV to attach expression values")
    parser.add_argument("--distance", type=int, default=100, help="Max gap (bp) to merge peaks together (default: 100)")
    parser.add_argument("--max-gene-distance", type=int, default=100_000,
                         help="Max bp distance to still assign a nearest gene when there's no overlap (default: 100000)")
    parser.add_argument("--include-alt-scaffolds", action="store_true",
                         help="Also use scaffold/patch/alt-haplotype contigs from the GTF (default: main chromosomes only)")
    parser.add_argument("--chrom-alias-map",
                         help="Optional 2-column TSV mapping non-standard chromosome/contig names "
                              "(e.g. 'GL000008.2') to their GENCODE/UCSC 'chr...' equivalents")
    parser.add_argument("--drop-unmapped-chroms", action="store_true",
                         help="Discard peaks on contigs that can't be normalized to 'chr...' style "
                              "(e.g. unplaced scaffold accessions with no alias-map entry)")
    parser.add_argument("--keep-per-sample-counts", action="store_true",
                         help="Keep one RNA-seq count column per sample/BAM instead of summing them "
                              "into a single combined column")
    parser.add_argument("--rnaseq-column-name", default="RNAseq_mapped_count",
                         help="Name of the combined RNA-seq count column (default: RNAseq_mapped_count)")
    parser.add_argument("--out", default="merged_regions.csv", help="Output CSV path")
    args = parser.parse_args()

    alias_map = load_chrom_alias_map(args.chrom_alias_map) if args.chrom_alias_map else None
    peaks = load_peaks(args.peaks, alias_map=alias_map, drop_unmapped_chroms=args.drop_unmapped_chroms)
    print(f"[info] loaded {len(peaks)} peaks from {len(args.peaks)} file(s)")

    regions = merge_peaks(peaks, args.distance)
    print(f"[info] merged into {len(regions)} regions (max gap = {args.distance} bp)")

    genes = load_gtf_genes(args.gtf, primary_chroms_only=not args.include_alt_scaffolds)
    print(f"[info] loaded {len(genes)} gene features from GTF")
    gene_index = build_gene_index(genes)

    regions = assign_genes(regions, gene_index, args.max_gene_distance)
    n_overlap = (regions["gene_relation"] == "overlap").sum()
    n_nearest = (regions["gene_relation"] == "nearest").sum()
    n_none = (regions["gene_relation"] == "none").sum()
    print(f"[info] gene assignment: {n_overlap} overlap, {n_nearest} nearest, {n_none} unassigned")

    if args.rnaseq:
        rnaseq, sample_cols = load_rnaseq(args.rnaseq)
        print(f"[info] loaded RNA-seq counts for {len(rnaseq)} genes, {len(sample_cols)} sample column(s)")
        regions = attach_expression(
            regions, rnaseq, sample_cols,
            combine_samples=not args.keep_per_sample_counts,
            column_name=args.rnaseq_column_name,
        )

    regions.to_csv(args.out, sep=",", index=False)
    print(f"[info] wrote {args.out}")


if __name__ == "__main__":
    main()