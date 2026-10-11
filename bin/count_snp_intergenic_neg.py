"""recalc_enrichment_intergenic.py — recompute SNP enrichment using a set of
randomly sampled intergenic background regions, instead of gene flanks.

    enrichment = gene_body_snp_density / intergenic_snp_density

Background-region sampling follows this method:

    10,000 intergenic regions, of length 150nt, were randomly selected from
    chr1-22, chrX and chrY (GENCODE vNN primary assembly, GRCh38) with
    BEDtools and custom Python scripts. Regions did not overlap one
    another or known gaps in GRCh38, were at least 10kbp away from the
    start or end coordinates of any gene in GENCODE vNN, and were
    distributed through the autosomes and sex chromosomes in proportion
    to the size of each chromosome.

This script implements that logic natively in Python (no BEDtools
dependency required) rather than shelling out to it — the merge/subtract/
random-sample operations are equivalent:

  1. Parse a genome-wide GENCODE/Ensembl/RefSeq gene annotation (GTF/GFF3
     "gene" features, or a plain BED of gene loci) — every annotated gene
     genome-wide, not just the genes in your input list.
  2. Pad every gene by --gene-buffer-bp (default 10,000) on each side and
     subtract it, along with any assembly gaps (--gaps-bed, optional), from
     chr1-22/X/Y to get the eligible intergenic space.
  3. Allocate --num-regions (default 10,000) total regions across those
     chromosomes in proportion to chromosome length, then randomly place
     non-overlapping --region-length (default 150) bp windows inside the
     eligible space of each chromosome.
  4. Count each source's SNPs in those windows to get one pooled,
     genome-wide intergenic SNP density per source.
  5. Recompute enrichment for every gene in the input CSV as
     snp_density_<source> / intergenic_density_<source>.

Input:
    The summary CSV produced by count_snp_density.py, i.e. one row per
    gene with columns:
      gene_name, gene_group, chromosome, start, end, gene_length,
      flank_length, snp_count_<src>, snp_density_<src>, flank_count_<src>,
      flank_density_<src>, enrichment_<src>   (for src in 1000genomes,
      gnomad, pangenome, dbsnp), plus locus_type, gene_type.

    A genome-wide gene annotation file (--gene-annotation), e.g. GENCODE's
    comprehensive gene annotation GTF (gencode.vNN.annotation.gtf.gz).

    Optionally, an assembly-gaps BED (--gaps-bed), e.g. the UCSC GRCh38
    "gap" table exported as BED3 (chrom, chromStart, chromEnd).

Output:
    The same rows, with two columns added per source:
      intergenic_density_<src>    — pooled genome-wide intergenic SNP rate
      enrichment_intergenic_<src> — snp_density_<src> / intergenic_density_<src>

This script reuses the fetch/count machinery from count_snp_density.py
(must be importable — keep it in the same directory):
  - cse.fetch_remote_vcf_regions() for 1000 Genomes and gnomAD ("gcs"
    source) — same `bcftools view -R <regions.bed>` pattern, pointed at
    the sampled intergenic regions.
  - cse.extract_positions() for the local pangenome VCF and dbSNP BigBed.

Usage:
    python recalc_enrichment_intergenic.py \
        --input-csv results/snp_enrichment_summary.csv \
        --gene-annotation data/downloads/gencode.v49.annotation.gtf.gz

    python recalc_enrichment_intergenic.py --gene-annotation genes.gtf.gz --gaps-bed gaps.bed
    python recalc_enrichment_intergenic.py --gene-annotation genes.gtf.gz --skip-fetch
"""

import argparse
import bisect
import csv
import gzip
import logging
from math import log2
import os
import random
import sys
import time
from collections import defaultdict

import count_snp_density as cse

logger = logging.getLogger("recalc_enrichment_intergenic")

# Primary-assembly chromosome lengths (bp), keyed the same way this script's
# window dicts are keyed: bare, uppercase ("1", "X", "Y").
GRCH38_CHROM_SIZES = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259,
    "6": 170805979, "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422,
    "11": 135086622, "12": 133275309, "13": 114364328, "14": 107043718, "15": 101991189,
    "16": 90338345, "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
    "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
}
GRCH37_CHROM_SIZES = {
    "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260,
    "6": 171115067, "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747,
    "11": 135006516, "12": 133851895, "13": 115169878, "14": 107349540, "15": 102531392,
    "16": 90354753, "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
    "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
}

# chr1-22, chrX, chrY only — matches the sampling method exactly (no chrM).
CANONICAL_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]

SOURCES = ("1000genomes", "gnomad", "pangenome", "dbsnp")


def setup_logging(verbose=False):
    # Configure count_snp_density's own logger (used inside the functions
    # we import from it) as well as this script's logger.
    cse.setup_logging(verbose)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
    logger.addHandler(handler)


def get_chrom_sizes(genome_build):
    return GRCH38_CHROM_SIZES if genome_build == "GRCh38" else GRCH37_CHROM_SIZES


def _chrom_key(chrom_raw):
    chrom = cse.normalize_chrom(chrom_raw, style="bare").upper()
    return "MT" if chrom == "M" else chrom


def read_summary_rows(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input CSV not found: {path}")
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    logger.info("Read %d gene row(s) from %s", len(rows), path)
    return rows, fieldnames


# ==========================================================================
# Parse genome-wide gene annotation and assembly gaps
# ==========================================================================

def _open_maybe_gz(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.lower().endswith(".gz") \
        else open(path, "rt", encoding="utf-8", errors="replace")


def parse_gene_annotation(path):
    """GTF/GFF3 'gene' feature rows, or a plain BED of gene loci. Returns
    EVERY annotated gene locus genome-wide (merged per chromosome) — not
    just the genes in the input CSV."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Gene annotation file not found: {path}")

    lower = path.lower()
    is_bed = ".bed" in lower.replace(".gz", "") and not (".gtf" in lower or ".gff" in lower)

    genes_by_chrom = defaultdict(list)
    n_lines = 0
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            fields = line.rstrip("\n").split("\t")

            if is_bed:
                if len(fields) < 3:
                    continue
                chrom_raw = fields[0]
                if not cse.is_standard_chromosome(chrom_raw):
                    continue
                try:
                    start, end = int(fields[1]) + 1, int(fields[2])  # BED is 0-based, half-open
                except ValueError:
                    continue
            else:
                if len(fields) < 5 or fields[2] != "gene":
                    continue
                chrom_raw = fields[0]
                if not cse.is_standard_chromosome(chrom_raw):
                    continue
                try:
                    start, end = int(fields[3]), int(fields[4])  # GTF/GFF3 is 1-based, inclusive
                except ValueError:
                    continue

            genes_by_chrom[_chrom_key(chrom_raw)].append((start, end))
            n_lines += 1

    if n_lines == 0:
        raise ValueError(
            f"No gene loci parsed from {path}. Expected a GTF/GFF3 with 'gene' feature rows "
            f"(field 3 == 'gene'), or a BED file of gene loci."
        )

    for chrom in genes_by_chrom:
        genes_by_chrom[chrom] = cse.merge_regions(genes_by_chrom[chrom])
    total = sum(len(v) for v in genes_by_chrom.values())
    logger.info("Parsed %d gene locus/loci across %d chromosome(s) from %s", total, len(genes_by_chrom), path)
    return genes_by_chrom


def parse_gaps_bed(path):
    """Assembly gaps (e.g. UCSC GRCh38 'gap' table exported as BED3:
    chrom, chromStart, chromEnd). Optional — returns {} if no path given."""
    if not path:
        return {}
    if not os.path.exists(path):
        raise FileNotFoundError(f"Gaps BED file not found: {path}")

    gaps_by_chrom = defaultdict(list)
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom_raw = fields[0]
            if not cse.is_standard_chromosome(chrom_raw):
                continue
            try:
                start, end = int(fields[1]) + 1, int(fields[2])  # BED is 0-based, half-open
            except ValueError:
                continue
            gaps_by_chrom[_chrom_key(chrom_raw)].append((start, end))

    for chrom in gaps_by_chrom:
        gaps_by_chrom[chrom] = cse.merge_regions(gaps_by_chrom[chrom])
    total = sum(len(v) for v in gaps_by_chrom.values())
    logger.info("Parsed %d gap region(s) across %d chromosome(s) from %s", total, len(gaps_by_chrom), path)
    return gaps_by_chrom


def compute_intergenic_regions(chrom_sizes, chroms, genes_by_chrom, gaps_by_chrom, buffer_bp):
    """Per chromosome: full length minus (genes padded by buffer_bp on each
    side) minus assembly gaps. This is the eligible space background
    regions are sampled from."""
    intergenic_by_chrom = {}
    for chrom in chroms:
        size = chrom_sizes[chrom]
        genes = genes_by_chrom.get(chrom, [])
        if buffer_bp > 0 and genes:
            genes = [(max(1, s - buffer_bp), e + buffer_bp) for s, e in genes]
        gaps = gaps_by_chrom.get(chrom, [])
        excluded = cse.merge_regions(genes + gaps)

        intergenic, cursor = [], 1
        for start, end in excluded:
            if start > cursor:
                intergenic.append((cursor, start - 1))
            cursor = max(cursor, end + 1)
        if cursor <= size:
            intergenic.append((cursor, size))

        intergenic_by_chrom[chrom] = intergenic
        excluded_bp = sum(e - s + 1 for s, e in excluded)
        logger.info("chr%s: %d bp excluded (genes ±%dbp + gaps), %d bp eligible intergenic (of %d bp total)",
                     chrom, excluded_bp, buffer_bp, size - excluded_bp, size)
    return intergenic_by_chrom


# ==========================================================================
# Allocate regions proportional to chromosome size, then sample them
# non-overlapping from the eligible intergenic space
# ==========================================================================

def allocate_regions_by_chrom_size(chrom_sizes, chroms, num_regions):
    """Largest-remainder allocation so the per-chromosome counts sum exactly
    to num_regions while staying proportional to chromosome length."""
    total_size = sum(chrom_sizes[c] for c in chroms)
    raw = {c: num_regions * chrom_sizes[c] / total_size for c in chroms}
    allocation = {c: int(raw[c]) for c in chroms}
    remainder = num_regions - sum(allocation.values())
    for c in sorted(chroms, key=lambda c: raw[c] - allocation[c], reverse=True)[:remainder]:
        allocation[c] += 1
    for c in chroms:
        logger.info("chr%s: allocated %d region(s) (%.2f%% of genome)", c, allocation[c], 100 * chrom_sizes[c] / total_size)
    return allocation


def sample_regions_in_intergenic(intergenic_by_chrom, allocation, region_length, seed):
    rng = random.Random(seed)
    windows_by_chrom = {}
    for chrom, intervals in sorted(intergenic_by_chrom.items()):
        target = allocation.get(chrom, 0)
        if target == 0:
            windows_by_chrom[chrom] = []
            continue

        usable = [(s, e) for s, e in intervals if e - s + 1 >= region_length]
        if not usable:
            logger.warning("chr%s: no eligible interval is >= region_length (%d bp); skipping.", chrom, region_length)
            windows_by_chrom[chrom] = []
            continue

        weights = [e - s + 1 - region_length + 1 for s, e in usable]  # valid start positions per interval
        chosen = []
        attempts, max_attempts = 0, max(target * 300, 5000)
        while len(chosen) < target and attempts < max_attempts:
            attempts += 1
            s, e = rng.choices(usable, weights=weights, k=1)[0]
            start = rng.randint(s, e - region_length + 1)
            end = start + region_length - 1
            if any(start <= w_end and end >= w_start for w_start, w_end in chosen):
                continue  # regions must not overlap one another
            chosen.append((start, end))

        if len(chosen) < target:
            logger.warning("chr%s: only sampled %d/%d intergenic region(s).", chrom, len(chosen), target)
        windows_by_chrom[chrom] = sorted(chosen)
    return windows_by_chrom


def build_window_index(background_windows):
    """{chrom: (sorted_starts, sorted_windows)} for bisect-based membership
    lookups against non-overlapping windows."""
    index = {}
    for chrom, windows in background_windows.items():
        sorted_windows = sorted(windows)
        index[chrom] = ([w[0] for w in sorted_windows], sorted_windows)
    return index


def position_in_windows(chrom, pos, window_index):
    entry = window_index.get(chrom)
    if not entry:
        return False
    starts, windows = entry
    i = bisect.bisect_right(starts, pos) - 1
    if i < 0:
        return False
    w_start, w_end = windows[i]
    return w_start <= pos <= w_end


def total_background_bp(background_windows, chroms=None):
    items = background_windows.items() if chroms is None else (
        (c, w) for c, w in background_windows.items() if c in chroms
    )
    return sum(end - start + 1 for _, windows in items for start, end in windows)


# ==========================================================================
# Count background SNPs per source
# ==========================================================================

def count_vcf_variants(path):
    if not path or not os.path.exists(path):
        return 0
    count = 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            count += 1
    return count


def fetch_and_count_1000genomes(background_windows, args):
    """bcftools -R already restricts the fetch to the sampled regions, so
    every variant returned counts toward the background total directly."""
    if not cse.check_bcftools():
        logger.warning("Skipping 1000 Genomes background — bcftools unavailable.")
        return (0, 0)

    template = cse.GRCH38_TEMPLATE if args.genome_build == "GRCh38" else cse.GRCH37_TEMPLATE
    os.makedirs(args.vcf_dir, exist_ok=True)
    total_count, counted_chroms = 0, set()

    for chrom, windows in sorted(background_windows.items()):
        if not windows:
            continue
        output_vcf = os.path.join(args.vcf_dir, f"background_1000genomes_chr{chrom}.vcf")
        if args.skip_fetch or (not args.force and os.path.exists(output_vcf) and os.path.getsize(output_vcf) > 0):
            if not os.path.exists(output_vcf):
                logger.warning("--skip-fetch set but %s is missing; skipping chr%s.", output_vcf, chrom)
                continue
            result = output_vcf
        else:
            remote_url = template.format(chrom=chrom)
            result = cse.fetch_remote_vcf_regions(remote_url, chrom, windows, output_vcf, "1000genomes-background")
        if not result:
            continue
        total_count += count_vcf_variants(result)
        counted_chroms.add(chrom)

    return (total_count, total_background_bp(background_windows, counted_chroms))


def fetch_and_count_gnomad_gcs(background_windows, args):
    if not cse.check_bcftools():
        logger.warning("Skipping gnomAD background — bcftools unavailable.")
        return (0, 0)

    os.makedirs(args.gnomad_dir, exist_ok=True)
    total_count, counted_chroms = 0, set()

    for chrom, windows in sorted(background_windows.items()):
        if not windows or chrom.lower() in cse.GNOMAD_GENOME_SKIP_CHROMS:
            continue
        chrom_style = f"chr{chrom}"
        output_vcf = os.path.join(args.gnomad_dir, f"background_gnomad_genomes_{chrom_style}.vcf")
        if args.skip_fetch or (not args.force and os.path.exists(output_vcf) and os.path.getsize(output_vcf) > 0):
            if not os.path.exists(output_vcf):
                logger.warning("--skip-fetch set but %s is missing; skipping chr%s.", output_vcf, chrom)
                continue
            result = output_vcf
        else:
            url = cse.GNOMAD_GCS_GENOME_TEMPLATE.format(chrom=chrom)
            result = cse.fetch_remote_vcf_regions(url, chrom_style, windows, output_vcf, "gnomad-background")
        if not result:
            continue
        total_count += count_vcf_variants(result)
        counted_chroms.add(chrom)

    return (total_count, total_background_bp(background_windows, counted_chroms))


def fetch_and_count_gnomad_api(background_windows, args):
    import requests
    session = requests.Session()
    total_count, counted_chroms = 0, set()

    for chrom, windows in sorted(background_windows.items()):
        if not windows or chrom.lower() in cse.GNOMAD_GENOME_SKIP_CHROMS:
            continue
        chrom_ok = True
        for start, end in windows:
            for chunk_start, chunk_end in cse.split_region(start, end, args.gnomad_max_region):
                variables = {"chrom": chrom, "start": chunk_start, "stop": chunk_end}
                ok = False
                for attempt in range(1, 6):
                    try:
                        response = session.post(cse.GNOMAD_API_URL,
                                                 json={"query": cse.GNOMAD_API_QUERY, "variables": variables},
                                                 timeout=120)
                    except Exception as exc:
                        logger.warning("gnomAD background network error %s:%d-%d (attempt %d/5): %s",
                                        chrom, chunk_start, chunk_end, attempt, exc)
                        time.sleep(5)
                        continue
                    time.sleep(1.0)
                    if response.status_code == 200:
                        data = response.json()
                        errors = data.get("errors")
                        if errors:
                            messages = "; ".join(e.get("message", str(e)) for e in errors)
                            logger.warning("gnomAD background errors %s:%d-%d (attempt %d/5): %s",
                                            chrom, chunk_start, chunk_end, attempt, messages)
                            time.sleep(5)
                            continue
                        variants = ((data.get("data") or {}).get("region") or {}).get("variants") or []
                        total_count += len(variants)
                        ok = True
                        break
                    elif response.status_code in (429, 500, 502, 503, 504):
                        time.sleep(5 * attempt)
                    else:
                        logger.error("gnomAD background HTTP %d for %s:%d-%d — not retrying.",
                                      response.status_code, chrom, chunk_start, chunk_end)
                        break
                if not ok:
                    logger.warning("Giving up on gnomAD background %s:%d-%d", chrom, chunk_start, chunk_end)
                    chrom_ok = False
        if chrom_ok:
            counted_chroms.add(chrom)

    return (total_count, total_background_bp(background_windows, counted_chroms))


def fetch_and_count_gnomad(background_windows, args):
    if args.gnomad_source == "gcs":
        return fetch_and_count_gnomad_gcs(background_windows, args)
    return fetch_and_count_gnomad_api(background_windows, args)


def count_local_vcf_background(path, background_windows):
    """Local pangenome VCF: assumed genome-wide, so read every position
    (cse.extract_positions' plain .vcf path) and filter to the sampled
    regions ourselves via bisect."""
    total_bp = total_background_bp(background_windows)
    if not path or not os.path.exists(path):
        logger.error("Pangenome VCF not found at %s — its intergenic density will be 0.", path)
        return (0, total_bp)

    positions_by_chrom = cse.extract_positions(path)
    window_index = build_window_index(background_windows)
    count = 0
    for chrom, positions in positions_by_chrom.items():
        if chrom not in window_index:
            continue
        for pos in positions:
            if position_in_windows(chrom, pos, window_index):
                count += 1
    return (count, total_bp)


def count_local_bigbed_background(path, background_windows):
    """Local dbSNP BigBed: use cse.extract_positions' random-access BigBed
    path, pointed at the sampled regions instead of gene+flank windows (it
    only reads the flank_start/flank_end keys, so a minimal fake window
    dict is enough)."""
    total_bp = total_background_bp(background_windows)
    if not path or not os.path.exists(path):
        logger.error("dbSNP BigBed not found at %s — its intergenic density will be 0.", path)
        return (0, total_bp)

    fake_windows_by_chrom = {
        chrom: [
            {"start": start, "end": end,
            "flank_start": start, "flank_end": end}
            for start, end in windows
        ]
        for chrom, windows in background_windows.items() if windows
    }
    positions_by_chrom = cse.extract_positions(path, fake_windows_by_chrom)
    count = sum(len(p) for p in positions_by_chrom.values())
    return (count, total_bp)


# ==========================================================================
# Recompute enrichment and write output
# ==========================================================================

PSEUDO_RATIO = 1 / 1000

def compute_enrichment_from_rates(gene_rate, background_rate):
    if background_rate == 0:
        return float("inf") if gene_rate > 0 else 0.0
    ratio = gene_rate / background_rate
    if ratio == 0:
        ratio = PSEUDO_RATIO  # log2(0.001) ≈ -9.97
    return log2(ratio)

def write_output(rows, fieldnames, source_stats, output_csv):
    new_cols = []
    for label in SOURCES:
        new_cols += [f"intergenic_density_{label}", f"enrichment_intergenic_{label}"]
    out_fieldnames = fieldnames + [c for c in new_cols if c not in fieldnames]

    rates = {}
    for label, (count, bp) in source_stats.items():
        rate = count / bp if bp > 0 else 0.0
        rates[label] = rate
        logger.info("%s intergenic background: %d variant(s) / %d bp = %.6g SNPs/bp", label, count, bp, rate)

    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
    with open(output_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_fieldnames)
        writer.writeheader()
        for row in rows:
            out_row = dict(row)
            for label in SOURCES:
                background_rate = rates.get(label, 0.0)
                try:
                    gene_rate = float(row.get(f"snp_density_{label}") or 0.0)
                except ValueError:
                    gene_rate = 0.0
                enrichment = compute_enrichment_from_rates(gene_rate, background_rate)
                out_row[f"intergenic_density_{label}"] = background_rate
                # Same convention as count_snp_density.py: write blank,
                # not "inf", so it reads back as NA in R/pandas.
                out_row[f"enrichment_intergenic_{label}"] = "" if enrichment == float("inf") else enrichment
            writer.writerow(out_row)
    logger.info("Wrote %d row(s) to %s", len(rows), output_csv)


# ==========================================================================
# Pipeline + CLI
# ==========================================================================

def run(args):
    setup_logging(args.verbose)
    rows, fieldnames = read_summary_rows(args.input_csv)

    chrom_sizes_full = get_chrom_sizes(args.genome_build)
    chroms = [c for c in CANONICAL_CHROMS if c in chrom_sizes_full]
    chrom_sizes = {c: chrom_sizes_full[c] for c in chroms}

    genes_by_chrom = parse_gene_annotation(args.gene_annotation)
    gaps_by_chrom = parse_gaps_bed(args.gaps_bed)

    intergenic_by_chrom = compute_intergenic_regions(chrom_sizes, chroms, genes_by_chrom, gaps_by_chrom,
                                                       buffer_bp=args.gene_buffer_bp)
    allocation = allocate_regions_by_chrom_size(chrom_sizes, chroms, args.num_regions)
    background_windows = sample_regions_in_intergenic(intergenic_by_chrom, allocation, args.region_length, args.seed)

    total_regions = sum(len(w) for w in background_windows.values())
    logger.info("Sampled %d/%d intergenic region(s) (%d bp each) across chr%s, total %d bp.",
                total_regions, args.num_regions, args.region_length, ", chr".join(chroms),
                total_background_bp(background_windows))

    source_stats = {
        "1000genomes": fetch_and_count_1000genomes(background_windows, args),
        "gnomad": fetch_and_count_gnomad(background_windows, args),
        "pangenome": count_local_vcf_background(args.pangenome_vcf, background_windows),
        "dbsnp": count_local_bigbed_background(args.dbsnp_bb, background_windows),
    }

    write_output(rows, fieldnames, source_stats, args.output_csv)
    logger.info("Done. Intergenic-based enrichment written to %s", args.output_csv)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Recompute SNP enrichment (gene body density / sampled intergenic-region density) "
                     "from a count_snp_density.py summary CSV, using randomly sampled 150nt intergenic "
                     "background regions distributed proportionally across chr1-22/X/Y."
    )
    parser.add_argument("-i", "--input-csv", default="results/snp_density_summary.csv",
                         help="Summary CSV produced by count_snp_density.py")
    parser.add_argument("-o", "--output-csv", default="results/snp_enrichment_summary_intergenic.csv")
    parser.add_argument("--gene-annotation", default="data/downloads/gencode.v50.annotation.gtf",
                         help="Genome-wide gene annotation: GTF/GFF3 with 'gene' feature rows "
                              "(e.g. GENCODE/Ensembl/RefSeq, .gtf/.gff3, optionally .gz), or a BED of gene loci. "
                              "Used to find EVERY annotated gene genome-wide, not just the genes in --input-csv.")
    parser.add_argument("--gaps-bed", default=None,
                         help="Optional BED3 of known assembly gaps (e.g. UCSC GRCh38 'gap' table exported as "
                              "chrom/chromStart/chromEnd). Sampled regions will avoid these too.")
    parser.add_argument("--gene-buffer-bp", type=int, default=10000,
                         help="Minimum distance (bp) a sampled region must keep from the start/end of any "
                              "annotated gene (default: 10000, i.e. 10kbp).")
    parser.add_argument("--num-regions", type=int, default=10000,
                         help="Total number of background regions to sample genome-wide, distributed across "
                              "chr1-22/X/Y in proportion to chromosome length (default: 10000).")
    parser.add_argument("--region-length", type=int, default=150,
                         help="Length in bp of each sampled background region (default: 150).")
    parser.add_argument("--genome-build", choices=["GRCh38", "GRCh37"], default="GRCh38")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for background region sampling")
    parser.add_argument("--vcf-dir", default="data/downloads/1000genomes_background_vcf",
                         help="Directory for fetched 1000 Genomes background VCFs")
    parser.add_argument("--gnomad-source", choices=["gcs", "api"], default="gcs")
    parser.add_argument("--gnomad-dir", default="data/downloads/gnomad_background_vcf",
                         help="Directory for fetched gnomAD background VCFs when --gnomad-source=gcs")
    parser.add_argument("--gnomad-max-region", type=int, default=cse.GNOMAD_API_MAX_REGION_SPAN,
                         help="Max bp span per GraphQL query when --gnomad-source=api")
    parser.add_argument("--pangenome-vcf", default="data/downloads/decomposed.vcf",
                         help="Local decomposed pangenome VCF (already downloaded)")
    parser.add_argument("--dbsnp-bb", default="data/downloads/dbSnp155.bb",
                         help="Local dbSNP BigBed file (already downloaded)")
    parser.add_argument("--skip-fetch", action="store_true",
                         help="Reuse whatever's already in --vcf-dir / --gnomad-dir instead of re-fetching")
    parser.add_argument("--force", action="store_true", help="Re-fetch even if output files already exist")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main():
    run(parse_args())


if __name__ == "__main__":
    main()