"""count_snp_density.py — fetch SNP data for a gene list from four sources
(1000 Genomes, gnomAD, a local pangenome VCF, and a local dbSNP BigBed) and
report, per gene and per source, the SNP count and SNP density (SNPs per bp)
within the gene body.

Usage:
    python count_snp_density.py --input-csv data/hgnc_genes_loc.csv
    python count_snp_density.py --skip-fetch   # reuse already-fetched files, just re-count

Data sources:
  - 1000 Genomes: remote variant-call VCFs (one per chromosome), fetched via
    `bcftools view -R <regions.bed> <remote url>` — a BED file rather than a
    single -r region string, because -r strings listing hundreds of
    individual gene regions can silently return zero matches. Requires
    bcftools (built with libcurl).
  - gnomAD: two fetch strategies are supported via --gnomad-source:
      * "gcs" (default): bcftools view -R against gnomAD's public per-chromosome
        genomes sites VCFs on Google Cloud Storage (exomes are not fetched).
      * "api": gnomAD's GraphQL API, queried per gene region (legacy path),
        chunked to stay under gnomAD's region-size limit, with GraphQL
        "errors" responses surfaced instead of silently treated as zero.
  - pangenome / dbSNP: local files you've already downloaded (a decomposed
    VCF and a BigBed respectively) — this script counts them but does not
    fetch them. Requires pyBigWig for the BigBed source.

Output: one row per gene with gene_length plus, for each source,
snp_count_<source> and snp_density_<source> (count / gene_length).
"""

import argparse
import csv
import glob
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict

import requests

logger = logging.getLogger("count_snp_density")

GRCH38_TEMPLATE = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/"
    "release/20190312_biallelic_SNV_and_INDEL/"
    "ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
)
GRCH37_TEMPLATE = (
    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
)

# gnomAD v4.1, GRCh38 genomes track only. Contigs are "chr"-prefixed, unlike
# the bare "1"/"2" contigs in the 1000 Genomes files above.
GNOMAD_GCS_GENOME_TEMPLATE = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/"
    "gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz"
)
# gnomAD v4 has no chrY genomes track and chrM is handled separately.
GNOMAD_GENOME_SKIP_CHROMS = {"y", "m", "mt"}

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
GNOMAD_API_QUERY = """
query RegionQuery($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4_genomes) { variant_id }
  }
}
"""
# The GraphQL API rejects region queries above a certain span (HTTP 200 with
# a populated "errors" array), so long genes are queried in chunks.
GNOMAD_API_MAX_REGION_SPAN = 300_000

# Real, primary-assembly chromosome names only: 1-22, X, Y, MT/M, optionally
# "chr"-prefixed. Rows like "1 Unlocalized Scaffold" are skipped.
_STANDARD_CHROM_RE = re.compile(r"^(chr)?([1-9]|1[0-9]|2[0-2]|x|y|mt|m)$", re.IGNORECASE)


# ==========================================================================
# Small helpers
# ==========================================================================

def setup_logging(verbose=False):
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
    logger.addHandler(handler)


def normalize_chrom(chrom, style="bare"):
    """Normalize a chromosome name to either 'bare' (1, X) or 'chr' (chr1, chrX) style."""
    bare = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return bare if style == "bare" else f"chr{bare}"


def is_standard_chromosome(chrom):
    return bool(_STANDARD_CHROM_RE.match(str(chrom).strip()))


def read_gene_ranges(csv_path):
    """Load gene ranges from a CSV with GeneName, Chromosome, Start, End
    columns (GeneGroup is used if present). Drops rows whose chromosome
    isn't a real contig."""
    logger.info("Reading gene ranges from %s", csv_path)
    with open(csv_path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    genes, skipped = [], []
    for row in rows:
        chrom = str(row["Chromosome"]).strip()
        gene = {
            "gene_name": row["GeneName"],
            "gene_group": row.get("GeneGroup") or "unknown",
            "chromosome": chrom,
            "start": int(row["Start"]),
            "end": int(row["End"]),
        }
        (genes if is_standard_chromosome(chrom) else skipped).append(gene)

    if skipped:
        sample = ", ".join(f'{g["gene_name"]} ({g["chromosome"]!r})' for g in skipped[:10])
        more = f" and {len(skipped) - 10} more" if len(skipped) > 10 else ""
        logger.warning("Skipping %d gene(s) with non-standard chromosome names: %s%s", len(skipped), sample, more)

    genes.sort(key=lambda g: (normalize_chrom(g["chromosome"]).upper(), g["start"]))
    logger.info("Loaded %d gene ranges", len(genes))
    return genes


def merge_regions(regions):
    """Merge overlapping/adjacent (start, end) regions."""
    merged = []
    for start, end in sorted(regions):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def split_region(start, end, max_span):
    """Split (start, end) into consecutive chunks no larger than max_span."""
    chunks = []
    cur = start
    while cur <= end:
        chunk_end = min(end, cur + max_span - 1)
        chunks.append((cur, chunk_end))
        cur = chunk_end + 1
    return chunks


def gene_regions_by_chrom(genes, chrom_style="bare"):
    """Group each gene's (start, end) region by chromosome."""
    by_chrom = defaultdict(list)
    for gene in genes:
        chrom = normalize_chrom(gene["chromosome"], style=chrom_style)
        by_chrom[chrom].append((gene["start"], gene["end"]))
    return by_chrom


# ==========================================================================
# Fetch: remote VCF regions via `bcftools view -R <regions.bed>`
# (shared by 1000 Genomes and the gnomAD "gcs" source)
# ==========================================================================

def check_bcftools():
    if shutil.which("bcftools") is None:
        logger.error(
            "bcftools not found on PATH. Install it (e.g. `conda install -c bioconda bcftools` "
            "or `apt-get install bcftools`) — it's required to fetch remote VCF regions "
            "efficiently without downloading full chromosome files."
        )
        return False
    return True


def fetch_remote_vcf_regions(remote_url, region_chrom, regions, output_vcf, label, retries=3, delay=5):
    """Fetch all `regions` (1-based, inclusive) on one chromosome from one
    remote, indexed VCF in a single bcftools call, using a BED file (-R)."""
    os.makedirs(os.path.dirname(output_vcf) or ".", exist_ok=True)
    bed_fd, bed_path = tempfile.mkstemp(suffix=".bed", prefix=f"{region_chrom}_")
    try:
        with os.fdopen(bed_fd, "w") as fh:
            for start, end in regions:  # BED is 0-based, half-open
                fh.write(f"{region_chrom}\t{start - 1}\t{end}\n")

        cmd = ["bcftools", "view", "-R", bed_path, "-Ov", "-o", output_vcf, remote_url]
        logger.info("Fetching %s %s (%d region(s)) from %s", label, region_chrom, len(regions), remote_url)

        for attempt in range(1, retries + 1):
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=600)
                if os.path.getsize(output_vcf) == 0:
                    logger.warning("bcftools wrote an empty file for %s %s (attempt %d/%d) — retrying.",
                                   label, region_chrom, attempt, retries)
                else:
                    logger.info("Wrote %s", output_vcf)
                    return output_vcf
            except subprocess.CalledProcessError as exc:
                logger.warning("bcftools failed for %s %s (attempt %d/%d): %s", label, region_chrom, attempt, retries,
                               exc.stderr.strip()[-500:] if exc.stderr else exc)
            except subprocess.TimeoutExpired:
                logger.warning("bcftools timed out for %s %s (attempt %d/%d)", label, region_chrom, attempt, retries)
            if attempt < retries:
                time.sleep(delay)
    finally:
        os.remove(bed_path)

    logger.error("Giving up on %s %s after %d attempts.", label, region_chrom, retries)
    return None


# ==========================================================================
# Fetch: 1000 Genomes VCFs, one per chromosome
# ==========================================================================

def fetch_1000genomes(genes, output_dir, genome_build, force=False):
    """Fetch 1000 Genomes calls for every gene body. Returns the list of VCF
    paths written (one per chromosome)."""
    if not check_bcftools():
        return []

    template = GRCH38_TEMPLATE if genome_build == "GRCh38" else GRCH37_TEMPLATE
    # Both releases name contigs bare ("1", "2", ..., "X"), not "chr1".
    by_chrom = gene_regions_by_chrom(genes, chrom_style="bare")

    written = []
    for chrom, regions in sorted(by_chrom.items()):
        output_vcf = os.path.join(output_dir, f"1000genomes_chr{chrom}.vcf")
        if not force and os.path.exists(output_vcf) and os.path.getsize(output_vcf) > 0:
            logger.info("Skipping chr%s — %s already exists.", chrom, output_vcf)
            written.append(output_vcf)
            continue
        remote_url = template.format(chrom=chrom)
        result = fetch_remote_vcf_regions(remote_url, chrom, merge_regions(regions), output_vcf, "1000genomes")
        if result:
            written.append(result)

    logger.info("1000 Genomes fetch done. Wrote %d VCF file(s) to %s", len(written), output_dir)
    return written


# ==========================================================================
# Fetch: gnomAD — "gcs" source (default)
# ==========================================================================

def fetch_gnomad_gcs(genes, output_dir, force=False):
    """Fetch gnomAD genomes sites VCFs (exomes excluded) for every gene body
    via bcftools. Returns the list of VCF paths written."""
    if not check_bcftools():
        return []

    by_chrom = gene_regions_by_chrom(genes, chrom_style="chr")

    written = []
    for chrom, regions in sorted(by_chrom.items()):
        bare = normalize_chrom(chrom, style="bare").lower()
        if bare in GNOMAD_GENOME_SKIP_CHROMS:
            logger.info("Skipping gnomAD genomes %s — no genomes track for this contig.", chrom)
            continue

        genome_vcf = os.path.join(output_dir, f"gnomad_genomes_{chrom}.vcf")
        if not force and os.path.exists(genome_vcf) and os.path.getsize(genome_vcf) > 0:
            logger.info("Skipping gnomAD genomes %s — %s already exists.", chrom, genome_vcf)
            written.append(genome_vcf)
            continue

        url = GNOMAD_GCS_GENOME_TEMPLATE.format(chrom=normalize_chrom(chrom, style="bare"))
        result = fetch_remote_vcf_regions(url, chrom, merge_regions(regions), genome_vcf, "gnomad-genomes")
        if result:
            written.append(result)

    logger.info("gnomAD genomes (gcs) fetch done. Wrote %d VCF file(s) to %s", len(written), output_dir)
    return written


# ==========================================================================
# Fetch: gnomAD — "api" source (legacy): GraphQL, per gene region
# ==========================================================================

def fetch_gnomad_api(genes, output_vcf, retries=5, delay=5, request_delay=1.0,
                     max_region_span=GNOMAD_API_MAX_REGION_SPAN, force=False):
    """Query gnomAD's GraphQL API for every gene body and write matches to a
    VCF. Oversized genes are chunked; GraphQL "errors" are logged and retried."""
    if not force and os.path.exists(output_vcf) and os.path.getsize(output_vcf) > 0:
        logger.info("Skipping gnomAD fetch — %s already exists.", output_vcf)
        return output_vcf

    session = requests.Session()
    rows = {}  # variant_id -> (chrom, pos); dict avoids duplicates from overlapping genes
    failed_genes = []

    for idx, gene in enumerate(genes, start=1):
        chrom = normalize_chrom(gene["chromosome"], style="bare")
        chunks = split_region(gene["start"], gene["end"], max_region_span)
        logger.info("gnomAD %d/%d: %s (%s:%d-%d, %d chunk(s))",
                    idx, len(genes), gene["gene_name"], chrom, gene["start"], gene["end"], len(chunks))

        gene_ok = True
        for chunk_start, chunk_end in chunks:
            variables = {"chrom": chrom, "start": chunk_start, "stop": chunk_end}
            chunk_ok = False

            for attempt in range(1, retries + 1):
                try:
                    response = session.post(GNOMAD_API_URL, json={"query": GNOMAD_API_QUERY, "variables": variables}, timeout=120)
                except requests.RequestException as exc:
                    logger.warning("gnomAD network error for %s (attempt %d/%d): %s", gene["gene_name"], attempt, retries, exc)
                    time.sleep(delay)
                    continue

                time.sleep(request_delay)
                if response.status_code == 200:
                    data = response.json()
                    errors = data.get("errors")
                    if errors:
                        messages = "; ".join(e.get("message", str(e)) for e in errors)
                        logger.warning("gnomAD returned errors for %s %s:%d-%d (attempt %d/%d): %s",
                                       gene["gene_name"], chrom, chunk_start, chunk_end, attempt, retries, messages)
                        time.sleep(delay)
                        continue
                    variants = ((data.get("data") or {}).get("region") or {}).get("variants") or []
                    for v in variants:
                        variant_id = v.get("variant_id", "")
                        parts = variant_id.split("-")
                        if len(parts) >= 2:
                            rows[variant_id] = (parts[0], parts[1])
                    chunk_ok = True
                    break
                elif response.status_code in (429, 500, 502, 503, 504):
                    logger.warning("gnomAD HTTP %d for %s (attempt %d/%d) — retrying.", response.status_code, gene["gene_name"], attempt, retries)
                    time.sleep(delay * attempt)
                else:
                    logger.error("gnomAD HTTP %d for %s — not retrying.", response.status_code, gene["gene_name"])
                    break

            if not chunk_ok:
                gene_ok = False

        if not gene_ok:
            failed_genes.append(gene["gene_name"])

    if failed_genes:
        sample = ", ".join(failed_genes[:10])
        more = f" and {len(failed_genes) - 10} more" if len(failed_genes) > 10 else ""
        logger.warning("gnomAD fetch had unresolved errors for %d gene(s): %s%s — their counts will be "
                       "artificially low. Re-run with --force to retry, or use --gnomad-source gcs instead.",
                       len(failed_genes), sample, more)

    os.makedirs(os.path.dirname(output_vcf) or ".", exist_ok=True)
    sorted_rows = sorted(((c, int(p), vid) for vid, (c, p) in rows.items()),
                         key=lambda r: (normalize_chrom(r[0]).upper(), r[1]))
    with open(output_vcf, "w", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, variant_id in sorted_rows:
            fh.write(f"{chrom}\t{pos}\t{variant_id}\t.\t.\t.\t.\t.\n")
    logger.info("Wrote %d gnomAD variant(s) to %s", len(sorted_rows), output_vcf)
    return output_vcf


# ==========================================================================
# Count: SNPs in each gene body
# ==========================================================================

def build_chrom_windows(genes):
    """Group each gene's body window by chromosome, sorted by start."""
    windows_by_chrom = defaultdict(list)
    for idx, gene in enumerate(genes):
        chrom = normalize_chrom(gene["chromosome"], style="bare").upper()
        windows_by_chrom[chrom].append({"idx": idx, "start": gene["start"], "end": gene["end"]})
    for windows in windows_by_chrom.values():
        windows.sort(key=lambda w: w["start"])
    return windows_by_chrom


def extract_positions(path, windows_by_chrom=None):
    """Read a source file exactly once and return {chrom: [positions]}.
    Positions are NOT assumed sorted — sweep_positions sorts them."""
    ext = os.path.splitext(path.lower())[1]
    positions_by_chrom = defaultdict(list)

    if ext == ".vcf":
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.split("\t", 2)
                if len(fields) < 2:
                    continue
                try:
                    pos = int(fields[1])
                except ValueError:
                    continue
                positions_by_chrom[normalize_chrom(fields[0]).upper()].append(pos)
        return positions_by_chrom

    if ext in (".bb", ".bigbed"):
        # dbSNP-scale BigBed files hold hundreds of millions of entries, so
        # only the merged gene-body windows are queried (random access)
        # rather than whole chromosomes.
        if windows_by_chrom is None:
            raise ValueError("BigBed sources require windows_by_chrom to query only gene regions")
        import pyBigWig
        bb = pyBigWig.open(path)
        try:
            bb_chrom_sizes = bb.chroms()
            for chrom, windows in windows_by_chrom.items():
                chrom_name = chrom if chrom in bb_chrom_sizes else normalize_chrom(chrom, style="chr")
                if chrom_name not in bb_chrom_sizes:
                    continue  # this BigBed has no data for this chromosome under either naming
                chrom_size = int(bb_chrom_sizes[chrom_name])
                for start, end in merge_regions([(w["start"], w["end"]) for w in windows]):
                    q_start = max(0, start - 1)  # BigBed is 0-based, half-open
                    q_end = min(chrom_size, end)
                    if q_end <= q_start:
                        continue
                    for entry in bb.entries(chrom_name, q_start, q_end) or []:
                        positions_by_chrom[chrom].append(int(entry[0]) + 1)  # -> 1-based
        finally:
            bb.close()
        return positions_by_chrom

    raise ValueError(f"Unsupported source file type: {path}")


def sweep_positions(positions_by_chrom, windows_by_chrom):
    """For each chromosome, sort its variant positions and sweep them against
    that chromosome's sorted gene windows, keeping a small 'active' list of
    genes the current position could still fall in (genes can overlap, so a
    variant may count toward more than one gene). Returns {gene_idx: count}."""
    results = defaultdict(int)
    for chrom, positions in positions_by_chrom.items():
        windows = windows_by_chrom.get(chrom)
        if not windows:
            continue
        positions.sort()

        next_idx, active = 0, []
        for pos in positions:
            while next_idx < len(windows) and windows[next_idx]["start"] <= pos:
                active.append(windows[next_idx])
                next_idx += 1
            active = [w for w in active if w["end"] >= pos]
            for w in active:
                results[w["idx"]] += 1
    return results


def count_source(path, windows_by_chrom):
    positions_by_chrom = extract_positions(path, windows_by_chrom)
    results = sweep_positions(positions_by_chrom, windows_by_chrom)
    total = sum(len(p) for p in positions_by_chrom.values())
    logger.info("%s: %d variant(s) read, %d gene(s) matched", path, total, len(results))
    return results


def write_summary(genes, counts_by_label, output_csv):
    """counts_by_label: {label: {gene_idx: count}, ...} — one entry per
    logical source. Writes one row per gene with each source's count and
    density (count / gene_length, SNPs per bp)."""
    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
    labels = list(counts_by_label.keys())

    fieldnames = ["gene_name", "gene_group", "chromosome", "start", "end", "gene_length"]
    for label in labels:
        fieldnames += [f"snp_count_{label}", f"snp_density_{label}"]

    with open(output_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for idx, gene in enumerate(genes):
            gene_length = gene["end"] - gene["start"] + 1
            row = {
                "gene_name": gene["gene_name"], "gene_group": gene["gene_group"], "chromosome": gene["chromosome"],
                "start": gene["start"], "end": gene["end"], "gene_length": gene_length,
            }
            for label in labels:
                count = counts_by_label[label].get(idx, 0)
                row[f"snp_count_{label}"] = count
                row[f"snp_density_{label}"] = count / gene_length if gene_length > 0 else 0.0
            writer.writerow(row)
    logger.info("Wrote %d rows to %s", len(genes), output_csv)


# ==========================================================================
# Pipeline + CLI
# ==========================================================================

def run_pipeline(input_csv, vcf_dir, genome_build, gnomad_source, gnomad_dir, gnomad_vcf,
                 gnomad_max_region, pangenome_vcf, dbsnp_bb, output_csv, skip_fetch=False, force=False):
    if not os.path.exists(input_csv):
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    genes = read_gene_ranges(input_csv)

    if not skip_fetch:
        written = fetch_1000genomes(genes, vcf_dir, genome_build, force=force)
        if not written:
            logger.error("1000 Genomes fetch produced no VCFs — aborting.")
            sys.exit(1)
        if gnomad_source == "gcs":
            gnomad_written = fetch_gnomad_gcs(genes, gnomad_dir, force=force)
            if not gnomad_written:
                logger.error("gnomAD (gcs) fetch produced no VCFs — aborting.")
                sys.exit(1)
        else:
            fetch_gnomad_api(genes, gnomad_vcf, max_region_span=gnomad_max_region, force=force)

    thousand_genomes_vcfs = sorted(glob.glob(os.path.join(vcf_dir, "*.vcf")))
    if not thousand_genomes_vcfs:
        logger.error("No 1000 Genomes VCFs found in %s. Run without --skip-fetch first.", vcf_dir)
        sys.exit(1)

    if gnomad_source == "gcs":
        gnomad_vcfs = sorted(glob.glob(os.path.join(gnomad_dir, "gnomad_genomes_*.vcf")))
        if not gnomad_vcfs:
            logger.error("No gnomAD genomes VCFs found in %s. Run without --skip-fetch first.", gnomad_dir)
            sys.exit(1)
    else:
        gnomad_vcfs = [gnomad_vcf]
        if not os.path.exists(gnomad_vcf):
            logger.error("gnomAD VCF not found at %s.", gnomad_vcf)
            sys.exit(1)

    for label, path in [("pangenome VCF", pangenome_vcf), ("dbSNP BigBed", dbsnp_bb)]:
        if not os.path.exists(path):
            logger.error("%s not found at %s.", label, path)
            sys.exit(1)

    windows_by_chrom = build_chrom_windows(genes)

    # Per-chromosome VCFs for 1000 Genomes (and gnomAD "gcs") share one
    # logical source, so their counts are summed under a single label.
    def pooled_counts(paths):
        pooled = defaultdict(int)
        for path in paths:
            for idx, count in count_source(path, windows_by_chrom).items():
                pooled[idx] += count
        return pooled

    counts_by_label = {
        "1000genomes": pooled_counts(thousand_genomes_vcfs),
        "gnomad": pooled_counts(gnomad_vcfs),
        "pangenome": count_source(pangenome_vcf, windows_by_chrom),
        "dbsnp": count_source(dbsnp_bb, windows_by_chrom),
    }

    write_summary(genes, counts_by_label, output_csv)
    logger.info("Pipeline complete. Summary written to %s", output_csv)


def main():
    parser = argparse.ArgumentParser(
        description="Fetch SNP data for gene regions from 1000 Genomes and gnomAD, count SNPs "
                    "from pangenome and dbSNP files, then report per-gene, per-source SNP "
                    "count and density within the gene body."
    )
    parser.add_argument("-i", "--input-csv", default="data/hgnc_genes_loc.csv", help="Gene range CSV (GeneName, Chromosome, Start, End, optional GeneGroup)")
    parser.add_argument("-o", "--output-csv", default="results/snp_density_summary.csv", help="Count/density summary CSV")
    parser.add_argument("--vcf-dir", default="data/downloads/1000genomes_vcf", help="Directory for fetched 1000 Genomes VCFs (one per chromosome)")
    parser.add_argument("--genome-build", choices=["GRCh38", "GRCh37"], default="GRCh38")
    parser.add_argument("--gnomad-source", choices=["gcs", "api"], default="gcs",
                        help="How to fetch gnomAD: 'gcs' (default) uses bcftools against gnomAD's public "
                             "genomes sites VCFs (exomes excluded) — robust, no rate limits. 'api' uses "
                             "gnomAD's GraphQL API, one request per gene region.")
    parser.add_argument("--gnomad-dir", default="data/downloads/gnomad_vcf", help="Directory for fetched gnomAD genomes VCFs when --gnomad-source=gcs (one per chromosome)")
    parser.add_argument("--gnomad-vcf", default="data/downloads/gnomad_region_data.vcf", help="Fetched gnomAD data when --gnomad-source=api, written as a single VCF")
    parser.add_argument("--gnomad-max-region", type=int, default=GNOMAD_API_MAX_REGION_SPAN,
                        help="Max bp span per GraphQL query when --gnomad-source=api; longer genes are chunked")
    parser.add_argument("--pangenome-vcf", default="data/downloads/decomposed.vcf", help="Local decomposed pangenome VCF (already downloaded; not fetched by this script)")
    parser.add_argument("--dbsnp-bb", default="data/downloads/dbSnp155.bb", help="Local dbSNP BigBed file (already downloaded; not fetched by this script)")
    parser.add_argument("--skip-fetch", action="store_true", help="Reuse whatever's already in --vcf-dir / --gnomad-dir / --gnomad-vcf instead of re-fetching (pangenome/dbSNP are always read locally)")
    parser.add_argument("--force", action="store_true", help="Re-fetch even if output files already exist")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    setup_logging(verbose=args.verbose)
    run_pipeline(
        args.input_csv, args.vcf_dir, args.genome_build,
        args.gnomad_source, args.gnomad_dir, args.gnomad_vcf, args.gnomad_max_region,
        args.pangenome_vcf, args.dbsnp_bb, args.output_csv,
        skip_fetch=args.skip_fetch, force=args.force,
    )


if __name__ == "__main__":
    main()