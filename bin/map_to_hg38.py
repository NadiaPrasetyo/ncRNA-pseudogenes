#!/usr/bin/env python3
"""
map_to_hg38.py

Simplified pipeline: give it a directory of .fastq.gz files, a genome fasta,
and a GTF -- everything else (indices, BAMs, peaks, counts) is derived and
cached automatically under --work-dir.

Samples are discovered directly from the fastq directory (no manifest files
needed). Each *.fastq.gz file's name (minus the extension) is the sample
name, e.g.:

    data/fastq/GSE20309_GSM509055_HeLa_TotalRNAseq.fastq.gz
        -> sample "GSE20309_GSM509055_HeLa_TotalRNAseq"

A sample is treated as RNA-seq if "rnaseq" (case-insensitive) appears
anywhere in its name; everything else is treated as ChIP-seq.

For every sample this script will:
  1. Align to the given genome
       - ChIP-seq samples -> Bowtie2 (ungapped, single-end)
       - RNA-seq samples  -> STAR (splice-aware)
  2. Sort + index the resulting BAM with samtools
  3. Annotate to genes using the given GTF
       - ChIP-seq: call peaks with genrich (using a matched Input/control
         sample where one is known), then assign each peak to its
         nearest gene with bedtools
       - RNA-seq: gene-level read counts with featureCounts

Requires on PATH: bowtie2, bowtie2-build, STAR, samtools, genrich, bedtools,
                  featureCounts (subread)
  conda install -c bioconda bowtie2 star samtools genrich bedtools subread

Usage:
    python3 map_to_hg38.py --genome-fasta GRCh38.p14.genome.fa \\
                           --gtf gencode.v50.chr_patch_hapl_scaff.annotation.gtf \\
                           --fastq-dir data/fastq
    python3 map_to_hg38.py --genome-fasta ... --gtf ... --fastq-dir ... chip   # ChIP-seq only
    python3 map_to_hg38.py --genome-fasta ... --gtf ... --fastq-dir ... rna    # RNA-seq only

Output (all under --work-dir, default "results/"):
    ref_work/                       derived indices + gene BED (built once, cached)
    bam/<sample>.sorted.bam(.bai)   per-sample alignments
    peaks/<sample>_peaks.annotated.tsv   peaks joined to nearest gene (ChIP-seq)
    counts/gene_counts.tsv          featureCounts gene-level matrix (RNA-seq)
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
import logging

# =====================================================================
# Matched Input/control sample per ChIP-seq sample, where known. Anything
# not listed here is peak-called without a control.
# =====================================================================
CHIP_CONTROL = {
    "GSE18184_GSM454595_IMR90hTert_RPC4": "GSE18184_GSM454598_IMR90hTert_Input",
    "GSE18184_GSM454596_IMR90hTert_Bdp1": "GSE18184_GSM454598_IMR90hTert_Input",
    "GSE18184_GSM454597_IMR90hTert_Brf1": "GSE18184_GSM454598_IMR90hTert_Input",
    "GSE18184_GSM454599_IMR90hTert_SNAP45": "GSE18184_GSM454598_IMR90hTert_Input",
    "GSE20309_GSM509049_HeLa_Brf1": "GSE20309_GSM509051_HeLa_Input",
    "GSE20309_GSM509050_HeLa_Brf2": "GSE20309_GSM509051_HeLa_Input",
    "GSE20309_GSM509052_HeLa_Pol3_Rep1": "GSE20309_GSM509051_HeLa_Input",
    "GSE20309_GSM509053_HeLa_Pol3_Rep2": "GSE20309_GSM509051_HeLa_Input",
    "GSE20309_GSM509056_HeLa_TFIIIC": "GSE20309_GSM509051_HeLa_Input",
}

# Samples that ARE inputs/controls themselves - never peak-call these.
IS_INPUT = {
    "GSE18184_GSM454598_IMR90hTert_Input",
    "GSE20309_GSM509051_HeLa_Input",
}


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

import os  # add to the top-level imports

def run(cmd, **kwargs):
    logging.info("+ " + " ".join(str(c) for c in cmd))
    try:
        subprocess.run(cmd, check=True, **kwargs)
    except subprocess.CalledProcessError as e:
        # Find whichever stream was redirected to a real file, so we can
        # show the user what actually went wrong.
        log_file = None
        for key in ("stderr", "stdout"):
            candidate = kwargs.get(key)
            if hasattr(candidate, "name"):
                log_file = candidate
                break
        if log_file is not None:
            try:
                log_file.flush()
                os.fsync(log_file.fileno())
            except (OSError, ValueError):
                pass
            try:
                with open(log_file.name) as f:
                    tail = f.readlines()[-20:]
            except OSError:
                tail = []
            if tail:
                logging.error(
                    f"Command failed (exit {e.returncode}). Last lines of {log_file.name}:\n"
                    + "".join(tail)
                )
            else:
                sig_note = f" (killed by signal {-e.returncode})" if e.returncode < 0 else ""
                logging.error(
                    f"Command failed with exit code {e.returncode}{sig_note} and wrote "
                    f"nothing to {log_file.name}."
                )
        raise

def check_tools(tools):
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        sys.exit(f"Missing required tool(s) on PATH: {', '.join(missing)}")

# =====================================================================
# Sample discovery (straight from the fastq dir, no manifests)
# =====================================================================
def discover_samples(fastq_dir: Path):
    """Return (chip_samples, rna_samples) lists, sorted, from *.fastq.gz files."""
    chip, rna = [], []
    for fq in sorted(fastq_dir.glob("*.fastq.gz")):
        sample = fq.name[: -len(".fastq.gz")]
        if "rnaseq" in sample.lower():
            rna.append(sample)
        else:
            chip.append(sample)
    return chip, rna


# =====================================================================
# Reference prep (local files only - nothing is downloaded)
# =====================================================================
def build_gene_bed(gtf_path: Path, gene_bed_path: Path) -> None:
    if gene_bed_path.exists():
        logging.info(f"Gene BED already present at {gene_bed_path}, skipping.")
        return

    logging.info(f"Deriving gene-level BED from {gtf_path.name} for peak annotation...")
    rows = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom, _, _, start, end, _, strand, _, attrs = fields

            gene_id = ""
            gene_name = ""
            for part in attrs.strip().split(";"):
                part = part.strip()
                if part.startswith("gene_id"):
                    gene_id = part.split('"')[1]
                elif part.startswith("gene_name"):
                    gene_name = part.split('"')[1]
            if not gene_name:
                gene_name = gene_id

            rows.append((chrom, int(start) - 1, int(end), gene_name, strand, gene_id))

    rows.sort(key=lambda r: (r[0], r[1]))
    gene_bed_path.parent.mkdir(parents=True, exist_ok=True)
    with open(gene_bed_path, "w") as out:
        for chrom, start, end, gene_name, strand, gene_id in rows:
            out.write(f"{chrom}\t{start}\t{end}\t{gene_name}\t.\t{strand}\t{gene_id}\n")
    logging.info(f"{len(rows)} gene records written to {gene_bed_path}")


def build_bowtie2_index(genome_fasta: Path, index_prefix: Path, threads: int) -> None:
    if Path(f"{index_prefix}.1.bt2").exists() or Path(f"{index_prefix}.1.bt2l").exists():
        logging.info("Bowtie2 index already built, skipping.")
        return
    logging.info("Building Bowtie2 index (this takes a while for a whole genome)...")
    index_prefix.parent.mkdir(parents=True, exist_ok=True)
    run(["bowtie2-build", "--threads", str(threads), str(genome_fasta), str(index_prefix)])


def build_star_index(genome_fasta: Path, gtf_path: Path, star_index_dir: Path, threads: int) -> None:
    if (star_index_dir / "SAindex").exists():
        logging.info("STAR index already built, skipping.")
        return
    logging.info("Building STAR index...")
    star_index_dir.mkdir(parents=True, exist_ok=True)
    run([
        "STAR", "--runMode", "genomeGenerate",
        "--runThreadN", str(threads),
        "--genomeDir", str(star_index_dir),
        "--genomeFastaFiles", str(genome_fasta),
        "--sjdbGTFfile", str(gtf_path),
        "--sjdbOverhang", "100",
    ])


# =====================================================================
# Alignment
# =====================================================================
def align_chip_sample(sample: str, fastq_dir: Path, bam_dir: Path, index_prefix: Path, threads: int) -> None:
    fq = fastq_dir / f"{sample}.fastq.gz"
    bam = bam_dir / f"{sample}.sorted.bam"

    if bam.exists():
        logging.info(f"{sample} already aligned, skipping.")
        return

    logging.info(f"Aligning (Bowtie2, single-end) {sample} ...")
    bt2_log = bam_dir / f"{sample}.bowtie2.log"
    bowtie2_cmd = ["bowtie2", "-x", str(index_prefix), "-U", str(fq),
                   "-p", str(threads), "--no-unal"]
    sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", str(bam), "-"]

    with open(bt2_log, "w") as errlog:
        bt2 = subprocess.Popen(bowtie2_cmd, stdout=subprocess.PIPE, stderr=errlog)
        run(sort_cmd, stdin=bt2.stdout)
        bt2.stdout.close()
        bt2.wait()
        if bt2.returncode != 0:
            sys.exit(f"bowtie2 failed for {sample}, see {bt2_log}")

    run(["samtools", "index", str(bam)])
    flagstat = subprocess.run(["samtools", "flagstat", str(bam)],
                               capture_output=True, text=True).stdout.splitlines()[0]
    logging.info(f"{sample} -> {bam} ({flagstat})")


def align_rna_sample(sample: str, fastq_dir: Path, bam_dir: Path, star_index_dir: Path, threads: int) -> Path:
    fq = fastq_dir / f"{sample}.fastq.gz"
    prefix = bam_dir / f"{sample}."
    bam = Path(f"{prefix}Aligned.sortedByCoord.out.bam")

    if bam.exists():
        logging.info(f"{sample} already aligned, skipping.")
        return bam

    logging.info(f"Aligning (STAR, splice-aware) {sample} ...")
    run([
        "STAR", "--runMode", "alignReads",
        "--runThreadN", str(threads),
        "--genomeDir", str(star_index_dir),
        "--readFilesIn", str(fq),
        "--readFilesCommand", "zcat",
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outFileNamePrefix", str(prefix),
        "--quantMode", "GeneCounts",
    ])
    run(["samtools", "index", str(bam)])
    logging.info(f"{sample} -> {bam}")
    return bam


# =====================================================================
# Annotation: ChIP-seq (peak calling + nearest-gene assignment)
# =====================================================================
def name_sort_bam(bam: Path, work_dir: Path, threads: int) -> Path:
    """Genrich requires reads sorted by name, not coordinate. Cached alongside peaks/."""
    namesorted = work_dir / f"{bam.stem}.namesorted.bam"
    if namesorted.exists():
        return namesorted
    run(["samtools", "sort", "-n", "-@", str(threads), "-o", str(namesorted), str(bam)])
    return namesorted


def sort_gene_bed(gene_bed: Path, ref_work_dir: Path) -> Path:
    """One-time sorted copy of genes.bed, matching the -k1,1 -k2,2n order
    bedtools expects. Cached so we don't re-sort it per sample."""
    sorted_bed = ref_work_dir / "genes.sorted.bed"
    if sorted_bed.exists():
        return sorted_bed
    with open(sorted_bed, "w") as out:
        run(["sort", "-k1,1", "-k2,2n", str(gene_bed)], stdout=out)
    return sorted_bed


def call_and_annotate_peaks(sample: str, bam_dir: Path, peak_dir: Path, gene_bed: Path, threads: int = 4) -> None:
    if sample in IS_INPUT:
        return  # never peak-call an input/control track

    bam = bam_dir / f"{sample}.sorted.bam"
    if not bam.exists():
        logging.info(f"WARNING: no BAM for {sample}, skipping peak calling.")
        return

    control = CHIP_CONTROL.get(sample)
    control_bam = bam_dir / f"{control}.sorted.bam" if control else None
    narrowpeak = peak_dir / f"{sample}_peaks.narrowPeak"

    if narrowpeak.exists():
        logging.info(f"Peaks already called for {sample}, skipping Genrich.")
    else:
        has_control = control_bam is not None and control_bam.exists()
        logging.info(f"Calling peaks for {sample} "
            f"{'vs control ' + control if has_control else '(no matched control)'}...")

        treatment_ns = name_sort_bam(bam, peak_dir, threads)
        genrich_log = peak_dir / f"{sample}.genrich.log"
        tmp_narrowpeak = peak_dir / f"{sample}_peaks.narrowPeak.tmp"
        cmd = ["Genrich", "-t", str(treatment_ns), "-o", str(tmp_narrowpeak), "-y"]
        if has_control:
            control_ns = name_sort_bam(control_bam, peak_dir, threads)
            cmd += ["-c", str(control_ns)]
        with open(genrich_log, "w") as errlog:
            run(cmd, stdout=errlog, stderr=subprocess.STDOUT)
        # Only becomes the "real" output once Genrich has exited 0, so a
        # crash can never look like a cached, already-finished peak call.
        os.replace(tmp_narrowpeak, narrowpeak)

    if not narrowpeak.exists():
        logging.info(f"WARNING: Genrich produced no narrowPeak for {sample}.")
        return

    # ---- nearest-gene annotation via bedtools closest ----
    annotated = peak_dir / f"{sample}_peaks.annotated.tsv"
    if annotated.exists():
        logging.info(f"{sample} already annotated, skipping.")
        return

    sorted_narrowpeak = peak_dir / f"{sample}_peaks.sorted.narrowPeak"
    with open(sorted_narrowpeak, "w") as out:
        run(["sort", "-k1,1", "-k2,2n", str(narrowpeak)], stdout=out)

    sorted_gene_bed = sort_gene_bed(gene_bed, gene_bed.parent)

    tmp_annotated = peak_dir / f"{sample}_peaks.annotated.tsv.tmp"
    with open(tmp_annotated, "w") as out:
        run([
            "bedtools", "closest",
            "-a", str(sorted_narrowpeak),
            "-b", str(sorted_gene_bed),
            "-d",
        ], stdout=out)
    os.replace(tmp_annotated, annotated)
    logging.info(f"{sample}: peaks annotated -> {annotated}")


# =====================================================================
# Annotation: RNA-seq (gene-level counts via featureCounts)
# =====================================================================
def quantify_rna_seq(rna_samples, bam_dir: Path, counts_dir: Path, gtf_path: Path, threads: int) -> None:
    if not rna_samples:
        logging.info("No RNA-seq samples found, skipping quantification.")
        return

    bams = []
    for sample in rna_samples:
        bam = bam_dir / f"{sample}.Aligned.sortedByCoord.out.bam"
        if bam.exists():
            bams.append(str(bam))

    if not bams:
        logging.info("No RNA-seq BAMs available yet, skipping featureCounts.")
        return

    logging.info(f"Running featureCounts on {len(bams)} RNA-seq sample(s)...")
    counts_dir.mkdir(parents=True, exist_ok=True)
    out_file = counts_dir / "gene_counts.tsv"
    run([
        "featureCounts", "-T", str(threads), "-a", str(gtf_path),
        "-o", str(out_file), "-g", "gene_id", "--extraAttributes", "gene_name",
        *bams,
    ])
    logging.info(f"Gene-level counts -> {out_file}")


# =====================================================================
# Driver
# =====================================================================
def run_chip(chip_samples, fastq_dir: Path, bam_dir: Path, peak_dir: Path, ref_work_dir: Path,
             genome_fasta: Path, gene_bed: Path, threads: int) -> None:
    if not chip_samples:
        logging.info("No ChIP-seq samples found, skipping.")
        return
    check_tools(["bowtie2", "bowtie2-build", "samtools", "Genrich", "bedtools", "sort"])

    index_prefix = ref_work_dir / "bowtie2_index" / genome_fasta.stem
    build_bowtie2_index(genome_fasta, index_prefix, threads)

    bam_dir.mkdir(parents=True, exist_ok=True)
    peak_dir.mkdir(parents=True, exist_ok=True)

    for sample in chip_samples:
        align_chip_sample(sample, fastq_dir, bam_dir, index_prefix, threads)
    for sample in chip_samples:
        try:
            call_and_annotate_peaks(sample, bam_dir, peak_dir, gene_bed, threads)
        except subprocess.CalledProcessError:
            logging.error(f"Peak calling failed for {sample}, skipping (see its .genrich.log). Continuing with remaining samples.")

def run_rna(rna_samples, fastq_dir: Path, bam_dir: Path, counts_dir: Path, ref_work_dir: Path,
            genome_fasta: Path, gtf_path: Path, threads: int) -> None:
    if not rna_samples:
        logging.info("No RNA-seq samples found, skipping.")
        return
    check_tools(["STAR", "samtools", "featureCounts"])
    star_index_dir = ref_work_dir / "star_index"
    build_star_index(genome_fasta, gtf_path, star_index_dir, threads)

    bam_dir.mkdir(parents=True, exist_ok=True)
    for sample in rna_samples:
        align_rna_sample(sample, fastq_dir, bam_dir, star_index_dir, threads)
    quantify_rna_seq(rna_samples, bam_dir, counts_dir, gtf_path, threads)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genome-fasta", required=True, type=Path,
                   help="Path to local hg38 genome fasta")
    p.add_argument("--gtf", required=True, type=Path,
                   help="Path to local GTF annotation")
    p.add_argument("--fastq-dir", default=Path("data/fastq"), type=Path,
                   help="Directory of *.fastq.gz files (sample name = filename minus .fastq.gz)")
    p.add_argument("--work-dir", default=Path("results"), type=Path,
                   help="All outputs (bam/, peaks/, counts/, ref_work/) are written here")
    p.add_argument("--threads", default=4, type=int)
    p.add_argument("targets", nargs="*", choices=["chip", "rna"], default=["chip", "rna"],
                   help="Which stage(s) to run (default: both)")
    return p.parse_args()


if __name__ == "__main__":
    setup_logging()
    args = parse_args()

    if not args.genome_fasta.exists():
        sys.exit(f"Genome fasta not found: {args.genome_fasta}")
    if not args.gtf.exists():
        sys.exit(f"GTF not found: {args.gtf}")
    if not args.fastq_dir.exists():
        sys.exit(f"Fastq dir not found: {args.fastq_dir}")

    ref_work_dir = args.work_dir / "ref_work"
    bam_dir = args.work_dir / "bam"
    peak_dir = args.work_dir / "peaks"
    counts_dir = args.work_dir / "counts"

    ref_work_dir.mkdir(parents=True, exist_ok=True)
    gene_bed = ref_work_dir / "genes.bed"
    build_gene_bed(args.gtf, gene_bed)

    chip_samples, rna_samples = discover_samples(args.fastq_dir)
    logging.info(f"Discovered {len(chip_samples)} ChIP-seq and {len(rna_samples)} RNA-seq sample(s) "
                 f"in {args.fastq_dir}")

    targets = args.targets or ["chip", "rna"]

    if "chip" in targets:
        run_chip(chip_samples, args.fastq_dir, bam_dir, peak_dir, ref_work_dir,
                  args.genome_fasta, gene_bed, args.threads)
    if "rna" in targets:
        run_rna(rna_samples, args.fastq_dir, bam_dir, counts_dir, ref_work_dir,
                args.genome_fasta, args.gtf, args.threads)

    logging.info("All alignment and annotation steps complete.")
    logging.info(f"  BAMs:    {bam_dir}/")
    logging.info(f"  Peaks:   {peak_dir}/*_peaks.annotated.tsv (ChIP-seq, nearest gene per peak)")
    logging.info(f"  Counts:  {counts_dir}/gene_counts.tsv (RNA-seq, gene-level)")