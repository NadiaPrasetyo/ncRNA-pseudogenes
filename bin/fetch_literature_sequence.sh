#!/bin/bash
# fetch_literature_seq_data.sh
#
# Fetches raw sequencing reads for three literature sources:
#   1. Canella et al.: GSE18184 (PRJNA119545, SRR035550-SRR035559) - Pol III machinery ChIP-seq, IMR90hTert
#   2. Oler et al.: GSE20309 (PRJNA125353, SRR036642-SRR036666) - Pol III ChIP-seq + RNA-seq, multiple cell lines
#   3. Moqtaderi et al.: ENCODE Yale/Snyder ChIP-seq raw data (hg18-era, UCSC hgdownload)
#
# Requires: sra-tools (prefetch, fasterq-dump), wget
#   conda install -c bioconda sra-tools
# Make sure sra-tools is up to date and configured (vdb-config --interactive)
# before running - see:
#   https://github.com/ncbi/sra-tools/wiki/Downloads
#   https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration
#
# Usage:
#   ./fetch_literature_seq_data.sh              # fetch all three sources
#   ./fetch_literature_seq_data.sh gse18184      # fetch just one source
#   ./fetch_literature_seq_data.sh gse20309 encode_yale
#
# Output: fastq/<sample_name>.fastq.gz for each sample, plus manifest_*.txt
# files listing sample names per assay type, consumed by align_hg38.sh.

set -euo pipefail

# =====================================================================
# Setup: output directory + logging
# =====================================================================
FASTQ_DIR="data/fastq"
mkdir -p "$FASTQ_DIR"

log () {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# =====================================================================
# Direct HTTP fetch helper (used by encode_yale - no SRA involved)
# =====================================================================
fetch_direct_fastq () {
    local sample_name="$1"
    local url="$2"
    local out_fastq="$FASTQ_DIR/${sample_name}.fastq.gz"

    if [[ -f "$out_fastq" ]]; then
        log "  $sample_name already exists ($out_fastq), skipping."
        return
    fi

    log "  Downloading $sample_name from $url"
    wget -q -O "$out_fastq" "$url"
    log "  -> $out_fastq"
}

# =====================================================================
# SRA download helper: prefetch + fasterq-dump per Run, then merge
# technical-replicate Runs into a single per-sample fastq.gz.
# =====================================================================
fetch_srr_and_merge_sra () {
    local sample_name="$1"
    shift
    local srrs=("$@")
    local out_fastq="$FASTQ_DIR/${sample_name}.fastq.gz"
    local tmp_dir="$FASTQ_DIR/.sra_tmp_${sample_name}"

    if [[ -f "$out_fastq" ]]; then
        log "  $sample_name already exists ($out_fastq), skipping."
        return
    fi

    log "  Fetching $sample_name from: ${srrs[*]}"
    mkdir -p "$tmp_dir"

    local tmp_files=()
    for srr in "${srrs[@]}"; do
        local se_fastq="$tmp_dir/${srr}.fastq"
        local pe_fastq="$tmp_dir/${srr}_1.fastq"

        # Already extracted from a previous (possibly killed) run of this
        # script - skip straight to using it, no need to re-prefetch/re-dump.
        if [[ -f "$se_fastq" ]]; then
            log "    $srr already extracted ($se_fastq), skipping download."
            tmp_files+=("$se_fastq")
            continue
        elif [[ -f "$pe_fastq" ]]; then
            log "    $srr already extracted ($pe_fastq), skipping download."
            log "    WARNING: $srr looks paired-end; using _1 only for $sample_name - verify this is correct."
            tmp_files+=("$pe_fastq")
            continue
        fi

        # prefetch itself resumes/skips if a complete .sra is already present
        # in $tmp_dir/$srr/, so this is safe to re-run after an interruption.
        log "    prefetch $srr"
        prefetch "$srr" -O "$tmp_dir"

        log "    fasterq-dump $srr"
        fasterq-dump --split-files -e 4 -O "$tmp_dir" "$tmp_dir/$srr/$srr.sra"

        if [[ -f "$tmp_dir/${srr}.fastq" ]]; then
            # single-end
            tmp_files+=("$tmp_dir/${srr}.fastq")
        elif [[ -f "$tmp_dir/${srr}_1.fastq" ]]; then
            # paired-end - all sources in this script are documented single-end;
            # flag it rather than silently dropping mate 2.
            log "    WARNING: $srr produced _1/_2 files (looks paired-end)."
            log "    WARNING: using _1 only for $sample_name - verify this is correct."
            tmp_files+=("$tmp_dir/${srr}_1.fastq")
        else
            log "    ERROR: no fastq output found for $srr in $tmp_dir"
            exit 1
        fi
    done

    log "    merging ${#tmp_files[@]} run(s) -> $out_fastq"
    cat "${tmp_files[@]}" | gzip > "$out_fastq"

    rm -rf "$tmp_dir"
    log "  -> $out_fastq"
}

# =====================================================================
# Source 1: Canella et al.: GSE18184 (PRJNA119545) - Pol III ChIP-seq, IMR90hTert
# =====================================================================
fetch_gse18184 () {
    log "=== Source 1: Canella et al.: GSE18184 (PRJNA119545) ==="

    # All ChIP-seq, single-end, IMR90hTert cell line. Grouped by GSM (technical
    # replicate runs merged per sample). Note SNAP45 runs are 76bp vs 35-36bp
    # for everything else in this series.
    declare -A CHIP_SAMPLES
    CHIP_SAMPLES[GSE18184_GSM454595_IMR90hTert_RPC4]="SRR035550 SRR035551"
    CHIP_SAMPLES[GSE18184_GSM454596_IMR90hTert_Bdp1]="SRR035552 SRR035553 SRR035554"
    CHIP_SAMPLES[GSE18184_GSM454597_IMR90hTert_Brf1]="SRR035555 SRR035556"
    CHIP_SAMPLES[GSE18184_GSM454598_IMR90hTert_Input]="SRR035557"
    CHIP_SAMPLES[GSE18184_GSM454599_IMR90hTert_SNAP45]="SRR035558 SRR035559"

    for sample in "${!CHIP_SAMPLES[@]}"; do
        fetch_srr_and_merge_sra "$sample" ${CHIP_SAMPLES[$sample]}
        echo "$sample" >> "$FASTQ_DIR/manifest_gse18184.txt"
    done

    log "GSE18184 fetch complete."
}

# =====================================================================
# Source 2: Oler et al.: GSE20309 (PRJNA125353) - Pol III ChIP-seq + RNA-seq
# =====================================================================
fetch_gse20309 () {
    log "=== Source 2: Oler et al.: GSE20309 (PRJNA125353) ==="

    # ChIP-seq samples (aligned with Bowtie2 in align_hg38.sh)
    declare -A CHIP_SAMPLES
    CHIP_SAMPLES[GSE20309_GSM509047_HEK293T_Pol3]="SRR036642 SRR036643"
    CHIP_SAMPLES[GSE20309_GSM509048_HFF_Pol3]="SRR036644 SRR036645"
    CHIP_SAMPLES[GSE20309_GSM509049_HeLa_Brf1]="SRR036646 SRR036647"
    CHIP_SAMPLES[GSE20309_GSM509050_HeLa_Brf2]="SRR036648 SRR036649"
    CHIP_SAMPLES[GSE20309_GSM509051_HeLa_Input]="SRR036650 SRR036651 SRR036652 SRR036653 SRR036654"
    CHIP_SAMPLES[GSE20309_GSM509052_HeLa_Pol3_Rep1]="SRR036655 SRR036656 SRR036657 SRR036658"
    CHIP_SAMPLES[GSE20309_GSM509053_HeLa_Pol3_Rep2]="SRR036659"
    CHIP_SAMPLES[GSE20309_GSM509054_HeLa_CappedRNAseq_Rep1]="SRR036660"
    CHIP_SAMPLES[GSE20309_GSM509056_HeLa_TFIIIC]="SRR036662 SRR036663"
    CHIP_SAMPLES[GSE20309_GSM509057_Jurkat_Pol3]="SRR036664 SRR036665 SRR036666"

    # RNA-seq samples (aligned with STAR in align_hg38.sh - splice-aware)
    declare -A RNA_SAMPLES
    RNA_SAMPLES[GSE20309_GSM509055_HeLa_TotalRNAseq]="SRR036661"

    for sample in "${!CHIP_SAMPLES[@]}"; do
        fetch_srr_and_merge_sra "$sample" ${CHIP_SAMPLES[$sample]}
        echo "$sample" >> "$FASTQ_DIR/manifest_chipseq.txt"
    done

    for sample in "${!RNA_SAMPLES[@]}"; do
        fetch_srr_and_merge_sra "$sample" ${RNA_SAMPLES[$sample]}
        echo "$sample" >> "$FASTQ_DIR/manifest_rnaseq.txt"
    done

    log "GSE20309 fetch complete."
}

# =====================================================================
# Source 3: Moqtaderi et al.: ENCODE Yale/Snyder ChIP-seq raw data (direct HTTP, no SRA)
# =====================================================================
fetch_encode_yale () {
    log "=== Source 3: Moqtaderi et al.: ENCODE Yale ChIP-seq (UCSC hgdownload) ==="

    # NOTE: these files are aligned to hg18 in UCSC's pre-made alignment tracks,
    # but the raw fastq itself is genome-agnostic - we align to hg38 ourselves
    # in align_hg38.sh. Some files carry a publication restriction date noted
    # in files.txt - check that before using data pre-dating its release.
    local base_url="http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeYaleChIPseq"
    local files_txt="$FASTQ_DIR/encode_yale_files.txt"

    if [[ ! -f "$files_txt" ]]; then
        log "Fetching files.txt manifest from UCSC..."
        wget -q -O "$files_txt" "$base_url/files.txt"
    fi
    log "$(wc -l < "$files_txt") entries listed in files.txt"

    # Pick out only the RawData fastq.gz entries for antibody=Pol3 (skip Alignments/Signal tracks, other antibodies)
    grep 'view=RawData' "$files_txt" | grep 'type=fastq' | grep -i 'antibody=Pol3' | cut -f1 \
        > "$FASTQ_DIR/encode_yale_rawdata_files.txt"
    log "$(wc -l < "$FASTQ_DIR/encode_yale_rawdata_files.txt") raw fastq files identified (antibody=Pol3 only)"

    while read -r fname; do
        sample_name="ENCODE_Yale_${fname%.fastq.gz}"
        fetch_direct_fastq "$sample_name" "$base_url/$fname"
        echo "$sample_name" >> "$FASTQ_DIR/manifest_encode_yale.txt"
    done < "$FASTQ_DIR/encode_yale_rawdata_files.txt"

    log "ENCODE Yale fetch complete."
    log "Reminder: check files.txt for any 'RESTRICTED until' dates before publishing results from these files."
}

# =====================================================================
# Main - run requested sources (default: all three)
# =====================================================================
requested=("$@")
[[ ${#requested[@]} -eq 0 ]] && requested=(gse18184 gse20309 encode_yale)

for source in "${requested[@]}"; do
    case "$source" in
        gse18184)    fetch_gse18184 ;;
        gse20309)    fetch_gse20309 ;;
        encode_yale) fetch_encode_yale ;;
        *)
            echo "Unknown source: $source (expected: gse18184, gse20309, encode_yale)" >&2
            exit 1
            ;;
    esac
done

log "All requested fetches complete. Run map_to_hg38.py next."