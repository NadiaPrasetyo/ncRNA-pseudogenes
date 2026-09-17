import requests
import re
import time
import csv
import pyBigWig

mane_gtf_path = 'data/downloads/MANE.GRCh38.v1.5.ensembl_genomic.gtf'  # Already unzipped. Fetch from https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/
bigbed_file_path = 'data/downloads/hgnc.bb'  # Final-resort backup. Fetch from https://hgdownload.soe.ucsc.edu/gbdb/hg38/hgnc/


# Regex helpers for pulling values out of the GTF attribute column, e.g.
# gene_id "ENSG00000141510"; transcript_id "ENST00000646891.2"; gene_name "BRAF"; tag "MANE_Select";
GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')
TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')


def _extract_attr(pattern, attributes):
    match = pattern.search(attributes)
    return match.group(1) if match else None


# Pick the best available NCBI RefSeq accession from an HGNC record: prefer the NM_/NR_
# entry from mane_select (it's the accession paired with the MANE Select transcript), and
# fall back to the first refseq_accession entry otherwise. Returns None if neither exists.
def get_best_ncbi_accession(hgnc_doc):
    mane_select = hgnc_doc.get('mane_select') or []
    mane_ncbi = next((v for v in mane_select if v.startswith(('NM_', 'NR_'))), None)
    if mane_ncbi:
        return mane_ncbi
    refseq = hgnc_doc.get('refseq_accession') or []
    return refseq[0] if refseq else None


# Parse the MANE GTF once and build two lookups keyed off the "transcript" feature rows:
#   by_transcript_id: versioned ENST id (e.g. "ENST00000646891.2") -> (chrom, start, end, gene_id, gene_name)
#   by_gene_name:     gene symbol -> (chrom, start, end, transcript_id, gene_id)   [fallback lookup]
# MANE GTFs contain one representative transcript per gene (MANE Select, or MANE Plus Clinical
# for a handful of genes), so this stays a clean mapping without disambiguating multiple rows.
# Note: GTF coordinates are 1-based, closed-interval (unlike BigBed's 0-based, half-open).
# Note: MANE only covers protein-coding genes -- it will not contain ncRNA genes/pseudogenes
# (e.g. RNU*, TRNA, RNY, VTRNA groups below, or RNU1-5P/RNU1-6P-style pseudogenes). Those fall
# through to the Ensembl REST API backup in fetch_ensembl_location() further down.
def load_mane_gene_locations(gtf_path):
    by_transcript_id = {}
    by_gene_name = {}

    with open(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            if feature != 'transcript':
                continue

            gene_name = _extract_attr(GENE_NAME_RE, attributes)
            gene_id = _extract_attr(GENE_ID_RE, attributes)
            transcript_id = _extract_attr(TRANSCRIPT_ID_RE, attributes)
            if not gene_name or not transcript_id:
                continue

            entry = (chrom, int(start), int(end), gene_id, gene_name)
            by_transcript_id[transcript_id] = entry
            by_gene_name[gene_name] = (chrom, int(start), int(end), transcript_id, gene_id)

    return by_transcript_id, by_gene_name


# Resolve a gene's location using, in priority order:
#   1. The MANE Select/Plus Clinical transcript ID from HGNC's mane_select field, matched
#      against the GTF's transcript rows (exact match, then version-stripped match).
#   2. A fallback lookup by gene symbol against the GTF's gene_name index.
# Also returns an Ensembl gene ID and NCBI accession, cross-checked against HGNC where possible.
def resolve_gene_location(hgnc_doc, by_transcript_id, by_gene_name):
    symbol = hgnc_doc.get('symbol')
    hgnc_ensembl_id = hgnc_doc.get('ensembl_gene_id')
    mane_select = hgnc_doc.get('mane_select') or []

    mane_transcript_id = next((v for v in mane_select if v.startswith('ENST')), None)
    ncbi_accession = get_best_ncbi_accession(hgnc_doc)

    location = None
    transcript_id_used = None
    gtf_gene_id = None

    if mane_transcript_id:
        location = by_transcript_id.get(mane_transcript_id)
        if not location:
            # Try matching without the version suffix in case of a version mismatch
            base_id = mane_transcript_id.split('.')[0]
            for tid, entry in by_transcript_id.items():
                if tid.split('.')[0] == base_id:
                    location = entry
                    break
        if location:
            transcript_id_used = mane_transcript_id

    if not location and symbol in by_gene_name:
        chrom, start, end, transcript_id, gtf_gene_id = by_gene_name[symbol]
        location = (chrom, start, end, gtf_gene_id, symbol)
        transcript_id_used = transcript_id

    if not location:
        return None

    chrom, start, end, gtf_gene_id, _ = location

    # Prefer the Ensembl gene ID pulled from the GTF entry itself, falling back to HGNC's.
    ensembl_id = gtf_gene_id or hgnc_ensembl_id or "not found"
    # Compare version-stripped, since the GTF's gene_id carries a ".N" version suffix
    # (e.g. "ENSG00000207389.1") that HGNC's ensembl_gene_id does not include.
    if gtf_gene_id and hgnc_ensembl_id and gtf_gene_id.split('.')[0] != hgnc_ensembl_id.split('.')[0]:
        print(f"  Note: Ensembl ID mismatch for {symbol} -- GTF: {gtf_gene_id}, HGNC: {hgnc_ensembl_id}")

    ncbi_accession = ncbi_accession or "not found"
    transcript_id_used = transcript_id_used or "not found"

    return chrom, start, end, transcript_id_used, ensembl_id, ncbi_accession


# Backup lookup for genes not covered by MANE (e.g. pseudogenes, ncRNAs) using Ensembl's
# REST API directly against the ensembl_gene_id HGNC gave us. Returns chrom/start/end plus
# whatever canonical transcript Ensembl reports, or None if the ID can't be resolved.
# Retries on timeouts/connection errors and honors 429 Retry-After, since a single slow
# response shouldn't be enough to fall all the way through to the BigBed backup.
def fetch_ensembl_location(ensembl_gene_id, max_retries=3, timeout=20):
    if not ensembl_gene_id:
        return None

    url = f'https://rest.ensembl.org/lookup/id/{ensembl_gene_id}?content-type=application/json'

    for attempt in range(1, max_retries + 1):
        try:
            response = requests.get(url, headers={'Content-Type': 'application/json'}, timeout=timeout)
        except requests.exceptions.RequestException as e:
            print(f"  Ensembl lookup attempt {attempt}/{max_retries} failed for {ensembl_gene_id}: {e}")
            if attempt < max_retries:
                time.sleep(attempt * 1.5)  # simple linear backoff between retries
                continue
            return None

        if response.status_code == 429:
            retry_after = float(response.headers.get('Retry-After', 1))
            print(f"  Ensembl rate limit hit for {ensembl_gene_id}; waiting {retry_after}s")
            time.sleep(retry_after)
            continue

        try:
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"  Ensembl lookup attempt {attempt}/{max_retries} failed for {ensembl_gene_id}: {e}")
            if attempt < max_retries:
                time.sleep(attempt * 1.5)
                continue
            return None

        try:
            data = response.json()
        except ValueError as e:
            print(f"  Error parsing Ensembl response for {ensembl_gene_id}: {e}")
            return None

        seq_region = data.get('seq_region_name')
        start = data.get('start')
        end = data.get('end')
        if seq_region is None or start is None or end is None:
            return None

        # Ensembl reports bare contig names (e.g. "1", "X"); normalize to match the GTF's "chr1" style.
        chrom = str(seq_region) if str(seq_region).startswith('chr') else f'chr{seq_region}'
        transcript_id = data.get('canonical_transcript', 'not found')

        return chrom, int(start), int(end), transcript_id

    return None


# Shared core for NCBI Datasets API calls: fetches a dataset_report URL with retry/backoff
# and 429 handling, then parses out chrom/start/end from the gene's current-assembly
# annotation and a representative transcript ID from the paired product report.
def _fetch_ncbi_datasets_report(url, label, max_retries=3, timeout=20):
    headers = {'Accept': 'application/json'}

    for attempt in range(1, max_retries + 1):
        try:
            response = requests.get(url, headers=headers, timeout=timeout)
        except requests.exceptions.RequestException as e:
            print(f"  NCBI Datasets lookup attempt {attempt}/{max_retries} failed for {label}: {e}")
            if attempt < max_retries:
                time.sleep(attempt * 1.5)
                continue
            return None

        if response.status_code == 429:
            retry_after = float(response.headers.get('Retry-After', 1))
            print(f"  NCBI Datasets rate limit hit for {label}; waiting {retry_after}s")
            time.sleep(retry_after)
            continue

        try:
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"  NCBI Datasets lookup attempt {attempt}/{max_retries} failed for {label}: {e}")
            if attempt < max_retries:
                time.sleep(attempt * 1.5)
                continue
            return None

        try:
            data = response.json()
        except ValueError as e:
            print(f"  Error parsing NCBI Datasets response for {label}: {e}")
            return None

        reports = data.get('reports') or []
        if not reports:
            return None

        gene = reports[0].get('gene') or {}
        product = reports[0].get('product') or {}

        chrom, start, end = None, None, None
        for annotation in gene.get('annotations') or []:
            for loc in annotation.get('genomic_locations') or []:
                genomic_range = loc.get('genomic_range') or {}
                begin = genomic_range.get('begin')
                g_end = genomic_range.get('end')
                seq_name = loc.get('sequence_name')
                if begin and g_end and seq_name:
                    chrom, start, end = seq_name, begin, g_end
                    break
            if chrom:
                break

        if chrom is None or start is None or end is None:
            return None

        chrom = str(chrom) if str(chrom).startswith('chr') else f'chr{chrom}'

        transcript_id = "not found"
        transcripts = product.get('transcripts') or []
        if transcripts:
            transcript_id = transcripts[0].get('ensembl_transcript') or transcripts[0].get('accession_version') or "not found"

        try:
            return chrom, int(start), int(end), transcript_id
        except (TypeError, ValueError):
            return None

    return None


# Third-tier backup, attempt A: look the gene up on NCBI's Datasets API by NCBI Gene ID
# (HGNC's entrez_id field). This is the more reliable of the two Datasets lookups since it
# doesn't depend on guessing an accession type -- any gene with an NCBI Gene ID resolves.
def fetch_ncbi_datasets_location_by_gene_id(entrez_id):
    if not entrez_id:
        return None
    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{entrez_id}/dataset_report'
    return _fetch_ncbi_datasets_report(url, f"gene ID {entrez_id}")


# Third-tier backup, attempt B: look the gene up on NCBI's Datasets API by transcript/protein
# accession. The /gene/accession/ endpoint only accepts transcript- or protein-level
# accessions (NM_/NR_/XM_/XR_/NP_/XP_) -- NOT genomic accessions like NG_ (RefSeqGene) or
# NC_ (chromosome), which come back with an empty report. Callers should only pass accessions
# of the supported types; see the prefix filter where this is called.
def fetch_ncbi_datasets_location_by_accession(ncbi_accession):
    if not ncbi_accession:
        return None
    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2/gene/accession/{ncbi_accession}/dataset_report'
    return _fetch_ncbi_datasets_report(url, ncbi_accession)


# Last-resort backup for genes that neither MANE, Ensembl, nor NCBI Datasets could resolve: search the
# original hgnc.bb BigBed file. This was the original data source; it's kept only as a final
# fallback since it's what was producing wrong coordinates for some genes in the first place --
# so anything it resolves is worth spot-checking rather than trusting outright.
# Note: BigBed coordinates are 0-based, half-open (unlike the GTF's 1-based, closed-interval),
# so the start is bumped by 1 here to keep all three sources consistent in the output CSV.
def fetch_bigbed_location(gene_symbol, bigbed_file):
    try:
        chromosomes = bigbed_file.chroms()
        for chrom in chromosomes:
            entries = bigbed_file.entries(chrom, 0, chromosomes[chrom])
            if not entries:
                continue
            for start, end, name in entries:
                if gene_symbol not in name:
                    continue
                parts = name.split()
                found_symbol = parts[0]
                if found_symbol != gene_symbol:
                    continue
                transcript_id = parts[1] if len(parts) > 1 else "not found"
                return chrom, start + 1, end, transcript_id
        return None
    except Exception as e:
        print(f"  Error fetching data from BigBed file for {gene_symbol}: {e}")
        return None


# Read gene symbols already present in the output CSV so re-runs don't redo work or
# duplicate rows. Returns an empty set if the file doesn't exist yet.
def load_existing_genes(output_file):
    existing = set()
    try:
        with open(output_file, 'r', newline='') as f:
            reader = csv.reader(f)
            next(reader, None)  # skip header
            for row in reader:
                if row:
                    existing.add(row[0])
    except FileNotFoundError:
        pass
    return existing


# Main function to get gene locations and print to the console
def get_gene_locations(gene_symbols, output_file, by_transcript_id, by_gene_name, bigbed_file, query='unknown'):

    if not isinstance(query, str) or not query.strip():
        query = 'unknown'

    existing_genes = load_existing_genes(output_file)

    # Open the output file to append results to a single CSV (create if missing)
    # Seek to start to check for existing header; write header if missing
    with open(output_file, 'a+') as file:
        file.seek(0)
        first_line = file.readline()
        if not first_line.startswith("GeneName"):
            file.write("GeneName,GeneGroup,Transcript_ID,EnsemblID,NCBI_Accession,Chromosome,Start,End\n")
        else:
            file.write("\n")  # Ensure we're at the end of the file for appending

        # Process each gene symbol
        for gene_symbol in gene_symbols:
            # Simple rule: it was returned for this query, so it's classified as this group.
            gene_group = query

            if gene_symbol in existing_genes:
                print(f"Skipping {gene_symbol} (already in output file)")
                continue

            print(f"Processing {gene_symbol} [{gene_group}]")

            # Fetch the full HGNC record for this symbol (search only returns minimal fields)
            hgnc_doc = fetch_hgnc_gene_details(gene_symbol)
            if not hgnc_doc:
                print(f"No HGNC record found for {gene_symbol}")
                print("------")
                continue

            resolved = resolve_gene_location(hgnc_doc, by_transcript_id, by_gene_name)

            # MANE only covers protein-coding genes, so pseudogenes/ncRNAs won't be in it.
            # Fall back to querying Ensembl directly by the ensembl_gene_id HGNC gave us.
            if not resolved:
                hgnc_ensembl_id = hgnc_doc.get('ensembl_gene_id')
                if hgnc_ensembl_id:
                    print(f"  No MANE match for {gene_symbol}; trying Ensembl lookup ({hgnc_ensembl_id})")
                    ensembl_location = fetch_ensembl_location(hgnc_ensembl_id)
                    time.sleep(1/10)  # brief pause to stay well under Ensembl's rate limit
                    if ensembl_location:
                        chrom, start, end, transcript_id = ensembl_location
                        ncbi_accession = get_best_ncbi_accession(hgnc_doc) or "not found"
                        resolved = (chrom, start, end, transcript_id, hgnc_ensembl_id, ncbi_accession)

            # Third resort: look the gene up on NCBI's Datasets API. Try by NCBI Gene ID first
            # (most reliable), then by accession -- but only accession types the /accession/
            # endpoint actually supports (NM_/NR_/XM_/XR_/NP_/XP_); genomic accessions like
            # NG_ (RefSeqGene) return an empty report, so skip those rather than wasting a call.
            if not resolved:
                entrez_id = hgnc_doc.get('entrez_id')
                if entrez_id:
                    print(f"  No Ensembl match for {gene_symbol}; trying NCBI Datasets lookup (gene ID {entrez_id})")
                    ncbi_location = fetch_ncbi_datasets_location_by_gene_id(entrez_id)
                    time.sleep(1/10)  # brief pause to stay well under NCBI's rate limit
                    if ncbi_location:
                        chrom, start, end, transcript_id = ncbi_location
                        hgnc_ensembl_id = hgnc_doc.get('ensembl_gene_id') or "not found"
                        ncbi_accession = get_best_ncbi_accession(hgnc_doc) or "not found"
                        resolved = (chrom, start, end, transcript_id, hgnc_ensembl_id, ncbi_accession)

                if not resolved:
                    ncbi_accession = get_best_ncbi_accession(hgnc_doc)
                    if ncbi_accession and ncbi_accession.startswith(('NM_', 'NR_', 'XM_', 'XR_', 'NP_', 'XP_')):
                        print(f"  No NCBI gene-ID match for {gene_symbol}; trying NCBI Datasets lookup by accession ({ncbi_accession})")
                        ncbi_location = fetch_ncbi_datasets_location_by_accession(ncbi_accession)
                        time.sleep(1/10)
                        if ncbi_location:
                            chrom, start, end, transcript_id = ncbi_location
                            hgnc_ensembl_id = hgnc_doc.get('ensembl_gene_id') or "not found"
                            resolved = (chrom, start, end, transcript_id, hgnc_ensembl_id, ncbi_accession)

            # Final resort: the original hgnc.bb BigBed file.
            if not resolved and bigbed_file is not None:
                print(f"  No NCBI Datasets match for {gene_symbol}; trying BigBed as last resort")
                bigbed_location = fetch_bigbed_location(gene_symbol, bigbed_file)
                if bigbed_location:
                    chrom, start, end, transcript_id = bigbed_location
                    hgnc_ensembl_id = hgnc_doc.get('ensembl_gene_id') or "not found"
                    ncbi_accession = get_best_ncbi_accession(hgnc_doc) or "not found"
                    resolved = (chrom, start, end, transcript_id, hgnc_ensembl_id, ncbi_accession)
                    print(f"  Note: {gene_symbol} location came from BigBed -- worth double-checking")

            if resolved:
                chrom, start, end, transcript_id, ensembl_id, ncbi_accession = resolved
                formatted_start = '{:,}'.format(start)
                formatted_end = '{:,}'.format(end)
                location_str = f"Genomic Sequence ({chrom}:{formatted_start}-{formatted_end})"

                print(f"{gene_symbol} [{gene_group}]: {location_str} (Transcript: {transcript_id}, Ensembl: {ensembl_id}, NCBI: {ncbi_accession})")
                file.write(f"{gene_symbol},{gene_group},{transcript_id},{ensembl_id},{ncbi_accession},{chrom},{start},{end}\n")
                existing_genes.add(gene_symbol)
            else:
                print(f"No location found for {gene_symbol}")

            print("------")

            # Optional: brief sleep to avoid rate limiting on the HGNC API calls above
            time.sleep(1/10)



# Fetch the full HGNC record for a single approved symbol (includes ensembl_gene_id,
# refseq_accession, mane_select, etc. -- fields the /search/ endpoint doesn't return)
def fetch_hgnc_gene_details(symbol):
    url = f'https://rest.genenames.org/fetch/symbol/{symbol}'
    headers = {'Accept': 'application/json'}
    try:
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Request error fetching {symbol}: {e}")
        return None

    try:
        data = response.json()
        docs = data['response']['docs']
        return docs[0] if docs else None
    except (ValueError, KeyError, IndexError) as e:
        print(f"Error parsing fetch response for {symbol}: {e}")
        return None


# Call function with the gene group and output file path
if __name__ == '__main__':
    input_txt = "results/outliers.txt"
    # parse the list of gene_symbols into a deduplicated array
    gene_symbols = []
    try:
        with open(input_txt, 'r') as f:
            gene_symbols = sorted({line.strip() for line in f if line.strip()})
    except FileNotFoundError:
        print(f"Input file not found: {input_txt}")

    query = 'outliers'

    # Parse the MANE GTF once up front instead of re-scanning it per gene
    print("Loading MANE GTF...")
    by_transcript_id, by_gene_name = load_mane_gene_locations(mane_gtf_path)
    print(f"Loaded {len(by_transcript_id)} transcripts from MANE GTF")

    # Open the BigBed file once for the final-resort backup
    try:
        bigbed_file = pyBigWig.open(bigbed_file_path)
    except Exception as e:
        print(f"Could not open BigBed backup file at {bigbed_file_path}: {e}")
        bigbed_file = None

    get_gene_locations(gene_symbols, 'results/outliers_hgnc_loc.csv', by_transcript_id, by_gene_name, bigbed_file, query=query)

    if bigbed_file is not None:
        bigbed_file.close()