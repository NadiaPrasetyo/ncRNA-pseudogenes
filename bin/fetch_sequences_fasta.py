"""
fetch_sequences_fasta.py
This script fetches sequences in FASTA format for a list of genomes based on their locations

Overview:
This script reads a CSV file containing gene information, including gene names, chromosome locations, and gene
types. It then fetches the corresponding sequences from a specified genome FASTA file and saves them in FASTA format to an output file.

Args:
    -i, --input_file: Path to the input CSV file containing gene data.
    -o, --output_file: Path to the output FASTA file to save the fetched sequences.
    -g, --genome: Path to the genome FASTA file (default: data/downloads/GRCh38.p14.genome.fa).
    -f, --flank_length: Length of flanking regions to include (default: 100).
    --verify: If set, verify the fetched sequences against the Ensembl REST API after writing them.
    --species: Species name to use when verifying against Ensembl (default: human).

Usage:
    python fetch_sequences_fasta.py -i path/to/input.csv -o path/to/output.fasta -g path/to/genome.fa -f 100
    python fetch_sequences_fasta.py -i path/to/input.csv -o path/to/output.fasta --verify

Output:
    The output will be a FASTA file containing the fetched sequences for the specified genes, including
    flanking regions as specified by the flank_length argument.

Requires:
    pyfaidx (pip install pyfaidx) - for fast, indexed random access to the genome FASTA file.
    pandas
    requests

Author: Nadia Prasetyo

"""

import argparse
import os
import sys
import re
import time
import requests
import pandas as pd
from pyfaidx import Fasta

ENSEMBL_REST_URL = "https://rest.ensembl.org"


def load_genes_list(input_file):
    """
    Load the genes list from the input CSV file.
    expected columns: 'gene_name', 'gene_group', 'chromosome', 'start', 'end', 'gene_type'

    Parameters:
    - input_file: Path to the input CSV file.

    Returns:
    - List of dictionaries containing gene information.
    """
    try:
        df = pd.read_csv(input_file)
        if 'gene_name' not in df.columns or 'chromosome' not in df.columns or 'start' not in df.columns or 'end' not in df.columns or 'gene_type' not in df.columns:
            sys.exit("Error: The CSV file must contain 'gene_name', 'chromosome', 'start', 'end', and 'gene_type' columns.")
        return df[['gene_name', 'chromosome', 'start', 'end', 'gene_type']].to_dict(orient='records')
    except Exception as e:
        sys.exit(f"Error reading input file '{input_file}': {e}")


def fetch_sequence_from_genome(genome_fasta, chromosome, start, end):
    """
    Fetch a sequence from an indexed local genome FASTA file.

    Coordinates are treated as 1-based, inclusive (the same convention used by
    the input CSV and by the Ensembl REST API), so this can be compared
    directly against verify_sequences()/fetch_ensembl_sequence() output.

    Parameters:
    - genome_fasta: An open pyfaidx.Fasta object for the genome.
    - chromosome: Chromosome name as it appears in the genome FASTA (e.g. 'chr1').
    - start: 1-based start position (inclusive).
    - end: 1-based end position (inclusive).

    Returns:
    - The sequence as a string.
    """
    if chromosome not in genome_fasta:
        raise KeyError(f"Chromosome '{chromosome}' not found in genome FASTA.")
    # pyfaidx's get_seq takes 1-based inclusive coordinates, matching Ensembl's convention.
    return str(genome_fasta.get_seq(chromosome, start, end))


def fetch_sequences(genes_list, genome_file, output_file, flank_length):
    """
    Fetch sequences in FASTA format for the given list of genes.

    Parameters:
    - genes_list: List of dictionaries containing gene information.
    - genome_file: Path to the genome FASTA file.
    - output_file: Path to the output FASTA file.
    - flank_length: Length of flanking regions to include.
    """
    try:
        # Open (and, if needed, build a .fai index for) the genome once, rather
        # than re-opening it for every gene.
        genome_fasta = Fasta(genome_file)

        with open(output_file, 'w') as fasta_out:
            for gene in genes_list:
                gene_name = gene['gene_name']
                chromosome = gene['chromosome']
                start = int(gene['start'])
                end = int(gene['end'])
                gene_type = gene['gene_type']

                # Adjust start and end positions based on flank_length.
                # Coordinates are 1-based, so the smallest valid start is 1 (not 0).
                adjusted_start = max(1, start - flank_length)
                adjusted_end = end + flank_length

                try:
                    sequence = fetch_sequence_from_genome(genome_fasta, chromosome, adjusted_start, adjusted_end)
                except KeyError as e:
                    print(f"Warning: skipping '{gene_name}' - {e}")
                    continue

                # Write to FASTA format
                fasta_out.write(f">{gene_name}|{gene_type}|{chromosome}:{adjusted_start}-{adjusted_end}\n")
                fasta_out.write(f"{sequence}\n")

        print(f"Sequences fetched and saved to '{output_file}' successfully.")

    except Exception as e:
        sys.exit(f"Error fetching sequences: {e}")


def normalize_chromosome(chromosome):
    """
    Convert a chromosome name from the CSV/genome-FASTA style (e.g. 'chr1', 'chrX',
    'chrM') into the style expected by the Ensembl REST API (e.g. '1', 'X', 'MT').

    Parameters:
    - chromosome: Chromosome name, with or without a 'chr' prefix.

    Returns:
    - Normalized chromosome name suitable for the Ensembl region string.
    """
    chrom = re.sub(r'^chr', '', chromosome, flags=re.IGNORECASE)
    if chrom.upper() in ('M', 'MT'):
        chrom = 'MT'
    return chrom


def fetch_ensembl_sequence(chromosome, start, end, species="human", strand=1):
    """
    Fetch a genomic sequence from the Ensembl REST API for comparison/verification.

    Parameters:
    - chromosome: Chromosome name (will be normalized, e.g. 'chr1' -> '1').
    - start: 1-based start position (inclusive).
    - end: 1-based end position (inclusive).
    - species: Species name/alias understood by Ensembl (default: 'human').
    - strand: 1 for forward strand, -1 for reverse strand.

    Returns:
    - The sequence as a string, as returned by Ensembl.
    """
    chrom = normalize_chromosome(chromosome)
    region = f"{chrom}:{start}..{end}:{strand}"
    url = f"{ENSEMBL_REST_URL}/sequence/region/{species}/{region}"
    headers = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    return response.json().get("seq", "")


def parse_fasta_records(fasta_file):
    """
    Yield (header, sequence) tuples from a FASTA file.
    """
    header = None
    seq_chunks = []
    with open(fasta_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def verify_sequences(output_file, species="human", request_delay=0.1):
    """
    Verify the sequences in a fetched FASTA file against the Ensembl REST API.

    Expects headers written by fetch_sequences(), in the form:
        gene_name | gene_type | chromosome:start-end

    Parameters:
    - output_file: Path to the FASTA file produced by fetch_sequences().
    - species: Species name/alias to query on Ensembl (default: 'human').
    - request_delay: Seconds to sleep between Ensembl requests, to stay within
      Ensembl's rate limits when verifying many sequences.

    Returns:
    - True if every sequence matched Ensembl, False if any mismatch or error occurred.
    """
    all_match = True
    location_pattern = re.compile(r'([\w.]+):(\d+)-(\d+)\s*$')

    for header, local_seq in parse_fasta_records(output_file):
        match = location_pattern.search(header)
        if not match:
            print(f"[SKIP] Could not parse location from header: {header}")
            continue

        chromosome, start, end = match.group(1), int(match.group(2)), int(match.group(3))

        try:
            ensembl_seq = fetch_ensembl_sequence(chromosome, start, end, species=species)
        except Exception as e:
            print(f"[ERROR] {header}: failed to fetch from Ensembl - {e}")
            all_match = False
            time.sleep(request_delay)
            continue

        if ensembl_seq.upper() == local_seq.upper():
            print(f"[MATCH] {header}")
        else:
            all_match = False
            print(f"[MISMATCH] {header}")
            print(f"          local length={len(local_seq)}, ensembl length={len(ensembl_seq)}")

        time.sleep(request_delay)  # be polite to the Ensembl REST API

    return all_match


def main():
    parser = argparse.ArgumentParser(description="Fetch sequences in FASTA format for a list of genomes based on their locations.")
    parser.add_argument("-i", "--input_file", required=True, help="Path to the input CSV file containing gene data.")
    parser.add_argument("-o", "--output_file", default="data/gene_sequences.fasta", help="Path to the output FASTA file to save the fetched sequences. default: data/gene_sequences.fasta.")
    parser.add_argument("-g", "--genome", default="data/downloads/GRCh38.p14.genome.fa", help="Path to the genome FASTA file (default: data/downloads/GRCh38.p14.genome.fa).")
    parser.add_argument("-f", "--flank_length", type=int, default=100, help="Length of flanking regions to include (default: 100).")
    parser.add_argument("--verify", action="store_true", help="Verify fetched sequences against the Ensembl REST API after writing them.")
    parser.add_argument("--species", default="human", help="Species name/alias to use when verifying against Ensembl (default: human).")
    args = parser.parse_args()

    # Ensure the input file exists
    if not os.path.isfile(args.input_file):
        sys.exit(f"Error: Input file '{args.input_file}' does not exist.")

    # Load the genes list from the input CSV file
    genes_list = load_genes_list(args.input_file)

    # Fetch sequences and save to output FASTA file
    fetch_sequences(genes_list, args.genome, args.output_file, args.flank_length)

    if args.verify:
        print("\nVerifying fetched sequences against Ensembl...")
        ok = verify_sequences(args.output_file, species=args.species)
        print("\nAll sequences matched Ensembl." if ok else "\nSome sequences did NOT match Ensembl - see log above.")


if __name__ == "__main__":
    main()