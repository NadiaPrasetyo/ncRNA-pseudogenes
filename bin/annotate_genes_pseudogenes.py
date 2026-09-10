#!/usr/bin/env python3
"""
annotate_genes.py

Reads a CSV of genes (with a 'gene_name' column) and annotates each row with
a gene_type classification: "pseudogene", "functional", or
"unknown", based on HGNC (https://www.genenames.org) locus_type data.

Classification logic (in order of precedence):
  1. Exact-symbol lookup against HGNC's /fetch/symbol endpoint. If found,
     the reported locus_type is used directly:
        - locus_type contains "pseudogene"      -> "pseudogene"
        - locus_type contains "RNA, transfer"    -> "functional_rna"
        - locus_type is any other approved type  -> "gene"
  2. If the symbol isn't found in HGNC directly (e.g. it's an older/alias
     symbol, or a member of a numbered RNA gene family such as RNU6-1,
     TRNA-Ala, RNY1 etc.), fall back to:
        a. Prefix match against a curated set of functional non-coding RNA
           gene-family prefixes (gene_groups) -> "functional"
        b. Wildcard search via search_hgnc_pseudogenes() -> "pseudogene"
        c. Wildcard search via search_hgnc_functional_genes() -> "functional"
  3. Otherwise -> "unknown"

Usage:
    python annotate_genes.py input.csv output.csv
    python annotate_genes.py input.csv output.csv --workers 4 --sleep 0.1
"""

import argparse
import csv
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests

HGNC_BASE = "https://rest.genenames.org"
HEADERS = {"Accept": "application/json"}
TIMEOUT = 10

# Curated prefixes for functional non-coding RNA gene families whose members
# are usually numbered (e.g. RNU6-1, TRNA-Ala-AGC-1-1, RNY1, VTRNA1-1, ...).
# HGNC's exact-symbol fetch often misses these because the CSV's gene_name
# may not match the canonical HGNC symbol formatting exactly.
GENE_GROUPS = {
    "RN7SK", "RN7SL", "RNU11", "RNU12", "RNU1", "RNU2", "RNU4ATAC",
    "RNU4", "RNU5", "RNU6ATAC", "RNU6", "RNU7", "RNY", "TRNA", "VTRNA",
}

# Simple in-memory cache so repeated gene names (or repeated runs within the
# same process) don't trigger duplicate HTTP calls.
_CACHE = {}


# --------------------------------------------------------------------------
# HGNC query helpers
# --------------------------------------------------------------------------

def hgnc_fetch_symbol(symbol):
    """
    Exact-match lookup of a gene symbol against HGNC's /fetch/symbol endpoint.
    Returns the locus_type string if found, otherwise None.
    """
    url = f"{HGNC_BASE}/fetch/symbol/{symbol}"
    try:
        response = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"[warn] HGNC fetch error for '{symbol}': {e}", file=sys.stderr)
        return None

    try:
        data = response.json()
        docs = data.get("response", {}).get("docs", [])
        if not docs:
            return None
        return docs[0].get("locus_type")
    except (ValueError, KeyError) as e:
        print(f"[warn] Error parsing HGNC fetch response for '{symbol}': {e}", file=sys.stderr)
        return None


def search_hgnc_pseudogenes(query):
    """Search HGNC for approved pseudogenes matching `query` (wildcard search)."""
    if not query or not isinstance(query, str):
        print("Invalid query: must be a non-empty string.", file=sys.stderr)
        return []
    url = (
        f"{HGNC_BASE}/search/*{query}*+AND+status:%22Approved%22"
        f"+AND+locus_type:%22pseudogene%22"
    )
    try:
        response = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Request error: {e}", file=sys.stderr)
        return []

    try:
        data = response.json()
        if data["response"]["numFound"] > 0:
            return [gene["symbol"] for gene in data["response"]["docs"]]
        return []
    except (ValueError, KeyError) as e:
        print(f"Error parsing response: {e}", file=sys.stderr)
        return []


def search_hgnc_functional_genes(query):
    """Search HGNC for approved functional (transfer RNA) genes matching `query`."""
    if not query or not isinstance(query, str):
        print("Invalid query: must be a non-empty string.", file=sys.stderr)
        return []
    url = (
        f"{HGNC_BASE}/search/*{query}*+AND+status:%22Approved%22"
        f"+AND+locus_type:%22RNA%2C%20transfer%22"
    )
    try:
        response = requests.get(url, headers=HEADERS, timeout=TIMEOUT)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Request error: {e}", file=sys.stderr)
        return []

    try:
        data = response.json()
        if data["response"]["numFound"] > 0:
            return [gene["symbol"] for gene in data["response"]["docs"]]
        return []
    except (ValueError, KeyError) as e:
        print(f"Error parsing response: {e}", file=sys.stderr)
        return []


# --------------------------------------------------------------------------
# Classification
# --------------------------------------------------------------------------

def _prefix_matches_gene_group(symbol):
    """Return True if symbol starts with one of the curated RNA family prefixes."""
    upper = symbol.upper()
    for prefix in GENE_GROUPS:
        # Match e.g. RNU6-1 -> RNU6, TRNA-Ala-AGC-1-1 -> TRNA, RNY1 -> RNY
        if upper.startswith(prefix):
            remainder = upper[len(prefix):]
            if remainder == "" or remainder[0] in "-0123456789":
                return True
    return False


def classify_gene(symbol, sleep=0.0):
    """
    Classify a single gene symbol as 'pseudogene', 'functional',
    or 'unknown' and also return the HGNC locus_type string, if recovered.
    Uses a module-level cache to avoid duplicate lookups.
    """
    if not symbol or not isinstance(symbol, str):
        return "unknown", "unknown"

    symbol = symbol.strip()
    if symbol in _CACHE:
        return _CACHE[symbol]

    locus_type = "unknown"
    result = "unknown"

    # 1. Exact HGNC lookup
    locus_type = hgnc_fetch_symbol(symbol)
    if sleep:
        time.sleep(sleep)

    if locus_type:
        lt = locus_type.lower()
        if "pseudogene" in lt:
            result = "pseudogene"
        else:
            result = "functional"
    else:
        # 2a. Curated RNA gene family prefix match
        if _prefix_matches_gene_group(symbol):
            result = "functional"
        else:
            # 2b. Wildcard pseudogene search
            pg_hits = search_hgnc_pseudogenes(symbol)
            if sleep:
                time.sleep(sleep)
            if symbol.upper() in {s.upper() for s in pg_hits}:
                result = "pseudogene"
            else:
                # 2c. Wildcard functional (tRNA) search
                fn_hits = search_hgnc_functional_genes(symbol)
                if sleep:
                    time.sleep(sleep)
                if symbol.upper() in {s.upper() for s in fn_hits}:
                    result = "functional"
                elif pg_hits or fn_hits:
                    # Query matched something but not an exact symbol; leave unknown
                    result = "unknown"

    _CACHE[symbol] = (result, locus_type)
    return result, locus_type


# --------------------------------------------------------------------------
# CSV processing
# --------------------------------------------------------------------------

def annotate_csv(input_path, output_path, sleep=0.0, workers=1):
    with open(input_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames

    # Guard against ragged rows: when a data row has MORE columns than the
    # header, DictReader stuffs the extras into a list under the key `None`
    # (its `restkey`). That's fine for reading, but DictWriter later chokes
    # on the `None` key since it's not in fieldnames. Warn and drop it.
    for i, row in enumerate(rows, start=2):  # start=2: header is line 1
        extra = row.pop(None, None)
        if extra:
            print(
                f"[warn] row {i} of {input_path} has {len(extra)} more column(s) "
                f"than the header ({', '.join(map(str, extra))}); ignoring the "
                f"extra value(s). Check this row for a stray comma or unquoted "
                f"field containing a comma.",
                file=sys.stderr,
            )
        # Rows with FEWER columns than the header get `None` filled in for
        # the missing ones (DictReader's `restval`). Normalize those to ""
        # so the later .strip() calls on gene_name/GeneName don't crash.
        for key, val in row.items():
            if val is None:
                row[key] = ""

    if "gene_name" not in fieldnames and "GeneName" not in fieldnames:
        raise ValueError("Input CSV must contain a 'gene_name' or 'GeneName' column.")

    unique_names = sorted({row.get("gene_name", row.get("GeneName", "")).strip() for row in rows})
    print(f"Found {len(unique_names)} unique gene names to classify...")

    classifications = {}
    locus_types = {}
    if workers > 1:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_name = {
                executor.submit(classify_gene, name, sleep): name for name in unique_names
            }
            for i, future in enumerate(as_completed(future_to_name), 1):
                name = future_to_name[future]
                try:
                    classifications[name], locus_types[name] = future.result()
                except Exception as e:
                    print(f"[warn] classification failed for '{name}': {e}", file=sys.stderr)
                    classifications[name] = "unknown"
                    locus_types[name] = "unknown"
                if i % 25 == 0 or i == len(unique_names):
                    print(f"  classified {i}/{len(unique_names)}")
    else:
        for i, name in enumerate(unique_names, 1):
            classifications[name], locus_types[name] = classify_gene(name, sleep=sleep)
            if i % 25 == 0 or i == len(unique_names):
                print(f"  classified {i}/{len(unique_names)}")

    out_fieldnames = list(fieldnames) + ["locus_type"] + ["gene_type"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames)
        writer.writeheader()
        for row in rows:
            name = row.get("gene_name", row.get("GeneName", "")).strip()
            row["locus_type"] = locus_types.get(name, "") or "unknown"
            row["gene_type"] = classifications.get(name, "unknown")
            writer.writerow(row)

    counts = {}
    for v in classifications.values():
        counts[v] = counts.get(v, 0) + 1
    print("Classification summary:")
    for k, v in sorted(counts.items()):
        print(f"  {k}: {v}")
    print(f"Annotated CSV written to: {output_path}")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Annotate genes in a CSV as gene/pseudogene/functional_rna using HGNC."
    )
    parser.add_argument("input_csv", help="Path to input CSV (must have a 'gene_name' column).")
    parser.add_argument("output_csv", help="Path to write the annotated output CSV.")
    parser.add_argument(
        "--sleep", type=float, default=0.1,
        help="Seconds to sleep between HGNC requests per worker (default: 0.1, be polite to the API)."
    )
    parser.add_argument(
        "--workers", type=int, default=1,
        help="Number of concurrent threads for HGNC lookups (default: 1)."
    )
    args = parser.parse_args()

    annotate_csv(args.input_csv, args.output_csv, sleep=args.sleep, workers=args.workers)


if __name__ == "__main__":
    main()