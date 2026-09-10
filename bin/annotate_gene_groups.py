#!/usr/bin/env python3
"""
annotate_gene_groups.py

Classify each `gene_name` in an input CSV into one of 15 small-RNA gene
family groups, using ONLY symbol lists actually returned by the HGNC
REST API (rest.genenames.org) -- no guessing from regex alone.

    Groups: RN7SK, RN7SL, RNU11, RNU12, RNU1, RNU2, RNU4ATAC, RNU4,
            RNU5, RNU6ATAC, RNU6, RNU7, RNY, TRNA, VTRNA

Method
------
1. For each of the 15 groups, query HGNC's /search endpoint with the
   group name as a wildcard substring (e.g. `*RNU1*`, `*RNU11*`, ...).
   This is a SUBSTRING search, so e.g. the "RNU1" query will also
   return RNU11, RNU12, RNU1-1, etc.

2. To resolve that overlap, every symbol HGNC returns (across ALL 15
   queries, pooled together) is routed to exactly one group using an
   ordered set of prefix rules that mirror real HGNC nomenclature
   (most specific prefixes -- RNU4ATAC, RNU6ATAC, RNU11, RNU12 --
   are checked before the shorter prefixes they'd otherwise collide
   with -- RNU4, RNU6, RNU1). The rules are only used to SORT symbols
   HGNC already returned as Approved -- they never invent a symbol
   HGNC didn't give us.

3. The input CSV's `gene_name` column is then cross-mapped against
   these HGNC-derived group sets by exact (case-insensitive) symbol
   match. Anything that isn't in any group's HGNC-derived set is
   reported as unmatched, both on stdout and in a separate CSV.

Optionally (--resolve-aliases), unmatched gene_names are given one
more chance: `search_hgnc_genes(gene_name)` is called directly on the
name itself, in case it's an alias/previous symbol for an approved
symbol that already sits in one of the group sets.

Usage
-----
    python hgnc_group_crossmap.py -i input.csv -o output.csv
    python hgnc_group_crossmap.py -i input.csv -o output.csv \
        --unmatched-file unmatched.csv --resolve-aliases
    python hgnc_group_crossmap.py -i input.csv -o output.csv \
        --group-cache hgnc_groups_cache.json
"""

import argparse
import csv
import json
import os
import re
import sys
import time

import requests

# --------------------------------------------------------------------------
# HGNC API helpers
# --------------------------------------------------------------------------


def search_hgnc_genes(query):
    if not query or not isinstance(query, str):
        print("Invalid query: must be a non-empty string.")
        return []
    url = f'https://rest.genenames.org/search/*{query}*+AND+status:%22Approved%22'
    headers = {'Accept': 'application/json'}
    try:
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Request error: {e}")
        return []

    try:
        data = response.json()
        if data['response']['numFound'] > 0:
            return [gene['symbol'] for gene in data['response']['docs']]
        else:
            print(f"No genes found for query: {query}")
            return []
    except (ValueError, KeyError) as e:
        print(f"Error parsing response: {e}")
        return []


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


# --------------------------------------------------------------------------
# Target groups + routing rules
#
# These regexes are ONLY used to route a symbol that HGNC already told us
# is Approved into exactly one of the 15 groups (needed because the
# substring searches overlap). They are never used to accept a symbol
# HGNC didn't return.
# --------------------------------------------------------------------------

GROUPS = [
    "RN7SK", "RN7SL", "RNU11", "RNU12", "RNU1", "RNU2", "RNU4ATAC",
    "RNU4", "RNU5", "RNU6ATAC", "RNU6", "RNU7", "RNY", "TRNA", "VTRNA",
]

# Order matters: most specific prefixes first.
ROUTING_RULES = [
    ("RNU4ATAC", re.compile(r'^RNU4ATAC', re.IGNORECASE)),
    ("RNU6ATAC", re.compile(r'^RNU6ATAC', re.IGNORECASE)),
    ("RNU11",    re.compile(r'^RNU11(?!\d)', re.IGNORECASE)),
    ("RNU12",    re.compile(r'^RNU12(?!\d)', re.IGNORECASE)),
    ("RNU1",     re.compile(r'^RNU1(?!\d)', re.IGNORECASE)),
    ("RNU2",     re.compile(r'^RNU2(?!\d)', re.IGNORECASE)),
    ("RNU4",     re.compile(r'^RNU4(?!\d)', re.IGNORECASE)),
    ("RNU5",     re.compile(r'^RNU5(?!\d)', re.IGNORECASE)),
    ("RNU6",     re.compile(r'^RNU6(?!\d)', re.IGNORECASE)),
    ("RNU7",     re.compile(r'^RNU7(?!\d)', re.IGNORECASE)),
    ("RN7SK",    re.compile(r'^RN7SK', re.IGNORECASE)),
    ("RN7SL",    re.compile(r'^RN7SL', re.IGNORECASE)),
    ("RNY",      re.compile(r'^RNY\d', re.IGNORECASE)),
    ("VTRNA",    re.compile(r'^VTRNA', re.IGNORECASE)),
    ("TRNA",     re.compile(r'^(TRNA(-|\d)|MT-T[A-Z])', re.IGNORECASE)),
]

UNMATCHED = "Unmatched"


def route_symbol(symbol):
    """Given an HGNC-approved symbol, return which of the 15 groups it
    belongs to, or None if it doesn't fit any of them."""
    if not symbol:
        return None
    s = symbol.strip()
    for group, pattern in ROUTING_RULES:
        if pattern.match(s):
            return group
    return None


# --------------------------------------------------------------------------
# Step 1 + 2: build authoritative group -> {symbols} sets from HGNC
# --------------------------------------------------------------------------


def build_group_symbol_sets(cache_path=None, api_delay=0.1):
    """
    Query HGNC once per group (substring search), pool all returned
    symbols, then route each pooled symbol into exactly one group.

    Returns: dict[group] -> set(symbols), plus a flat symbol->group map.
    """
    cache = {}
    if cache_path and os.path.exists(cache_path):
        try:
            with open(cache_path) as f:
                cache = json.load(f)
        except (ValueError, OSError):
            print(f"Warning: could not read group cache {cache_path}, refetching.")
            cache = {}

    pooled_symbols = set()

    for group in GROUPS:
        if group in cache:
            symbols = cache[group]
        else:
            print(f"Querying HGNC for group: {group}")
            symbols = search_hgnc_genes(group)
            cache[group] = symbols
            if api_delay:
                time.sleep(api_delay)
        pooled_symbols.update(symbols)

    if cache_path:
        try:
            with open(cache_path, "w") as f:
                json.dump(cache, f, indent=2, sort_keys=True)
        except OSError as e:
            print(f"Warning: could not write group cache {cache_path}: {e}")

    group_sets = {g: set() for g in GROUPS}
    symbol_to_group = {}
    for symbol in pooled_symbols:
        group = route_symbol(symbol)
        if group:
            group_sets[group].add(symbol)
            symbol_to_group[symbol.upper()] = group

    print("\nHGNC-derived group sizes:")
    for g in GROUPS:
        print(f"  {g:10s}: {len(group_sets[g])}")

    return group_sets, symbol_to_group


# --------------------------------------------------------------------------
# Step 2 + 3: cross-map the CSV's gene_name column, report unmatched
# --------------------------------------------------------------------------


def crossmap_csv(input_path, output_path, symbol_to_group,
                  unmatched_path=None, resolve_aliases=False, api_delay=0.1):
    with open(input_path, newline="", encoding="utf-8-sig") as f_in:
        reader = csv.DictReader(f_in)
        if reader.fieldnames is None or "gene_name" not in reader.fieldnames:
            print("Error: input CSV must have a 'gene_name' column.")
            sys.exit(1)

        fieldnames = list(reader.fieldnames) + ["gene_group"]
        rows_out = []
        unmatched_rows = []
        alias_cache = {}
        counts = {g: 0 for g in GROUPS}
        counts[UNMATCHED] = 0

        for row in reader:
            gene_name = (row.get("gene_name") or "").strip()
            group = symbol_to_group.get(gene_name.upper())

            if not group and resolve_aliases and gene_name:
                if gene_name in alias_cache:
                    candidates = alias_cache[gene_name]
                else:
                    candidates = search_hgnc_genes(gene_name)
                    alias_cache[gene_name] = candidates
                    if api_delay:
                        time.sleep(api_delay)
                for cand in candidates:
                    g = symbol_to_group.get(cand.upper())
                    if g:
                        group = g
                        break

            row["gene_group"] = group or UNMATCHED
            counts[row["gene_group"]] = counts.get(row["gene_group"], 0) + 1
            rows_out.append(row)
            if not group:
                unmatched_rows.append(row)

    with open(output_path, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)

    if unmatched_path and unmatched_rows:
        with open(unmatched_path, "w", newline="", encoding="utf-8") as f_un:
            writer = csv.DictWriter(f_un, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(unmatched_rows)

    print("\nCross-map summary:")
    total = len(rows_out)
    for g in GROUPS + [UNMATCHED]:
        n = counts.get(g, 0)
        if n:
            print(f"  {g:10s}: {n}")
    print(f"  {'TOTAL':10s}: {total}")

    if unmatched_rows:
        print(f"\n{len(unmatched_rows)} gene_name(s) did NOT fit any HGNC group:")
        for row in unmatched_rows:
            print(f"  - {row['gene_name']}")
        if unmatched_path:
            print(f"(also written to {unmatched_path})")
    else:
        print("\nEvery gene_name matched an HGNC group.")

    return output_path


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Cross-map a gene_name CSV column against HGNC-derived "
                     "small RNA gene family groups."
    )
    parser.add_argument("-i", "--input", required=True, help="Input CSV path")
    parser.add_argument("-o", "--output", required=True, help="Output CSV path")
    parser.add_argument(
        "--unmatched-file", default=None,
        help="Optional path to write rows whose gene_name matched no HGNC group"
    )
    parser.add_argument(
        "--group-cache", default="hgnc_groups_cache.json",
        help="JSON cache of the per-group HGNC search results "
             "(default: hgnc_groups_cache.json)"
    )
    parser.add_argument(
        "--resolve-aliases", action="store_true",
        help="For gene_names that don't directly match, also query HGNC "
             "search on the name itself to resolve aliases/previous symbols"
    )
    parser.add_argument(
        "--api-delay", type=float, default=0.1,
        help="Seconds to sleep between HGNC API calls (default: 0.1)"
    )
    args = parser.parse_args()

    group_sets, symbol_to_group = build_group_symbol_sets(
        cache_path=args.group_cache, api_delay=args.api_delay
    )

    crossmap_csv(
        input_path=args.input,
        output_path=args.output,
        symbol_to_group=symbol_to_group,
        unmatched_path=args.unmatched_file,
        resolve_aliases=args.resolve_aliases,
        api_delay=args.api_delay,
    )


if __name__ == "__main__":
    main()