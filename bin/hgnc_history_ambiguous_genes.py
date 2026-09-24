#!/usr/bin/env python3
"""
hgnc_history_ambigous_genes.py

Compile the list of HGNC gene symbols that may have been misannotated according to the predicted
functional probability and fetch the hgnc symbol change history for them.

Ambiguous genes are detected from the `Gene_Type` field vs `functional_probability`, and can
optionally be supplemented with specific genes of interest via --goi (e.g. --goi RNU2-2).

Output a nicely formatted csv containing:
gene_symbol, current_annotation, predicted_annotation, functional_probability, last_modified_in_hgnc
"""

import argparse
import csv
import sys
import time

import pandas as pd
import requests

PATH = "https://rest.genenames.org/fetch/symbol/"
TIMEOUT = 0.5  # polite delay between HGNC REST requests (seconds)
REQUEST_TIMEOUT = 10  # seconds, per-request HTTP timeout

HEADERS = {"Accept": "application/json"}


def fetch_hgnc_data(symbol):
    """
    Query the HGNC REST API for a gene symbol and return its change history.

    Returns a dict with:
        - hgnc_current_symbol: the symbol currently recognised by HGNC
        - prev_symbols: list of previous/withdrawn symbols
        - date_modified: last modification date recorded by HGNC

    Returns None if the symbol could not be found or the request failed.
    """
    url = f"{PATH}{symbol}"
    try:
        response = requests.get(url, headers=HEADERS, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
    except requests.RequestException as exc:
        print(f"Warning: could not fetch HGNC data for '{symbol}': {exc}", file=sys.stderr)
        return None

    try:
        payload = response.json()
        docs = payload.get("response", {}).get("docs", [])
    except ValueError as exc:
        print(f"Warning: invalid JSON returned for '{symbol}': {exc}", file=sys.stderr)
        return None

    if not docs:
        print(f"Warning: no HGNC record found for '{symbol}'", file=sys.stderr)
        return None

    doc = docs[0]
    return {
        "hgnc_current_symbol": doc.get("symbol", symbol),
        "prev_symbols": doc.get("prev_symbol", []),
        "date_modified": doc.get("date_modified", "N/A"),
    }


def get_ambiguous_genes(input_csv, goi=None):
    """
    Read the input csv and identify ambiguous / potentially misannotated genes.

    A gene is considered ambiguous if:
      - its `Gene_Type` already marks it as ambiguous (e.g. genes such as RNU2-2
        that are typed as "Ambiguous" in the input), OR
      - it is typed as a functional/protein-coding gene (i.e. not a pseudogene)
        but functional_probability < 0.5, OR
      - it is typed as a pseudogene but functional_probability > 0.5

    `goi` (genes of interest) is an optional iterable of gene symbols (e.g. from
    --goi) that are force-included in the result regardless of the automatic
    criteria above. If a requested gene is present in the input csv, its actual
    Gene_Type / functional_probability are used; if it is not present, it is
    still included (with "N/A" placeholders) so its HGNC history can be fetched.

    Returns a DataFrame with columns:
        gene_symbol, current_annotation, predicted_annotation, functional_probability
    """
    df = pd.read_csv(input_csv)

    required_cols = {"Gene", "Gene_Type", "functional_probability"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Input csv is missing required column(s): {sorted(missing)}")

    gene_type_norm = df["Gene_Type"].astype(str).str.strip().str.lower()
    prob = df["functional_probability"]

    is_ambiguous_type = gene_type_norm.str.contains("ambig", na=False)
    is_pseudogene = gene_type_norm.str.contains("pseudo", na=False)
    is_functional_but_low = ~is_pseudogene & ~is_ambiguous_type & (prob < 0.5)
    is_pseudogene_but_high = is_pseudogene & (prob > 0.5)

    mask = is_ambiguous_type | is_functional_but_low | is_pseudogene_but_high
    ambiguous = df.loc[mask].copy()

    ambiguous["predicted_annotation"] = ambiguous["functional_probability"].apply(
        lambda p: "Functional" if p >= 0.5 else "Pseudogene"
    )

    result = ambiguous.rename(columns={"Gene": "gene_symbol", "Gene_Type": "current_annotation"})[
        ["gene_symbol", "current_annotation", "predicted_annotation", "functional_probability"]
    ]

    if goi:
        goi_set = {g.strip() for g in goi if g.strip()}
        already_included = set(result["gene_symbol"])
        extra_needed = goi_set - already_included

        if extra_needed:
            df_by_gene = df.set_index("Gene")
            extra_rows = []
            for symbol in sorted(extra_needed):
                if symbol in df_by_gene.index:
                    row = df_by_gene.loc[symbol]
                    if isinstance(row, pd.DataFrame):  # duplicate Gene entries
                        row = row.iloc[0]
                    prob_val = row["functional_probability"]
                    current = row["Gene_Type"]
                    predicted = "Functional" if prob_val >= 0.5 else "Pseudogene"
                else:
                    print(
                        f"Warning: gene of interest '{symbol}' not found in input csv; "
                        "will still fetch its HGNC history.",
                        file=sys.stderr,
                    )
                    prob_val, current, predicted = "N/A", "N/A", "N/A"

                extra_rows.append(
                    {
                        "gene_symbol": symbol,
                        "current_annotation": current,
                        "predicted_annotation": predicted,
                        "functional_probability": prob_val,
                    }
                )
            result = pd.concat([result, pd.DataFrame(extra_rows)], ignore_index=True)

    result = result.drop_duplicates(subset="gene_symbol").reset_index(drop=True)
    return result


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compile HGNC gene symbols that may be misannotated (based on predicted "
            "functional probability) and fetch their HGNC symbol change history."
        )
    )
    parser.add_argument("input_csv", help="Path to the input csv of gene predictions")
    parser.add_argument(
        "-o",
        "--output",
        default="results/ambiguous_genes_hgnc_history.csv",
        help="Path to write the output csv (default: %(default)s)",
    )
    parser.add_argument(
        "--goi",
        nargs="+",
        default=None,
        metavar="SYMBOL",
        help=(
            "Optional gene(s) of interest to force-include in the ambiguous set, "
            "in addition to those detected automatically (e.g. --goi RNU2-2)"
        ),
    )
    args = parser.parse_args()

    ambiguous_df = get_ambiguous_genes(args.input_csv, goi=args.goi)

    if ambiguous_df.empty:
        print("No ambiguous genes found.")
        return

    rows = []
    total = len(ambiguous_df)
    for i, row in ambiguous_df.iterrows():
        symbol = row["gene_symbol"]
        print(f"[{i + 1}/{total}] Fetching HGNC history for {symbol}...", file=sys.stderr)

        hgnc_data = fetch_hgnc_data(symbol)
        last_modified = hgnc_data["date_modified"] if hgnc_data else "N/A"
        prev_symbols = hgnc_data["prev_symbols"] if hgnc_data else "N/A"

        rows.append(
            {
                "gene_symbol": symbol,
                "current_annotation": row["current_annotation"],
                "predicted_annotation": row["predicted_annotation"],
                "functional_probability": row["functional_probability"],
                "last_modified_in_hgnc": last_modified,
                "prev_symbols": prev_symbols
            }
        )

        if i < total - 1:
            time.sleep(TIMEOUT)

    fieldnames = [
        "gene_symbol",
        "current_annotation",
        "predicted_annotation",
        "functional_probability",
        "last_modified_in_hgnc",
        "prev_symbols"
    ]
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} ambiguous gene(s) to {args.output}")


if __name__ == "__main__":
    main()