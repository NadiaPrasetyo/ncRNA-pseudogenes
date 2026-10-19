#!/usr/bin/env python3
"""
hgnc_history_ambigous_genes.py

Compile the list of ambiguous HGNC gene symbols and fetch the hgnc symbol change history for them.

By default the ambiguous set is the curated LDA_OUTLIERS list. It can optionally be supplemented with:
  --include-low-functional   genes typed as functional but with functional_probability < 0.5
  --include-high-pseudogene  genes typed as pseudogene but with functional_probability > 0.5
  --goi SYMBOL [SYMBOL ...]  specific genes of interest (e.g. --goi RNU2-2)

Output a nicely formatted csv containing:
gene_symbol, current_annotation, predicted_annotation, functional_probability,
selection_reason, last_modified_in_hgnc, prev_symbols
"""

import argparse
import csv
import os
import sys
import time

import pandas as pd
import requests

PATH = "https://rest.genenames.org/fetch/symbol/"
TIMEOUT = 0.5  # polite delay between HGNC REST requests (seconds)
REQUEST_TIMEOUT = 10  # seconds, per-request HTTP timeout

HEADERS = {"Accept": "application/json"}

LDA_OUTLIERS = {"RN7SL3", "VTRNA1-3", "RNU1-4", "RNU1-3", "VTRNA1-1", "VTRNA2-1", "RN7SKP1", "RNU6-1171P", "RNU6-184P",
"RNU6-109P", "RNU6-1076P", "RNU6-1083P", "RN7SL275P", "RN7SL560P", "RNU1-27P", "RNU1-28P", "RNU2-4P", "RNU2-63P", "RNU4-63P",
"RNU5A-8P", "RNU5B-5P", "RNU5E-4P", "RNU5E-6P", "RNU5F-5P", "RNU6-1194P", "RNU6-1334P", "RNU6-202P", "RNU6-693P", "RNU6-781P",
"RNU4ATAC10P", "RNU4ATAC11P", "RNU4ATAC13P", "RNU4ATAC15P", "RNU4ATAC18P", "RNU4ATAC2P", "RNU4ATAC4P", "RNU4ATAC5P", "RNU4ATAC9P",
"RNU6ATAC10P", "RNU6ATAC11P", "RNU6ATAC15P", "RNU6ATAC16P", "RNU6ATAC18P", "RNU6ATAC20P", "RNU6ATAC21P", "RNU6ATAC23P",
"RNU6ATAC25P", "RNU6ATAC26P", "RNU6ATAC27P", "RNU6ATAC29P", "RNU6ATAC2P", "RNU6ATAC30P", "RNU6ATAC31P", "RNU6ATAC33P",
"RNU6ATAC34P", "RNU6ATAC37P", "RNU6ATAC38P", "RNU6ATAC3P", "RNU6ATAC41P", "RNU6ATAC4P", "RNU6ATAC5P", "RNU6ATAC6P",
"RNU6ATAC7P", "RNU6ATAC8P", "RNU11-3P", "RNU11-4P", "RNU11-6P", "RNU12-2P", "RN7SL471P", "RN7SKP253", "RN7SKP70"}

REQUIRED_COLS = {"Gene", "Gene_Type", "functional_probability"}
OUT_COLS = ["gene_symbol", "current_annotation", "predicted_annotation", "functional_probability", "selection_reason"]


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


def predicted_label(prob):
    """Map a functional probability to a predicted annotation."""
    if pd.isna(prob):
        return "N/A"
    return "Functional" if prob >= 0.5 else "Pseudogene"


def get_ambiguous_genes(input_csv, include_low_functional=False, include_high_pseudogene=False, goi=None):
    """
    Read the input csv and identify the ambiguous genes.

    The ambiguous set is, by default, the LDA_OUTLIERS list. Optionally add:
      - include_low_functional: genes typed as functional (not pseudogene) with
        functional_probability < 0.5
      - include_high_pseudogene: genes typed as pseudogene with
        functional_probability > 0.5
      - goi: iterable of gene symbols that are force-included regardless of the
        criteria above.

    Any LDA outlier or gene of interest that is absent from the input csv is
    still included (with "N/A" placeholders) so its HGNC history can be fetched.

    Returns a DataFrame with columns:
        gene_symbol, current_annotation, predicted_annotation,
        functional_probability, selection_reason
    """
    df = pd.read_csv(input_csv)

    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        raise ValueError(f"Input csv is missing required column(s): {sorted(missing)}")

    gene_type_norm = df["Gene_Type"].astype(str).str.strip().str.lower()
    prob = df["functional_probability"]
    is_pseudogene = gene_type_norm.str.contains("pseudo", na=False)

    is_outlier = df["Gene"].isin(LDA_OUTLIERS)
    is_low_functional = (
        ~is_outlier & ~is_pseudogene & (prob < 0.5) if include_low_functional
        else pd.Series(False, index=df.index)
    )
    is_high_pseudogene = (
        ~is_outlier & is_pseudogene & (prob > 0.5) if include_high_pseudogene
        else pd.Series(False, index=df.index)
    )

    reason = pd.Series("", index=df.index)
    reason[is_low_functional] = "functional_low_probability"
    reason[is_high_pseudogene] = "pseudogene_high_probability"
    reason[is_outlier] = "LDA_outlier"

    selected = df.loc[is_outlier | is_low_functional | is_high_pseudogene].copy()
    selected["selection_reason"] = reason.loc[selected.index]
    selected["predicted_annotation"] = selected["functional_probability"].apply(predicted_label)

    result = selected.rename(columns={"Gene": "gene_symbol", "Gene_Type": "current_annotation"})[OUT_COLS]

    # Genes that must be present but weren't selected from the csv:
    #   - LDA outliers missing from the csv
    #   - genes of interest not already selected
    goi_set = {g.strip() for g in goi if g.strip()} if goi else set()
    already_included = set(result["gene_symbol"])
    forced = [(s, "LDA_outlier") for s in sorted(LDA_OUTLIERS - already_included)]
    forced += [(s, "goi") for s in sorted(goi_set - already_included - LDA_OUTLIERS)]

    if forced:
        df_by_gene = df.drop_duplicates(subset="Gene").set_index("Gene")
        extra_rows = []
        for symbol, why in forced:
            if symbol in df_by_gene.index:
                row = df_by_gene.loc[symbol]
                prob_val = row["functional_probability"]
                current = row["Gene_Type"]
                predicted = predicted_label(prob_val)
            else:
                print(
                    f"Warning: '{symbol}' ({why}) not found in input csv; "
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
                    "selection_reason": why,
                }
            )
        result = pd.concat([result, pd.DataFrame(extra_rows)], ignore_index=True)

    result = result.drop_duplicates(subset="gene_symbol").reset_index(drop=True)
    return result


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compile ambiguous HGNC gene symbols (LDA outliers by default) and "
            "fetch their HGNC symbol change history."
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
        "--include-low-functional",
        action="store_true",
        help="Also include genes typed as functional but with functional_probability < 0.5",
    )
    parser.add_argument(
        "--include-high-pseudogene",
        action="store_true",
        help="Also include genes typed as pseudogene but with functional_probability > 0.5",
    )
    parser.add_argument(
        "--goi",
        nargs="+",
        default=None,
        metavar="SYMBOL",
        help=(
            "Optional gene(s) of interest to force-include in the ambiguous set "
            "(e.g. --goi RNU2-2)"
        ),
    )
    args = parser.parse_args()

    ambiguous_df = get_ambiguous_genes(
        args.input_csv,
        include_low_functional=args.include_low_functional,
        include_high_pseudogene=args.include_high_pseudogene,
        goi=args.goi,
    )

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
                "selection_reason": row["selection_reason"],
                "last_modified_in_hgnc": last_modified,
                "prev_symbols": prev_symbols,
            }
        )

        if i < total - 1:
            time.sleep(TIMEOUT)

    fieldnames = [
        "gene_symbol",
        "current_annotation",
        "predicted_annotation",
        "functional_probability",
        "selection_reason",
        "last_modified_in_hgnc",
        "prev_symbols",
    ]

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} ambiguous gene(s) to {args.output}")


if __name__ == "__main__":
    main()