"""
fishers_exact_test.py

Performs Fisher's exact test on a 2x2 contingency table:
has_polIII_motif vs no_polIII_motif
high SNP enrichment vs not high SNP enrichment

"High SNP enrichment" is defined with the same empirical-null FDR approach as the
R function fdr_threshold(): pseudogenes form the null distribution, each functional
gene gets an empirical p-value (1 + #null >= value) / (n_null + 1), p-values are
BH-adjusted, and the threshold is the smallest enrichment value with q <= FDR.

Author: Nadia Prasetyo
"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import argparse


def bh_adjust(p):
    """Benjamini-Hochberg adjustment (equivalent to R's p.adjust(method = "BH"))."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    # enforce monotonicity from the largest p-value downwards
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.minimum(ranked, 1.0)
    return q


def fdr_threshold(x, cls, null_label="pseudogene", func_label="functional", fdr=0.05):
    """Python port of the R fdr_threshold() function."""
    x = pd.Series(x).reset_index(drop=True)
    cls = pd.Series(cls).reset_index(drop=True)

    null = x[(cls == null_label) & x.notna()].to_numpy()
    func = (x.notna() & (cls == func_label)).to_numpy()

    p = np.full(len(x), np.nan)
    q = np.full(len(x), np.nan)
    p[func] = [(1 + np.sum(null >= v)) / (len(null) + 1) for v in x[func]]
    q[func] = bh_adjust(p[func])

    sig = ~np.isnan(q) & (q <= fdr)
    threshold = x[sig].min() if sig.any() else np.nan
    return {
        "p": p, "q": q, "sig": sig, "func": func,
        "threshold": threshold,
        "n_null": len(null), "n_func": int(func.sum()),
    }


def setup_2x2_table(df):
    # Create a contingency table
    contingency_table = pd.crosstab(df['has_polIII_motif'], df['high_SNP_enrichment'])

    # Make sure all four cells exist even if a category is empty
    contingency_table = contingency_table.reindex(index=[True, False], columns=[True, False], fill_value=0)

    # Give the row/column axes clear, unambiguous labels
    contingency_table.index = contingency_table.index.map(
        {True: 'has_polIII_motif', False: 'no_polIII_motif'}
    )
    contingency_table.columns = contingency_table.columns.map(
        {True: 'high_SNP_enrichment', False: 'not_high_SNP_enrichment'}
    )
    contingency_table.index.name = None
    contingency_table.columns.name = None

    print("Contingency table created:")
    return contingency_table


def process_and_count_input_file(input_file, class_column, null_label, func_label, fdr):
    # input csv has the following columns: "gene_name","expected_promoter",
    # "contains_type1_promoter","contains_type2_promoter","contains_type3_promoter",
    # "gene_group","locus_type","gene_type","enrichment_intergenic_dbsnp"
    df = pd.read_csv(input_file)
    print(f"Input file {input_file} read successfully. Number of rows: {len(df)}")

    # Create a new column 'has_polIII_motif' based on the presence of any polIII motif based on their counts:
    # type1_motif_counts,type2_motif_counts,type3_motif_counts
    df['has_polIII_motif'] = (df['type1_motif_counts'] + df['type2_motif_counts'] + df['type3_motif_counts']) > 0

    # High SNP enrichment threshold: FDR against the pseudogene null
    res = fdr_threshold(df['enrichment_intergenic_dbsnp'], df[class_column],
                        null_label=null_label, func_label=func_label, fdr=fdr)
    print(f"FDR = {fdr}: n_null ({null_label}) = {res['n_null']}, "
          f"n_func ({func_label}) = {res['n_func']}, "
          f"n_sig = {int(res['sig'].sum())}, threshold = {res['threshold']}")

    if np.isnan(res['threshold']):
        raise SystemExit("No functional genes pass the FDR threshold; cannot define high SNP enrichment.")
    if res['n_func'] == 0:
        raise SystemExit(f"No rows with {class_column} == '{func_label}'. Check --class_column / --functional_label.")

    # Apply the FDR-derived threshold to every gene (functional and pseudogene).
    # For functional genes, enrichment >= threshold is identical to q <= fdr.
    # Rows with missing enrichment can't be classified, so they are dropped.
    n_missing = int(df['enrichment_intergenic_dbsnp'].isna().sum())
    if n_missing:
        print(f"Dropping {n_missing} rows with missing enrichment_intergenic_dbsnp")
    df = df[df['enrichment_intergenic_dbsnp'].notna()].copy()
    df['high_SNP_enrichment'] = df['enrichment_intergenic_dbsnp'] >= res['threshold']
    print(f"Genes in contingency table: {len(df)}")
    return df


def perform_fishers_exact_test(contingency_table):
    # Perform Fisher's exact test
    odds_ratio, p_value = fisher_exact(contingency_table)
    print(f"Fisher's exact test results: Odds Ratio = {odds_ratio}, P-value = {p_value}")
    return odds_ratio, p_value


def main():
    parser = argparse.ArgumentParser(description="Perform Fisher's exact test on a 2x2 contingency table.")
    parser.add_argument('-i', '--input_file', help='Path to the input CSV file.')
    parser.add_argument('-o', '--output_file', default='results/fishers_exact_test_results.csv',
                        help='Path to the output CSV file for results.')
    parser.add_argument('--fdr', type=float, default=0.05, help='FDR level for the high-SNP threshold.')
    parser.add_argument('--class_column', default='Gene_Type',
                        help='Column that distinguishes pseudogenes (null) from functional genes.')
    parser.add_argument('--null_label', default='Pseudogene', help='Value in class_column for the null set.')
    parser.add_argument('--functional_label', default='Functional',
                        help='Value in class_column for functional genes.')
    args = parser.parse_args()

    # Process the input file and create the contingency table
    df = process_and_count_input_file(args.input_file, args.class_column,
                                      args.null_label, args.functional_label, args.fdr)
    contingency_table = setup_2x2_table(df)
    print(contingency_table)

    # Perform Fisher's exact test
    odds_ratio, p_value = perform_fishers_exact_test(contingency_table)

    # save the results to a CSV file
    results_df = pd.DataFrame({
        'Odds Ratio': [odds_ratio],
        'P-value': [p_value]
    })

    # save the contingency table to a CSV file
    contingency_table.to_csv('results/contingency_table.csv', index_label='polIII_motif_status')
    results_df.to_csv(args.output_file, index=False)
    print(f"Results saved to {args.output_file} and contingency table saved to results/contingency_table.csv")


if __name__ == "__main__":
    main()