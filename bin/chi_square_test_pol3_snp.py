"""
chi_square_test_pol3_snp.py

Performs three Chi-square tests of independence, each on its own
contingency table:

    1. SNP enrichment category         vs  PolIII motif status
    2. gene_type                       vs  PolIII motif status
    3. SNP enrichment category         vs  gene_type

For every test the script reports the Chi-square statistic, p-value,
degrees of freedom, and the Chi-square critical value at alpha = 0.05
(so the statistic can be compared against the critical value directly,
in addition to using the p-value).

All three contingency tables are written into a single CSV file, each
preceded by a '#' comment line identifying which table it is.
All three test summaries are written into a single results file
(CSV by default; use --output_format txt for a plain-text summary).

Author: Nadia Prasetyo
"""
import os
import argparse
import pandas as pd
from scipy import stats


def categorize_snp_enrichment(df):
    # Create a new column to divide SNP enrichment into high, medium, and low
    # categories based on quantiles.
    quantiles = df['enrichment_intergenic_dbsnp'].quantile([0.25, 0.75])
    low_threshold = quantiles[0.25]   # snp enrichment values below this threshold are considered low
    high_threshold = quantiles[0.75]  # snp enrichment values above this threshold are considered high
    # snp enrichment in between these thresholds is considered medium
    df['snp_enrichment_category'] = pd.cut(
        df['enrichment_intergenic_dbsnp'],
        bins=[-float('inf'), low_threshold, high_threshold, float('inf')],
        labels=['low', 'medium', 'high']
    )
    return df


def determine_polIII_motif(row):
    # Determine the observed_polIII_motif category based on the presence of polIII motifs
    # if contains_type1_promoter is the only True, then observed_polIII_motif = type1
    # if contains_type2_promoter is the only True, then observed_polIII_motif = type2
    # if contains_type3_promoter is the only True, then observed_polIII_motif = type3
    # if more than one of the contains_type*_promoter columns are True, then we have
    # type1_and_type2, type1_and_type3, type2_and_type3, or type1_and_type2_and_type3
    if row['contains_type1_promoter'] and not row['contains_type2_promoter'] and not row['contains_type3_promoter']:
        return 'type1'
    elif not row['contains_type1_promoter'] and row['contains_type2_promoter'] and not row['contains_type3_promoter']:
        return 'type2'
    elif not row['contains_type1_promoter'] and not row['contains_type2_promoter'] and row['contains_type3_promoter']:
        return 'type3'
    elif row['contains_type1_promoter'] and row['contains_type2_promoter'] and not row['contains_type3_promoter']:
        return 'type1_and_type2'
    elif row['contains_type1_promoter'] and not row['contains_type2_promoter'] and row['contains_type3_promoter']:
        return 'type1_and_type3'
    elif not row['contains_type1_promoter'] and row['contains_type2_promoter'] and row['contains_type3_promoter']:
        return 'type2_and_type3'
    elif row['contains_type1_promoter'] and row['contains_type2_promoter'] and row['contains_type3_promoter']:
        return 'type1_and_type2_and_type3'
    else:
        return 'no_polIII_motif'


def process_and_count_input_file(input_file):
    # Read the input CSV file into a DataFrame
    # input csv has the following columns: "gene_name","expected_promoter",
    # "contains_type1_promoter","contains_type2_promoter","contains_type3_promoter",
    # "gene_group","locus_type","gene_type","enrichment_intergenic_dbsnp"
    df = pd.read_csv(input_file)
    print(f"Input file {input_file} read successfully. Number of rows: {len(df)}")

    # create a new column of category: observed_polIII_motif
    df['observed_polIII_motif'] = df.apply(determine_polIII_motif, axis=1)

    # Create a new column 'snp_enrichment_category' based on the enrichment_intergenic_dbsnp value
    df = categorize_snp_enrichment(df)
    return df


def setup_contingency_table(df, row_col, col_col, row_order=None, col_order=None):
    """
    Build a contingency table (crosstab) of counts between two categorical
    columns of df. If row_order / col_order are provided, the table is
    reindexed to that fixed category order (missing categories are filled
    with 0); otherwise the categories found in the data are used, sorted.
    """
    contingency_table = pd.crosstab(df[row_col], df[col_col])

    if row_order is not None:
        contingency_table = contingency_table.reindex(index=row_order, fill_value=0)
    else:
        contingency_table = contingency_table.sort_index(axis=0)

    if col_order is not None:
        contingency_table = contingency_table.reindex(columns=col_order, fill_value=0)
    else:
        contingency_table = contingency_table.sort_index(axis=1)

    return contingency_table


def perform_chi_square_test(contingency_table, alpha=0.05):
    """
    Run a Chi-square test of independence on a contingency table.

    Empty rows/columns produce zero expected frequencies and cause
    scipy.stats.chi2_contingency to fail, so they are excluded from the
    test while the complete table (with the empty rows/columns) is still
    kept for reporting purposes.

    Returns a dict with the Chi-square statistic, p-value, degrees of
    freedom, and the Chi-square critical value at the given alpha
    (default 0.05).
    """
    test_table = contingency_table.loc[
        contingency_table.sum(axis=1).gt(0),
        contingency_table.sum(axis=0).gt(0),
    ]
    if test_table.shape[0] < 2 or test_table.shape[1] < 2:
        raise ValueError(
            "At least two non-empty rows and columns are required for a "
            "chi-square test."
        )

    chi2, p_value, dof, expected = stats.chi2_contingency(test_table)
    critical_value = stats.chi2.ppf(1 - alpha, dof)
    significant = chi2 > critical_value

    print(
        f"Chi-square test results: Chi-square = {chi2}, P-value = {p_value}, "
        f"dof = {dof}, critical value (alpha={alpha}) = {critical_value}"
    )

    return {
        'Chi-square': chi2,
        'P-value': p_value,
        'Degrees_of_freedom': dof,
        'Critical_value_0.05': critical_value,
        'Significant_at_0.05': significant,
    }


def write_contingency_tables_csv(tables, output_path):
    """
    Write a list of (title, index_label, dataframe) contingency tables into
    a single CSV file, each preceded by a '#' comment line naming the table.
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w', newline='') as f:
        for i, (title, index_label, table) in enumerate(tables):
            f.write(f"# {title}\n")
            table.to_csv(f, index_label=index_label)
            if i != len(tables) - 1:
                f.write("\n")
    print(f"Contingency tables saved to {output_path}")


def write_results_summary(results, output_path, output_format):
    """
    Write the summary of all Chi-square test results (Chi-square, p-value,
    degrees of freedom, and 0.05 critical value) to a CSV or txt file.
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    results_df = pd.DataFrame([
        {'Test': name, **stats_dict} for name, stats_dict in results.items()
    ])

    if output_format == 'csv':
        results_df.to_csv(output_path, index=False)
    else:  # txt
        with open(output_path, 'w') as f:
            for name, stats_dict in results.items():
                f.write(f"# {name}\n")
                f.write(f"Chi-square           : {stats_dict['Chi-square']}\n")
                f.write(f"P-value              : {stats_dict['P-value']}\n")
                f.write(f"Degrees of freedom   : {stats_dict['Degrees_of_freedom']}\n")
                f.write(f"Critical value (0.05): {stats_dict['Critical_value_0.05']}\n")
                f.write(f"Significant at 0.05  : {stats_dict['Significant_at_0.05']}\n")
                f.write("\n")
    print(f"Results summary saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Perform three Chi-square tests: (1) SNP enrichment category vs "
                    "PolIII motif status, (2) gene_type vs PolIII motif status, and "
                    "(3) SNP enrichment category vs gene_type."
    )
    parser.add_argument('-i', '--input_file', help='Path to the input CSV file.', required=True)
    parser.add_argument('-o', '--output_file', default='results/chi_square_test_results.csv',
                         help='Path to the output file for the summary of results (Chi-square, '
                              'P-value, degrees of freedom, critical value).')
    parser.add_argument('-c', '--contingency_file', default='results/chi_squared_contingency_tables.csv',
                         help='Path to the output CSV file for all contingency tables.')
    parser.add_argument('-f', '--output_format', choices=['csv', 'txt'], default='csv',
                         help="Format for the results summary file: 'csv' (default) or 'txt'.")
    parser.add_argument('-a', '--alpha', type=float, default=0.05,
                         help='Significance level used for the critical value (default: 0.05).')
    args = parser.parse_args()

    # Process the input file and derive the categorical columns needed
    df = process_and_count_input_file(args.input_file)

    promoter_order = [
        'type1', 'type2', 'type3',
        'type1_and_type2', 'type1_and_type3',
        'type2_and_type3', 'type1_and_type2_and_type3',
        'no_polIII_motif'
    ]
    snp_order = ['low', 'medium', 'high']

    # --- Test 1: SNP enrichment category vs PolIII motif status ---
    table1 = setup_contingency_table(
        df, 'snp_enrichment_category', 'observed_polIII_motif',
        row_order=snp_order, col_order=promoter_order
    )
    results1 = perform_chi_square_test(table1, alpha=args.alpha)

    # --- Test 2: gene_type vs PolIII motif status ---
    table2 = setup_contingency_table(
        df, 'gene_type', 'observed_polIII_motif',
        row_order=None, col_order=promoter_order
    )
    results2 = perform_chi_square_test(table2, alpha=args.alpha)

    # --- Test 3: SNP enrichment category vs gene_type ---
    table3 = setup_contingency_table(
        df, 'snp_enrichment_category', 'gene_type',
        row_order=snp_order, col_order=None
    )
    results3 = perform_chi_square_test(table3, alpha=args.alpha)

    # Save all three contingency tables into a single CSV, each labeled with a '#' comment
    write_contingency_tables_csv(
        [
            ("Contingency table: snp_enrichment_category (rows) vs observed_polIII_motif (columns)",
             'snp_enrichment_category', table1),
            ("Contingency table: gene_type (rows) vs observed_polIII_motif (columns)",
             'gene_type', table2),
            ("Contingency table: snp_enrichment_category (rows) vs gene_type (columns)",
             'snp_enrichment_category', table3),
        ],
        args.contingency_file
    )

    # Save all three test summaries into a single results file
    results = {
        'snp_enrichment_category_vs_polIII_status': results1,
        'gene_type_vs_polIII_status': results2,
        'snp_enrichment_category_vs_gene_type': results3,
    }
    write_results_summary(results, args.output_file, args.output_format)


if __name__ == "__main__":
    main()