"""
chi_square_test_pol3_snp.py

Performs Chi's square test on a 2x2 contingency table: 
has_polIII_motif vs no_polIII_motif 
high SNP enrichment vs not high SNP enrichment

Author: Nadia Prasetyo
"""
import sys
import pandas as pd
from scipy import stats
import argparse

def categorize_snp_enrichment(df):
    # Create a new column to divide SNP enrichment into high and medium, and low categories based on quantiles
    quantiles = df['enrichment_intergenic_dbsnp'].quantile([0.25, 0.75])
    low_threshold = quantiles[0.25] # snp enrichment values below this threshold are considered low
    high_threshold = quantiles[0.75] # snp enrichment values above this threshold are considered high
    #snp enrichment in between these thresholds are considered medium
    df['snp_enrichment_category'] = pd.cut(df['enrichment_intergenic_dbsnp'], bins=[-float('inf'), low_threshold, high_threshold, float('inf')], labels=['low', 'medium', 'high'])
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
    # input csv has the following columns: "gene_name","expected_promoter","
    # contains_type1_promoter","contains_type2_promoter","contains_type3_promoter",
    # "gene_group","locus_type","gene_type","enrichment_intergenic_dbsnp"
    df = pd.read_csv(input_file)
    print(f"Input file {input_file} read successfully. Number of rows: {len(df)}")
    # create a new column of category: observed_polIII_motif
    df['observed_polIII_motif'] = df.apply(determine_polIII_motif, axis=1)

    # Create a new column 'snp_enrichment_category' based on the enrichment_intergenic_dbsnp value
    df = categorize_snp_enrichment(df)
    return df

def setup_contingency_table(df):
    # Create a cross table containing counts for each promoter type and SNP
    # enrichment category.
    promoter_types = [
        'type1', 'type2', 'type3',
        'type1_and_type2', 'type1_and_type3',
        'type2_and_type3', 'type1_and_type2_and_type3',
        'no_polIII_motif'
    ]
    snp_categories = ['low', 'medium', 'high']

    contingency_table = pd.crosstab(
        df['observed_polIII_motif'],
        df['snp_enrichment_category']
    ).reindex(index=promoter_types, columns=snp_categories, fill_value=0)


    return contingency_table

def perform_chi_square_test(contingency_table):
    # Empty promoter categories produce zero expected frequencies and cause
    # scipy.stats.chi2_contingency to fail. Exclude them from the test while
    # retaining the complete table for reporting.
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
    print(f"Chi-square test results: Chi-square = {chi2}, P-value = {p_value}")
    return chi2, p_value

def main():
    parser = argparse.ArgumentParser(description="Perform Chi-square test on a 2x2 contingency table.")
    parser.add_argument('-i', '--input_file', help='Path to the input CSV file.')
    parser.add_argument('-o', '--output_file', default='results/chi_square_test_results.csv', help='Path to the output CSV file for results.')
    args = parser.parse_args()

    # Process the input file and create the contingency table
    df = process_and_count_input_file(args.input_file)
    contingency_table = setup_contingency_table(df)

    # Perform Chi-square test
    chi2, p_value = perform_chi_square_test(contingency_table)

    # save the results to a CSV file
    results_df = pd.DataFrame({
        'Chi-square': [chi2],
        'P-value': [p_value]
    })

    # save the contigency table to a CSV file
    contingency_table.to_csv('results/chi_squared_contingency_table.csv', index_label='polIII_motif_status')
    results_df.to_csv(args.output_file, index=False)
    print(f"Results saved to {args.output_file} and contingency table saved to results/contingency_table.csv")

if __name__ == "__main__":
    main()
    
