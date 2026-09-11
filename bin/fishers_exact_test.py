"""
fishers_exact_test.py

Performs Fisher's exact test on a 2x2 contingency table: 
has_polIII_motif vs no_polIII_motif 
high SNP enrichment vs not high SNP enrichment

Author: Nadia Prasetyo
"""
import sys
import pandas as pd
from scipy.stats import fisher_exact
import argparse

def setup_2x2_table(df):
    # Create a contingency table
    contingency_table = pd.crosstab(df['has_polIII_motif'], df['high_SNP_enrichment'])
    
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

def process_and_count_input_file(input_file):
    # Read the input CSV file into a DataFrame
    # input csv has the following columns: "gene_name","expected_promoter","
    # contains_type1_promoter","contains_type2_promoter","contains_type3_promoter",
    # "gene_group","locus_type","gene_type","enrichment_intergenic_dbsnp"
    df = pd.read_csv(input_file)
    print(f"Input file {input_file} read successfully. Number of rows: {len(df)}")
    # Create a new column 'has_polIII_motif' based on the presence of any polIII motif
    df['has_polIII_motif'] = df[['contains_type1_promoter', 'contains_type2_promoter', 'contains_type3_promoter']].any(axis=1)
    # Create a new column 'high_SNP_enrichment' based on the enrichment_intergenic_dbsnp value
    # get the top 75% of enrichment_intergenic_dbsnp values
    threshold = df['enrichment_intergenic_dbsnp'].quantile(0.75)
    df['high_SNP_enrichment'] = df['enrichment_intergenic_dbsnp'] > threshold
    return df

def perform_fishers_exact_test(contingency_table):
    # Perform Fisher's exact test
    odds_ratio, p_value = fisher_exact(contingency_table)
    print(f"Fisher's exact test results: Odds Ratio = {odds_ratio}, P-value = {p_value}")
    return odds_ratio, p_value

def main():
    parser = argparse.ArgumentParser(description="Perform Fisher's exact test on a 2x2 contingency table.")
    parser.add_argument('-i', '--input_file', help='Path to the input CSV file.')
    parser.add_argument('-o', '--output_file', default='results/fishers_exact_test_results.csv', help='Path to the output CSV file for results.')
    args = parser.parse_args()

    # Process the input file and create the contingency table
    df = process_and_count_input_file(args.input_file)
    contingency_table = setup_2x2_table(df)

    # Perform Fisher's exact test
    odds_ratio, p_value = perform_fishers_exact_test(contingency_table)

    # save the results to a CSV file
    results_df = pd.DataFrame({
        'Odds Ratio': [odds_ratio],
        'P-value': [p_value]
    })

    # save the contigency table to a CSV file
    contingency_table.to_csv('results/contingency_table.csv', index_label='polIII_motif_status')
    results_df.to_csv(args.output_file, index=False)
    print(f"Results saved to {args.output_file} and contingency table saved to results/contingency_table.csv")

if __name__ == "__main__":
    main()
    
