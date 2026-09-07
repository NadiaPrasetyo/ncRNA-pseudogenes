import pandas
import csv
import argparse

#e.g:
#QTGAL,TRNA,chr17,82942149,83051770,109622,10962200,4005,0.03653463720785974,23514,0.00214500738902775,17.03240622607808,37979,0.34645417890569413,238094,0.02171954534673697,15.951262946567324,2244,0.02047034354417909,68649,0.006262337851891044,3.2688021675479617,45687,0.4167685318640419,2265411,0.20665660177701556,2.016720144821403,0.021371370351502803,1.7095130825474032,0.2118689081338088,1.6352289817196117,0.007101333333333334,2.882605643660217,0.2879,1.4476156021675648,gene with protein product,functional
def remove_protein_coding_genes(input_file, output_file):
    # Read the input CSV file
    df = pandas.read_csv(input_file)

    # Filter out rows where the 'locus_type' column is 'gene with protein product'
    filtered_df = df[df['locus_type'] != 'gene with protein product']
    print(f"Filtered out {len(df) - len(filtered_df)} protein-coding genes.")

    # Write the filtered DataFrame to the output CSV file
    filtered_df.to_csv(output_file, index=False, quoting=csv.QUOTE_NONNUMERIC)
    print(f"Filtered data saved to '{output_file}' successfully.")

def filter_by_gene_name(input_file, output_file, gene_name_csv):
    # Read the input CSV file
    df = pandas.read_csv(input_file)
    gene_names_df = pandas.read_csv(gene_name_csv)

    # Filter rows by gene name (case-insensitive)
    filtered_df = df[df['gene_name'].str.lower().isin(gene_names_df['Gene'].str.lower())]
    print(f"Filtered to {len(filtered_df)} rows with gene names from '{gene_name_csv}'.")

    # Write the filtered DataFrame to the output CSV file
    filtered_df.to_csv(output_file, index=False, quoting=csv.QUOTE_NONNUMERIC)
    print(f"Filtered data saved to '{output_file}' successfully.")

def main():
    parser = argparse.ArgumentParser(description='Remove protein-coding genes from a CSV file.')
    parser.add_argument('input_file', help='Path to the input CSV file')
    parser.add_argument('output_file', help='Path to the output CSV file')
    parser.add_argument('--filter-gene-name', help='Optional: Filter by gene name (case-insensitive) to only include specific genes in this csv', default=None)

    args = parser.parse_args()

    if args.filter_gene_name:
        print(f"Filtering by gene name using '{args.filter_gene_name}'...")
        filter_by_gene_name(args.input_file, args.output_file, args.filter_gene_name)


    remove_protein_coding_genes(args.input_file, args.output_file)

if __name__ == '__main__':
    main()