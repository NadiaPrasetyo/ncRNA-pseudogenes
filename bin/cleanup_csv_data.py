import csv
import os
import sys

# Path to the input file (update this path as needed)
input_file = "results/combined_gene_data.csv"  # Input CSV file with the specified fields
# Ensure the input file exists
if not os.path.isfile(input_file):
    sys.exit(f"Error: Input file '{input_file}' does not exist.")

# Function to clean the data
def clean_gene_data():
    print(f"Cleaning data in {input_file}...")

    valid_rows = []  # List to store valid rows
    invalid_rows = []  # List to store invalid rows (for logging)

    try:
        # Open the CSV file for reading
        with open(input_file, 'r', newline='', encoding='utf-8') as infile:
            reader = csv.DictReader(infile)
            if set(reader.fieldnames) != {"Gene", "PhastCons30_median", "PhyloP100_median", 
                                          "PhyloP447_median", "GTEX_max", "ENCODE_max", 
                                          "Gene_Type", "Gene_group"}:
                sys.exit("Error: The CSV file does not have the expected fields.")

            for row in reader:
                gene = row["Gene"]
                gene_group = row["Gene_group"]

                # Apply the rule: remove rows where gene has no hyphen ONLY IF gene group is TRNA
                if gene_group == "TRNA" and '-' not in gene:
                    invalid_rows.append(row)  # Log invalid rows
                    continue

                # If the row is valid, add it to the list
                valid_rows.append(row)

        # Overwrite the input file with only the valid rows
        with open(input_file, 'w', newline='', encoding='utf-8') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
            writer.writeheader()
            writer.writerows(valid_rows)

        print(f"Data cleaned successfully. Invalid rows have been removed.")
        if invalid_rows:
            print(f"\nInvalid rows removed (these had no hyphen in the gene and belonged to 'TRNA'):")
            for row in invalid_rows:
                print(f"Gene: {row['Gene']}, Full Row: {row}")

    except Exception as e:
        print(f"An unexpected error occurred while cleaning data: {e}")


# Main execution logic
if __name__ == "__main__":
    clean_gene_data()
