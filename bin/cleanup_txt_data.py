import re
import os
import sys

# Regular expression to extract chromosome, start, and end information
# allows for X and Y chromosomes in addition to numeric chromosomes
pattern = re.compile(r"(\d+|X|Y):([\d,]+)-([\d,]+)")

# Path to the input file
input_file = "data/RNY_data.txt"  # Update this to the path of your actual file

# Ensure the input file exists
if not os.path.isfile(input_file):
    sys.exit(f"Error: Input file '{input_file}' does not exist.")

# Temporary list to store valid lines
valid_lines = []
invalid_lines = []  # List to store invalid lines for later printing

try:
    with open(input_file, 'r') as infile:
        raw_data = infile.read()

        # Split the raw data based on '------' to get each gene's information
        genes_data = raw_data.split("------")

        for gene_data in genes_data:
            # Extract the gene symbol
            gene_match = re.search(r"Processing (.+) with Transcript ID", gene_data)
            if gene_match:
                current_gene = gene_match.group(1)
            else:
                print(f"Warning: Gene symbol not found in data segment:\n{gene_data}\n")
                invalid_lines.append(gene_data)  # Store for invalid lines
                continue

            # Extract the location string
            location_match = re.search(pattern, gene_data)
            if location_match:
                chrom = f"chr{location_match.group(1)}"
                try:
                    start = int(location_match.group(2).replace(',', ''))
                    end = int(location_match.group(3).replace(',', ''))
                except ValueError as e:
                    print(f"Error parsing start/end positions for gene {current_gene}: {e}")
                    invalid_lines.append(gene_data)  # Store for invalid lines
                    continue

                # If both gene and location are valid, store the gene data
                valid_lines.append(gene_data)
            else:
                print(f"Warning: Location pattern not found for gene {current_gene} in data segment.")
                invalid_lines.append(gene_data)  # Store for invalid lines

    # If invalid lines exist, print them
    if invalid_lines:
        print("\nInvalid lines detected (these will be removed):\n")
        for invalid_line in invalid_lines:
            print(invalid_line)

    # Check if any changes have been made
    if valid_lines == genes_data:
        print("Nothing changed.")
    else:
        # Overwrite the original file with only the valid lines
        with open(input_file, 'w') as outfile:
            outfile.write("------".join(valid_lines))  # Join valid lines and write them back
        print(f"\nData cleaned successfully. Invalid lines have been removed.")

except Exception as e:
    print(f"An unexpected error occurred: {e}")
