import re
import pyBigWig
import csv
import os
import sys

# Define the file path to your bigWig file
bw_file = "data/hg38.phastCons30way.bw"
input_file = "data/RNU6_data.txt"  # Update with the path to your file containing the gene data
output_file = "data/RNU6_cons.csv"  # Define the output file for CSV format
temp_output_file = f"{output_file}.tmp"  # Temporary file for output

# Regular expression to extract chromosome, start, and end information
# allows for X and Y chromosomes in addition to numeric chromosomes
pattern = re.compile(r"(\d+|X|Y):([\d,]+)-([\d,]+)")

# Ensure input files exist
if not os.path.isfile(bw_file):
    sys.exit(f"Error: BigWig file '{bw_file}' does not exist.")
if not os.path.isfile(input_file):
    sys.exit(f"Error: Input file '{input_file}' does not exist.")

try:
    # Open the bigWig file
    bw = pyBigWig.open(bw_file)
    if not bw.isBigWig():
        sys.exit(f"Error: File '{bw_file}' is not a valid BigWig file.")

    # Open the input file and the output CSV file
    with open(input_file, 'r') as infile, open(temp_output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        raw_data = infile.read()

        # Split the raw data based on '------' to get each gene's information
        genes_data = raw_data.split("------")

        # Write the CSV header
        writer.writerow(['Gene', 'Chromosome', 'Position', 'Score'])

        current_gene = None
        for gene_data in genes_data:
            # Extract the gene symbol
            gene_match = re.search(r"Processing (.+) with Transcript ID", gene_data)
            if gene_match:
                current_gene = gene_match.group(1)
            else:
                print(f"Warning: Gene symbol not found in data segment:\n{gene_data}\n")
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
                    continue

                # Get individual scores for each position in the range
                for i in range(start, end + 1):
                    try:
                        scores = bw.values(chrom, i, i + 1)
                        if scores:  # Ensure that scores are returned for the position
                            score = scores[0]  # Extract the score for the current position
                            writer.writerow([current_gene, chrom, i, score])
                        else:
                            print(f"Warning: No score data for {chrom}:{i}")
                    except RuntimeError as e:
                        print(f"Error fetching scores for {chrom}:{i}: {e}")
                        continue
            else:
                print(f"Warning: Location pattern not found for gene {current_gene} in data segment.")
        
    # Rename the temporary output file to the final output file after successful completion
    os.rename(temp_output_file, output_file)
    print(f"Data successfully saved to {output_file}")

except Exception as e:
    print(f"An unexpected error occurred: {e}")
    if os.path.exists(temp_output_file):
        os.remove(temp_output_file)
finally:
    # Close the bigWig file if it was successfully opened
    if 'bw' in locals() and not bw.close:
        bw.close()
