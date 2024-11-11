import re
import pyBigWig
import csv

# Define the file path to your bigWig file
bw_file = "data/hg38.phastCons30way.bw"
input_file = "data/RNU5_data.txt"  # Update with the path to your file containing the gene data
output_file = "data/RNU5_cons.csv"  # Define the output file for CSV format

# Open the bigWig file
bw = pyBigWig.open(bw_file)

# Regular expression to extract chromosome, start, and end information
pattern = re.compile(r"Genomic Sequence \((\d+):(\d+),(\d+)\)")

# Open the input file and the output CSV file
with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    writer = csv.writer(outfile)
    raw_data = infile.read()
    
    # Split the raw data based on '------' to get each gene's information
    genes_data = raw_data.split("------")

    # Regular expression to extract chromosome, start, and end positions
    pattern = r"(\d+):([\d,]+)-([\d,]+)"
    
    # Write the CSV header
    writer.writerow(['Gene', 'Chromosome', 'Position', 'Score'])
    
    current_gene = None
    for gene_data in genes_data:
        # Extract the gene symbol
        gene_match = re.search(r"Processing (.+) with Transcript ID", gene_data)
        if gene_match:
            current_gene = gene_match.group(1)
        
        # Extract the location string
        location_match = re.search(pattern, gene_data)
        if location_match:
            chrom = f"chr{location_match.group(1)}"
            start = int(location_match.group(2).replace(',', ''))
            end = int(location_match.group(3).replace(',', ''))
            
            # Get individual scores for each position in the range
            for i in range(start, end + 1):
                scores = bw.values(chrom, i, i + 1)
                if scores:  # Ensure that scores are returned for the position
                    score = scores[0]  # Extract the score for the current position
                    # Write the gene information to the CSV file for each position and score
                    writer.writerow([current_gene, chrom, i, score])
            

# Close the bigWig file
bw.close()

print(f"Data successfully saved to {output_file}")
