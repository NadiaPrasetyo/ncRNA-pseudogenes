import requests
import pyBigWig
import time

# Function to fetch chromosomal locations based on the provided BigBed data
def fetch_chromosomal_location(gene_symbol, bigbed_file):
    try:
        # Get the list of chromosomes in the BigBed file
        chromosomes = bigbed_file.chroms()

        # Search through each chromosome's entries to find the gene symbol
        for chrom in chromosomes:
            entries = bigbed_file.entries(chrom, 0, chromosomes[chrom])
            for start, end, name in entries:
                # The gene symbol is embedded within a tab-separated string in the name field
                if gene_symbol in name:
                    # Extract the gene symbol and transcript ID from the name field
                    # Assuming the gene symbol is the first part of the name (before any spaces or tabs)
                    parts = name.split()
                    gene_name = parts[0]  # Assuming gene symbol is the first element
                    
                    # Extract transcript ID from the name if available (adjust based on your format)
                    # In the case that the transcript ID is not part of the name, we return "not found"
                    transcript_id = parts[1] if len(parts) > 1 else "not found"
                    
                    return chrom, start, end, gene_name, transcript_id
        return None
    except Exception as e:
        print(f"Error fetching data from BigBed file for {gene_symbol}: {e}")
        return None

# Generate UCSC Genome Browser link
def generate_ucsc_link(chrom, start, end):
    return f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={chrom}%3A{start}-{end}"

# Main function to get gene locations and print to the console
def get_gene_locations(query, output_file, bigbed_file_path):
    if not query or not isinstance(query, str):
        print("Invalid query: must be a non-empty string.")
        return
    
    # Open the output file to write the results
    with open(output_file, 'w') as file:
        # Fetch gene symbols based on the query
        gene_symbols = search_hgnc_genes(query)
        if not gene_symbols:
            return
        
        # Open the BigBed file
        bigbed_file = pyBigWig.open(bigbed_file_path)

        # Process each gene symbol
        for gene_symbol in gene_symbols:
            # Print processing message for each gene
            print(f"Processing {gene_symbol} with Transcript ID: not found")
            file.write(f"Processing {gene_symbol} with Transcript ID: not found\n")

            # Fetch the genomic location for the gene
            location_data = fetch_chromosomal_location(gene_symbol, bigbed_file)
            if location_data:
                chrom, start, end, gene_name, transcript_id = location_data
                formatted_start = '{:,}'.format(start)
                formatted_end = '{:,}'.format(end)
                location_str = f"Genomic Sequence ({chrom}:{formatted_start}-{formatted_end})"
                
                # Generate UCSC link
                ucsc_link = generate_ucsc_link(chrom, start, end)
                
                # Print and write the output to the file
                print(f"{gene_symbol} ({gene_name}): {location_str}")
                print(f"UCSC Genome Browser link: {ucsc_link}")
                file.write(f"{gene_name} ({transcript_id}): {location_str}\n")
                file.write(f"UCSC Genome Browser link: {ucsc_link}\n")
            else:
                print(f"No location found for {gene_symbol}")
                file.write(f"No location found for {gene_symbol}\n")
            
            # Print separator between genes
            print("------")
            file.write("------\n")
            
            # Optional: brief sleep to avoid rate limiting
            time.sleep(1/10)

# Example function to search HGNC genes
def search_hgnc_genes(query):
    if not query or not isinstance(query, str):
        print("Invalid query: must be a non-empty string.")
        return []
    
    url = f'https://rest.genenames.org/search/*{query}*+AND+(locus_type:%22RNA%2C%20transfer%22+OR+locus_type:%22pseudogene%22)+AND+status:%22Approved%22'
    headers = {'Accept': 'application/json'}
    try:
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"Request error: {e}")
        return []
    
    try:
        data = response.json()
        if data['response']['numFound'] > 0:
            return [gene['symbol'] for gene in data['response']['docs']]
        else:
            print(f"No genes found for query: {query}")
            return []
    except (ValueError, KeyError) as e:
        print(f"Error parsing response: {e}")
        return []

# Define gene group and BigBed file path
gene_group = 'TRNA'
bigbed_file_path = 'data/hgnc.bb'  # Specify the path to your BigBed file

# Call function with the gene group and output file path
get_gene_locations(gene_group, f'data/{gene_group}_data_temp.txt', bigbed_file_path)
