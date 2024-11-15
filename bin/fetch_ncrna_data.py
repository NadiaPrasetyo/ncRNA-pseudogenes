import requests
import biomart
import time

gene_group = 'TRNA'

# Function to search HGNC database for gene symbols starting with a query term.
# This function performs an HTTP request to the HGNC database and retrieves gene symbols
# that start with the query term.
def search_hgnc_genes(query):
    # Check if query is valid (non-empty string)
    if not query or not isinstance(query, str):
        print("Invalid query: must be a non-empty string.")
        return []  # Return an empty list if the query is invalid

    # URL for HGNC REST API search with a wildcard on gene symbol that is APPROVED
    url = f'https://rest.genenames.org/search/*{query}*+AND+status:%22Approved%22'
    headers = {'Accept': 'application/json'}  # Set Accept header for JSON response

    try:
        # Send GET request to the HGNC API with a timeout to avoid hanging requests
        response = requests.get(url, headers=headers, timeout=10)
        # Raise an exception if the response status is not 200 OK
        response.raise_for_status()
    except requests.exceptions.Timeout:
        print(f"Request timed out while searching HGNC for {query}.")
        return []  # Return empty list on timeout
    except requests.exceptions.RequestException as e:
        print(f"Request error while searching HGNC for {query}: {e}")
        return []  # Return empty list on request error

    try:
        # Parse JSON response
        data = response.json()
        # Check if there are any gene symbols found
        if data['response']['numFound'] > 0:
            # Extract gene symbols from response data
            print(f"Found {data['response']['numFound']} genes for query: {query}")
            return [gene['symbol'] for gene in data['response']['docs']]
        else:
            print(f"No genes found for query: {query}")
            return []  # Return empty list if no genes found
    except (ValueError, KeyError) as e:
        # Handle cases where JSON parsing fails or the expected keys are missing
        print(f"Error parsing response from HGNC for {query}: {e}")
        return []  # Return empty list on error


# Function to get Ensembl transcript IDs for a given gene symbol.
# This function queries the Ensembl BioMart service to fetch transcript IDs
# associated with a particular gene symbol.
def fetch_ensembl_transcript_ids(gene_symbol):
    # Initialize connection to Ensembl BioMart server
    server = biomart.BiomartServer('http://www.ensembl.org/biomart')
    dataset = server.datasets['hsapiens_gene_ensembl']  # Dataset for human genes
    attributes = ['ensembl_transcript_id', 'external_gene_name']  # Attributes to fetch (transcript ID)

    try:
        # Query Ensembl BioMart for the specified gene symbol
        response = dataset.search({
            'attributes': attributes,
            'filters': {'hgnc_symbol': gene_symbol}
        })
    except Exception as e:
        # Handle any errors when querying Ensembl (e.g., network issues or invalid query)
        print(f"Error fetching data from Ensembl for {gene_symbol}: {e}")
        return []  # Return an empty list if error occurs

    transcript_ids = []  # List to store the transcript IDs

    try:
        # Decode the response content and split by line to extract transcript IDs
        data = response.content.decode('utf-8')
        for line in data.split('\n'):
            if line:  # Ignore empty lines
                ensembl_transcript_id, external_gene_name = line.split('\t')
                transcript_ids.append(ensembl_transcript_id)  # Append transcript ID to list
    except (AttributeError, ValueError) as e:
        # Handle errors when processing the response (e.g., malformed data)
        print(f"Error processing response for {gene_symbol}: {e}")
    
    return transcript_ids  # Return the list of transcript IDs


# Function to get genomic coordinates from Ensembl REST API for a given transcript ID.
# This function queries Ensembl's REST API to get the genomic location (chromosome, start, end)
# for a given transcript ID.
def fetch_ensembl_genome_location(transcript_id):
    # URL for Ensembl's REST API to get genomic data for a transcript ID
    url = f"https://rest.ensembl.org/lookup/id/{transcript_id}?expand=1"
    headers = {'Content-Type': 'application/json'}  # Set Content-Type header for JSON response

    try:
        # Send GET request to Ensembl API to get genomic location data
        response = requests.get(url, headers=headers, timeout=10)
        # Raise an exception if the response status is not 200 OK
        response.raise_for_status()
    except requests.exceptions.Timeout:
        print(f"Request timed out while fetching location for transcript {transcript_id}.")
        return None  # Return None if request times out
    except requests.exceptions.RequestException as e:
        # Catch all request exceptions (e.g., network errors, invalid status codes)
        print(f"Request error fetching location for transcript {transcript_id}: {e}")
        return None  # Return None if there is any error

    try:
        # Parse JSON response from Ensembl API
        data = response.json()
        # Check if the required genomic location data is present
        if 'seq_region_name' in data and 'start' in data and 'end' in data:
            chrom = data['seq_region_name']
            start = data['start']
            end = data['end']
            return chrom, start, end  # Return chromosome, start, and end positions
        else:
            print(f"No genomic location found for Transcript ID {transcript_id}")
    except (ValueError, KeyError) as e:
        # Handle cases where the response format is incorrect or missing expected fields
        print(f"Error parsing location data for Transcript ID {transcript_id}: {e}")
    
    return None  # Return None if no genomic location is found or there is an error


# Function to generate UCSC Genome Browser link for a given genomic location.
# This function takes chromosome, start, and end positions and generates a link
# to view the genomic location on the UCSC Genome Browser.
def generate_ucsc_link(chrom, start, end):
    # Format the UCSC link using the provided chromosome and genomic coordinates
    return f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={chrom}%3A{start}-{end}"


# Main function to get gene locations based on a query.
# This function integrates all the above functions to search for genes, fetch transcript IDs,
# and fetch genomic locations for each transcript.
# Modified Main function to get gene locations based on a query and write to file
# Modified Main function to get gene locations based on a query or a list of gene symbols
def get_gene_locations(query=None, gene_symbols=None, output_file='data/output.txt'):
    # Validate input query or gene_symbols list
    if query and gene_symbols:
        print("Error: Please provide either a query or a list of gene symbols, not both.")
        return
    
    if not query and not gene_symbols:
        print("Error: Please provide either a query or a list of gene symbols.")
        return

    # If gene symbols are provided directly, use them
    if gene_symbols:
        gene_symbols_to_process = gene_symbols
    else:
        # Otherwise, fetch gene symbols from HGNC based on the query
        gene_symbols_to_process = search_hgnc_genes(query)
        if not gene_symbols_to_process:
            return  # Exit if no gene symbols were found

    # Open the output file in write mode
    with open(output_file, 'w') as file:
        #count = 0  # To control which genes to process
        
        # Loop through each gene symbol to process
        for gene_symbol in gene_symbols_to_process:
            # count += 1
            
            # # Skip genes before a specific count (if required)
            # if count < 116:
            #     continue
            
            # #Extract the numeric part of the gene symbol for additional filtering
            # substr_gene_symbol = gene_symbol[5:-1]
                    
            # #Skip symbols with no numeric part or if the numeric part is less than the specified amount
            # if substr_gene_symbol == "" or int(substr_gene_symbol) < 5:
            #     continue
        
            # Fetch Ensembl transcript IDs for the gene symbol
            transcript_ids = fetch_ensembl_transcript_ids(gene_symbol)
            if not transcript_ids:
                print(f"No Ensembl Transcript IDs found for {gene_symbol}")
                file.write(f"No Ensembl Transcript IDs found for {gene_symbol}\n")
                continue  # Skip if no transcript IDs found

            # Loop through each transcript ID to fetch genomic location
            for transcript_id in transcript_ids:
                msg = f"Processing {gene_symbol} with Transcript ID: {transcript_id}"
                print(msg)
                file.write(msg + '\n')

                # Fetch genomic coordinates for the transcript
                location_data = fetch_ensembl_genome_location(transcript_id)
                
                # If location data is found, generate and print the UCSC link
                if location_data:
                    chrom, start, end = location_data
                    formatted_start = '{:,}'.format(start)  # Format start coordinate
                    formatted_end = '{:,}'.format(end)  # Format end coordinate
                    location_str = f"Genomic Sequence ({chrom}:{formatted_start}-{formatted_end})"
                    
                    # Generate the UCSC link for the location
                    ucsc_link = generate_ucsc_link(chrom, start, end)
                    
                    # Print and write the gene information and UCSC link
                    output_str = f"{gene_symbol} ({transcript_id}): {location_str}"
                    link_str = f"UCSC Genome Browser link: {ucsc_link}"
                    print(output_str)
                    print(link_str)
                    file.write(output_str + '\n')
                    file.write(link_str + '\n')
                else:
                    no_loc_msg = f"No location found for {gene_symbol} ({transcript_id})"
                    print(no_loc_msg)
                    file.write(no_loc_msg + '\n')
                
                # Sleep briefly to respect Ensembl's API rate limiting
                time.sleep(1/10)
                print("------")
                file.write("------\n")

# Example usage
# If you have a list of gene symbols to look up:
gene_symbols_list = ['MT-TA', 'MT-TC', 'MT-TD', 'MT-TE', 'MT-TF', 'MT-TG', 'MT-TH', 'MT-TI', 'MT-TK', 'MT-TL1', 'MT-TL2', 'MT-TM', 'MT-TN', 'MT-TP', 'MT-TQ', 'MT-TS1', 'MT-TS2', 'MT-TT', 'MT-TV', 'MT-TW', 'MT-TY', 'NMTRL-TAA1-1', 'NMTRL-TAA4-1', 'NMTRQ-TTG3-1', 'NMTRQ-TTG5-1', 'NMTRQ-TTG14-1', 'NMTRS-TGA1-1', 'TRA-AGC9-2', 'TRA-AGC12-2', 'TRA-AGC13-1', 'TRA-AGC13-3', 'TRA-AGC16-1', 'TRA-AGC17-1', 'TRA-AGC18-1', 'TRA-AGC18-2', 'TRA-AGC19-1', 'TRA-AGC20-1', 'TRA-AGC21-1', 'TRA-AGC22-1', 'TRA-TGC8-1', 'TRC-GCA24-1', 'TRD-GTC4-1', 'TRD-GTC5-1', 'TRD-GTC6-1', 'TRD-GTC7-1', 'TRD-GTC8-1', 'TRD-GTC9-1', 'TRE-CTC3-1', 'TRE-CTC5-1', 'TRE-CTC6-1', 'TRE-CTC8-1', 'TRE-CTC17-1', 'TRE-TTC5-1', 'TRE-TTC8-2', 'TRE-TTC16-1', 'TRG-GCC5-1', 'TRG-GCC6-1', 'TRH-GTG2-1', 'TRH-GTG3-1', 'TRK-CTT10-1', 'TRK-CTT11-1', 'TRK-TTT11-1', 'TRK-TTT16-1', 'TRL-AAG5-1', 'TRL-AAG8-1', 'TRL-CAG3-1', 'TRN-ATT1-1', 'TRN-ATT1-2', 'TRN-GTT3-2', 'TRN-GTT11-1', 'TRN-GTT11-2', 'TRN-GTT12-1', 'TRN-GTT13-1', 'TRN-GTT14-1', 'TRN-GTT15-1', 'TRN-GTT15-2', 'TRN-GTT16-1', 'TRN-GTT16-2', 'TRN-GTT16-3', 'TRN-GTT16-4', 'TRN-GTT17-1', 'TRN-GTT18-1', 'TRN-GTT19-1', 'TRN-GTT19-2', 'TRN-GTT20-1', 'TRP-AGG3-1', 'TRQ-CTG8-1', 'TRQ-CTG8-2', 'TRQ-CTG8-3', 'TRQ-CTG10-1', 'TRQ-CTG12-1', 'TRQ-CTG14-1', 'TRQ-CTG15-1', 'TRQ-CTG18-1', 'TRR-CCT5-1', 'TRR-TCG6-1', 'TRS-AGA6-1', 'TRSUP-CTA1-1', 'TRSUP-TTA1-1', 'TRSUP-TTA2-1', 'TRT-AGT7-1', 'TRT-CGT6-1', 'TRU-TCA3-1', 'TRV-AAC7-1', 'TRV-CAC7-1', 'TRV-CAC9-1', 'TRV-CAC10-1', 'TRV-CAC12-1', 'TRW-CCA6-1', 'TRX-CAT2-1', 'TRY-GTA9-1', 'TRY-GTA10-1', 'TRE-TTC8-1', 'TRE-TTC11-1', 'TRE-TTC12-1', 'TRE-TTC13-1', 'TRF-GAA7-1', 'TRG-CCC8-1', 'TRK-CTT15-1', 'TRK-TTT12-1', 'TRL-AAG6-1', 'TRMT10BP1', 'TRMT112P8', 'TRN-GTT16-5', 'TRN-GTT21-1', 'TRQ-CTG17-1', 'TRR-CCT6-1', 'TRX-CAT3-1']
get_gene_locations(gene_symbols=gene_symbols_list, output_file=f'data/{gene_group}_partial_data.txt')

# If you want to use a query to search HGNC:
# get_gene_locations(query=gene_group, output_file=f'data/{gene_group}_data.txt')
