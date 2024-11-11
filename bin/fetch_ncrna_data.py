import requests
import biomart
import time

# Function to search HGNC database for gene symbols starting with a query term.
# This function performs an HTTP request to the HGNC database and retrieves gene symbols
# that start with the query term.
def search_hgnc_genes(query):
    # Check if query is valid (non-empty string)
    if not query or not isinstance(query, str):
        print("Invalid query: must be a non-empty string.")
        return []  # Return an empty list if the query is invalid

    # URL for HGNC REST API search with a wildcard on gene symbol
    url = f'https://rest.genenames.org/search/symbol:{query}*'  
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
def get_gene_locations(query):
    # Validate input query
    if not query or not isinstance(query, str):
        print("Invalid query: must be a non-empty string.")
        return
    
    # Fetch gene symbols matching the query from HGNC database
    gene_symbols = search_hgnc_genes(query)
    if not gene_symbols:
        return  # Exit if no gene symbols were found

    # Loop through each gene symbol to process
    for gene_symbol in gene_symbols:
        # Extract the numeric part of the gene symbol for additional filtering
        substr_gene_symbol = gene_symbol[5:-1]
        
        # Skip symbols with no numeric part or if the numeric part is less than the specified amount
        if substr_gene_symbol == "" or int(substr_gene_symbol) < 1137:
            continue
        
        # Fetch Ensembl transcript IDs for the gene symbol
        transcript_ids = fetch_ensembl_transcript_ids(gene_symbol)
        if not transcript_ids:
            print(f"No Ensembl Transcript IDs found for {gene_symbol}")
            continue  # Skip if no transcript IDs found

        # Loop through each transcript ID to fetch genomic location
        for transcript_id in transcript_ids:
            print(f"Processing {gene_symbol} with Transcript ID: {transcript_id}")

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
                
                # Print the gene information and UCSC link
                print(f"{gene_symbol} ({transcript_id}): {location_str}")
                print(f"UCSC Genome Browser link: {ucsc_link}")
            else:
                print(f"No location found for {gene_symbol} ({transcript_id})")
            
            # Sleep briefly to respect Ensembl's API rate limiting
            time.sleep(1/10)
            print("------")


# Example usage
# get_gene_locations('RNU2')  # Uncomment to run with specific query
# get_gene_locations('RNU1-')
# get_gene_locations('RNU4-')
# get_gene_locations('RNU5')
get_gene_locations('RNU6')
