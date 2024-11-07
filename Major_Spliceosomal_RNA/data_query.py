import requests
import biomart
import time


# Function to search HGNC database for gene symbols starting with a query term
def search_hgnc_genes(query):
    url = f'https://rest.genenames.org/search/symbol:{query}*'  # Using wildcard to match any genes starting with the query
    headers = {'Accept': 'application/json'}
    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        data = response.json()
        if data['response']['numFound'] > 0:
            # Extract gene symbols from the search result
            return [gene['symbol'] for gene in data['response']['docs']]
        else:
            print(f"No genes found for query: {query}")
            return []
    else:
        print(f"Error searching HGNC for {query}. HTTP Status: {response.status_code}")
        return []

# Function to get Ensembl transcript IDs for a given gene symbol
def fetch_ensembl_transcript_ids(gene_symbol):
    server = biomart.BiomartServer('http://www.ensembl.org/biomart')
    dataset = server.datasets['hsapiens_gene_ensembl']
    attributes = ['ensembl_transcript_id', 'external_gene_name']
    response = dataset.search({
        'attributes': attributes,
        'filters': {
            'hgnc_symbol': gene_symbol
        }
    })
    
    transcript_ids = []
    try:
        data = response.content.decode('utf-8')
        for line in data.split('\n'):
            if line:
                ensembl_transcript_id, external_gene_name = line.split('\t')
                transcript_ids.append(ensembl_transcript_id)
    except AttributeError:
        print("Error decoding response or no data found.")
        
    return transcript_ids

# Function to get genomic coordinates from Ensembl REST API
def fetch_ensembl_genome_location(transcript_id):
    url = f"https://rest.ensembl.org/lookup/id/{transcript_id}?expand=1"
    headers = {'Content-Type': 'application/json'}

    response = requests.get(url, headers=headers)
    if response.ok:
        data = response.json()
        if 'seq_region_name' in data and 'start' in data and 'end' in data:
            chrom = data['seq_region_name']
            start = data['start']
            end = data['end']
            return chrom, start, end
        else:
            print(f"No genomic location found in Ensembl for Transcript ID {transcript_id}")
    else:
        print(f"Error fetching from Ensembl: {response.status_code}")
    return None

# Function to generate UCSC Genome Browser link
def generate_ucsc_link(chrom, start, end):
    return f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={chrom}%3A{start}-{end}"

# Main function to get gene locations
def get_gene_locations(query):
    # Fetch genes with the given query term
    gene_symbols = search_hgnc_genes(query)
    if not gene_symbols:
        return
    

    # Loop through each gene symbol
    for gene_symbol in gene_symbols:
        # # if the number in the gene symbol is more than 99, we will process it otherwise we will skip it
        # substr_gene_symbol = gene_symbol[5:-1]
        # if (substr_gene_symbol == ""):
        #     continue
        # if int(substr_gene_symbol) < 100:
        #     continue        
        
        transcript_ids = fetch_ensembl_transcript_ids(gene_symbol)
        if not transcript_ids:
            print(f"No Ensembl Transcript IDs found for {gene_symbol}")
            continue

        for transcript_id in transcript_ids:
            print(f"Processing {gene_symbol} with Transcript ID: {transcript_id}")
            
            # Fetch location using Ensembl API
            location_data = fetch_ensembl_genome_location(transcript_id)
            
            if location_data:
                chrom, start, end = location_data
                formatted_start = '{:,}'.format(start)
                formatted_end = '{:,}'.format(end)
                location_str = f"Genomic Sequence ({chrom}:{formatted_start}-{formatted_end})"
                
                # Generate UCSC Genome Browser link
                ucsc_link = generate_ucsc_link(chrom, start, end)
                
                print(f"{gene_symbol} ({transcript_id}): {location_str}")
                print(f"UCSC Genome Browser link: {ucsc_link}")
            else:
                print(f"No location found for {gene_symbol} ({transcript_id})")
            time.sleep(1)  # To respect Ensembl API rate limiting
            print("------")

# Example usage
# get_gene_locations('RNU2')
# get_gene_locations('RNU1-')
get_gene_locations('RNU4-')

