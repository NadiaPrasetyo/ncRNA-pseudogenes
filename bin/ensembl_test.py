import biomart

# Connect to the Ensembl BioMart server
server = biomart.BiomartServer('http://www.ensembl.org/biomart')

# Select the dataset you want to query (e.g., Ensembl Genes for Homo sapiens)
dataset = server.datasets['hsapiens_gene_ensembl']

# Define the attributes you want to fetch: Ensembl gene ID, HGNC symbol, and gene biotype
attributes = [
    'ensembl_transcript_id',
    'external_gene_name',
    'gene_biotype'
]

# Define a list of gene symbols you want to search for
filter_attribute = 'hgnc_symbol'
filter_values = [f'RNU2-{i}' + ('P' if i > 1 else '') for i in range(1, 70)]  # Generates RNU2-1, RNU2-2P, ..., RNU2-69P

# Perform the BioMart query with multiple filter values
response = dataset.search({
    'attributes': attributes,
    'filters': {
        filter_attribute: filter_values  # Pass the list of gene symbols
    }
})

# Try-except block to handle if the data is not found in the response
try:
    data = response.content.decode('utf-8')
except AttributeError:
    data = ""

# Create a dictionary to map gene symbols to Ensembl IDs
genesymbol_to_ensembl = {}
for line in data.split('\n'):
    if line:
        ensembl_transcript_id, external_gene_name, gene_biotype = line.split('\t')
        genesymbol_to_ensembl[external_gene_name] = ensembl_transcript_id

# Print the dictionary with gene symbols mapped to Ensembl IDs
print(genesymbol_to_ensembl)

# Optionally, print the raw response for debugging
print(response)
