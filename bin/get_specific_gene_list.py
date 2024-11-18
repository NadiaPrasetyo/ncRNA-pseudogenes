import requests
import time

# Define gene group
gene_group = 'TRNA'

# Example function to search HGNC genes
def search_hgnc_pseudogenes(query):
    if not query or not isinstance(query, str):
        print("Invalid query: must be a non-empty string.")
        return []
    # URL with pseudogene query
    url = f'https://rest.genenames.org/search/*{query}*+AND+status:%22Approved%22+AND+locus_type:%22pseudogene%22'
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
    
    # Example function to search HGNC genes
def search_hgnc_functional_genes(query):
    if not query or not isinstance(query, str):
        print("Invalid query: must be a non-empty string.")
        return []
    # URL with functional transfer RNA query
    url = f'https://rest.genenames.org/search/*{query}*+AND+status:%22Approved%22+AND+locus_type:%22RNA%2C%20transfer%22'
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

# # Call function with the gene group
# gene_symbols = search_hgnc_pseudogenes(gene_group)
# print(f"Found {len(gene_symbols)} genes for query '{gene_group}':")
# print(gene_symbols)

#Call function with the gene group
gene_symbols = search_hgnc_functional_genes(gene_group)
print(f"Found {len(gene_symbols)} genes for query '{gene_group}':")


# Separate genes with a 'P' and those without a 'P'
genes_with_p = [gene for gene in gene_symbols if 'P' in gene]
genes_without_p = [gene for gene in gene_symbols if 'P' not in gene]

# Print results
print(f"Genes with 'P': {genes_with_p}")
print(f"Genes without 'P': {genes_without_p}")
