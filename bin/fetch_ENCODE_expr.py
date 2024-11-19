import csv

# Step 1: Extract Gene Symbols from the TXT file
def extract_gene_symbols(txt_file):
    gene_symbols = []
    with open(txt_file, 'r') as file:
        for line in file:
            if line.startswith("Processing"):
                # Extract the gene symbol from the line
                gene_symbol = line.split()[1]
                gene_symbols.append(gene_symbol)
    return set(gene_symbols)  # Use set to ensure uniqueness

# Step 2: Fetch Max TPM and FPKM for each Gene Symbol
def fetch_max_expression_data(tsv_file, gene_symbols):
    result_data = []

    with open(tsv_file, 'r') as file:
        # Skip the first line of metadata
        next(file)
        
        # Now use DictReader to read the actual data
        reader = csv.DictReader(file, delimiter='\t')

        # Strip spaces from column headers if needed
        reader.fieldnames = [field.strip() for field in reader.fieldnames]
        gene_data = {gene: {'TPM': [], 'FPKM': []} for gene in gene_symbols}

        # Store TPM and FPKM values for each gene symbol
        for row in reader:
            gene_symbol = row['Gene symbol']  # Adjust after checking headers
            if gene_symbol in gene_symbols:
                gene_data[gene_symbol]['TPM'].append(float(row['TPM']))
                gene_data[gene_symbol]['FPKM'].append(float(row['FPKM']))
        
        # Find max TPM and FPKM for each gene
        for gene, values in gene_data.items():
            if values['TPM']:
                max_tpm = max(values['TPM'])
            else:
                max_tpm = None
            if values['FPKM']:
                max_fpkm = max(values['FPKM'])
            else:
                max_fpkm = None
            
            result_data.append({'Gene': gene, 'Max_TPM': max_tpm, 'Max_FPKM': max_fpkm})

    return result_data

# Step 3: Save Results to CSV
def save_to_csv(data, output_csv):
    with open(output_csv, 'w', newline='') as file:
        fieldnames = ['Gene', 'Max_TPM', 'Max_FPKM']
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        
        # Write the header
        writer.writeheader()
        
        # Write the data rows
        writer.writerows(data)

# Main Execution
if __name__ == "__main__":
    
    gene_group = "VTRNA"  # Replace with your gene group
    
    # Input file paths
    txt_file = f"data/{gene_group}_data.txt"  # Replace with the path to your TXT file
    tsv_file = "data/rna_expression_report_2024_11_19_2h_13m.tsv" 
    output_csv = f"data/ENCODE-expr_summary/{gene_group}_expr.csv"
    
    # Process the files
    gene_symbols = extract_gene_symbols(txt_file)
    max_expression_data = fetch_max_expression_data(tsv_file, gene_symbols)
    save_to_csv(max_expression_data, output_csv)
    
    print(f"Results saved to {output_csv}")
