# script that takes a tsv file that contains a list of HGNC pseudogenes and categorizes them into various types of ncRNA groups:
# long non-coding RNA (lncRNA), microRNA (miRNA), small nucleolar RNA (snoRNA), small nuclear RNA (snRNA), and other types of ncRNAs. The script reads the input file, processes the data, and outputs a categorized list of pseudogenes.

import argparse
import pandas as pd
import csv
import os

# download link: blob:https://www.genenames.org/1e51d098-b9c5-4d50-a070-112a7b8fa8dd
def categorize_pseudogenes(symbol, name):
    # Ribosomal: 5S ribosomal, 18S ribosomal, 28S ribosomal, 5.8S ribosomal, MT-RNR2
    # coding RNA: protein, domain, ribonuclease, BICD, Rhox homeobox family, septin 14, double homeobox A, vimentin, CXXC,
    # RPS27A, ATP binding cassette, prolyl-tRNA synthetase, apolipoprotein, RNA exonuclease, 
    # RNA polymerase, CD99, nuclear RNA export factor, RNA binding, RNA guanylyltransferase, LSM6, POLR1H,
    # ribosomal RNA processing, elongation factor,small nuclear RNA activating complex polypeptide,FAU, 
    # ZRSR2, WAP, ATE1, CDIP, MRM3, RBMX, TARBP2, REXO1, mediator complex subunit, hydrolase, leucine,
    # testis expressed, RRN3, UBTF, SYF2, RBMX2, LSM3, SSU72, cutA, LSM1, SRP14, open reading frame, PRKRA,
    # helicase, BMS1
    # Y RNA: RNY
    # small nucleolar RNA: small Cajal body-specific RNA, small nucleolar RNA, C/D box
    # vault RNA
    # tRNA
    # small nuclear RNA: U1, U2, U4, U5, U6, U7, U11, U12, U4atac, U6atac
    # long non-coding RNA: long independently transcribed non-coding RNA, 7SL, brain cytoplasmic, RN7SK
    # miscelaneous RNA: RNA component of mitochondrial RNA processing endoribonuclease

    # first, we check if the name and symbol contains any of the key words and sentence fragments that are associated with the coding RNA
    if any(keyword in name or keyword in symbol for keyword in ['coding', 'protein', 'domain', 'ribonuclease', 'BICD', 'Rhox', 'septin', 'double homeobox',
     'vimentin', 'CXXC', 'RPS27A', 'ATP binding cassette', 'prolyl-tRNA synthetase', 'apolipoprotein', 'RNA exonuclease', 'RNA polymerase', 'CD99', 
     'nuclear RNA export factor', 'RNA binding', 'RNA guanylyltransferase', 'LSM6', 'POLR1H', 'ribosomal RNA processing', 'elongation factor', 
     'small nuclear RNA activating complex polypeptide', 'FAU', 'ZRSR2', 'WAP', 'ATE1', 'CDIP', 'MRM3', 'RBMX', 'TARBP2', 'REXO1', 'mediator complex subunit', 
     'hydrolase', 'leucine rich repeat', 'testis expressed', 'RRN3', 'UBTF', 'SYF2', 'RBMX2', 'LSM3', 'SSU72', 'cutA', 'LSM1', 'SRP14', 'open reading frame', 
     'PRKRA', 'helicase', 'BMS1', 'leucine']):
        return 'coding RNA'

    # then we check if the name and symbol contains any of the key words and sentence fragments that are associated with the ribosomal RNA
    elif any(keyword in name or keyword in symbol for keyword in ['ribosomal', 'MT-RNR2']):
        return 'ribosomal RNA'

    # then we check if the name and symbol contains any of the key words and sentence fragments that are associated with the Y RNA
    elif any(keyword in name or keyword in symbol for keyword in ['RNY']):
        return 'Y RNA'

    # then we check if the name and symbol contains any of the key words and sentence fragments that are associated with the small nucleolar RNA
    elif any(keyword in name or keyword in symbol for keyword in ['small Cajal body-specific RNA', 'small nucleolar RNA', 'C/D box']):
        return 'small nucleolar RNA'

    # then we check if the name and symbol contains any of the key words and sentence fragments that are associated with the vault RNA
    elif any(keyword in name or keyword in symbol for keyword in ['vault']):
        return 'vault RNA'

    # then we check if the name and symbol contains any of the key words and sentence fragments that are associated with the tRNA
    elif any(keyword in name or keyword in symbol for keyword in ['tRNA']):
        return 'tRNA'

    # then we check if the name and symbol contains any of the key words and sentence fragments that are associated with the small nuclear RNA
    elif any(keyword in name or keyword in symbol for keyword in ['U1', 'U2', 'U4', 'U5', 'U6', 'U7', 'U11', 'U12', 'U4atac', 'U6atac']):
        return 'small nuclear RNA'

    # then we check if the name and symbol contains any of the key words and sentence fragments that are associated with the long non-coding RNA
    elif any(keyword in name or keyword in symbol for keyword in ['long independently transcribed non-coding RNA', '7SL', 'brain cytoplasmic', 'RN7SK']):
        return 'long non-coding RNA'

    # then we check if the name and symbol contains any of the key words and sentence fragments that are associated with the miscelaneous RNA
    elif any(keyword in name or keyword in symbol for keyword in ['RNA component of mitochondrial RNA processing endoribonuclease']):
        return 'miscelaneous RNA'

def separate_ribosomal_rna(symbol, name):
    # separate into the subgroups:
    # Ribosomal: 5S ribosomal, 18S ribosomal, 28S ribosomal, 5.8S ribosomal, MT-RNR2
    if any(keyword in name or keyword in symbol for keyword in ['5S ribosomal']):
        return '5S ribosomal RNA'
    elif any(keyword in name or keyword in symbol for keyword in ['18S ribosomal']):
        return '18S ribosomal RNA'
    elif any(keyword in name or keyword in symbol for keyword in ['28S ribosomal']):
        return '28S ribosomal RNA'
    elif any(keyword in name or keyword in symbol for keyword in ['5.8S ribosomal']):
        return '5.8S ribosomal RNA'
    elif any(keyword in name or keyword in symbol for keyword in ['MT-RNR2']):
        return 'MT-RNR2 ribosomal RNA'
    else:
        return 'other ribosomal RNA'

def separate_small_nucleolar_rna(symbol, name):
    if any(keyword in name or keyword in symbol for keyword in ['small Cajal body-specific RNA']):
        return 'small Cajal body-specific RNA'
    if any(keyword in name or keyword in symbol for keyword in ['C/D box']):
        return 'C/D box small nucleolar RNA'
    if any(keyword in name or keyword in symbol for keyword in ['H/ACA box']):
        return 'H/ACA box small nucleolar RNA'
    else:
        return 'other small nucleolar RNA'

def separate_small_nuclear_rna(symbol, name):
    if any(keyword in name or keyword in symbol for keyword in ['U1']):
        return 'U1 small nuclear RNA'
    if any(keyword in name or keyword in symbol for keyword in ['U2']):
        return 'U2 small nuclear RNA'
    if any(keyword in name or keyword in symbol for keyword in ['U4']):
        return 'U4 small nuclear RNA'
    if any(keyword in name or keyword in symbol for keyword in ['U5']):
        return 'U5 small nuclear RNA'
    if any(keyword in name or keyword in symbol for keyword in ['U6']):
        return 'U6 small nuclear RNA'
    if any(keyword in name or keyword in symbol for keyword in ['U7']):
        return 'U7 small nuclear RNA'
    if any(keyword in name or keyword in symbol for keyword in ['U11']):
        return 'U11 small nuclear RNA'
    if any(keyword in name or keyword in symbol for keyword in ['U12']):
        return 'U12 small nuclear RNA'
    if any(keyword in name or keyword in symbol for keyword in ['U4atac']):
        return 'U4atac small nuclear RNA'
    if any(keyword in name or keyword in symbol for keyword in ['U6atac']):
        return 'U6atac small nuclear RNA'
    else:
        return 'other small nuclear RNA'

def separate_long_non_coding_rna(symbol, name):
    if any(keyword in name or keyword in symbol for keyword in ['long independently transcribed non-coding RNA']):
        return 'long independently transcribed non-coding RNA'
    if any(keyword in name or keyword in symbol for keyword in ['7SL']):
        return '7SL long non-coding RNA'
    if any(keyword in name or keyword in symbol for keyword in ['brain cytoplasmic']):
        return 'brain cytoplasmic long non-coding RNA'
    if any(keyword in name or keyword in symbol for keyword in ['RN7SK']):
        return '7SK long non-coding RNA'
    else:
        return 'other long non-coding RNA'

def summarize_and_separate_groups(output_file):
    print(f"Summarizing and separating ncRNA subtypes from: {output_file}")
    summary = {}
    # separate into the subgroups: and write it to a the output file
    with open(output_file, 'r') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        fieldnames = reader.fieldnames + ['ncRNA_subtype']
        rows = []
        for row in reader:
            symbol = row['Symbol']
            name = row['Name']
            ncRNA_type = row['ncRNA_type']
            if ncRNA_type == 'ribosomal RNA':
                ncRNA_subtype = separate_ribosomal_rna(symbol, name)
            elif ncRNA_type == 'small nucleolar RNA':
                ncRNA_subtype = separate_small_nucleolar_rna(symbol, name)
            elif ncRNA_type == 'small nuclear RNA':
                ncRNA_subtype = separate_small_nuclear_rna(symbol, name)
            elif ncRNA_type == 'long non-coding RNA':
                ncRNA_subtype = separate_long_non_coding_rna(symbol, name)
            else:
                ncRNA_subtype = ncRNA_type  # for other types, we just use the main type as the subtype
            row['ncRNA_subtype'] = ncRNA_subtype
            rows.append(row)

            # summarize the counts of each ncRNA type and subtype
            if ncRNA_type not in summary:
                summary[ncRNA_type] = {}
            if ncRNA_subtype not in summary[ncRNA_type]:
                summary[ncRNA_type][ncRNA_subtype] = 0
            summary[ncRNA_type][ncRNA_subtype] += 1

    print(f"Writing categorized pseudogenes with subtypes to: {output_file}")
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    # write the summary to a new file
    with open('data/HGNC_pseudogenes_ncRNA_summary.tsv', 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(['ncRNA_type', 'ncRNA_subtype', 'count'])
        for ncRNA_type, subtypes in summary.items():
            print(f"Summary for {ncRNA_type}, subtypes:")
            for ncRNA_subtype, count in subtypes.items():
                print(f"  {ncRNA_subtype}: {count}")
                writer.writerow([ncRNA_type, ncRNA_subtype, count])

    print(f"Summary of ncRNA types and subtypes saved to: data/HGNC_pseudogenes_ncRNA_summary.tsv")

def main():
    args = argparse.ArgumentParser(description='Categorize HGNC pseudogenes into various types of ncRNA groups.')
    args.add_argument('input_file', type=str, help='Path to the input TSV file containing HGNC pseudogenes.')
    args.add_argument('-o', '--output_file', type=str, default='data/categorized_pseudogenes.tsv', help='Path to the output TSV file to save categorized pseudogenes. default: data/categorized_pseudogenes.tsv')

    args = args.parse_args()

    # check if the input file exists and create the output directory if it doesn't exist
    input_file = args.input_file
    output_file = args.output_file
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' does not exist.")
    if not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))
    
    print(f"Processing input file: {input_file}")
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        fieldnames = reader.fieldnames + ['ncRNA_type']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for row in reader:
            symbol = row['Symbol']
            name = row['Name']
            ncRNA_type = categorize_pseudogenes(symbol, name)
            row['ncRNA_type'] = ncRNA_type
            if ncRNA_type == 'coding RNA':
                print(f"Categorized {symbol} ({name}) as {ncRNA_type}, skipping writing to output file.")
            else:
                writer.writerow(row)

    print(f"Categorized pseudogenes saved to: {output_file}")

    summarize_and_separate_groups(output_file)


if __name__ == '__main__':
    main()