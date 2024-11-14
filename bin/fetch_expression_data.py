import os
import pyBigWig
import csv
import re

# Script to fetch expression data for specified genes from RNAseq data
# Need to download all data files and place in a directory named 'GTEX-RNAseq' in the 'data' folder
# RNAseq data files are in BigWig format and contain expression values for different tissues
# Gene data file should be in the 'data' folder and contain gene names and genomic locations

# List of filenames to process
rna_seq_files = [
    "data/GTEX-RNAseq/GTEX-1C475-0726-SM-73KVL.Esophagus_Muscularis.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1C475-1826-SM-73KWA.Skin_Sun_Exposed_Lower_leg.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1GMR3-0626-SM-9WYT3.Artery_Coronary.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1H1E6-0826-SM-9WG83.Pancreas.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1HFI6-0011-R7b-SM-CM2SS.Brain_Putamen_basal_ganglia.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1HGF4-0011-R5b-SM-CM2ST.Brain_Caudate_basal_ganglia.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1HSGN-0726-SM-A9G2F.Thyroid.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1HSKV-0011-R1b-SM-CMKH7.Brain_Hippocampus.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1I1GU-1226-SM-A9SKT.Esophagus_Gastroesophageal_Junction.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1IDJC-1326-SM-CL53H.Colon_Transverse.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1IDJU-1026-SM-AHZ2U.Vagina.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1JKYN-1026-SM-CGQG4.Testis.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1JN76-0626-SM-CKZOQ.Skin_Not_Sun_Exposed_Suprapubic.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1KXAM-1926-SM-D3LAG.Colon_Sigmoid.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1LG7Z-0005-SM-DKPQ6.Whole_Blood.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1MA7W-1526-SM-DHXKS.Uterus.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1PIEJ-1526-SM-E6CP8.Small_Intestine_Terminal_Ileum.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-11NSD-1126-SM-5N9BQ.Esophagus_Mucosa.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-13OVI-1126-SM-5KLZF.Kidney_Cortex.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-13S86-0326-SM-5SI6K.Heart_Atrial_Appendage.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-13X6J-0011-R11a-SM-5P9HE.Brain_Cerebellar_Hemisphere.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-14BIN-0011-R6a-SM-5S2RH.Brain_Nucleus_accumbens_basal_ganglia.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-14BMU-0626-SM-73KZ6.Adipose_Visceral_Omentum.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-14DAR-1026-SM-73KV3.Prostate.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-14PKU-0526-SM-6871A.Spleen.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-14PN4-0011-R3b-SM-686ZU.Brain_Anterior_cingulate_cortex_BA24.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-117XS-0008-SM-5Q5DQ.Cells_Cultured_fibroblasts.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-145MH-2926-SM-5Q5D2.Brain_Cerebellum.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-1122O-0003-SM-5Q5DL.Cells_EBV-transformed_lymphocytes.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-NFK9-0326-SM-3MJGV.Adipose_Subcutaneous.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-NFK9-0626-SM-2HMIV.Muscle_Skeletal.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-NFK9-0926-SM-2HMJU.Heart_Left_Ventricle.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-NFK9-1526-SM-3LK7B.Stomach.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-OHPK-2326-SM-3MJH2.Fallopian_Tube.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-S3XE-1226-SM-4AD4L.Bladder.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-S341-1126-SM-4AD6T.Cervix_Ectocervix.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-T5JC-0011-R4A-SM-32PLT.Brain_Amygdala.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-T5JC-0011-R8A-SM-32PLM.Brain_Hypothalamus.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-T5JC-0011-R10A-SM-32PM2.Brain_Frontal_Cortex_BA9.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-T5JC-1626-SM-EZ6KW.Kidney_Medulla.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-TML8-1626-SM-32QOO.Nerve_Tibial.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-UTHO-3026-SM-3GAFB.Brain_Cortex.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-WYVS-0426-SM-4ONDL.Artery_Aorta.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-XPT6-2226-SM-4B66R.Artery_Tibial.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-Y5LM-0126-SM-4VBRL.Adrenal_Gland.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-Y5LM-0426-SM-4VBRO.Liver.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-Y5LM-1826-SM-4VDT9.Minor_Salivary_Gland.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-Y5V5-0826-SM-4VBQD.Lung.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-Y111-2926-SM-4TT25.Pituitary.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-YFC4-0011-R9a-SM-4SOK4.Brain_Spinal_cord_cervical_c-1.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-Z93S-0011-R2a-SM-4RGNG.Brain_Substantia_nigra.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-ZPIC-1326-SM-DO91Y.Cervix_Endocervix.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-ZT9W-2026-SM-51MRA.Breast_Mammary_Tissue.RNAseq.bw",
    "data/GTEX-RNAseq/GTEX-ZVT2-0326-SM-5E44G.Ovary.RNAseq.bw"
]

gene_group = "RNU6ATAC" #Update on the gene group you're working on

#Alter these file paths for each specified gene groups
# Path to input gene file and output CSV file
gene_file_path = f"data/{gene_group}_data.txt"
output_csv = f"data/{gene_group}_expr.csv"

# Regular expressions for extracting gene name and genomic location
gene_pattern = re.compile(r"Processing (.+?) with Transcript ID")
location_pattern = re.compile(r"(\d+|X|Y):([\d,]+)-([\d,]+)")

# Step 1: Parse genes from the input file
genes = []
try:
    with open(gene_file_path, 'r') as gene_file:
        raw_data = gene_file.read()
        gene_data_sections = raw_data.split("------")  # Splitting sections by '------'
        
        for section in gene_data_sections:
            # Extract gene name
            gene_match = gene_pattern.search(section)
            if not gene_match:
                print(f"Warning: No gene found in section:\n{section}\n")
                continue
            gene_name = gene_match.group(1)
            
            # Extract chromosome, start, and end positions
            location_match = location_pattern.search(section)
            if not location_match:
                print(f"Warning: No location found for gene {gene_name} in section:\n{section}\n")
                continue
            chromosome = f"chr{location_match.group(1)}"
            start = int(location_match.group(2).replace(',', ''))
            end = int(location_match.group(3).replace(',', ''))
            
            # Append parsed gene details to the genes list
            genes.append({"name": gene_name, "chromosome": chromosome, "start": start, "end": end})
except FileNotFoundError:
    print(f"Gene file not found: {gene_file_path}")
except Exception as e:
    print(f"Error reading gene file: {e}")

# Initialize dictionary to store max expression values for each gene
gene_max_info = {gene["name"]: {"max_value": float('-inf'), "location": None} for gene in genes}

# Step 2: Process each RNAseq file to find max expression for each gene
for file_path in rna_seq_files:
    try:
        with pyBigWig.open(file_path) as bw:
            # Extract tissue/body location from the filename (e.g., "Esophagus_Muscularis")
            body_location = os.path.basename(file_path).split('.')[1].replace('_', ' ')
            
            for gene in genes:
                chrom = gene["chromosome"]
                start = gene["start"]
                end = gene["end"]
                
                # Fetch expression values in the specified genomic range
                values = bw.values(chrom, start, end)
                valid_values = [value for value in values if value is not None]
                
                # Update max expression if new value is higher
                if valid_values:
                    max_value = max(valid_values)
                    if max_value > gene_max_info[gene["name"]]["max_value"]:
                        gene_max_info[gene["name"]]["max_value"] = max_value
                        gene_max_info[gene["name"]]["location"] = body_location
                
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except pyBigWig.open.Error as e:
        print(f"Error opening {file_path}: {e}")
    except Exception as e:
        print(f"Error processing {file_path}: {e}")

# Step 3: Write results to CSV
try:
    with open(output_csv, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(["Gene", "Max Expression", "Location"])  # CSV header
        
        for gene, info in gene_max_info.items():
            csv_writer.writerow([gene, info["max_value"], info["location"]])
    print(f"Results written to {output_csv}")
except IOError as e:
    print(f"Error writing to CSV file {output_csv}: {e}")