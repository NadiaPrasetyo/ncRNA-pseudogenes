import pyBigWig

# List of RNAseq BigWig file paths
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

# Initialize a set to store unique chromosomes
unique_chromosomes = set()

# Loop through each file to extract chromosomes
for file_path in rna_seq_files:
    try:
        with pyBigWig.open(file_path) as bw:
            # Get chromosomes from the BigWig file
            chromosomes = bw.chroms().keys()
            unique_chromosomes.update(chromosomes)
            print(f"File: {file_path}")
            print(f"Chromosomes: {chromosomes}\n")
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")

# Print all unique chromosomes across all files
print("Unique Chromosomes Across All BigWig Files:")
print(sorted(unique_chromosomes))
