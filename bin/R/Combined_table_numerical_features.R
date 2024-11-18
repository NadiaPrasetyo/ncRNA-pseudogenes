# Set the directory where your files are located
data_dir <- "../../data/"

# List of gene groups
gene_groups <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6", "RNU4ATAC", 
                 "RNU6ATAC", "RNU11", "RNU12", "VTRNA", "RNY", "TRNA", 
                 "RN7SL", "RNU7", "RN7SK")

# Define the list of exceptions (genes that are pseudogenes but do not end with "P",from TRNA gene groups)
exception_genes <- c('TRA-AGC23-1', 'TRA-TGC9-1', 'TRC-ACA1-1', 'TRC-GCA25-1', 'TRE-CTC7-1', 
                     'TRE-CTC16-1', 'TRE-TTC6-1', 'TRE-TTC7-1', 'TRE-TTC8-1', 'TRE-TTC9-1',
                     'TRE-TTC10-1', 'TRE-TTC11-1', 'TRE-TTC12-1', 'TRE-TTC13-1', 'TRF-GAA7-1',
                     'TRF-GAA8-1', 'TRF-GAA9-1', 'TRF-GAA10-1', 'TRF-GAA11-1', 'TRF-GAA12-1',
                     'TRG-CCC8-1', 'TRG-TCC5-1', 'TRG-TCC6-1', 'TRI-AAT10-1', 'TRI-AAT11-1', 
                     'TRK-CTT12-1', 'TRK-CTT13-1', 'TRK-CTT14-1', 'TRK-CTT15-1', 'TRK-CTT16-1',
                     'TRK-TTT10-1', 'TRK-TTT12-1', 'TRK-TTT13-1', 'TRK-TTT15-1', 'TRL-AAG6-1', 
                     'TRL-AAG7-1', 'TRL-TAA5-1', 'TRL-TAG4-1', 'TRN-GTT16-5', 'TRN-GTT21-1',
                     'TRN-GTT22-1', 'TRN-GTT23-1', 'TRQ-CTG9-1', 'TRQ-CTG11-1', 'TRQ-CTG13-1', 
                     'TRQ-CTG16-1', 'TRQ-CTG17-1', 'TRQ-TTG5-1', 'TRQ-TTG6-1', 'TRQ-TTG10-1',
                     'TRR-CCT6-1', 'TRR-CCT7-1', 'TRR-CCT9-1', 'TRS-ACT1-1', 'TRUND-NNN3-1',
                     'TRUND-NNN4-1', 'TRUND-NNN6-1', 'TRUND-NNN7-1', 'TRUND-NNN8-1', 'TRUND-NNN9-1', 
                     'TRUND-NNN10-1', 'TRV-CAC11-1', 'TRV-CAC11-2', 'TRX-CAT3-1', 'TRY-GTA11-1',
                     'TRY-GTA12-1')

# Initialize an empty data frame to store the combined data
combined_data <- data.frame(Gene_group = character(),
                            PhastCons30_median = numeric(),
                            PhyloP100_median = numeric(),
                            PhyloP447_median = numeric(),
                            GTEX_max = numeric(),
                            Gene_type = character())

# Loop through each file, and add the median conservation and max expression into the combined_data frame
for (gene in gene_groups) {
  # Construct file names based on the gene group
  phastCons_file <- file.path(data_dir, paste0("phastCons30_summary/", gene, "_phastCons30_summary_metrics.csv"))
  phyloP100_file <- file.path(data_dir, paste0("phyloP100_summary/", gene, "_phyloP100_summary_metrics.csv"))
  phyloP447_file <- file.path(data_dir, paste0("phyloP447_summary/", gene, "_phyloP447_summary_metrics.csv"))
  expr_file <- file.path(data_dir, paste0(gene, "_expr.csv"))
  
  # Read the files
  phastCons_data <- read.csv(phastCons_file)
  phyloP100_data <- read.csv(phyloP100_file)
  phyloP447_data <- read.csv(phyloP447_file)
  expr_data <- read.csv(expr_file)
  
  # Filter or select only the median conservation columns
  # For phastCons, phyloP100, and phyloP447, we'll assume columns are named "median" for median conservation
  phastCons_data <- phastCons_data[, c("Gene", "Median_Conservation")]  
  phyloP100_data <- phyloP100_data[, c("Gene", "Median_Conservation")]  
  phyloP447_data <- phyloP447_data[, c("Gene", "Median_Conservation")] 
  expr_data <- expr_data[, c("Gene", "Max.Expression")]  
  
  # Combine all data frames for this gene group
  # Initialize gene_combined_data with phastCons_data
  gene_combined_data <- phastCons_data
  gene_combined_data <- merge(gene_combined_data, phyloP100_data, by = "Gene", all = TRUE)
  gene_combined_data <- merge(gene_combined_data, phyloP447_data, by = "Gene", all = TRUE)
  gene_combined_data <- merge(gene_combined_data, expr_data, by = "Gene", all = TRUE)
  
  # Rename columns to match the initial ones
  colnames(gene_combined_data) <- c("Gene", "PhastCons30_median", "PhyloP100_median", "PhyloP447_median", "GTEX_max")
  
  # Create a new column 'Gene_Type' to differentiate functional genes and pseudogenes
  gene_combined_data <- gene_combined_data %>%
    mutate(Gene_Type = ifelse(grepl("P", Gene) & !Gene %in% exception_genes, "Pseudogene", "Functional"))
  
  # Append the combined gene data to the main combined data frame
  gene_combined_data$Gene_group <- gene  # Add the gene group name
  combined_data <- rbind(combined_data, gene_combined_data)
}

# Write the combined data to a CSV file
output_file <- "../../results/combined_gene_data.csv"  # Specify the file name and path
write.csv(combined_data, file = output_file, row.names = FALSE)

# View a message confirming the file has been written
cat("CSV file has been written to:", output_file, "\n")
