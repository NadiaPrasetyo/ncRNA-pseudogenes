# Loading the required packages
library(dplyr)
library(readr)

# Define the gene groups and conservation types
gene_groups <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6", "RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12", "VTRNA",
                 "RNY","TRNA","RN7SL","RNU7","RN7SK")
conservation_types <- c("phastCons30", "phyloP100", "phyloP447")

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

# Create an empty data frame to store the test results with sample counts (without Gene_Type column)
test_results <- data.frame(Conservation_Type = character(),
                           Gene_group = character(),
                           Statistic = numeric(),
                           p_value = numeric(),
                           Functional_Sample_Count = integer(),
                           Pseudogene_Sample_Count = integer(),
                           stringsAsFactors = FALSE)

# Loop through each conservation type (folder)
for (conservation in conservation_types) {
  
  # Create empty vectors for pooled functional and pseudogene data and sample counts
  pooled_functional <- c()
  pooled_pseudogenes <- c()
  pooled_functional_count <- 0
  pooled_pseudogene_count <- 0
  
  # Loop through each gene group (e.g., RNU1, RNU2, etc.)
  for (gene in gene_groups) {
    
    # Construct the file path for the current conservation type and gene group
    file_path <- paste0("../../data/", conservation, "_summary/", gene, "_", conservation, "_summary_metrics.csv")
    
    # Check if the file exists
    if (!file.exists(file_path)) {
      cat("File does not exist for gene", gene, "in conservation type", conservation, "\n")
      
      # Record the result with "Insufficient data" message in Statistic and p_value columns
      test_results <- rbind(test_results, data.frame(Conservation_Type = conservation,
                                                     Gene_group = gene,
                                                     Statistic = NA,
                                                     p_value = NA,
                                                     Functional_Sample_Count = 0,
                                                     Pseudogene_Sample_Count = 0,
                                                     stringsAsFactors = FALSE))
      
      next  # Skip to the next gene
    }
    
    # Read the data for the current gene group and conservation type
    data <- read_csv(file_path)
    
    # Create a new column 'Gene_Type' to differentiate functional genes and pseudogenes (no need to add this column to results)
    data <- data %>%
      mutate(Gene_Type = ifelse(grepl("P", Gene) & !Gene %in% exception_genes, "Pseudogene", "Functional"))    
    
    # Separate functional genes and pseudogenes
    functional_genes <- data$Median_Conservation[data$Gene_Type == "Functional"]
    pseudogenes <- data$Median_Conservation[data$Gene_Type == "Pseudogene"]
    
    # Update the sample counts
    functional_count <- length(functional_genes)
    pseudogene_count <- length(pseudogenes)
    
    # Append to pooled lists
    pooled_functional <- c(pooled_functional, functional_genes)
    pooled_pseudogenes <- c(pooled_pseudogenes, pseudogenes)
    pooled_functional_count <- pooled_functional_count + functional_count
    pooled_pseudogene_count <- pooled_pseudogene_count + pseudogene_count
    
    # Perform KS test only if both groups have data
    if (functional_count > 0 && pseudogene_count > 0) {
      ks_result <- ks.test(functional_genes, pseudogenes)
      test_results <- rbind(test_results, data.frame(Conservation_Type = conservation,
                                                     Gene_group = gene,
                                                     Statistic = ks_result$statistic,
                                                     p_value = ks_result$p.value,
                                                     Functional_Sample_Count = functional_count,
                                                     Pseudogene_Sample_Count = pseudogene_count,
                                                     stringsAsFactors = FALSE))
    } else {
      # If either group is empty, handle it here (e.g., no KS test performed)
      test_results <- rbind(test_results, data.frame(Conservation_Type = conservation,
                                                     Gene_group = gene,
                                                     Statistic = NA,
                                                     p_value = NA,
                                                     Functional_Sample_Count = functional_count,
                                                     Pseudogene_Sample_Count = pseudogene_count,
                                                     stringsAsFactors = FALSE))
      cat("Skipping KS test for gene", gene, "due to insufficient data\n")
    }
    
    # Print a message indicating that the gene has been processed
    cat("Processed:", gene, "in conservation type", conservation, "\n")
  }
  
  # Perform the KS test on pooled data (only if both groups have data)
  if (pooled_functional_count > 0 && pooled_pseudogene_count > 0) {
    ks_result <- ks.test(pooled_functional, pooled_pseudogenes)
    test_results <- rbind(test_results, data.frame(Conservation_Type = conservation,
                                                   Gene_group = "Pooled",
                                                   Statistic = ks_result$statistic,
                                                   p_value = ks_result$p.value,
                                                   Functional_Sample_Count = pooled_functional_count,
                                                   Pseudogene_Sample_Count = pooled_pseudogene_count,
                                                   stringsAsFactors = FALSE))
  } else {
    # If either group is empty, handle it here for pooled data
    test_results <- rbind(test_results, data.frame(Conservation_Type = conservation,
                                                   Gene_group = "Pooled",
                                                   Statistic = NA,
                                                   p_value = NA,
                                                   Functional_Sample_Count = pooled_functional_count,
                                                   Pseudogene_Sample_Count = pooled_pseudogene_count,
                                                   stringsAsFactors = FALSE))
    cat("Skipping KS test for pooled data due to insufficient data\n")
  }
}

# Write the results into a CSV file
write.csv(test_results, "../../results/ks_test_conservation_results.csv", row.names = FALSE)

# Test random data Kolmogorov-Smirnov test (example)
x1 <- rnorm(100)
ks.test(x1, "pnorm")