# Load the required packages
library(dplyr)
library(readr)

# Define the gene groups
gene_groups <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6", "RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12", "VTRNA",
                 "RNY","TRNA","RN7SL","RNU7","RN7SK")

expression_types <- c("ENCODE", "GTEX")

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

#Define list of 2nd exceptions: have P but are functional genes
exception_genes_2 <- c('MT-TP', 'NMTRP-TGG1-1', 'TRP-AGG1-1', 'TRP-AGG2-1', 'TRP-AGG2-2', 'TRP-AGG2-3', 
                       'TRP-AGG2-4', 'TRP-AGG2-5', 'TRP-AGG2-6', 'TRP-AGG2-7', 'TRP-AGG2-8', 'TRP-AGG3-1', 
                       'TRP-AGG5-1', 'TRP-CGG1-1', 'TRP-CGG1-2', 'TRP-CGG1-3', 'TRP-CGG2-1', 'TRP-GGG1-1',
                       'TRP-TGG1-1', 'TRP-TGG2-1', 'TRP-TGG3-1', 'TRP-TGG3-2', 'TRP-TGG3-3', 'TRP-TGG3-4', 
                       'TRP-TGG3-5', 'TRSUP-CTA1-1', 'TRSUP-CTA2-1', 'TRSUP-CTA3-1', 'TRSUP-TTA1-1', 'TRSUP-TTA2-1',
                       'TRSUP-TTA3-1')

# Create an empty data frame to store the test results
test_results <- data.frame(Expression_Type = character(),
                           Gene_group = character(),
                           Statistic = numeric(),
                           p_value = numeric(),
                           Functional_Sample_Count = integer(),
                           Pseudogene_Sample_Count = integer(),
                           stringsAsFactors = FALSE)

# Loop through each Expression type (folder)
for (expression in expression_types) {
  
  # Create empty vectors for pooled functional and pseudogene data and sample counts
  all_functional_genes <- c()
  all_pseudogenes <- c()
  
  # Loop through each gene group (e.g., RNU1, RNU2, etc.)
  for (gene in gene_groups) {
    # Construct the file path for the current expression type and gene group
    file_path <- paste0("../../data/", expression, "-expr_summary/", gene, "_expr.csv")
    
    # Check if the file exists
    if (!file.exists(file_path)) {
      cat("File does not exist for gene", gene, "in expression type", expression, "\n")
      
      # Record the result with "Insufficient data" message in Statistic and p_value columns
      test_results <- rbind(test_results, data.frame(Expression_Type = expression,
                                                     Gene_group = gene,
                                                     Statistic = NA,
                                                     p_value = NA,
                                                     Functional_Sample_Count = 0,
                                                     Pseudogene_Sample_Count = 0,
                                                     stringsAsFactors = FALSE))
      
      next  # Skip to the next gene
    }
    
    # Read the data for the current gene group
    data <- read_csv(file_path)
    
    # Create a new column 'Gene_Type' to differentiate functional genes and pseudogenes
    data <- data %>%
      mutate(Gene_Type = case_when(
        # First check for exception_genes_2, where genes with "P" are functional
        Gene %in% exception_genes_2 ~ "Functional",
        
        # Then check for exception_genes, where genes do not end with "P"
        Gene %in% exception_genes ~ "Pseudogene",
        
        # For genes that has "P" but not in exception_genes_2, classify as Pseudogene
        grepl("P", Gene) ~ "Pseudogene",
        
        # All other genes are considered Functional
        TRUE ~ "Functional"
      )) 
    
    # Separate functional genes and pseudogenes
    if(expression == 'GTEX'){
      functional_genes <- data$`Max Expression`[data$Gene_Type == "Functional"]
      pseudogenes <- data$`Max Expression`[data$Gene_Type == "Pseudogene"]
    }else{
      functional_genes <- data$`Max_FPKM`[data$Gene_Type == "Functional"]
      pseudogenes <- data$`Max_FPKM`[data$Gene_Type == "Pseudogene"]
    }
    
    # Add to the pooled data for later testing
    all_functional_genes <- c(all_functional_genes, functional_genes)
    all_pseudogenes <- c(all_pseudogenes, pseudogenes)
    
    # Count the sample sizes for functional genes and pseudogenes
    functional_n <- length(functional_genes)
    pseudogene_n <- length(pseudogenes)
    
    # Always run the KS test for each gene group, regardless of the sample sizes
    ks_result <- tryCatch({
      ks.test(functional_genes, pseudogenes)
    }, error = function(e) {
      # If there's an error (e.g., because of too few data points), return NULL
      NULL
    })
    
    # If KS test was successful, store the results, otherwise, store NA
    if (!is.null(ks_result)) {
      test_results <- rbind(test_results, data.frame(Expression_Type = expression,
                                                     Gene_group = gene,
                                                     Statistic = ks_result$statistic,
                                                     p_value = ks_result$p.value,
                                                     Functional_Sample_Count = functional_n,
                                                     Pseudogene_Sample_Count = pseudogene_n,
                                                     stringsAsFactors = FALSE))
    } else {
      test_results <- rbind(test_results, data.frame(Expression_Type = expression,
                                                     Gene_group = gene,
                                                     Statistic = NA,
                                                     p_value = NA,
                                                     Functional_Sample_Count = functional_n,
                                                     Pseudogene_Sample_Count = pseudogene_n,
                                                     stringsAsFactors = FALSE))
    }
    
    # Print a message indicating that the gene has been processed
    cat("Processed:", gene, "\n")
  }
  
  # Perform the KS test on the pooled functional and pseudogene data
  pooled_functional_n <- length(all_functional_genes)
  pooled_pseudogene_n <- length(all_pseudogenes)
  
  # Always run the KS test on pooled data, regardless of sample size
  pooled_ks_result <- tryCatch({
    ks.test(all_functional_genes, all_pseudogenes)
  }, error = function(e) {
    NULL
  })
  
  # If KS test on pooled data was successful, store the results, otherwise, store NA
  if (!is.null(pooled_ks_result)) {
    test_results <- rbind(test_results, data.frame(Expression_Type = expression,
                                                   Gene_group = "Pooled",
                                                   Statistic = pooled_ks_result$statistic,
                                                   p_value = pooled_ks_result$p.value,
                                                   Functional_Sample_Count = pooled_functional_n,
                                                   Pseudogene_Sample_Count = pooled_pseudogene_n,
                                                   stringsAsFactors = FALSE))
  } else {
    test_results <- rbind(test_results, data.frame(Expression_Type = expression,
                                                   Gene_group = "Pooled",
                                                   Statistic = NA,
                                                   p_value = NA,
                                                   Functional_Sample_Count = pooled_functional_n,
                                                   Pseudogene_Sample_Count = pooled_pseudogene_n,
                                                   stringsAsFactors = FALSE))
  }
}

# Write the results into a CSV file
write.csv(test_results, "../../results/ks_test_expression_results.csv", row.names = FALSE)

# Test random data Kolmogorov-Smirnov test
x1 <- rnorm(100)
ks.test(x1, "pnorm")
