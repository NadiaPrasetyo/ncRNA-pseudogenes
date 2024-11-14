# Loading the required packages
library(dplyr)
library(readr)

# Define the gene groups and conservation types
gene_groups <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6", "RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12", "VTRNA")
conservation_types <- c("phastCons30", "phyloP100", "phyloP447")

# Create an empty data frame to store the test results with sample counts
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
    
    # Create a new column 'Gene_Type' to differentiate functional genes and pseudogenes
    data <- data %>%
      mutate(Gene_Type = ifelse(grepl("P$", Gene), "Pseudogene", "Functional"))
    
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
    
    # Perform KS test (even with insufficient data, we will just record NA values)
    ks_result <- ks.test(functional_genes, pseudogenes)
    test_results <- rbind(test_results, data.frame(Conservation_Type = conservation,
                                                   Gene_group = gene,
                                                   Statistic = ks_result$statistic,
                                                   p_value = ks_result$p.value,
                                                   Functional_Sample_Count = functional_count,
                                                   Pseudogene_Sample_Count = pseudogene_count,
                                                   stringsAsFactors = FALSE))
    
    # Print a message indicating that the gene has been processed
    cat("Processed:", gene, "in conservation type", conservation, "\n")
  }
  
  # Perform the KS test on pooled data (always attempt KS, even if sample size is low)
  ks_result <- ks.test(pooled_functional, pooled_pseudogenes)
  
  # Store the pooled test result with sample counts
  test_results <- rbind(test_results, data.frame(Conservation_Type = conservation,
                                                 Gene_group = "Pooled",
                                                 Statistic = ks_result$statistic,
                                                 p_value = ks_result$p.value,
                                                 Functional_Sample_Count = pooled_functional_count,
                                                 Pseudogene_Sample_Count = pooled_pseudogene_count,
                                                 stringsAsFactors = FALSE))
}

# Write the results into a CSV file
write.csv(test_results, "../../results/ks_test_conservation_results.csv", row.names = FALSE)

# Test random data Kolmogorov-Smirnov test (example)
x1 <- rnorm(100)
ks.test(x1, "pnorm")
