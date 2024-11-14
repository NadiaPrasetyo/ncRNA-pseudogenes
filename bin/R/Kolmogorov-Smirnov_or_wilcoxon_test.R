
# loading the required package
library(dplyr)
library("dgof")
library(readr)

# Loading the required package
library(dplyr)
library("dgof")
library(readr)

# Define the gene groups and conservation types
gene_groups <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6")
conservation_types <- c("phastCons30", "phyloP100", "phyloP447")

# Create an empty list to store the test results
test_results <- data.frame(Conservation_Type = character(),
                           Gene_group = character(),
                           Test = character(),
                           Statistic = numeric(),
                           p_value = numeric(),
                           stringsAsFactors = FALSE)

# Loop through each conservation type (folder)
for (conservation in conservation_types) {
  
  # Create empty vectors for pooled functional and pseudogene data
  pooled_functional <- c()
  pooled_pseudogenes <- c()
  
  # Loop through each gene group (e.g., RNU1, RNU2, etc.)
  for (gene in gene_groups) {
    
    # Construct the file path for the current conservation type and gene group
    file_path <- paste0("../../data/", conservation, "_summary/", gene, "_", conservation, "_summary_metrics.csv")
    
    # Check if the file exists
    if (!file.exists(file_path)) {
      cat("Skipping: File does not exist for gene", gene, "in conservation type", conservation, "\n")
      next
    }
    
    # Read the data for the current gene group and conservation type
    data <- read_csv(file_path)
    
    # Create a new column 'Gene_Type' to differentiate functional genes and pseudogenes
    data <- data %>%
      mutate(Gene_Type = ifelse(grepl("P$", Gene), "Pseudogene", "Functional"))
    
    # Separate functional genes and pseudogenes
    functional_genes <- data$Median_Conservation[data$Gene_Type == "Functional"]
    pseudogenes <- data$Median_Conservation[data$Gene_Type == "Pseudogene"]
    
    # Append to pooled lists
    pooled_functional <- c(pooled_functional, functional_genes)
    pooled_pseudogenes <- c(pooled_pseudogenes, pseudogenes)
    
    # Print a message indicating that the gene has been processed
    cat("Processed:", gene, "in conservation type", conservation, "\n")
  }
  
  # Perform the KS test on pooled data if there's enough data in both groups
  if(length(pooled_functional) >= 2 && length(pooled_pseudogenes) >= 2) {
    cat("Running KS test on pooled data for conservation type", conservation, "\n")
    ks_result <- ks.test(pooled_functional, pooled_pseudogenes)
    
    # Store the pooled test result
    test_results <- rbind(test_results, data.frame(Conservation_Type = conservation,
                                                   Gene_group = "Pooled",
                                                   Test = "KS",
                                                   Statistic = ks_result$statistic,
                                                   p_value = ks_result$p.value,
                                                   stringsAsFactors = FALSE))
  } else {
    cat("Skipping pooled KS test for conservation type", conservation, "due to insufficient data\n")
  }
}

# Write the results into a CSV file
write.csv(test_results, "../../results/ks_test_conservation_pooled.csv", row.names = FALSE)


#test random data kolmogorov smirnov test
x1 <- rnorm(100)
ks.test(x1, "pnorm")

