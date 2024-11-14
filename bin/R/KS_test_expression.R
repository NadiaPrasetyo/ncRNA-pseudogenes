# Load the required packages
library(dplyr)
library(readr)

# Define the gene groups
gene_groups <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6", "RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12", "VTRNA")

# Create an empty data frame to store the test results
test_results <- data.frame(Gene_group = character(),
                           Test = character(),
                           Statistic = numeric(),
                           p_value = numeric(),
                           Functional_n = integer(),
                           Pseudogene_n = integer(),
                           stringsAsFactors = FALSE)

# Create empty vectors to collect pooled data for functional and pseudogenes across all genes
all_functional_genes <- c()
all_pseudogenes <- c()

# Loop through each gene group (e.g., RNU1, RNU2, etc.)
for (gene in gene_groups) {
  
  # Read the data for the current gene group
  data <- read_csv(paste0("../../data/", gene, "_expr.csv"))
  
  # Create a new column 'Gene_Type' to differentiate functional genes and pseudogenes
  data <- data %>%
    mutate(Gene_Type = ifelse(grepl("P$", Gene), "Pseudogene", "Functional"))
  
  # Separate functional genes and pseudogenes
  functional_genes <- data$`Max Expression`[data$Gene_Type == "Functional"]
  pseudogenes <- data$`Max Expression`[data$Gene_Type == "Pseudogene"]
  
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
    test_results <- rbind(test_results, data.frame(Gene_group = gene,
                                                   Test = "KS",
                                                   Statistic = ks_result$statistic,
                                                   p_value = ks_result$p.value,
                                                   Functional_n = functional_n,
                                                   Pseudogene_n = pseudogene_n,
                                                   stringsAsFactors = FALSE))
  } else {
    test_results <- rbind(test_results, data.frame(Gene_group = gene,
                                                   Test = "KS",
                                                   Statistic = NA,
                                                   p_value = NA,
                                                   Functional_n = functional_n,
                                                   Pseudogene_n = pseudogene_n,
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
  test_results <- rbind(test_results, data.frame(Gene_group = "All_Genes_Pooled",
                                                 Test = "KS",
                                                 Statistic = pooled_ks_result$statistic,
                                                 p_value = pooled_ks_result$p.value,
                                                 Functional_n = pooled_functional_n,
                                                 Pseudogene_n = pooled_pseudogene_n,
                                                 stringsAsFactors = FALSE))
} else {
  test_results <- rbind(test_results, data.frame(Gene_group = "All_Genes_Pooled",
                                                 Test = "KS",
                                                 Statistic = NA,
                                                 p_value = NA,
                                                 Functional_n = pooled_functional_n,
                                                 Pseudogene_n = pooled_pseudogene_n,
                                                 stringsAsFactors = FALSE))
}

# Write the results into a CSV file
write.csv(test_results, "../../results/ks_test_expression_results.csv", row.names = FALSE)

# Test random data Kolmogorov-Smirnov test
x1 <- rnorm(100)
ks.test(x1, "pnorm")
