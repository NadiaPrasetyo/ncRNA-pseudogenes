
# loading the required package
library(dplyr)
library("dgof")
library(readr)

# Define the gene groups and conservation types
gene_groups <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6")

# Create an empty list to store the test results
test_results <- data.frame(Gene_group = character(),
                         Test = character(),
                         Statistic = numeric(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE)

# Loop through each gene group (e.g., RNU1, RNU2, etc.)
for (gene in gene_groups) {

  # Read the data for the current gene group and conservation type
  data <- read_csv(paste0("../../data/", gene, "_expr.csv"))
  
  # Create a new column 'Gene_Type' to differentiate functional genes and pseudogenes
  data <- data %>%
    mutate(Gene_Type = ifelse(grepl("P$", Gene), "Pseudogene", "Functional"))
  
  # Separate functional genes and pseudogenes
  functional_genes <- data$`Max Expression`[data$Gene_Type == "Functional"]
  pseudogenes <- data$`Max Expression`[data$Gene_Type == "Pseudogene"]
  
  # Always run the KS test for RNU2
  if (gene == "RNU2") {
    ks_result <- ks.test(functional_genes, pseudogenes)
    test_results <- rbind(test_results, data.frame(Gene_group = gene,
                                                   Test = "KS",
                                                   Statistic = ks_result$statistic,
                                                   p_value = ks_result$p.value,
                                                   stringsAsFactors = FALSE))
  } else {
    # Check if there is enough data for either functional or pseudogenes (at least 2 samples in each group)
    if(length(functional_genes) >= 2 && length(pseudogenes) >= 2) {
      # Perform the Kolmogorov-Smirnov test for other genes with sufficient data
      ks_result <- ks.test(functional_genes, pseudogenes)
      test_results <- rbind(test_results, data.frame(Gene_group = gene,
                                                     Test = "KS",
                                                     Statistic = ks_result$statistic,
                                                     p_value = ks_result$p.value,
                                                     stringsAsFactors = FALSE))
    } else {
      cat("Skipping:", gene,  "due to insufficient data\n")
    }
  }
  
  # Print a message indicating that the gene has been processed
  cat("Processed:", gene, "\n")
}


# Write the results into a CSV file
write.csv(test_results, "../../results/ks_test_expression_results.csv", row.names = FALSE)


#test random data kolmogorov smirnov test
x1 <- rnorm(100)
ks.test(x1, "pnorm")

