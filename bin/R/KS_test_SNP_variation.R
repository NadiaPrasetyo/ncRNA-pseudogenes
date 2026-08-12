# Loading the required packages
library(dplyr)
library(readr)

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------

# Path to the single input CSV containing all genes / all SNP databases
input_file <- "../../data/snp_enrichment_genes_pseudogenes.csv"

# Path for the output results CSV
output_file <- "../../results/ks_test_snp_enrichment_results.csv"

# Define the gene groups (same convention as the conservation analysis)
gene_groups <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6", "RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12", "VTRNA",
                 "RNY", "TRNA", "RN7SL", "RNU7", "RN7SK")

# Define the SNP databases (each has its own enrichment_<db> column in the input file)
snp_databases <- c("1000genomes", "gnomad", "pangenome", "dbsnp")

# ---------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

all_data <- read_csv(input_file, show_col_types = FALSE)

# Standardize Gene_Type using the gene_type column already provided in the input
# (handles minor variations in capitalization/whitespace, e.g. "Pseudogene", "pseudogene ")
all_data <- all_data %>%
  mutate(Gene_Type = ifelse(grepl("pseudo", tolower(trimws(gene_type))), "Pseudogene", "Functional"))

# ---------------------------------------------------------------------------
# RUN KS TESTS
# ---------------------------------------------------------------------------

# Create an empty data frame to store the test results with sample counts
test_results <- data.frame(SNP_Database = character(),
                           Gene_group = character(),
                           Statistic = numeric(),
                           p_value = numeric(),
                           Functional_Sample_Count = integer(),
                           Pseudogene_Sample_Count = integer(),
                           stringsAsFactors = FALSE)

# Loop through each SNP database (e.g. 1000genomes, gnomad, pangenome, dbsnp)
for (snp_db in snp_databases) {
  
  enrichment_col <- paste0("enrichment_", snp_db)
  
  if (!enrichment_col %in% colnames(all_data)) {
    cat("Column", enrichment_col, "not found in input data, skipping SNP database", snp_db, "\n")
    next
  }
  
  # Create empty vectors for pooled functional and pseudogene data and sample counts
  pooled_functional <- c()
  pooled_pseudogenes <- c()
  pooled_functional_count <- 0
  pooled_pseudogene_count <- 0
  
  # Loop through each gene group (e.g., RNU1, RNU2, etc.)
  for (gene in gene_groups) {
    
    # Subset the data for the current gene group
    data <- all_data %>% filter(gene_group == gene)
    
    if (nrow(data) == 0) {
      cat("No data found for gene group", gene, "in SNP database", snp_db, "\n")
      
      test_results <- rbind(test_results, data.frame(SNP_Database = snp_db,
                                                     Gene_group = gene,
                                                     Statistic = NA,
                                                     p_value = NA,
                                                     Functional_Sample_Count = 0,
                                                     Pseudogene_Sample_Count = 0,
                                                     stringsAsFactors = FALSE))
      
      next  # Skip to the next gene
    }
    
    # Separate functional genes and pseudogenes, using the enrichment values for this database
    functional_genes <- data[[enrichment_col]][data$Gene_Type == "Functional"]
    pseudogenes <- data[[enrichment_col]][data$Gene_Type == "Pseudogene"]
    
    # Drop any NA/NaN/Inf enrichment values before testing (enrichment ratios can be
    # Inf when flank_density is 0, which breaks ks.test's internal ECDF calculation)
    functional_genes <- functional_genes[is.finite(functional_genes)]
    pseudogenes <- pseudogenes[is.finite(pseudogenes)]
    
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
      test_results <- rbind(test_results, data.frame(SNP_Database = snp_db,
                                                     Gene_group = gene,
                                                     Statistic = ks_result$statistic,
                                                     p_value = ks_result$p.value,
                                                     Functional_Sample_Count = functional_count,
                                                     Pseudogene_Sample_Count = pseudogene_count,
                                                     stringsAsFactors = FALSE))
    } else {
      # If either group is empty, handle it here (e.g., no KS test performed)
      test_results <- rbind(test_results, data.frame(SNP_Database = snp_db,
                                                     Gene_group = gene,
                                                     Statistic = NA,
                                                     p_value = NA,
                                                     Functional_Sample_Count = functional_count,
                                                     Pseudogene_Sample_Count = pseudogene_count,
                                                     stringsAsFactors = FALSE))
      cat("Skipping KS test for gene", gene, "in SNP database", snp_db, "due to insufficient data\n")
    }
    
    # Print a message indicating that the gene has been processed
    cat("Processed:", gene, "in SNP database", snp_db, "\n")
  }
  
  # Perform the KS test on pooled data (only if both groups have data)
  if (pooled_functional_count > 0 && pooled_pseudogene_count > 0) {
    ks_result <- ks.test(pooled_functional, pooled_pseudogenes)
    test_results <- rbind(test_results, data.frame(SNP_Database = snp_db,
                                                   Gene_group = "Pooled",
                                                   Statistic = ks_result$statistic,
                                                   p_value = ks_result$p.value,
                                                   Functional_Sample_Count = pooled_functional_count,
                                                   Pseudogene_Sample_Count = pooled_pseudogene_count,
                                                   stringsAsFactors = FALSE))
  } else {
    # If either group is empty, handle it here for pooled data
    test_results <- rbind(test_results, data.frame(SNP_Database = snp_db,
                                                   Gene_group = "Pooled",
                                                   Statistic = NA,
                                                   p_value = NA,
                                                   Functional_Sample_Count = pooled_functional_count,
                                                   Pseudogene_Sample_Count = pooled_pseudogene_count,
                                                   stringsAsFactors = FALSE))
    cat("Skipping KS test for pooled data in SNP database", snp_db, "due to insufficient data\n")
  }
}

# ---------------------------------------------------------------------------
# WRITE RESULTS
# ---------------------------------------------------------------------------

write.csv(test_results, output_file, row.names = FALSE)

# Test random data Kolmogorov-Smirnov test (example)
x1 <- rnorm(100)
ks.test(x1, "pnorm")