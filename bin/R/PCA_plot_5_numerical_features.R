# Step 1: Load necessary libraries
library(ggplot2)
library(ggfortify)

# Step 2: Load the data
data <- read.csv("../../results/combined_gene_data.csv")  # Replace with your file path

# Remove rows with NA values
data <- data[complete.cases(data), ]

# Step 3: Select numerical columns for PCA
numerical_features <- data[, c("PhastCons30_median", "PhyloP100_median", "PhyloP447_median", "GTEX_max", "ENCODE_max")]

# Step 4: Standardize the data
numerical_features_scaled <- scale(numerical_features)

# Step 5: Perform PCA
pca_result <- prcomp(numerical_features_scaled, scale. = TRUE)

# Step 6: Extract PCA scores (principal components) for each observation
pca_scores <- pca_result$x

# Step 7: Combine PCA scores with the original data (e.g., Gene_Type) for easier interpretation
pca_data <- cbind(data, pca_scores)

# Step 8: Print the PCA points into a table (or save as a CSV)
write.csv(pca_data, "../../results/PCA_points_table.csv", row.names = FALSE)

# Optionally, display the first few rows of the table
head(pca_data)

# Step 9: Plot data
plot <- autoplot(pca_result, data = data, colour = "Gene_Type", loadings = TRUE, 
                 loadings.colour = 'cornflowerblue', loadings.label = TRUE, 
                 loadings.label.colour = "cornflowerblue", loadings.label.size  = 3)

# Step 10: Save the plot
ggsave("../../results/PCA_plot.pdf", plot = plot)

# Variance explained by each PC
variance_explained <- (pca_result$sdev)^2

# Proportion of variance explained (converted to percentage)
prop_variance_explained <- (variance_explained / sum(variance_explained)) * 100

# Create a data frame for display
results <- data.frame(
  PC = paste0("PC", 1:length(variance_explained)),
  Variance = variance_explained,
  Proportion = prop_variance_explained,
  Cumulative = cumsum(prop_variance_explained)
)

# Display results
print(results)

# Scree plot with labels on bars
plot <- ggplot(data.frame(PC = seq_along(variance_explained),
                          Variance = prop_variance_explained),
               aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = round(Variance, 2)), vjust = -0.5, size = 3) +
  geom_line(aes(y = cumsum(prop_variance_explained), group = 1), color = "red") +
  geom_point(aes(y = cumsum(prop_variance_explained)), color = "red") +
  labs(title = "Scree Plot",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal()

ggsave("../../results/PCA_plot_variance_summary.pdf", plot = plot)

# Loadings (contributions of variables to PCs)
loadings <- pca_result$rotation

print(loadings)

# Extract loadings for PC1 and PC2
loadings_df <- data.frame(
  Variable = rownames(loadings),
  PC1 = loadings[, "PC1"],
  PC2 = loadings[, "PC2"]
)

# Loadings plot
plot <- ggplot(loadings_df, aes(x = PC1, y = PC2, label = Variable)) +
  geom_point(size = 3, color = "blue") +
  geom_text(vjust = 1.5, size = 4) +  # Labels for variables
  labs(
    title = "Loadings Plot for PC1 and PC2",
    x = "PC1 Loadings",
    y = "PC2 Loadings"
  ) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed") +  # Add horizontal axis
  geom_vline(xintercept = 0, linetype = "dashed")    # Add vertical axis

ggsave("../../results/PCA_plot_loadings.pdf", plot = plot)
