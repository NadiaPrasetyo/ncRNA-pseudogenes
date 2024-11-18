# Load necessary libraries
library(ggplot2)

# Step 1: Load the data
data <- read.csv("../../results/combined_gene_data.csv")  # Replace with your file path

# Remove rows with NA values
data <- data[complete.cases(data), ]

# Step 2: Select numerical columns for PCA
numerical_features <- data[, c("PhastCons30_median", "PhyloP100_median", "PhyloP447_median", "GTEX_max")]

# Step 3: Standardize the data
numerical_features_scaled <- scale(numerical_features)

# Step 4: Perform PCA
pca_result <- prcomp(numerical_features_scaled, center = TRUE, scale. = TRUE)

# Step 5: Extract PCA components for plotting
pca_data <- data.frame(pca_result$x)
pca_data$Gene_Type <- data$Gene_Type  # Add the gene type for coloring in the plot

# Step 6: Plot the PCA
plot<-ggplot(pca_data, aes(x = PC1, y = PC2, color = Gene_Type)) +
  geom_point(size = 3) +
  labs(title = "PCA of Gene types", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()

# Step 7: Save the plot
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
plot<-ggplot(data.frame(PC = seq_along(variance_explained),
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

# Display loadings for each PC
print(loadings)



