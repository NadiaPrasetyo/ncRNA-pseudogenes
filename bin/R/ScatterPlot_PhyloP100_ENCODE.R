# Load necessary libraries
library(ggplot2)
library(readr)  # For reading the CSV file

# Load the dataset from the CSV file
data <- read_csv("../../results/combined_gene_data.csv")

# Create the scatter plot with log scale on the ENCODE_max (y-axis) and alpha set to 0.8
plot <- ggplot(data, aes(x = PhyloP100_median, y = ENCODE_max, color = Gene_Type)) +
  geom_point(alpha = 0.8) +  # Add points with transparency (alpha = 0.8)
  labs(title = "Scatter Plot of PhyloP_median vs ENCODE_max",
       x = "PhyloP_median",
       y = "ENCODE_max (Log Scale)") +
  theme_minimal() +  # Clean minimal theme
  scale_color_manual(values = c("Pseudogene" = "cornflowerblue", "Functional" = "firebrick")) +  # Custom colors for Gene_Type
  scale_y_log10()  # Apply log scale to the y-axis

# Save the plot as a PDF
ggsave("../../results/Scatterplot_phyloP100_vs_ENCODE.pdf", plot = plot)
