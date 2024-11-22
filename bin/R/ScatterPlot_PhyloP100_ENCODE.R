# Load necessary libraries
library(ggplot2)
library(readr)  # For reading the CSV file

# Load the dataset from the CSV file
data <- read_csv("../../results/combined_gene_data.csv")

data <- data %>%
  group_by(Gene_group) %>%
  mutate(ENCODE_max_scaled = scale(ENCODE_max)) %>%
  ungroup()

# Create the scatter plot with labels for each point
plot <- ggplot(data, aes(x = PhyloP100_median, y = ENCODE_max_scaled, color = Gene_Type, label = Gene)) +
  geom_point(alpha = 0.8, size = 3) +  # Add points with transparency (alpha = 0.8)
  geom_text(aes(label = Gene), angle = 90, vjust = -1 , hjust = 0.5, size = 5, check_overlap = FALSE, family = "Times") +  # Add text labels
  labs(title = "Scatter Plot of PhyloP_median vs ENCODE_max",
       x = "PhyloP_median",
       y = "ENCODE_max (Scaled)") +
  theme_minimal() +  # Clean minimal theme
  scale_color_manual(values = c("Pseudogene" = "cornflowerblue", "Functional" = "firebrick")) +  # Custom colors for Gene_Type
  theme(
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    plot.title = element_text(face = "bold", size = 24),  # Larger title font size
    axis.title.x = element_text(size = 22),              # Larger x-axis label font size
    axis.title.y = element_text(size = 22),              # Larger y-axis label font size
    axis.text.x = element_text(size = 24),                # Larger X axis text
    axis.text.y = element_text(size = 22),                # Larger y-axis tick font size
    text = element_text(family = "serif")  # Set font family to serif (Times New Roman)
  )
#scale_y_log10()  # Apply log scale to the y-axis

# Save the plot as a PDF
ggsave("../../results/Scatterplot_phyloP100_vs_ENCODE_with_labels.pdf", plot = plot, width = 10)
