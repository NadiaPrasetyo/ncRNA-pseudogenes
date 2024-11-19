# Load required libraries
library(ggplot2)
library(ggrepel)  # For non-overlapping labels
library(dplyr)    # For data manipulation
library(stringr)

# List of datasets and corresponding colors
datasets <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6","RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12", "VTRNA", "RNY","TRNA","RN7SL","RNU7","RN7SK")
colors <- c("cornflowerblue", "firebrick2", "goldenrod1", "mediumseagreen", "mediumorchid2", "maroon1", "cyan2", "sienna", "snow3", "green4", "dodgerblue4", "darkorchid4", "lightsalmon4", "red3", "darkolivegreen")  # Different colors for each dataset

# Loop through each dataset
for (i in 1:length(datasets)) {
  
  # Load the data from the CSV file for the current dataset
  data <- read.csv(paste0("../../data/ENCODE-expr/", datasets[i], "_expr_data.csv"))
  
  # Calculate 10% of the count for minimum y value
  minimum_y <- nrow(data) * 0.05
  
  # Calculate thresholds for high max expression (2 standard deviations above the mean)
  mean_expr <- mean(data$Max_FPKM, na.rm = TRUE)
  sd_expr <- sd(data$Max_FPKM, na.rm = TRUE)
  high_expr_threshold <- mean_expr + 2 * sd_expr
  
  # Filter data for genes with high max expression
  label_data <- data %>%
    filter(Max_FPKM > high_expr_threshold)
  
  # Create a new column for styling: check if 'Gene' does not contain a "P"
  label_data <- label_data %>%
    mutate(is_bold = !str_detect(Gene.Symbol, "P"),
           label_color = if_else(is_bold, "Gene.Symbol", "Pseudogene"),
           fontface = if_else(is_bold, "bold", "plain"))
  
  # Create the density data for overlay
  density_data <- density(data$Max_FPKM, na.rm = TRUE)
  density_median_df <- data.frame(x = density_data$x, y = density_data$y)
  
  data <- data %>%
    filter(is.finite(Max_FPKM))
  
  # Create the density plot for Max Expression without labels
  max_expr_plot_no_labels <- ggplot(data, aes(x = Max_FPKM)) +
    geom_histogram(fill = colors[i], color = colors[i], alpha = 0.8) +
    labs(title = paste("Distribution of Max Expression Across Genes -", datasets[i]),
         x = "Max Expression",
         y = "Density") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Create the density plot for Max Expression with non-overlapping gene labels
  max_expr_plot_labeled <- ggplot(data, aes(x = Max_FPKM)) +
    geom_histogram(fill = colors[i], color = colors[i], alpha = 0.8) +
    labs(title = paste("Distribution of Max Expression Across Genes -", datasets[i]),
         x = "Max Expression",
         y = "Density") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_line(data = density_median_df, aes(x = x, y = y), color = colors[i]) +  # Add density curve
    geom_text_repel(data = label_data, 
                    aes(x = Max_FPKM, y = pmax(minimum_y, 0), label = Gene.Symbol,
                        color = label_color, fontface = fontface),
                    nudge_y = 0,
                    box.padding = 0.5,
                    direction = "y",
                    max.overlaps = 200,
                    segment.size = 0,
                    segment.color = NA,
                    alpha = 0.7,
                    force = 1,
                    nudge_x = 0)
  
  # Save the unlabeled plot as a PDF with dataset-specific filenames
  ggsave(paste0("../../results/rough-plots/unlabeled/", datasets[i], "_Max_ENCODE_expr_Distribution.pdf"), 
         plot = max_expr_plot_no_labels, width = 10, height = 6)
  
  # Save the labeled plot as a PDF with dataset-specific filenames
  ggsave(paste0("../../results/rough-plots/", datasets[i], "_Max_ENCODE_expr_Distribution_Labeled.pdf"), 
         plot = max_expr_plot_labeled, width = 10, height = 6)
}