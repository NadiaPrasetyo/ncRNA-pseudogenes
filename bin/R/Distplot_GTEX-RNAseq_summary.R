# Load required libraries
library(ggplot2)
library(ggrepel)  # For non-overlapping labels
library(dplyr)    # For data manipulation
library(stringr)

# List of datasets and corresponding colors
datasets <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6")
colors <- c("cornflowerblue", "firebrick2", "goldenrod1", "mediumorchid2", "mediumseagreen")  # Different colors for each dataset

# Loop through each dataset
for (i in 1:length(datasets)) {
  
  # Load the data from the CSV file for the current dataset
  # Assuming your files are named "RNU1_expr.csv", "RNU2_expr.csv", etc.
  data <- read.csv(paste0("../../data/", datasets[i], "_expr.csv"))
  
  # Calculate thresholds for high max expression (2 standard deviations above the mean)
  mean_expr <- mean(data$Max.Expression, na.rm = TRUE)
  sd_expr <- sd(data$Max.Expression, na.rm = TRUE)
  high_expr_threshold <- mean_expr + 2 * sd_expr
  
  # Filter data for genes with high max expression
  label_data <- data %>%
    filter(Max.Expression > high_expr_threshold)
  
  # Create a new column for styling: check if 'Gene' does not contain a "P"
  label_data <- label_data %>%
    mutate(is_bold = !str_detect(Gene, "P"),  # TRUE if 'Gene' does not contain "P"
           label_color = if_else(is_bold, "Gene", "Pseudogene"),  # Set color to blue if not "P"
           fontface = if_else(is_bold, "bold", "plain"))  # Set fontface to bold if not "P"
  
  # Get the min and max values of Max.Expression for setting x-axis limits
  min_expr <- min(data$Max.Expression, na.rm = TRUE)
  max_expr <- max(data$Max.Expression, na.rm = TRUE)
  
  # Set custom extensions for the x-axis to ensure the tails are fully visible
  x_axis_extension <- (max_expr - min_expr) * 0.02  # Extend by 10% beyond the data range
  
  # Create the density plot for Max Expression without labels
  max_expr_plot_no_labels <- ggplot(data, aes(x = Max.Expression)) +
    geom_density(fill = colors[i], color = colors[i], alpha = 0.8) +  # Density plot with dataset-specific color
    labs(title = paste("Distribution of Max Expression Across Genes -", datasets[i]),
         x = "Max Expression",
         y = "Density") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_continuous(limits = c(min_expr - x_axis_extension, max_expr + x_axis_extension))  # Extend x-axis to show tails
  
  # Create the density plot for Max Expression with non-overlapping gene labels
  max_expr_plot_labeled <- ggplot(data, aes(x = Max.Expression)) +
    geom_density(fill = colors[i], color = colors[i], alpha = 0.8) +  # Density plot with dataset-specific color
    labs(title = paste("Distribution of Max Expression Across Genes -", datasets[i]),
         x = "Max Expression",
         y = "Density") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text_repel(data = label_data, aes(x = Max.Expression, y = 0, label = Gene,
                                           color = label_color, fontface = fontface),
                    nudge_y = 0.05,    # Nudge labels slightly above the points
                    direction = "y",   # Stack vertically
                    max.overlaps = 200) +  # Increase the number of allowed overlaps
    scale_x_continuous(limits = c(min_expr - x_axis_extension, max_expr + x_axis_extension))  # Extend x-axis to show the full tail
  
  # Save the unlabeled plot as a PDF with dataset-specific filenames
  ggsave(paste0("../../results/rough-plots/", datasets[i], "_Max_Expression_Distribution.pdf"), 
         plot = max_expr_plot_no_labels, width = 10, height = 6)
  
  # Save the labeled plot as a PDF with dataset-specific filenames
  ggsave(paste0("../../results/rough-plots/", datasets[i], "_Max_Expression_Distribution_Labeled.pdf"), 
         plot = max_expr_plot_labeled, width = 10, height = 6)
}
