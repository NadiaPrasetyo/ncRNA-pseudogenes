# Load necessary libraries
library(ggplot2)
library(ggrepel)  # For non-overlapping labels
library(dplyr)    # For data manipulation
library(stringr)

# Define the list of datasets and corresponding colors
datasets <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6")
colors <- c("lightblue", "lightpink", "lightyellow", "plum2", "lightgreen")  # Different colors for each dataset

# Loop through each dataset
for (i in 1:length(datasets)) {
  
  # Read the CSV data for the current dataset
  data <- read.csv(paste0("../../data/phastCons30_summary/", datasets[i], "_summary_metrics.csv"))
  
  # Calculate thresholds for high values
  max_threshold <- mean(data$Max_Conservation, na.rm = TRUE) + 2 * sd(data$Max_Conservation, na.rm = TRUE)
  median_threshold <- mean(data$Median_Conservation, na.rm = TRUE) + 2 * sd(data$Median_Conservation, na.rm = TRUE)
  
  # Filter data for genes with high max or median conservation
  label_data <- data %>%
    filter(Max_Conservation > max_threshold | Median_Conservation > median_threshold)
  
  # Create a new column for styling (check if 'Gene' does not contain a "p")
  label_data <- label_data %>%
    mutate(is_bold = !str_detect(Gene, "P"),  # TRUE if 'Gene' does not contain "p"
           label_color = if_else(is_bold, "Gene", "Pseudogene"),  # Set color to blue if not "p"
           fontface = if_else(is_bold, "bold", "plain"))  # Set fontface to bold if not "p"
  
  # Calculate the density for the data
  density_median <- density(data$Median_Conservation, na.rm = TRUE)
  density_max <- density(data$Max_Conservation, na.rm = TRUE)
  
  # Convert density data to data frames
  density_median_df <- data.frame(x = density_median$x, y = density_median$y)
  density_max_df <- data.frame(x = density_max$x, y = density_max$y)
  
  # Plot Median Conservation with non-overlapping gene labels for high values
  median_plot <- ggplot(data, aes(x = Median_Conservation)) +
    geom_density(fill = colors[i], alpha = 0.8) +
    labs(title = paste("Distribution of Median Conservation for", datasets[i]),
         x = "Median Conservation",
         y = "Density") +
    theme_minimal() +
    geom_line(data = density_median_df, aes(x = x, y = y), color = colors[i]) +  # Add density curve
    geom_text_repel(data = label_data, aes(x = Median_Conservation, y = 0, label = Gene,
                                           color = label_color, fontface = fontface),
                    nudge_y = 0.05,   # Nudge labels slightly above the points
                    direction = "y",   # Stack vertically
                    max.overlaps = 200)  # Increase the number of allowed overlaps
  
  # Save the median plot with the appropriate file name
  ggsave(paste0("../../results/rough-plots/", datasets[i], "_Median_Conservation_Distribution_Labeled.pdf"), 
         plot = median_plot, width = 20, height = 6)
  
  # Plot Median Conservation without gene labels
  median_plot_no_labels <- ggplot(data, aes(x = Median_Conservation)) +
    geom_density(fill = colors[i], alpha = 0.8) +
    labs(title = paste("Distribution of Median Conservation for", datasets[i]),
         x = "Median Conservation",
         y = "Density") +
    theme_minimal() +
    geom_line(data = density_median_df, aes(x = x, y = y), color = colors[i])  # Add density curve
  
  # Save the median plot without labels
  ggsave(paste0("../../results/rough-plots/", datasets[i], "_Median_Conservation_Distribution.pdf"), 
         plot = median_plot_no_labels, width = 20, height = 6)
  
  # Plot Max Conservation with non-overlapping gene labels for high values
  max_plot <- ggplot(data, aes(x = Max_Conservation)) +
    geom_density(fill = colors[i], alpha = 0.8) +
    labs(title = paste("Distribution of Max Conservation for", datasets[i]),
         x = "Max Conservation",
         y = "Density") +
    theme_minimal() +
    geom_line(data = density_max_df, aes(x = x, y = y), color = colors[i]) +  # Add density curve
    geom_text_repel(data = label_data, aes(x = Max_Conservation, y = 0, label = Gene,
                                           color = label_color, fontface = fontface),
                    nudge_y = 0.05,    # Nudge labels slightly above the points
                    direction = "y",   # Stack vertically
                    max.overlaps = 200)  # Increase the number of allowed overlaps
  
  # Save the max plot with the appropriate file name
  ggsave(paste0("../../results/rough-plots/", datasets[i], "_Max_Conservation_Distribution_Labeled.pdf"), 
         plot = max_plot, width = 20, height = 6)
  
  # Plot Max Conservation without gene labels
  max_plot_no_labels <- ggplot(data, aes(x = Max_Conservation)) +
    geom_density(fill = colors[i], alpha = 0.8) +
    labs(title = paste("Distribution of Max Conservation for", datasets[i]),
         x = "Max Conservation",
         y = "Density") +
    theme_minimal() +
    geom_line(data = density_max_df, aes(x = x, y = y), color = colors[i])  # Add density curve
  
  # Save the max plot without labels
  ggsave(paste0("../../results/rough-plots/", datasets[i], "_Max_Conservation_Distribution.pdf"), 
         plot = max_plot_no_labels, width = 20, height = 6)
  
}
