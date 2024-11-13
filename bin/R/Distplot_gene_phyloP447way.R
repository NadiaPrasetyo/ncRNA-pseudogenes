# Load necessary libraries
library(ggplot2)
library(ggrepel)  # For non-overlapping labels
library(dplyr)    # For data manipulation
library(stringr)

# Define the list of datasets and corresponding colors
datasets <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6")
colors <- c("lightskyblue", "salmon", "khaki1", "orchid1", "seagreen1")  # Different colors for each dataset

# Loop through each dataset
for (i in 1:length(datasets)) {
  
  # Read the CSV data for the current dataset
  data <- read.csv(paste0("../../data/phyloP447_summary/", datasets[i], "_phyloP447_summary_metrics.csv"))

  #calculate 10% of the count
  minimum_y <- nrow(data) * 0.1
  
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
    geom_histogram(fill = colors[i], alpha = 0.8) +
    labs(title = paste("Distribution of Median Conservation for", datasets[i]),
         x = "Median Conservation",
         y = "Density") +
    theme_minimal() +
    geom_line(data = density_median_df, aes(x = x, y = y), color = colors[i]) +  # Add density curve
    geom_text_repel(data = label_data, 
                    aes(x = Median_Conservation, y = pmax(minimum_y, 0), label = Gene,  # Ensure y >= 0.1
                        color = label_color, fontface = fontface),
                    nudge_y = 0,   # Nudge labels slightly above the points
                    box.padding = 0.5,      # Add padding around each label
                    direction = "y",   # Stack vertically
                    max.overlaps = 200,  # Allow overlap of up to 200 labels
                    segment.size = 0,    # Remove the connecting lines between the labels
                    segment.color = NA,  # Ensure no segment lines are drawn
                    alpha = 0.7,
                    force = 1,          # Adjust the repelling strength for better stacking
                    nudge_x = 0)        # Avoid horizontal nudging so labels stay aligned vertically
  
  # Save the median plot with the appropriate file name
  ggsave(paste0("../../results/rough-plots/", datasets[i], "_phyloP447_Median_Conservation_Distribution_Labeled.pdf"), 
         plot = median_plot, width = 15, height = 6)
  
  # Plot Median Conservation without gene labels
  median_plot_no_labels <- ggplot(data, aes(x = Median_Conservation)) +
    geom_histogram(fill = colors[i], alpha = 0.8) +
    labs(title = paste("Distribution of Median Conservation for", datasets[i]),
         x = "Median Conservation",
         y = "Density") +
    theme_minimal() +
    geom_line(data = density_median_df, aes(x = x, y = y), color = colors[i])  # Add density curve
  
  # Save the median plot without labels
  ggsave(paste0("../../results/rough-plots/unlabeled/", datasets[i], "_phyloP447_Median_Conservation_Distribution.pdf"), 
         plot = median_plot_no_labels, width = 15, height = 6)
  
  # Plot Max Conservation with non-overlapping gene labels for high values
  max_plot <- ggplot(data, aes(x = Max_Conservation)) +
    geom_histogram(fill = colors[i], alpha = 0.8) +
    labs(title = paste("Distribution of Max Conservation for", datasets[i]),
         x = "Max Conservation",
         y = "Density") +
    theme_minimal() +
    geom_line(data = density_max_df, aes(x = x, y = y), color = colors[i]) +  # Add density curve
    geom_text_repel(data = label_data, 
                    aes(x = Max_Conservation, y = pmax(minimum_y, 0), label = Gene,  # Ensure y >= 0.1
                        color = label_color, fontface = fontface),
                    nudge_y = 0,   # Nudge labels slightly above the points
                    box.padding = 0.5,      # Add padding around each label
                    direction = "y",   # Stack vertically
                    max.overlaps = 200,  # Allow overlap of up to 200 labels
                    segment.size = 0,    # Remove the connecting lines between the labels
                    segment.color = NA,  # Ensure no segment lines are drawn
                    alpha = 0.7,
                    force = 1,          # Adjust the repelling strength for better stacking
                    nudge_x = 0)        # Avoid horizontal nudging so labels stay aligned vertically
  
  # Save the max plot with the appropriate file name
  ggsave(paste0("../../results/rough-plots/", datasets[i], "_phyloP447_Max_Conservation_Distribution_Labeled.pdf"), 
         plot = max_plot, width = 15, height = 6)
  
  # Plot Max Conservation without gene labels
  max_plot_no_labels <- ggplot(data, aes(x = Max_Conservation)) +
    geom_histogram(fill = colors[i], alpha = 0.8) +
    labs(title = paste("Distribution of Max Conservation for", datasets[i]),
         x = "Max Conservation",
         y = "Density") +
    theme_minimal() +
    geom_line(data = density_max_df, aes(x = x, y = y), color = colors[i])  # Add density curve
  
  # Save the max plot without labels
  ggsave(paste0("../../results/rough-plots/unlabeled/", datasets[i], "_phyloP447_Max_Conservation_Distribution.pdf"), 
         plot = max_plot_no_labels, width = 15, height = 6)
  
}
