# Load necessary libraries
library(ggplot2)
library(ggrepel)  # For non-overlapping labels
library(dplyr)    # For data manipulation

# Read the CSV data into R
data <- read.csv("../../results/RNU1_summary_metrics.csv")

# Calculate thresholds for high values
max_threshold <- mean(data$Max_Conservation, na.rm = TRUE) + 2 * sd(data$Max_Conservation, na.rm = TRUE)
median_threshold <- mean(data$Median_Conservation, na.rm = TRUE) + 2 * sd(data$Median_Conservation, na.rm = TRUE)

# Filter data for genes with high max or median conservation
label_data <- data %>%
  filter(Max_Conservation > max_threshold | Median_Conservation > median_threshold)

# Calculate the density for the data
density_median <- density(data$Median_Conservation, na.rm = TRUE)
density_max <- density(data$Max_Conservation, na.rm = TRUE)

# Convert density data to data frames
density_median_df <- data.frame(x = density_median$x, y = density_median$y)
density_max_df <- data.frame(x = density_max$x, y = density_max$y)

#Plot Median conservation distrbution with no labels
median_plot <- ggplot(data, aes(x = Median_Conservation)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Distribution of Median Conservation",
       x = "Median Conservation",
       y = "Density") +
  theme_minimal()

# Save the median plot as a PDF
ggsave("../../results/Median_Conservation_Distribution.pdf", plot = median_plot, width = 10, height = 6)

# Plot Median Conservation with non-overlapping gene labels for high values
median_plot <- ggplot(data, aes(x = Median_Conservation)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Distribution of Median Conservation",
       x = "Median Conservation",
       y = "Density") +
  theme_minimal() +
  geom_line(data = density_median_df, aes(x = x, y = y), color = "blue") +  # Add density curve
  geom_text_repel(data = label_data, aes(x = Median_Conservation, y = 0, label = Gene),
                  nudge_y = 0.05,   # Nudge labels slightly above the points
                  direction = "y",   # Stack vertically
                  color = "black")

# Save the median plot as a PDF
ggsave("../../results/Median_Conservation_Distribution_Labeled.pdf", plot = median_plot, width = 10, height = 6)

#Plot Max conservation with no labels
max_plot <- ggplot(data, aes(x = Max_Conservation)) +
  geom_density(fill = "red", alpha = 0.5) +
  labs(title = "Distribution of Max Conservation",
       x = "Max Conservation",
       y = "Density") +
  theme_minimal() 

# Save the max plot as a PDF
ggsave("../../results/Max_Conservation_Distribution.pdf", plot = max_plot, width = 10, height = 6)

# Plot Max Conservation with non-overlapping gene labels for high values
max_plot <- ggplot(data, aes(x = Max_Conservation)) +
  geom_density(fill = "red", alpha = 0.5) +
  labs(title = "Distribution of Max Conservation",
       x = "Max Conservation",
       y = "Density") +
  theme_minimal() +
  geom_line(data = density_max_df, aes(x = x, y = y), color = "red") +  # Add density curve
  geom_text_repel(data = label_data, aes(x = Max_Conservation, y = 0, label = Gene),
                  nudge_y = 0.05,    # Nudge labels slightly above the points
                  direction = "y",   # Stack vertically
                  color = "black")

# Save the max plot as a PDF
ggsave("../../results/Max_Conservation_Distribution_Labeled.pdf", plot = max_plot, width = 10, height = 6)
