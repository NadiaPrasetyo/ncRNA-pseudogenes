# Load necessary libraries
library(ggplot2)
library(dplyr)

# Data-loading step
data <- read.csv("../../results/combined_gene_data.csv")

# Define custom colors for each Gene_group and Gene_Type combination
custom_colors <- c(
  "RNU1.Functional" = "#1f77b4",
  "RNU1.Pseudogene" = "#aec7e8",  
  "RNU2.Functional" = "#d62728",
  "RNU2.Pseudogene" = "#ff9896",  
  "RNU4.Functional" = "#ffc23c",
  "RNU4.Pseudogene" = "#fff26c",  
  "RNU5.Functional" = "#2ca02c",
  "RNU5.Pseudogene" = "#98df8a",  
  "RNU6.Functional" = "#9467bd",
  "RNU6.Pseudogene" = "#c5b0d5",
  "RNU4ATAC.Functional" = "#8c564b",
  "RNU4ATAC.Pseudogene" = "#c5b0a0",
  "RNU6ATAC.Functional" = "#e377c2",
  "RNU6ATAC.Pseudogene" = "#f7b6d0",
  "RNU11.Functional" = "#7f7f7f",
  "RNU11.Pseudogene" = "#bfbfbf",
  "RNU12.Functional" = "#bcbd22",
  "RNU12.Pseudogene" = "#d7e300",
  "VTRNA.Functional" = "#17becf",
  "VTRNA.Pseudogene" = "#a7d8e7",
  "RNY.Functional" = "#ffb6c1",
  "RNY.Pseudogene" = "#ff8ca0",
  "TRNA.Functional" = "#8c8c8c",
  "TRNA.Pseudogene" = "#b2b2b2",
  "RN7SL.Functional" = "#f4a300",
  "RN7SL.Pseudogene" = "#f9c100",
  "RNU7.Functional" = "#e9d3ff",
  "RNU7.Pseudogene" = "#d6a9ff",
  "RN7SK.Functional" = "#1f77b4",
  "RN7SK.Pseudogene" = "#88c9f2"
)

# Clean the data: Remove rows with non-finite PhyloP100_median
clean_data <- data %>%
  filter(is.finite(PhyloP100_median))

# Normalize data by Gene_group and Gene_Type
normalized_data <- clean_data %>%
  group_by(Gene_group) %>%
  mutate(
    control_median = ifelse(
      any(Gene_Type == "Pseudogene"),
      median(PhyloP100_median[Gene_Type == "Pseudogene"]),
      NA
    ),
    control_mad = ifelse(
      any(Gene_Type == "Pseudogene"),
      mad(PhyloP100_median[Gene_Type == "Pseudogene"]),
      1  # Use mad = 1 to avoid division by zero
    ),
    Z_score = ifelse(
      !is.na(control_median),
      (PhyloP100_median - control_median) / control_mad,
      NA
    )
  ) %>%
  ungroup() %>%
  mutate(Z_score = ifelse(is.na(Z_score), 0, Z_score))

# Create Gene_Type_combined for all gene groups (not just the combined ones)
normalized_data <- normalized_data %>%
  mutate(
    Gene_Type_combined = paste(Gene_group, Gene_Type, sep = ".")  # Combine Gene_group and Gene_Type for color mapping
  )

# Define the combined groups for RNU1, RNU2, RNU4, RNU5, and RNU6
combined_groups <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6")

# Create combined data for the selected Gene_groups
combined_data <- normalized_data %>%
  filter(Gene_group %in% combined_groups) %>%
  mutate(
    Gene_group = factor(Gene_group, levels = combined_groups),
    Gene_Type_label = ifelse(Gene_Type == "Pseudogene", paste(Gene_group, "(P)"), paste(Gene_group, "(F)")),
    Gene_Type_combined = paste(Gene_group, Gene_Type, sep = ".")  # Combine Gene_group and Gene_Type for color mapping
  )

# Define function to calculate outlier boundaries for jitter
calculate_outliers <- function(data) {
  Q1 <- quantile(data$Z_score, 0.25, na.rm = TRUE)
  Q3 <- quantile(data$Z_score, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  return(c(lower_bound, upper_bound))
}

# Calculate outliers for the combined data
outlier_bounds <- calculate_outliers(combined_data)

# Filter jitter points to only include those outside the boxplot range
combined_data_jitter <- combined_data %>%
  filter(Z_score < outlier_bounds[1] | Z_score > outlier_bounds[2])

# Plot combined data with jitter points
combined_plot <- ggplot(combined_data, aes(x = Gene_Type_label, y = Z_score, color = Gene_Type_combined)) +
  geom_jitter(width = 0.2, height = 0, size = 1.5, alpha = 0.7) +  # Use jitter instead of violin and boxplot
  stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +  # Highlight median
  theme_minimal() +
  labs(
    title = "PhyloP100 Z-scores for RNU1, RNU2, RNU4, RNU5, RNU6",
    x = "Gene Type",
    y = "Z-score"
  ) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14)
  )

# Save the combined jitter plot
ggsave(filename = "../../results/Z_scores_combined_RNU1_RNU2_RNU4_RNU5_RNU6_jitter.pdf", plot = combined_plot, width = 8, height = 6)

# Filter data for RNU12
rnu12_data <- normalized_data %>%
  filter(Gene_group == "RNU12")

rnu12_data$Z_score = rnu12_data$PhyloP100_median / rnu12_data$PhyloP100_median[rnu12_data$Gene == "RNU12-2P"]

# Plot for RNU12 with only the two singular points
rnu12_plot <- ggplot(rnu12_data, aes(x = Gene_Type, y = Z_score, color = Gene_Type_combined)) +
  geom_point(size = 4, alpha = 0.7) +  # Singular data points
  scale_color_manual(values = custom_colors) +  # Custom colors for RNU12
  labs(
    title = "PhyloP100 Z-scores for RNU12",
    x = "Gene Type",
    y = "Z-score"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14)
  ) +
  expand_limits(y = c(0, 3))

# Save the plot for RNU12
ggsave(filename = "../../results/Z_scores_RNU12_jitter.pdf", plot = rnu12_plot, width = 8, height = 6)

# Create and save individual plots for other Gene_groups
other_gene_groups <- setdiff(unique(normalized_data$Gene_group), combined_groups)
# Exclude "RNU12" directly from the vector
other_gene_groups <- other_gene_groups[other_gene_groups != "RNU12"]

for (group in other_gene_groups) {
  group_data <- normalized_data %>%
    filter(Gene_group == group) %>%
    filter(is.finite(Z_score))  # Ensure no NA values in Z_scores
  
  # Calculate outliers for the group
  outlier_bounds <- calculate_outliers(group_data)
  
  # Filter jitter points to only include those outside the boxplot range
  group_data_jitter <- group_data %>%
    filter(Z_score < outlier_bounds[1] | Z_score > outlier_bounds[2])
  
  # Create and save individual jitter plot
  p <- ggplot(group_data, aes(x = Gene_Type, y = Z_score, color = interaction(Gene_group, Gene_Type))) +
    geom_jitter(width = 0.2, height = 0, size = 1.5, alpha = 0.7) +  # Use jitter
    stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +  # Highlight median
    theme_minimal() +
    labs(
      title = paste("PhyloP100 Z-scores for Gene Group:", group),
      x = "Gene Type",
      y = "Z-score"
    ) +
    scale_color_manual(values = custom_colors) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14)
    )
  
  ggsave(filename = paste0("../../results/Z_scores_", group, "_jitter.pdf"), plot = p, width = 8, height = 6)
}