# Load necessary libraries
library(ggplot2)
library(dplyr)

# Data-loading step
data <- read.csv("../../results/combined_gene_data.csv")

# Define custom colors for each Gene_group and Gene_Type combination
custom_colors <- c(
  "RNU1.Functional" = "#174a76",  # Darker, more saturated blue
  "RNU1.Pseudogene" = "#7a99c4",  # More saturated, lighter blue  
  "RNU2.Functional" = "#9f1f1f",  # Deeper red
  "RNU2.Pseudogene" = "#e54d4a",  # Richer red-pink  
  "RNU4.Functional" = "#d26b2c",  # Darker, more saturated yellow-orange
  "RNU4.Pseudogene" = "#e4c039",  # Darker, richer yellow  
  "RNU5.Functional" = "#1b6f1b",  # Darker green
  "RNU5.Pseudogene" = "#73b761",  # Richer, more saturated green  
  "RNU6.Functional" = "#5b3686",  # Darker, more vibrant purple
  "RNU6.Pseudogene" = "#8e68b2",  # Darker, more saturated lavender
  "RNU4ATAC.Functional" = "#6b3b32",  # Darker, warmer brown
  "RNU4ATAC.Pseudogene" = "#9e7b66",  # More saturated light brown
  "RNU6ATAC.Functional" = "#b04e90",  # Darker pinkish purple
  "RNU6ATAC.Pseudogene" = "#e07fbb",  # Richer, more saturated pink
  "RNU11.Functional" = "#4c4c4c",  # Darker grey
  "RNU11.Pseudogene" = "#8f8f8f",  # Darker, richer grey
  "RNU12.Functional" = "#8a7e1f",  # Darker, more saturated olive
  "RNU12.Pseudogene" = "#9f9e00",  # Darker, richer yellow-green
  "VTRNA.Functional" = "#137f8b",  # Deeper, more saturated teal
  "VTRNA.Pseudogene" = "#6ba8bb",  # Richer, darker teal-blue
  "RNY.Functional" = "#e5788a",  # Darker, deeper pink
  "RNY.Pseudogene" = "#f35f7f",  # Richer, more saturated pink-red
  "TRNA.Functional" = "#585858",  # Darker, richer grey
  "TRNA.Pseudogene" = "#8e8e8e",  # More saturated grey
  "RN7SL.Functional" = "#c57800",  # Darker, richer amber-orange
  "RN7SL.Pseudogene" = "#d78f00",  # Darker, more saturated gold
  "RNU7.Functional" = "#b48df5",  # Darker lavender purple
  "RNU7.Pseudogene" = "#9d6ee0",  # Richer, darker purple
  "RN7SK.Functional" = "#174a76",  # Same as RNU1 Functional
  "RN7SK.Pseudogene" = "#4da7d8"   # Darker, more saturated light blue
)

# Clean the data: Remove rows with non-finite ENCODE_max
clean_data <- data %>%
  filter(is.finite(ENCODE_max))

# Normalize data by Gene_group and Gene_Type
normalized_data <- clean_data %>%
  group_by(Gene_group) %>%
  mutate(
    control_median = ifelse(
      any(Gene_Type == "Pseudogene"),
      median(ENCODE_max[Gene_Type == "Pseudogene"], na.rm = TRUE),
      NA
    ),
    control_mad = ifelse(
      any(Gene_Type == "Pseudogene"),
      mad(ENCODE_max[Gene_Type == "Pseudogene"], na.rm = TRUE),
      1  # Default MAD to 1 to avoid division by zero
    ),
    control_mad = ifelse(control_mad == 0, 1e-6, control_mad),  # Avoid MAD of 0 by using a small constant
    Z_score = ifelse(
      !is.na(control_median),
      (ENCODE_max - control_median) / control_mad,
      NA
    )
  ) %>%
  ungroup() %>%
  mutate(Z_score = ifelse(is.finite(Z_score), Z_score, 0))  # Replace remaining non-finite values with 0


# Apply log transformation to Z-score
normalized_data <- normalized_data %>%
  mutate(
    Z_score_log = ifelse(
      Z_score > 0, log10(Z_score + 1),  # Log-transform positive Z-scores
      0  # Keep Z-scores <= 0 as 0
    )
  )

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

# Plot combined data with jitter points
combined_plot <- ggplot(combined_data, aes(x = Gene_Type_label, y = Z_score_log, color = Gene_Type_combined)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +  # Use jitter instead of violin and boxplot
  stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +  # Highlight median
  theme_minimal() +
  labs(
    title = "Gene Expression (ENCODE) of Major Spliceosomal RNAs",
    x = "Gene Type",
    y = "Z-score"
  ) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 22),  # Larger title font size
    axis.title.x = element_text(size = 20),              # Larger x-axis label font size
    axis.title.y = element_text(size = 20),              # Larger y-axis label font size
    axis.text.x = element_text(size = 18, angle = 25, hjust = 1),    #Diagonal X axis text
    axis.text.y = element_text(size = 18),                # Larger y-axis tick font size
    text = element_text(family = "serif")  # Set font family to serif (Times New Roman)
  )

# Save the combined jitter plot
ggsave(filename = "../../results/ENCODE_Z_combined_RNU1_RNU2_RNU4_RNU5_RNU6_jitter.pdf", plot = combined_plot, width = 12, height = 7)

# Define groups
combined_groups_2 <- c("RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12")
remaining_ncRNAs <- setdiff(unique(normalized_data$Gene_group), c(combined_groups, combined_groups_2))

# --- Plot Combined Group 2 ---
# Filter combined data for the second group
combined_data_2 <- normalized_data %>%
  filter(Gene_group %in% combined_groups_2) %>%
  mutate(
    Gene_group = factor(Gene_group, levels = combined_groups_2),
    Gene_Type_label = ifelse(Gene_Type == "Pseudogene", paste(Gene_group, "(P)"), paste(Gene_group, "(F)")),
    Gene_Type_combined = paste(Gene_group, Gene_Type, sep = ".")  # Combine Gene_group and Gene_Type for color mapping
  )

# Ensure combined_data_2 has valid ENCODE_max values
combined_data_2 <- combined_data_2 %>%
  mutate(
    Z_score = ifelse(
      Gene == "RNU12-2P",
      ENCODE_max / ENCODE_max[Gene == "RNU12-2P"],  # Reference normalization
      Z_score  # Keep the existing Z_score otherwise
    ),
    # Log transformation for Z_score, using log1p to handle zero values safely
    Z_score_log = ifelse(Z_score > 0, log(Z_score), 0)
  )

# Plot combined group 2
combined_plot_2 <- ggplot(combined_data_2, aes(x = Gene_Type_label, y = Z_score_log, color = Gene_Type_combined)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +  # Use jitter 
  stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +  # Highlight median
  theme_minimal() +
  labs(
    title = "Gene Expression (ENCODE) of Minor Spliceosomal RNAs",
    x = "Gene Type",
    y = "Z-score"
  ) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 22),  # Larger title font size
    axis.title.x = element_text(size = 20),              # Larger x-axis label font size
    axis.title.y = element_text(size = 20),              # Larger y-axis label font size
    axis.text.x = element_text(size = 18, angle = 25, hjust = 1),    #Diagonal X axis text
    axis.text.y = element_text(size = 18),                # Larger y-axis tick font size
    text = element_text(family = "serif")  # Set font family to serif (Times New Roman)
  )

# Save the combined plot for group 2
ggsave(filename = "../../results/ENCODE_Z_combined_RNU4ATAC_RNU6ATAC_RNU11_RNU12_jitter.pdf", plot = combined_plot_2, width = 12, height = 7)

# --- Plot Remaining ncRNAs ---
# Filter data for remaining ncRNAs
remaining_data <- normalized_data %>%
  filter(Gene_group %in% remaining_ncRNAs) %>%
  mutate(
    Gene_Type_label = ifelse(Gene_Type == "Pseudogene", paste(Gene_group, "(P)"), paste(Gene_group, "(F)")),
    Gene_Type_combined = paste(Gene_group, Gene_Type, sep = ".")  # Combine Gene_group and Gene_Type for color mapping
  )


# Plot remaining ncRNAs
remaining_plot <- ggplot(remaining_data, aes(x = Gene_Type_label, y = Z_score_log, color = Gene_Type_combined)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +  # Use jitter instead of violin and boxplot
  stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +  # Highlight median
  theme_minimal() +
  labs(
    title = "Gene Expression (ENCODE) of other sncRNAs",
    x = "Gene Type",
    y = "Z-score"
  ) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 22),  # Larger title font size
    axis.title.x = element_text(size = 20),              # Larger x-axis label font size
    axis.title.y = element_text(size = 20),              # Larger y-axis label font size
    axis.text.x = element_text(size = 18, angle = 25, hjust = 1),    #Diagonal X axis text
    axis.text.y = element_text(size = 18),                # Larger y-axis tick font size
    text = element_text(family = "serif")  # Set font family to serif (Times New Roman)
  )

# Save the plot for remaining ncRNAs
ggsave(filename = "../../results/ENCODE_Z_remaining_ncRNAs_jitter.pdf", plot = remaining_plot, width = 12, height = 7)
