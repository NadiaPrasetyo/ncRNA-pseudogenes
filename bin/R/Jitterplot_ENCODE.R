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
# Step 1: Compute control_median, control_mad, and control_mean_abs_dev outside of mutate
control_stats <- clean_data %>%
  filter(Gene_Type == "Pseudogene") %>%
  group_by(Gene_group) %>%
  summarise(
    control_median = median(ENCODE_max, na.rm = TRUE),
    control_mad = mad(ENCODE_max, constant = 1.4826, na.rm = TRUE),
    control_mean_abs_dev = mean(abs(ENCODE_max - median(ENCODE_max, na.rm = TRUE)), na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Join the computed control statistics back to the original data
normalized_data <- clean_data %>%
  left_join(control_stats, by = "Gene_group") %>%
  group_by(Gene_group) %>%
  mutate(
    # Avoid division by zero or undefined MAD
    Z_score = case_when(
      !is.na(control_mad) & control_mad > 1e-6 ~ 
        (ENCODE_max - control_median) / (1.4826 * control_mad),
      !is.na(control_mean_abs_dev) & control_mad <= 1e-6 ~ 
        (ENCODE_max - control_median) / (1.2533 * control_mean_abs_dev),
      TRUE ~ NA_real_  # Assign NA for cases where Z-score can't be computed
    )
  ) %>%
  ungroup() %>%
  mutate(Z_score = ifelse(is.finite(Z_score), Z_score, 1))  # Replace non-finite Z-scores with 1


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
combined_plot <- ggplot(combined_data, aes(x = Gene_Type_label, y = Z_score, color = Gene_Type_combined)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +  # Use jitter instead of violin and boxplot
  stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +  # Highlight median
  theme_minimal() +
  labs(
    title = "Gene Expression (ENCODE) of Major Spliceosomal RNAs",
    x = "Gene Type",
    y = "Z-score (Log10 scale)"
  ) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  scale_y_continuous(trans = "log10", labels = scientific)+
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 24),  # Larger title font size
    axis.title.x = element_text(size = 22),              # Larger x-axis label font size
    axis.title.y = element_text(size = 22),              # Larger y-axis label font size
    axis.text.x = element_text(size = 24, angle = 25, hjust = 1),    #Diagonal X axis text
    axis.text.y = element_text(size = 22),                # Larger y-axis tick font size
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
    )
  )

# Plot combined group 2
combined_plot_2 <- ggplot(combined_data_2, aes(x = Gene_Type_label, y = Z_score, color = Gene_Type_combined)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +  # Use jitter 
  stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +  # Highlight median
  theme_minimal() +
  labs(
    title = "Gene Expression (ENCODE) of Minor Spliceosomal RNAs",
    x = "Gene Type",
    y = "Z-score (Log10 scale)"
  ) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  scale_y_continuous(trans = "log10", labels = scientific)+
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 24),  # Larger title font size
    axis.title.x = element_text(size = 22),              # Larger x-axis label font size
    axis.title.y = element_text(size = 22),              # Larger y-axis label font size
    axis.text.x = element_text(size = 24, angle = 25, hjust = 1),    #Diagonal X axis text
    axis.text.y = element_text(size = 22),                # Larger y-axis tick font size
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
remaining_plot <- ggplot(remaining_data, aes(x = Gene_Type_label, y = Z_score, color = Gene_Type_combined)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +  # Use jitter instead of violin and boxplot
  stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +  # Highlight median
  theme_minimal() +
  labs(
    title = "Gene Expression (ENCODE) of other sncRNAs",
    x = "Gene Type",
    y = "Z-score (Log10 scale)"
  ) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  scale_y_continuous(trans = "log10", labels = scientific)+
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 24),  # Larger title font size
    axis.title.x = element_text(size = 22),              # Larger x-axis label font size
    axis.title.y = element_text(size = 22),              # Larger y-axis label font size
    axis.text.x = element_text(size = 24, angle = 25, hjust = 1),    #Diagonal X axis text
    axis.text.y = element_text(size = 22),                # Larger y-axis tick font size
    text = element_text(family = "serif")  # Set font family to serif (Times New Roman)
  )

# Save the plot for remaining ncRNAs
ggsave(filename = "../../results/ENCODE_Z_remaining_ncRNAs_jitter.pdf", plot = remaining_plot, width = 12, height = 7)


# Combine all pseudogenes and functional genes together
combined_gene_type_data <- normalized_data %>%
  filter(Gene_group!="TRNA") %>%
  mutate(
    Gene_Type_combined = ifelse(Gene_Type == "Pseudogene", "Pseudogene", "Functional")
  )

# Plot combined pseudogenes and functional genes
combined_gene_type_plot <- ggplot(combined_gene_type_data, aes(x = Gene_Type_combined, y = Z_score, color = Gene_Type_combined)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +  # Use jitter instead of violin and boxplot
  stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +  # Highlight median
  theme_minimal() +
  labs(
    title = "Gene Expression (ENCODE) of Pseudogenes vs Functional Genes",
    x = "Gene Type",
    y = "Z-score (Log10 scale)"
  ) +
  scale_color_manual(values = c("Functional" = "firebrick", "Pseudogene" = "cornflowerblue")) +  # Custom colors for each Gene_Type
  scale_y_continuous(trans = "log10", labels = scales::scientific)+
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 24),  # Larger title font size
    axis.title.x = element_text(size = 22),              # Larger x-axis label font size
    axis.title.y = element_text(size = 22),              # Larger y-axis label font size
    axis.text.x = element_text(size = 24),    # Diagonal X axis text
    axis.text.y = element_text(size = 22),                # Larger y-axis tick font size
    text = element_text(family = "serif")  # Set font family to serif (Times New Roman)
  )

# Save the plot
ggsave(filename = "../../results/ENCODE_Z_combined_Pseudogene_vs_Functional_jitter.pdf", plot = combined_gene_type_plot, width = 12, height = 7)
