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
# First, compute control_median and control_mad outside of mutate
control_stats <- clean_data %>%
  filter(Gene_Type == "Pseudogene") %>%
  group_by(Gene_group) %>%
  summarise(
    control_median = median(PhyloP100_median, na.rm = TRUE),
    control_mad = mad(PhyloP100_median, constant = 1.4826, na.rm = TRUE),
    .groups = "drop"
  )

# Now, join the computed control stats back to the original data
normalized_data <- clean_data %>%
  left_join(control_stats, by = "Gene_group") %>%
  group_by(Gene_group) %>%
  mutate(
    control_mean_abs_dev = ifelse(
      !is.na(control_median), 
      mean(abs(PhyloP100_median[Gene_Type == "Pseudogene"] - control_median), na.rm = TRUE),
      NA
    ),
    # Avoid division by zero or undefined MAD
    Z_score = case_when(
      !is.na(control_mad) & control_mad > 1e-6 ~ 
        (PhyloP100_median - control_median) / (1.4826 * control_mad),
      !is.na(control_mean_abs_dev) & control_mad <= 1e-6 ~ 
        (PhyloP100_median - control_median) / (1.2533 * control_mean_abs_dev),
      TRUE ~ NA_real_  # Assign NA for cases where Z-score can't be computed
    )
  ) %>%
  ungroup()


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
    title = "Vertebrate Conservation (PhyloP100) of Major Spliceosomal RNAs",
    x = "Gene Type",
    y = "Z-score"
  ) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
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
ggsave(filename = "../../results/Z_scores_combined_RNU1_RNU2_RNU4_RNU5_RNU6_jitter.pdf", plot = combined_plot, width = 12, height = 7)

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

combined_data_2$Z_score = combined_data_2$PhyloP100_median / combined_data_2$PhyloP100_median[combined_data_2$Gene == "RNU12-2P"]

# Plot combined group 2
combined_plot_2 <- ggplot(combined_data_2, aes(x = Gene_Type_label, y = Z_score, color = Gene_Type_combined)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +  # Use jitter 
  stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +  # Highlight median
  theme_minimal() +
  labs(
    title = "Vertebrate Conservation (PhyloP100) of Minor Spliceosomal RNAs",
    x = "Gene Type",
    y = "Z-score"
  ) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
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
ggsave(filename = "../../results/Z_scores_combined_RNU4ATAC_RNU6ATAC_RNU11_RNU12_jitter.pdf", plot = combined_plot_2, width = 12, height = 7)

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
    title = "Vertebrate Conservation (PhyloP100) of other sncRNAs",
    x = "Gene Type",
    y = "Z-score"
  ) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
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
ggsave(filename = "../../results/Z_scores_remaining_ncRNAs_jitter.pdf", plot = remaining_plot, width = 12, height = 7)

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
    title = "Vertebrate Conservation (PhyloP100) of Pseudogenes vs Functional Genes",
    x = "Gene Type",
    y = "Z-score"
  ) +
  scale_color_manual(values = c("Functional" = "firebrick", "Pseudogene" = "cornflowerblue")) +  # Custom colors for each Gene_Type
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 24),  # Larger title font size
    axis.title.x = element_text(size = 22),              # Larger x-axis label font size
    axis.title.y = element_text(size = 22),              # Larger y-axis label font size
    axis.text.x = element_text(size = 24),    # Diagonal X axis text
    axis.text.y = element_text(size = 22),                # Larger y-axis tick font size
    text = element_text(family = "serif")  # Set font family to serif (Times New Roman)
  )

median_values <- combined_gene_type_data %>%
  group_by(Gene_Type_combined) %>%
  summarize(Median_Z = median(Z_score, na.rm = TRUE))

# Save the plot
ggsave(filename = "../../results/PhyloP100_Z_combined_Pseudogene_vs_Functional_jitter.pdf", plot = combined_gene_type_plot, width = 12, height = 7)
