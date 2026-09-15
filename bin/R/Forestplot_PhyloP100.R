# Load necessary libraries
library(ggplot2)
library(dplyr)

# Data-loading step
data <- read.csv("../../results/combined_gene_data_with_snp.csv")

# Define ONE vibrant base color per gene group, then auto-derive a pastel
# ("washed out") version of the same hue for that group's pseudogenes.
# This way Functional = vibrant, Pseudogene = pastel, for every group,
# and color + shape together fully encode Gene_group x Gene_Type.
base_colors <- c(
  "RNU1"      = "#1f77b4",
  "RNU2"      = "#d62728",
  "RNU4"      = "#e8a800",
  "RNU5"      = "#2ca02c",
  "RNU6"      = "#9467bd",
  "RNU4ATAC"  = "#8c564b",
  "RNU6ATAC"  = "#e377c2",
  "RNU11"     = "#555555",
  "RNU12"     = "#8c8d00",
  "VTRNA"     = "#17becf",
  "RNY"       = "#e0507a",
  "TRNA"      = "#7f5aa2",
  "RN7SL"     = "#f4a300",
  "RNU7"      = "#6a3fbf",
  "RN7SK"     = "#3b8f6b"
)

# Blend a hex color toward white to get a pastel version (amount = 0-1,
# higher = more washed out)
make_pastel <- function(hex, amount = 0.6) {
  rgb_val <- grDevices::col2rgb(hex)
  pastel_val <- rgb_val + (255 - rgb_val) * amount
  grDevices::rgb(pastel_val[1, ], pastel_val[2, ], pastel_val[3, ], maxColorValue = 255)
}

pastel_colors <- setNames(make_pastel(base_colors), names(base_colors))

# Build the final color lookup keyed by "GeneGroup.Functional" / "GeneGroup.Pseudogene"
custom_colors <- c(
  setNames(base_colors, paste0(names(base_colors), ".Functional")),
  setNames(pastel_colors, paste0(names(pastel_colors), ".Pseudogene"))
)

# Clean the data: Remove rows with non-finite PhyloP100_median
clean_data <- data %>%
  filter(is.finite(PhyloP100_median))

# -----------------------------------------------------------------------
# 4. POOLED Z-SCORE: use ALL pseudogenes (across every Gene_group) as a
#    single global "negative control" distribution, rather than computing
#    a separate median/MAD per Gene_group.
# -----------------------------------------------------------------------
global_control_median <- median(
  clean_data$PhyloP100_median[clean_data$Gene_Type == "Pseudogene"],
  na.rm = TRUE
)

global_control_mad <- mad(
  clean_data$PhyloP100_median[clean_data$Gene_Type == "Pseudogene"],
  constant = 1.4826,
  na.rm = TRUE
)

# Fallback mean absolute deviation, in case the global MAD is ~0
global_control_mean_abs_dev <- mean(
  abs(clean_data$PhyloP100_median[clean_data$Gene_Type == "Pseudogene"] - global_control_median),
  na.rm = TRUE
)

if (!is.na(global_control_mad) && global_control_mad > 1e-6) {
  normalized_data <- clean_data %>%
    mutate(Z_score = (PhyloP100_median - global_control_median) / global_control_mad)
} else {
  normalized_data <- clean_data %>%
    mutate(Z_score = (PhyloP100_median - global_control_median) / (1.2533 * global_control_mean_abs_dev))
}

# -----------------------------------------------------------------------
# 3. Keep the SAME ordering of Gene_groups that was used previously
#    (combined_groups, then combined_groups_2, then everything else),
#    but now put every group into ONE long plot.
# -----------------------------------------------------------------------
combined_groups   <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6")
combined_groups_2 <- c("RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12")
remaining_ncRNAs  <- setdiff(unique(normalized_data$Gene_group), c(combined_groups, combined_groups_2))

gene_group_order <- c(combined_groups, combined_groups_2, remaining_ncRNAs)

normalized_data <- normalized_data %>%
  filter(Gene_group %in% gene_group_order) %>%
  mutate(
    Gene_group = factor(Gene_group, levels = gene_group_order),
    Gene_Type  = factor(Gene_Type, levels = c("Functional", "Pseudogene")),
    Gene_Type_combined = paste(Gene_group, Gene_Type, sep = ".")
  )

# -----------------------------------------------------------------------
# 1 & 2. Single plot, forest-plot style (Z-score on x-axis, gene GROUP on
#    y-axis - no "(F)"/"(P)" suffix). Functional vs Pseudogene is encoded
#    ONLY via color (vibrant vs pastel) and shape (solid vs hollow). Both
#    types are drawn on the same row per Gene_group, with a small vertical
#    dodge so the two point clouds don't sit exactly on top of each other.
# -----------------------------------------------------------------------
group_levels_rev <- rev(gene_group_order)   # so the first group ends up at the top of the plot
normalized_data <- normalized_data %>%
  mutate(
    Gene_group_f = factor(Gene_group, levels = group_levels_rev),
    y_num    = as.numeric(Gene_group_f),
    y_offset = ifelse(Gene_Type == "Functional", 0.15, -0.15),
    y_pos    = y_num + y_offset
  )

# Precompute one median Z-score per Gene_group x Gene_Type, placed at the
# SAME fixed y_offset used for that type's dodge, so we get exactly one
# diamond per row-half.
median_data <- normalized_data %>%
  group_by(Gene_group, Gene_Type, Gene_Type_combined) %>%
  summarize(
    Median_Z = median(Z_score, na.rm = TRUE),
    y_pos = first(y_num) + ifelse(first(Gene_Type) == "Functional", 0.15, -0.15),
    .groups = "drop"
  )

# -----------------------------------------------------------------------
# Some Gene_group x Gene_Type combos have too few points (1-2) to draw a
# meaningful violin. Split those out so they're shown as plain points
# instead, using the same solid/hollow shape convention.
# -----------------------------------------------------------------------
group_counts <- normalized_data %>%
  count(Gene_Type_combined, name = "n_points")

thin_groups <- group_counts %>% filter(n_points <= 2) %>% pull(Gene_Type_combined)

violin_data <- normalized_data %>% filter(!(Gene_Type_combined %in% thin_groups))
thin_data   <- normalized_data %>% filter(Gene_Type_combined %in% thin_groups)

# -----------------------------------------------------------------------
# 5. Outlier detection (also exported to CSV further below).
#    - Pseudogenes with Z > 2  -> "significant" (more conserved than expected)
#    - Pseudogenes with Z > 3  -> "extreme"
#    - Functional genes with Z < 0 -> "functional_outlier"
#    These outliers are the ONLY individual points drawn as jitter, on top
#    of the violins. Solid = functional, hollow = pseudogene, same
#    vibrant/pastel color scheme.
# -----------------------------------------------------------------------
pseudogene_outliers <- normalized_data %>%
  filter(Gene_Type == "Pseudogene", Z_score > 2) %>%
  mutate(Outlier_category = ifelse(Z_score > 3, "extreme", "significant"))

functional_outliers <- normalized_data %>%
  filter(Gene_Type == "Functional", Z_score < 0) %>%
  mutate(Outlier_category = "functional_outlier")

outlier_points <- bind_rows(pseudogene_outliers, functional_outliers) %>%
  arrange(Gene_group, Gene_Type, desc(Z_score))

forest_plot <- ggplot() +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 3, linetype = "dotted", color = "grey60", linewidth = 0.4) +
  geom_violin(
    data = violin_data,
    aes(x = Z_score, y = y_pos, group = Gene_Type_combined,
        fill = Gene_Type_combined, color = Gene_Type_combined),
    orientation = "y",
    scale = "width",
    width = 0.28,
    alpha = 0.55,
    trim = TRUE,
    linewidth = 0.4
  ) +
  geom_jitter(
    data = thin_data,
    aes(x = Z_score, y = y_pos, color = Gene_Type_combined, shape = Gene_Type),
    height = 0.06, width = 0, size = 2.8, alpha = 0.9, stroke = 0.8
  ) +
  geom_jitter(
    data = outlier_points,
    aes(x = Z_score, y = y_pos, color = Gene_Type_combined, shape = Gene_Type),
    height = 0.06, width = 0, size = 2.8, alpha = 0.9, stroke = 0.8
  ) +
  geom_point(
    data = median_data,
    aes(x = Median_Z, y = y_pos),
    shape = 23, size = 2, fill = "white", color = "black"
  ) +
  scale_shape_manual(values = c("Functional" = 16, "Pseudogene" = 1)) +  # solid vs hollow
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(
    breaks = seq_along(group_levels_rev),
    labels = group_levels_rev,
    expand = expansion(add = 0.6)
  ) +
  theme_minimal() +
  labs(
    title = "Vertebrate Conservation (PhyloP100) of Functional Genes vs Pseudogenes",
    x = "Z-score (relative to pooled pseudogene distribution)",
    y = NULL
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 13),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    #panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    text = element_text(family = "serif")
  )

# Save the single long forest-style plot (tall, since every gene group is
# its own row, with both types dodged within that row)
n_rows <- length(gene_group_order)
ggsave(
  filename = "../../results/phylop100_Z_scores_all_groups_forest_plot.pdf",
  plot = forest_plot,
  width = 10,
  height = max(7, n_rows * 0.7)
)

# Export the same outlier set used in the plot to CSV
write.csv(
  outlier_points,
  file = "../../results/phylop100_Z_score_outliers.csv",
  row.names = FALSE
)