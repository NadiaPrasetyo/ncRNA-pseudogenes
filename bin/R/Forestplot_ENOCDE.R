# Load necessary libraries
library(ggplot2)
library(dplyr)
library(scales)

# Data-loading step
data <- read.csv("../../results/combined_gene_data_with_snp.csv")

# Define ONE vibrant base color per gene group, then auto-derive a pastel
# ("washed out") version of the same hue for that group's pseudogenes.
# Functional = vibrant, Pseudogene = pastel, for every group; color + shape
# together fully encode Gene_group x Gene_Type.
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
  "RN7SK"     = "#3b8f6b",
  "Pooled"    = "#000000"   # NEW: distinct near-black for the pooled summary row
)

make_pastel <- function(hex, amount = 0.6) {
  rgb_val <- grDevices::col2rgb(hex)
  pastel_val <- rgb_val + (255 - rgb_val) * amount
  grDevices::rgb(pastel_val[1, ], pastel_val[2, ], pastel_val[3, ], maxColorValue = 255)
}

pastel_colors <- setNames(make_pastel(base_colors), names(base_colors))

custom_colors <- c(
  setNames(base_colors, paste0(names(base_colors), ".Functional")),
  setNames(pastel_colors, paste0(names(pastel_colors), ".Pseudogene"))
)

# Clean the data: Remove rows with non-finite ENCODE_max
clean_data <- data %>%
  filter(is.finite(ENCODE_max))

# -----------------------------------------------------------------------
# POOLED Z-SCORE: use ALL pseudogenes (across every Gene_group) as a single
# global "negative control" distribution, rather than a separate
# median/MAD per Gene_group. Same logic as the PhyloP100 forest plot.
# (Unchanged.)
# -----------------------------------------------------------------------
global_control_median <- median(
  clean_data$ENCODE_max[clean_data$Gene_Type == "Pseudogene"],
  na.rm = TRUE
)

global_control_mad <- mad(
  clean_data$ENCODE_max[clean_data$Gene_Type == "Pseudogene"],
  constant = 1.4826,
  na.rm = TRUE
)

# Fallback mean absolute deviation, in case the global MAD is ~0
global_control_mean_abs_dev <- mean(
  abs(clean_data$ENCODE_max[clean_data$Gene_Type == "Pseudogene"] - global_control_median),
  na.rm = TRUE
)

if (!is.na(global_control_mad) && global_control_mad > 1e-6) {
  normalized_data <- clean_data %>%
    mutate(Z_score = (ENCODE_max - global_control_median) / global_control_mad)
} else {
  normalized_data <- clean_data %>%
    mutate(Z_score = (ENCODE_max - global_control_median) / (1.2533 * global_control_mean_abs_dev))
}

# -----------------------------------------------------------------------
# Keep the SAME ordering of Gene_groups used in the PhyloP100 forest plot.
# A "Pooled" pseudo-group is appended LAST so it renders at the BOTTOM of
# the forest plot.
# -----------------------------------------------------------------------
combined_groups   <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6")
combined_groups_2 <- c("RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12")
remaining_ncRNAs  <- setdiff(unique(normalized_data$Gene_group), c(combined_groups, combined_groups_2))

gene_group_order <- c(combined_groups, combined_groups_2, remaining_ncRNAs, "Pooled")

# NEW: build the "Pooled" rows by duplicating every existing row (all
# Gene_groups) with Gene_group relabeled to "Pooled". This reuses the
# SAME Z_score values already computed above -- nothing is recalculated
# or rescaled -- it just lets the Functional/Pseudogene violins be drawn
# once more using every point across the whole dataset.
pooled_data <- normalized_data %>%
  mutate(Gene_group = "Pooled")

normalized_data <- bind_rows(normalized_data, pooled_data) %>%
  filter(Gene_group %in% gene_group_order) %>%
  mutate(
    Gene_group = factor(Gene_group, levels = gene_group_order),
    Gene_Type  = factor(Gene_Type, levels = c("Functional", "Pseudogene")),
    Gene_Type_combined = paste(Gene_group, Gene_Type, sep = ".")
  )

# -----------------------------------------------------------------------
# Forest-plot style layout: Z-score on x-axis, gene GROUP on y-axis (no
# "(F)"/"(P)" suffix). Functional vs Pseudogene encoded via color
# (vibrant vs pastel) and shape (solid vs hollow). Both types dodged
# within the same row per Gene_group.
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
# SAME fixed y_offset used for that type's dodge.
median_data <- normalized_data %>%
  group_by(Gene_group, Gene_Type, Gene_Type_combined) %>%
  summarize(
    Median_Z = median(Z_score, na.rm = TRUE),
    y_pos = first(y_num) + ifelse(first(Gene_Type) == "Functional", 0.15, -0.15),
    .groups = "drop"
  )

# -----------------------------------------------------------------------
# Some Gene_group x Gene_Type combos have too few points (1-2) to draw a
# meaningful violin. Show those as plain points instead.
# -----------------------------------------------------------------------
group_counts <- normalized_data %>%
  count(Gene_Type_combined, name = "n_points")

thin_groups <- group_counts %>% filter(n_points <= 2) %>% pull(Gene_Type_combined)

violin_data <- normalized_data %>% filter(!(Gene_Type_combined %in% thin_groups))
thin_data   <- normalized_data %>% filter(Gene_Type_combined %in% thin_groups)

# -----------------------------------------------------------------------
# Outlier detection (also exported to CSV further below).
#    - Pseudogenes with Z > 2  -> "significant" (more expressed than expected)
#    - Pseudogenes with Z > 3  -> "extreme"
#    - Functional genes with Z < 0 -> "functional_outlier"
#    These outliers (plus any thin groups) are the only individual points
#    drawn on top of the violins. Solid = functional, hollow = pseudogene.
#    NOTE: the "Pooled" pseudo-group is excluded from the outlier export
#    since its points are duplicates of the per-group rows above it and
#    would otherwise double-count every outlier.
# -----------------------------------------------------------------------
pseudogene_outliers <- normalized_data %>%
  filter(Gene_Type == "Pseudogene", Z_score > 2, Gene_group != "Pooled") %>%
  mutate(Outlier_category = ifelse(Z_score > 3, "extreme", "significant"))

functional_outliers <- normalized_data %>%
  filter(Gene_Type == "Functional", Z_score < 0, Gene_group != "Pooled") %>%
  mutate(Outlier_category = "functional_outlier")

outlier_points <- bind_rows(pseudogene_outliers, functional_outliers) %>%
  arrange(Gene_group, Gene_Type, desc(Z_score))

# For the Pooled row itself, jitter is still drawn (below) using the
# ordinary thin/outlier logic per Gene_Type_combined, so it visually
# matches all the other rows -- but the outlier CSV above stays specific
# to individual gene groups.
pooled_jitter <- normalized_data %>%
  filter(
    Gene_group == "Pooled",
    (Gene_Type == "Pseudogene" & Z_score > 2) |
      (Gene_Type == "Functional" & Z_score > 100)
  )

# -----------------------------------------------------------------------
# Log-scaled x-axis: Z-scores can be negative (and are often clustered
# near zero for pseudogenes), so a plain log10 scale would drop/NaN those
# values. Use a pseudo-log transform instead - it behaves like log10 for
# large |Z| but stays linear (and defined) near zero, including negatives.
# (Unchanged.)
# -----------------------------------------------------------------------
z_breaks <- c(-10, -3, -1, 0, 1, 2, 3, 10, 100, 1000, 10000, 100000)

forest_plot <- ggplot() +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 3, linetype = "dotted", color = "grey60", linewidth = 0.4) +
  # Separator line between the per-group rows and the Pooled summary row
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "grey70", linewidth = 0.4) +
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
  geom_jitter(
    data = pooled_jitter,
    aes(x = Z_score, y = y_pos, color = Gene_Type_combined, shape = Gene_Type),
    height = 0.06, width = 0, size = 2.8, alpha = 0.8, stroke = 0.8
  ) +
  geom_point(
    data = median_data,
    aes(x = Median_Z, y = y_pos),
    shape = 23, size = 2, fill = "white", color = "black"
  ) +
  scale_shape_manual(values = c("Functional" = 16, "Pseudogene" = 1)) +  # solid vs hollow
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = z_breaks,
    labels = z_breaks
  ) +
  scale_y_continuous(
    breaks = seq_along(group_levels_rev),
    labels = group_levels_rev,
    expand = expansion(add = 0.6)
  ) +
  theme_minimal() +
  labs(
    title = "Gene Expression (ENCODE) of Functional Genes vs Pseudogenes",
    #subtitle = "Z-scores computed against the pooled pseudogene median/MAD (vibrant = functional, pastel = pseudogene); points shown are outliers or groups with n\u22642; x-axis is pseudo-log scaled",
    x = "Z-score (relative to pooled pseudogene distribution, pseudo-log scale)",
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

# Save the single long forest-style plot. +1 row to account for the added
# "Pooled" row.
n_rows <- length(gene_group_order)
ggsave(
  filename = "../../results/ENCODE_Z_scores_all_groups_forest_plot.pdf",
  plot = forest_plot,
  width = 10,
  height = max(7, n_rows * 0.7)
)

# Export the outlier set used in the plot to CSV (per-group only,
# "Pooled" excluded to avoid duplicate rows -- see note above)
write.csv(
  outlier_points,
  file = "../../results/ENCODE_Z_score_outliers.csv",
  row.names = FALSE
)