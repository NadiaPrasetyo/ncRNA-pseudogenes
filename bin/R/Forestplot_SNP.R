# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Data-loading step
data <- read.csv("../../results/combined_gene_data_with_snp.csv")
sources <- c("gnomad", "1000genomes", "pangenome", "dbsnp")

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
# NOTE: no Z-score here. "enrichment_intergenic" is already normalized
# against its own flanking control region, so it's plotted directly on the
# x-axis rather than converted to a Z-score first. The x-axis itself uses a
# log2 (pseudo-log) transform, see the "X-AXIS SCALING" note below.
metric_col  <- "enrichment_intergenic"
snp_sources <- c("1000genomes", "gnomad", "pangenome", "dbsnp")

# ---------------------------------------------------------------------------
# Vibrant (Functional) / pastel (Pseudogene) color scheme, same convention
# used in the PhyloP100 and ENCODE forest plots.
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Standardize key columns
# ---------------------------------------------------------------------------
# UPDATED: the new CSV already provides "Gene_Type" and "Gene_group" in
# title case (rather than lowercase "gene_type"/"gene_group"), so there's
# no separate Gene_group alias to create anymore - just normalize the
# Functional/Pseudogene labels in place.
data <- data %>%
  mutate(
    Gene_Type = ifelse(tolower(Gene_Type) == "pseudogene", "Pseudogene", "Functional")
  )

# ---------------------------------------------------------------------------
# Reshape from wide (one set of columns per source) to long format
# ---------------------------------------------------------------------------
# UPDATED: the new CSV only carries "snp_density_*" and
# "enrichment_intergenic_*" per source (no snp_count/flank_count/
# flank_density/plain "enrichment" columns), so the pivot only needs to
# match those two value prefixes.
long_data <- data %>%
  pivot_longer(
    cols = matches("^(snp_density|enrichment_intergenic)_"),
    names_to = c(".value", "source"),
    names_pattern = paste0("(snp_density|enrichment_intergenic)_(",
                           paste(snp_sources, collapse = "|"), ")")
  )

# Keep zero/negative enrichment values (real data points) - they're handled
# at plot time via a pseudo-log2 x-axis, which stays defined through zero
# (and negative values, if any slip through) while behaving like a log2
# scale everywhere away from zero. See "X-AXIS SCALING" below.
clean_data <- long_data %>%
  filter(is.finite(.data[[metric_col]]))

normalized_data <- clean_data %>%
  mutate(
    Gene_Type_combined = paste(Gene_group, Gene_Type, sep = "."),
    source = factor(source, levels = snp_sources)
  ) %>%
  filter(source %in% sources) %>%
  mutate(source = factor(source, levels = sources))

# ---------------------------------------------------------------------------
# Same gene-group ordering used in the PhyloP100/ENCODE forest plots. A
# "Pooled" pseudo-group is appended LAST so it renders at the BOTTOM of
# each forest plot.
# ---------------------------------------------------------------------------
combined_groups   <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6")
combined_groups_2 <- c("RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12")
remaining_ncRNAs  <- setdiff(unique(normalized_data$Gene_group), c(combined_groups, combined_groups_2))

gene_group_order <- c(combined_groups, combined_groups_2, remaining_ncRNAs, "Pooled")

# NEW: build the "Pooled" rows by duplicating every existing row (all
# Gene_groups, all sources) with Gene_group relabeled to "Pooled". This
# reuses the SAME enrichment_intergenic values already in the data --
# nothing is recalculated or rescaled -- it just lets each per-source
# forest plot draw one more Functional/Pseudogene comparison pooling
# every gene group together. Since downstream steps (IQR outlier bounds,
# medians, thin-group detection) all group by source x Gene_Type_combined,
# "Pooled" is treated as its own group and gets its own IQR/median
# computed from every pooled point, same as any other group.
pooled_data <- normalized_data %>%
  mutate(Gene_group = "Pooled")

normalized_data <- bind_rows(normalized_data, pooled_data) %>%
  filter(Gene_group %in% gene_group_order) %>%
  mutate(
    Gene_group = factor(Gene_group, levels = gene_group_order),
    Gene_Type  = factor(Gene_Type, levels = c("Functional", "Pseudogene")),
    Gene_Type_combined = paste(Gene_group, Gene_Type, sep = ".")
  )

group_levels_rev <- rev(gene_group_order)   # first group ends up at the top of the plot
normalized_data <- normalized_data %>%
  mutate(
    Gene_group_f = factor(Gene_group, levels = group_levels_rev),
    y_num    = as.numeric(Gene_group_f),
    y_offset = ifelse(Gene_Type == "Functional", 0.15, -0.15),
    y_pos    = y_num + y_offset
  )

# ---------------------------------------------------------------------------
# Outlier definition (no Z-score available here, so we fall back to the
# standard boxplot rule): within each Gene_group x Gene_Type x source, a
# point is an outlier if it falls beyond 1.5x the IQR from the nearest
# quartile. These, plus any thin (n<=2) groups, are the only individual
# points drawn on top of the violins - solid = functional, hollow =
# pseudogene, matching the convention from the other forest plots.
# ---------------------------------------------------------------------------
iqr_bounds <- normalized_data %>%
  group_by(source, Gene_Type_combined) %>%
  summarise(
    q1  = quantile(.data[[metric_col]], 0.25, na.rm = TRUE),
    q3  = quantile(.data[[metric_col]], 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    .groups = "drop"
  ) %>%
  mutate(
    lower_bound = q1 - 1.5 * iqr,
    upper_bound = q3 + 1.5 * iqr
  )

normalized_data <- normalized_data %>%
  left_join(iqr_bounds, by = c("source", "Gene_Type_combined")) %>%
  mutate(
    is_outlier = .data[[metric_col]] < lower_bound | .data[[metric_col]] > upper_bound
  )

# Precompute one median enrichment value per source x Gene_group x
# Gene_Type, placed at the same fixed y_offset used for that type's dodge.
median_data <- normalized_data %>%
  group_by(source, Gene_group, Gene_Type, Gene_Type_combined) %>%
  summarize(
    Median_val = median(.data[[metric_col]], na.rm = TRUE),
    y_pos = first(y_num) + ifelse(first(Gene_Type) == "Functional", 0.15, -0.15),
    .groups = "drop"
  )

# Groups with too few points (1-2) per source to draw a meaningful violin
group_counts <- normalized_data %>%
  count(source, Gene_Type_combined, name = "n_points")

thin_lookup <- group_counts %>% filter(n_points <= 2)

z_theme <- theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 13),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    text = element_text(family = "serif")
  )

n_rows <- length(gene_group_order)
all_outliers <- list()

# ---------------------------------------------------------------------------
# X-AXIS SCALING (log2 / pseudo-log)
# ---------------------------------------------------------------------------
# UPDATED: switched from a linear x-axis to a base-2 pseudo-log transform
# so the 0-1 range (depletion) gets as much visual room as the >1 range
# (enrichment) - e.g. the gap between 0.25 and 0.5 now reads the same
# width as the gap between 2 and 4. A true log2 is undefined at 0 (and
# negative values), so we use scales::pseudo_log_trans(base = 2), which
# behaves like log2 away from zero but stays smoothly defined through
# zero for any 0/negative enrichment points.
#
# UPDATED AGAIN: view is now centered on 1 (no enrichment/depletion), with
# ticks at powers of 2 out to 4x in either direction (0.25, 0.5, 1, 2, 4).
# In log2 space these are symmetric (-2, -1, 0, 1, 2), so centering the
# view on this range puts the "1" gridline in the middle of the plot.
# x_view_range sets that symmetric window via coord_cartesian() below -
# it only zooms the display; points/violins beyond it are still part of
# the underlying data (medians, IQR outlier flags, etc. are unaffected).
#
# UPDATED AGAIN: with the default sigma = 1, pseudo_log_trans is still in
# its "linear near zero" regime at x = 0.25-0.5 (sigma sets the width of
# that linear region), so those ticks were visibly compressed relative to
# 2-4 even though they're the same log2 distance from 1. Shrinking sigma
# to 0.001 pushes the linear region down near true zero, so everything at
# 0.25 and above behaves like an honest log2 scale (evenly spaced ticks),
# while 0 (and any stray negative values) still transform safely instead
# of producing -Inf/NaN the way a plain log2 scale would.
x_trans      <- scales::pseudo_log_trans(base = 2, sigma = 0.1)
x_breaks     <- c(0.25, 0.5, 1, 2, 4)
x_labels     <- x_breaks
x_view_range <- c(0, 6)

# ---------------------------------------------------------------------------
# One long forest-style plot per SNP source
# ---------------------------------------------------------------------------
for (src in sources) {
  
  src_data <- normalized_data %>% filter(source == src)
  
  thin_groups_src <- thin_lookup %>% filter(source == src) %>% pull(Gene_Type_combined)
  
  violin_data <- src_data %>% filter(!(Gene_Type_combined %in% thin_groups_src) & !is_outlier)
  thin_data   <- src_data %>% filter(Gene_Type_combined %in% thin_groups_src)
  outlier_points <- src_data %>% filter(!(Gene_Type_combined %in% thin_groups_src) & is_outlier)
  
  all_outliers[[src]] <- outlier_points
  
  median_data_src <- median_data %>% filter(source == src)
  
  if (src == "dbsnp"){
    src <- "dbSNP155"
    print("changing dbsnp to dbSNP155")
  }
  
  forest_plot <- ggplot() +
    geom_vline(xintercept = 1, linetype = "solid", color = "grey40", linewidth = 0.4) +  # 1 = no enrichment/depletion
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.4) +
    # Separator line between the per-group rows and the Pooled summary row
    geom_hline(yintercept = 1.5, linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_violin(
      data = violin_data,
      aes(x = .data[[metric_col]], y = y_pos, group = Gene_Type_combined,
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
      aes(x = .data[[metric_col]], y = y_pos, color = Gene_Type_combined, shape = Gene_Type),
      height = 0.06, width = 0, size = 2.8, alpha = 0.9, stroke = 0.8
    ) +
    geom_jitter(
      data = outlier_points,
      aes(x = .data[[metric_col]], y = y_pos, color = Gene_Type_combined, shape = Gene_Type),
      height = 0.06, width = 0, size = 2.8, alpha = 0.9, stroke = 0.8
    ) +
    geom_point(
      data = median_data_src,
      aes(x = Median_val, y = y_pos),
      shape = 23, size = 2, fill = "white", color = "black"
    ) +
    scale_shape_manual(values = c("Functional" = 16, "Pseudogene" = 1)) +  # solid vs hollow
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors) +
    scale_x_continuous(
      trans  = x_trans,
      breaks = x_breaks,
      labels = x_labels
    ) +
    coord_cartesian(xlim = x_view_range) +
    scale_y_continuous(
      breaks = seq_along(group_levels_rev),
      labels = group_levels_rev,
      expand = expansion(add = 0.6)
    ) +
    labs(
      title = paste0("SNP Enrichment (", src, ") of Functional Genes vs Pseudogenes"),
      #subtitle = "Vibrant = functional, pastel = pseudogene; points shown are IQR outliers or groups with n\u22642; solid vline at 1 = no enrichment/depletion",
      x = "SNP Enrichment (pseudo-log2 scale)",
      y = NULL
    ) +
    z_theme
  
  ggsave(
    filename = paste0("../../results/SNP_Enrichment_forest_plot_", src, ".pdf"),
    plot = forest_plot,
    width = 10,
    height = max(7, n_rows * 0.7)
  )
}

# ---------------------------------------------------------------------------
# Export outliers to a single CSV, restricted to dbSNP only. The "Pooled"
# pseudo-group is excluded here since its points are duplicates of the
# per-group rows above it and would otherwise double-count every outlier.
# ---------------------------------------------------------------------------
outliers_df <- bind_rows(all_outliers) %>%
  filter(source == "dbsnp") %>%
  filter(Gene_group != "Pooled") %>%
  select(-lower_bound, -upper_bound) %>%
  arrange(Gene_group, Gene_Type, desc(.data[[metric_col]]))

write.csv(
  outliers_df,
  file = "../../results/SNP_Enrichment_outliers_dbsnp.csv",
  row.names = FALSE
)