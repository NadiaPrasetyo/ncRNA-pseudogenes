# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Data-loading step
data <- read.csv("../../data/snp_intergenic_genes_pseudogenes.csv")

# ---------------------------------------------------------------------------
# CONFIG - adjust these if your data differs from these assumptions
# ---------------------------------------------------------------------------
# Metric to plot per source, used directly (no Z-score) since "enrichment"
# is already normalized against the flanking control region.
metric_col <- "enrichment_intergenic"

# The four SNP sources present as wide-format column suffixes
snp_sources <- c("1000genomes", "gnomad", "pangenome", "dbsnp")

# ---------------------------------------------------------------------------
# Define custom colors for each Gene_group and Gene_Type combination
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Standardize key columns to match the original script's naming/casing
# ---------------------------------------------------------------------------
data <- data %>%
  mutate(
    Gene_group = gene_group,
    Gene_Type = ifelse(tolower(gene_type) == "pseudogene", "Pseudogene", "Functional")
  )

# ---------------------------------------------------------------------------
# Reshape from wide (one set of columns per source) to long format
# (.value, source) split via names_pattern
# ---------------------------------------------------------------------------
long_data <- data %>%
  pivot_longer(
    cols = matches("^(snp_count|snp_density|flank_count|flank_density|enrichment|enrichment_intergenic)_"),
    names_to = c(".value", "source"),
    names_pattern = paste0("(snp_count|snp_density|flank_count|flank_density|enrichment|enrichment_intergenic)_(",
                           paste(snp_sources, collapse = "|"), ")")
  )

# Clean the data: Remove rows with non-finite values for the chosen metric
clean_data <- long_data %>%
  filter(is.finite(.data[[metric_col]]))

# ---------------------------------------------------------------------------
# No Z-score normalization: "enrichment_intergenic" is already normalized
# against the flanking control region, so it's plotted directly.
# ---------------------------------------------------------------------------
normalized_data <- clean_data %>%
  mutate(
    Gene_Type_combined = paste(Gene_group, Gene_Type, sep = "."),
    source = factor(source, levels = snp_sources)
  )

# ---------------------------------------------------------------------------
# Common theme block, matching the original script's styling
# ---------------------------------------------------------------------------
z_theme <- theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 24),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18, angle = 25, hjust = 1),
    axis.text.y = element_text(size = 18),
    strip.text = element_text(size = 18, face = "bold"),
    plot.caption = element_text(size = 12, face = "italic", hjust = 0),
    text = element_text(family = "serif")
  )

# Define groups (same groupings as the original script)
combined_groups <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6")
combined_groups_2 <- c("RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12")
remaining_ncRNAs <- setdiff(unique(normalized_data$Gene_group), c(combined_groups, combined_groups_2))

# Non-faceted axis text angle (used since each PDF now covers a single source)
single_source_theme <- z_theme +
  theme(axis.text.x = element_text(size = 24, angle = 25, hjust = 1))

# ---------------------------------------------------------------------------
# Shared y-axis view per plot type, computed ACROSS all sources, so the
# four "Major Spliceosomal" PDFs (one per source) share identical axis
# limits/ticks and are directly comparable, and likewise for each other
# plot type.
#
# A handful of extreme outliers (enrichment_intergenic > 400 in some groups)
# were stretching the axis so badly that the bulk of the data (mostly
# 0.5-2) collapsed into a thin band near the bottom. Instead of using the
# full min/max range, the axis is now "zoomed" to a quantile-based cap that
# comfortably covers the typical values, via coord_cartesian() rather than
# scale_y_continuous(limits = ...). This only changes what's VISIBLE -
# every point is still included in the underlying data, jitter, and the
# median/mean summary stats; points above the cap are simply clipped off
# the top of the panel, and the plot caption reports how many were clipped
# and what the true max was, so nothing is hidden silently.
#
# CAP_QUANTILE controls how tight the zoom is - e.g. 0.97 sets the visible
# upper bound just above the 97th percentile of that group's data. Lower it
# (e.g. 0.90) for an even tighter zoom, or raise it if too many points are
# getting clipped.
# ---------------------------------------------------------------------------
CAP_QUANTILE <- 0.97

make_y_scale <- function(values, n_breaks = 6, cap_quantile = CAP_QUANTILE) {
  values <- values[is.finite(values)]
  rng <- range(values, na.rm = TRUE)
  cap <- quantile(values, cap_quantile, na.rm = TRUE)
  # Headroom above the cap so points near it aren't drawn on the panel edge
  upper <- cap * 1.15
  lower <- min(0, rng[1])
  pad <- (upper - lower) * 0.05
  ylim <- c(lower - pad, upper)
  
  n_clipped <- sum(values > upper)
  max_value <- rng[2]
  
  list(
    ylim = ylim,
    breaks = pretty(c(lower, upper), n = n_breaks),
    n_clipped = n_clipped,
    max_value = max_value
  )
}

# Builds a caption describing any points clipped off the top of the panel,
# or NULL if nothing was clipped (so labs(caption = NULL) omits it cleanly).
clip_caption <- function(scale_info) {
  if (scale_info$n_clipped == 0) return(NULL)
  sprintf(
    "%d point(s) above the axis limit not shown (max enrichment_intergenic = %.1f)",
    scale_info$n_clipped, scale_info$max_value
  )
}

y_scale_1 <- make_y_scale(normalized_data$enrichment_intergenic[normalized_data$Gene_group %in% combined_groups])
y_scale_2 <- make_y_scale(normalized_data$enrichment_intergenic[normalized_data$Gene_group %in% combined_groups_2])
y_scale_3 <- make_y_scale(normalized_data$enrichment_intergenic[normalized_data$Gene_group %in% remaining_ncRNAs])
y_scale_combined <- make_y_scale(normalized_data$enrichment_intergenic[normalized_data$Gene_group != "TRNA"])

# ---------------------------------------------------------------------------
# SUMMARY TABLE: average enrichment_intergenic per source x Gene_group
# (Functional and Pseudogene rows pooled together), plus a "pooled" row per
# source that pools ALL Gene_groups together. Columns: source, Gene_group,
# n, mean, median, sd. Example rows:
#   pangenome, RNU1,   ...
#   pangenome, pooled, ...
# ---------------------------------------------------------------------------
enrichment_by_group <- normalized_data %>%
  group_by(source, Gene_group) %>%
  summarise(
    n = n(),
    mean_enrichment = mean(enrichment_intergenic, na.rm = TRUE),
    median_enrichment = median(enrichment_intergenic, na.rm = TRUE),
    sd_enrichment = sd(enrichment_intergenic, na.rm = TRUE),
    .groups = "drop"
  )

enrichment_pooled <- normalized_data %>%
  group_by(source) %>%
  summarise(
    n = n(),
    mean_enrichment = mean(enrichment_intergenic, na.rm = TRUE),
    median_enrichment = median(enrichment_intergenic, na.rm = TRUE),
    sd_enrichment = sd(enrichment_intergenic, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Gene_group = "pooled", .after = source)

enrichment_summary <- bind_rows(enrichment_by_group, enrichment_pooled) %>%
  arrange(source, Gene_group != "pooled", Gene_group)

print(enrichment_summary, n = Inf)

write.csv(
  enrichment_summary,
  file = "../../results/SNP_Enrichment_summary_by_source_group.csv",
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# Loop over each SNP source and generate the four plots as separate PDFs
# ---------------------------------------------------------------------------
for (src in snp_sources) {
  
  source_data <- normalized_data %>% filter(source == src)
  
  # --- Plot 1: Major spliceosomal RNAs ---
  data_1 <- source_data %>%
    filter(Gene_group %in% combined_groups) %>%
    mutate(
      Gene_group = factor(Gene_group, levels = combined_groups),
      Gene_Type_label = ifelse(Gene_Type == "Pseudogene", paste(Gene_group, "(P)"), paste(Gene_group, "(F)"))
    )
  
  plot_1 <- ggplot(data_1, aes(x = Gene_Type_label, y = enrichment_intergenic, color = Gene_Type_combined)) +
    geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
    stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +
    scale_color_manual(values = custom_colors) +
    scale_y_continuous(breaks = y_scale_1$breaks) +
    coord_cartesian(ylim = y_scale_1$ylim) +
    labs(
      title = paste("SNP Enrichment of Major Spliceosomal RNAs -", src),
      x = "Gene Type",
      y = "SNP Enrichment (intergenic)",
      caption = clip_caption(y_scale_1)
    ) +
    single_source_theme
  
  ggsave(filename = paste0("../../results/SNP_Enrichment_", src, "_RNU1_RNU2_RNU4_RNU5_RNU6.pdf"),
         plot = plot_1, width = 12, height = 7)
  
  # --- Plot 2: Minor spliceosomal RNAs ---
  data_2 <- source_data %>%
    filter(Gene_group %in% combined_groups_2) %>%
    mutate(
      Gene_group = factor(Gene_group, levels = combined_groups_2),
      Gene_Type_label = ifelse(Gene_Type == "Pseudogene", paste(Gene_group, "(P)"), paste(Gene_group, "(F)"))
    )
  
  plot_2 <- ggplot(data_2, aes(x = Gene_Type_label, y = enrichment_intergenic, color = Gene_Type_combined)) +
    geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
    stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +
    scale_color_manual(values = custom_colors) +
    scale_y_continuous(breaks = y_scale_2$breaks) +
    coord_cartesian(ylim = y_scale_2$ylim) +
    labs(
      title = paste("SNP Enrichment of Minor Spliceosomal RNAs -", src),
      x = "Gene Type",
      y = "SNP Enrichment (intergenic)",
      caption = clip_caption(y_scale_2)
    ) +
    single_source_theme
  
  ggsave(filename = paste0("../../results/SNP_Enrichment_", src, "_RNU4ATAC_RNU6ATAC_RNU11_RNU12.pdf"),
         plot = plot_2, width = 12, height = 7)
  
  # --- Plot 3: Remaining ncRNAs ---
  data_3 <- source_data %>%
    filter(Gene_group %in% remaining_ncRNAs) %>%
    mutate(
      Gene_Type_label = ifelse(Gene_Type == "Pseudogene", paste(Gene_group, "(P)"), paste(Gene_group, "(F)"))
    )
  
  plot_3 <- ggplot(data_3, aes(x = Gene_Type_label, y = enrichment_intergenic, color = Gene_Type_combined)) +
    geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
    stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +
    scale_color_manual(values = custom_colors) +
    scale_y_continuous(breaks = y_scale_3$breaks) +
    coord_cartesian(ylim = y_scale_3$ylim) +
    labs(
      title = paste("SNP Enrichment of Other sncRNAs -", src),
      x = "Gene Type",
      y = "SNP Enrichment (intergenic)",
      caption = clip_caption(y_scale_3)
    ) +
    single_source_theme
  
  ggsave(filename = paste0("../../results/SNP_Enrichment_", src, "_remaining_ncRNAs.pdf"),
         plot = plot_3, width = 12, height = 7)
  
  # --- Plot 4: Combined Functional vs Pseudogene ---
  combined_gene_type_data <- source_data %>%
    filter(Gene_group != "TRNA") %>%
    mutate(
      Gene_Type_combined_simple = ifelse(Gene_Type == "Pseudogene", "Pseudogene", "Functional")
    )
  
  combined_gene_type_plot <- ggplot(combined_gene_type_data,
                                    aes(x = Gene_Type_combined_simple, y = enrichment_intergenic, color = Gene_Type_combined_simple)) +
    geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
    stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +
    scale_color_manual(values = c("Functional" = "firebrick", "Pseudogene" = "cornflowerblue")) +
    scale_y_continuous(breaks = y_scale_combined$breaks) +
    coord_cartesian(ylim = y_scale_combined$ylim) +
    labs(
      title = paste("SNP Enrichment: Pseudogenes vs Functional Genes -", src),
      x = "Gene Type",
      y = "SNP Enrichment (intergenic)",
      caption = clip_caption(y_scale_combined)
    ) +
    single_source_theme +
    theme(axis.text.x = element_text(size = 24, angle = 0, hjust = 0.5))
  
  ggsave(filename = paste0("../../results/SNP_Enrichment_", src, "_Pseudogene_vs_Functional.pdf"),
         plot = combined_gene_type_plot, width = 12, height = 7)
}

# ---------------------------------------------------------------------------
# Plot 5: Distribution by SOURCE, pooled across all gene groups and NOT
# split by Functional/Pseudogene - i.e. "what does the overall
# enrichment_intergenic distribution look like for 1000genomes vs gnomad vs
# pangenome vs dbsnp, and how do their means/medians compare?"
#
# Each source's jittered points show the spread; a triangle marks the mean
# and a diamond marks the median so the two can be compared directly (they
# diverge more for sources with heavier outlier tails).
# ---------------------------------------------------------------------------
source_colors <- c(
  "1000genomes" = "#1f77b4",
  "gnomad"      = "#d62728",
  "pangenome"   = "#2ca02c",
  "dbsnp"       = "#9467bd"
)

pooled_by_source <- normalized_data  # all Gene_groups, Functional + Pseudogene combined

source_summary_stats <- pooled_by_source %>%
  group_by(source) %>%
  summarise(
    Mean = mean(enrichment_intergenic, na.rm = TRUE),
    Median = median(enrichment_intergenic, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(Mean, Median), names_to = "stat_type", values_to = "value")

y_scale_source <- make_y_scale(pooled_by_source$enrichment_intergenic)

source_distribution_plot <- ggplot(pooled_by_source, aes(x = source, y = enrichment_intergenic)) +
  geom_jitter(aes(color = source), width = 0.25, height = 0, size = 2, alpha = 0.35) +
  geom_point(
    data = source_summary_stats,
    aes(x = source, y = value, shape = stat_type),
    size = 5, color = "black", fill = "white", stroke = 1.2
  ) +
  scale_color_manual(values = source_colors, guide = "none") +
  scale_shape_manual(name = "Statistic", values = c("Mean" = 24, "Median" = 23)) +
  scale_y_continuous(breaks = y_scale_source$breaks) +
  coord_cartesian(ylim = y_scale_source$ylim) +
  labs(
    title = "SNP Enrichment Distribution by Source (all gene groups pooled)",
    x = "SNP Source",
    y = "SNP Enrichment (intergenic)",
    caption = clip_caption(y_scale_source)
  ) +
  z_theme +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5)
  )

ggsave(filename = "../../results/SNP_Enrichment_distribution_by_source.pdf",
       plot = source_distribution_plot, width = 12, height = 7)