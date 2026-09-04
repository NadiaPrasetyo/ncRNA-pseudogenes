# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Data-loading step
data <- read.csv("../../data/snp_intergenic_genes_pseudogenes.csv")
sources <- "gnomad"

# Filter out chrY and chrMT genes: incomplete/non-representative coverage
chrom_col_candidates <- c("chr", "chrom", "chromosome", "Chr", "Chrom", "Chromosome", "CHR")
chrom_col <- intersect(chrom_col_candidates, names(data))

if (length(chrom_col) == 0) {
  stop("Could not find a chromosome column to filter on. Columns present: ",
       paste(names(data), collapse = ", "))
}
chrom_col <- chrom_col[1]

n_before <- nrow(data)
excluded_chroms <- data %>% filter(grepl("^(chr)?(Y|MT|M)$", .data[[chrom_col]], ignore.case = TRUE))
data <- data %>% filter(!grepl("^(chr)?(Y|MT|M)$", .data[[chrom_col]], ignore.case = TRUE))

message(sprintf("Filtered out %d of %d rows on chrY/chrMT (column: '%s'); %d rows remain.",
                nrow(excluded_chroms), n_before, chrom_col, nrow(data)))

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

# Clean the data: Remove rows with non-finite values for the chosen metric.
# NOTE: this deliberately keeps rows where enrichment_intergenic is 0 (or
# negative) - those are real data points, just ones a log axis can't
# display at their true position. They get handled per-plot below.
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
# IMPORTANT: no global flooring of zero/negative values happens here anymore.
#
# Previously, every enrichment_intergenic <= 0 was floored to a fixed
# LOG_FLOOR = 0.02 *before* any axis limits were computed. That's fragile:
# the axis lower bound is computed separately per plot (from real positive
# values), so a fixed 0.02 could land above the true axis floor, below it,
# or in some arbitrary spot in the middle of the panel depending on the
# data - it was never guaranteed to actually be "at the bottom".
#
# Now: enrichment_intergenic is left untouched here (zeros/negatives stay
# as-is), and every stats table (enrichment_summary, source_summary_stats)
# is computed on these true raw values. Only at the point of *plotting*
# each panel do we pin non-positive values to that panel's own computed
# axis floor (see `floor_for_plot()` below) - so they always sit exactly
# at the bottom of whichever axis they end up on, no matter what that
# bottom number turns out to be.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Shared y-axis view per plot type, computed ACROSS all sources, so the
# four "Major Spliceosomal" PDFs (one per source) share identical axis
# limits/ticks and are directly comparable, and likewise for each other
# plot type.
#
# The axis is now a TRUE log10 scale (via scale_y_continuous(trans =
# "log10"))
#
# The padding multiplier below (PAD_FACTOR) was widened from 1.15 to 1.3 so
# There's more headroom/footroom around the real data on the log axis. this is
# the "spread it out" fix requested: extra vertical log-space around the
# points keeps them from being squashed against the panel edges.
#
# CAP_QUANTILE controls how tight the zoom is - e.g. 0.97 sets the visible
# upper bound just above the 97th percentile of that group's data. Lower it
# (e.g. 0.90) for an even tighter zoom, or raise it if too many points are
# getting clipped.
# ---------------------------------------------------------------------------
CAP_QUANTILE <- 1
PAD_FACTOR <- 1.3

# Physical PDF height in inches, used by every ggsave() call below. Bump
# this up if labels still feel tight - taller panels give the (sparser)
# breaks defined below more vertical room between them.
PLOT_HEIGHT <- 7

# Which multiples-of-each-decade get a labeled tick. 9 ticks/decade (the
# old 1:9 pattern) is what caused the overlap: on a log axis 2-3-4...9 all
# bunch up near the top of each decade, so their labels collide regardless
# of plot height. c(1, 2, 5) gives 3 evenly-legible ticks per decade
# instead. Use 1 alone (just powers of 10) for an even sparser axis.
BREAK_MULTIPLES <- c(1, 2, 5)


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
# make_y_scale(): computes axis bounds/breaks from the REAL (strictly
# positive) values only - zero/negative values never influence where the
# axis lower bound sits. The resulting `floor_value` (== the axis lower
# bound, `ylim[1]`) is what non-positive points get pinned to at plot time,
# so they always land exactly at the bottom of the panel - whatever that
# bottom number happens to be for this particular group of data. Clipped-
# point counts for the caption are handled separately by
# panel_clip_info(), computed per-panel rather than baked in here (see
# note below) so a single-source panel's caption reflects that panel's
# own data, not whatever pooled-across-sources data built the shared axis.
# ---------------------------------------------------------------------------
make_y_scale <- function(values, cap_quantile = CAP_QUANTILE,
                         display_floor = 0.01, pad_factor = PAD_FACTOR,
                         break_multiples = BREAK_MULTIPLES) {
  values <- values[is.finite(values)]
  real_values <- values[values > 0]
  
  cap <- quantile(real_values, cap_quantile, na.rm = TRUE)
  upper <- cap * pad_factor
  
  lower <- max(display_floor, min(real_values, na.rm = TRUE) / pad_factor)
  
  mag_lo <- floor(log10(lower))
  mag_hi <- ceiling(log10(upper))
  
  # Labeled breaks: sparse (e.g. 1, 2, 5 per decade) so labels don't overlap.
  candidate_breaks <- unlist(lapply(mag_lo:mag_hi, function(m) break_multiples * 10^m))
  breaks <- sort(unique(candidate_breaks[candidate_breaks >= lower & candidate_breaks <= upper]))
  
  # Minor breaks: the full 1:9 grid, unlabeled - keeps the visual sense of
  # log spacing (via faint gridlines) without cluttering the axis text.
  minor_candidate <- unlist(lapply(mag_lo:mag_hi, function(m) (1:9) * 10^m))
  minor_breaks <- sort(unique(minor_candidate[minor_candidate >= lower & minor_candidate <= upper]))
  
  list(
    ylim = c(lower, upper),
    breaks = breaks,
    minor_breaks = minor_breaks,
    floor_value = lower     # <-- where non-positive points get pinned
  )
}

# NOTE: make_y_scale() no longer tracks clipped-point counts itself. Several
# of the shared scales above (y_scale_1, y_scale_2, y_scale_3,
# y_scale_combined) are intentionally computed by pooling values ACROSS ALL
# FOUR raw sources (so per-source PDFs share identical axis limits), but
# each individual PDF only ever plots ONE source's subset of that data. If
# the clipped/pinned-point COUNT were baked into the shared scale, the
# caption on a single-source panel would report how many points were
# pinned across all sources combined - not how many are actually pinned in
# the panel you're looking at, which is misleading whenever the sources
# don't have the same number of zero/near-zero rows. `panel_clip_info()`
# below fixes that by counting clipped points fresh from whatever data
# frame is actually being plotted.
panel_clip_info <- function(df, scale_info, col = "enrichment_intergenic") {
  values <- df[[col]]
  values <- values[is.finite(values)]
  real_values <- values[values > 0]
  list(
    n_clipped_low = sum(values <= 0),
    n_clipped_high = sum(real_values > scale_info$ylim[2]),
    max_value = if (length(real_values) > 0) max(real_values) else NA_real_
  )
}

# Shared log10 y-scale layer, built once and reused on every plot so the
# transform/label formatting stays consistent.
log_y_scale <- function(scale_info) {
  scale_y_continuous(
    trans = "log10",
    breaks = scale_info$breaks,
    minor_breaks = scale_info$minor_breaks,
    labels = scales::label_number(accuracy = 0.01, drop0trailing = TRUE)
  )
}

# Pins zero/negative values in a data frame's metric column to this panel's
# own computed axis floor (scale_info$floor_value), so they always render
# exactly at the bottom of the axis no matter what that number is. Positive
# values pass through unchanged. Returns the data frame with an added
# `<col>_plot` column meant for use as the plotting aesthetic - the
# original column is left intact for any downstream stats.
floor_for_plot <- function(df, scale_info, col = "enrichment_intergenic") {
  plot_col <- paste0(col, "_plot")
  df[[plot_col]] <- ifelse(is.finite(df[[col]]) & df[[col]] <= 0,
                           scale_info$floor_value,
                           df[[col]])
  df
}

# Builds a caption describing any points clipped off the top of the panel,
# or pinned to the bottom because they were zero/negative. `clip_info` must
# come from panel_clip_info() run on the SAME data frame being plotted, so
# the counts always match what's actually visible in that panel - not some
# other (e.g. all-sources-pooled) group the axis limits were derived from.
# Returns NULL if nothing needed adjusting (so labs(caption = NULL) omits
# it cleanly).
clip_caption <- function(clip_info, scale_info) {
  parts <- c()
  if (clip_info$n_clipped_high > 0) {
    parts <- c(parts, sprintf("%d point(s) above the axis limit not shown (max = %.1f)",
                              clip_info$n_clipped_high, clip_info$max_value))
  }
  if (clip_info$n_clipped_low > 0) {
    parts <- c(parts, sprintf("%d point(s) at/below zero pinned to the axis floor (%.3g), not shown at true scale",
                              clip_info$n_clipped_low, scale_info$floor_value))
  }
  if (length(parts) == 0) return(NULL)
  paste(parts, collapse = "; ")
}

# Scales computed from TRUE raw values (zeros/negatives excluded from the
# bound calculation, but counted for the caption).
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
#
# These are computed on the TRUE raw enrichment_intergenic values (zeros
# included as zero) - not on any plot-floor substitute - so the numbers
# reflect the real data regardless of how any individual panel's axis ends
# up being drawn.
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

# Stats computed on true raw values (real zeros included, not floored).
source_summary_stats <- pooled_by_source %>%
  group_by(source) %>%
  summarise(
    Mean = mean(enrichment_intergenic, na.rm = TRUE),
    Median = median(enrichment_intergenic, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(Mean, Median), names_to = "stat_type", values_to = "value")

y_scale_source <- make_y_scale(pooled_by_source$enrichment_intergenic)
clip_info_source <- panel_clip_info(pooled_by_source, y_scale_source)

# Pin zero/negative points to this panel's own axis floor for display only.
pooled_by_source <- floor_for_plot(pooled_by_source, y_scale_source)

source_distribution_plot <- ggplot(pooled_by_source, aes(x = source, y = enrichment_intergenic_plot)) +
  geom_jitter(aes(color = source), width = 0.25, height = 0, size = 2, alpha = 0.35) +
  geom_point(
    data = source_summary_stats,
    aes(x = source, y = value, shape = stat_type),
    size = 5, color = "black", fill = "white", stroke = 1.2
  ) +
  scale_color_manual(values = source_colors, guide = "none") +
  scale_shape_manual(name = "Statistic", values = c("Mean" = 24, "Median" = 23)) +
  log_y_scale(y_scale_source) +
  coord_cartesian(ylim = y_scale_source$ylim) +
  labs(
    title = "SNP Enrichment Distribution by Source (all gene groups pooled)",
    x = "SNP Source",
    y = "SNP Enrichment (intergenic, log10 scale)",
    caption = clip_caption(clip_info_source, y_scale_source)
  ) +
  z_theme +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5)
  )

ggsave(filename = "../../results/SNP_Enrichment_distribution_by_source.pdf",
       plot = source_distribution_plot, width = 12, height = PLOT_HEIGHT)


# ---------------------------------------------------------------------------
# Restrict everything downstream to the source(s) defined in `sources`
# (e.g. "gnomad"). snp_sources is still used above for the pivot_longer
# regex since the raw data has columns for all four sources, but nothing
# past this point will compute or plot the others.
# ---------------------------------------------------------------------------
normalized_data <- normalized_data %>%
  filter(source %in% sources) %>%
  mutate(source = factor(source, levels = sources))

# ---------------------------------------------------------------------------
# Loop over each SNP source and generate the four plots as separate PDFs
# ---------------------------------------------------------------------------
for (src in sources) {
  
  source_data <- normalized_data %>% filter(source == src)
  
  # --- Plot 1: Major spliceosomal RNAs ---
  data_1 <- source_data %>%
    filter(Gene_group %in% combined_groups) %>%
    mutate(
      Gene_group = factor(Gene_group, levels = combined_groups),
      Gene_Type_label = ifelse(Gene_Type == "Pseudogene", paste(Gene_group, "(P)"), paste(Gene_group, "(F)"))
    )
  clip_info_1 <- panel_clip_info(data_1, y_scale_1)
  data_1 <- data_1 %>% floor_for_plot(y_scale_1)
  
  plot_1 <- ggplot(data_1, aes(x = Gene_Type_label, y = enrichment_intergenic_plot, color = Gene_Type_combined)) +
    geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
    stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +
    scale_color_manual(values = custom_colors) +
    log_y_scale(y_scale_1) +
    coord_cartesian(ylim = y_scale_1$ylim) +
    labs(
      title = paste("SNP Enrichment of Major Spliceosomal RNAs -", src),
      x = "Gene Type",
      y = "SNP Enrichment (intergenic, log10 scale)",
      caption = clip_caption(clip_info_1, y_scale_1)
    ) +
    single_source_theme
  
  ggsave(filename = paste0("../../results/SNP_Enrichment_", src, "_RNU1_RNU2_RNU4_RNU5_RNU6.pdf"),
         plot = plot_1, width = 12, height = PLOT_HEIGHT)
  
  # --- Plot 2: Minor spliceosomal RNAs ---
  data_2 <- source_data %>%
    filter(Gene_group %in% combined_groups_2) %>%
    mutate(
      Gene_group = factor(Gene_group, levels = combined_groups_2),
      Gene_Type_label = ifelse(Gene_Type == "Pseudogene", paste(Gene_group, "(P)"), paste(Gene_group, "(F)"))
    )
  clip_info_2 <- panel_clip_info(data_2, y_scale_2)
  data_2 <- data_2 %>% floor_for_plot(y_scale_2)
  
  plot_2 <- ggplot(data_2, aes(x = Gene_Type_label, y = enrichment_intergenic_plot, color = Gene_Type_combined)) +
    geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
    stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +
    scale_color_manual(values = custom_colors) +
    log_y_scale(y_scale_2) +
    coord_cartesian(ylim = y_scale_2$ylim) +
    labs(
      title = paste("SNP Enrichment of Minor Spliceosomal RNAs -", src),
      x = "Gene Type",
      y = "SNP Enrichment (intergenic, log10 scale)",
      caption = clip_caption(clip_info_2, y_scale_2)
    ) +
    single_source_theme
  
  ggsave(filename = paste0("../../results/SNP_Enrichment_", src, "_RNU4ATAC_RNU6ATAC_RNU11_RNU12.pdf"),
         plot = plot_2, width = 12, height = PLOT_HEIGHT)
  
  # --- Plot 3: Remaining ncRNAs ---
  data_3 <- source_data %>%
    filter(Gene_group %in% remaining_ncRNAs) %>%
    mutate(
      Gene_Type_label = ifelse(Gene_Type == "Pseudogene", paste(Gene_group, "(P)"), paste(Gene_group, "(F)"))
    )
  clip_info_3 <- panel_clip_info(data_3, y_scale_3)
  data_3 <- data_3 %>% floor_for_plot(y_scale_3)
  
  plot_3 <- ggplot(data_3, aes(x = Gene_Type_label, y = enrichment_intergenic_plot, color = Gene_Type_combined)) +
    geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
    stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +
    scale_color_manual(values = custom_colors) +
    log_y_scale(y_scale_3) +
    coord_cartesian(ylim = y_scale_3$ylim) +
    labs(
      title = paste("SNP Enrichment of Other sncRNAs -", src),
      x = "Gene Type",
      y = "SNP Enrichment (intergenic, log10 scale)",
      caption = clip_caption(clip_info_3, y_scale_3)
    ) +
    single_source_theme
  
  ggsave(filename = paste0("../../results/SNP_Enrichment_", src, "_remaining_ncRNAs.pdf"),
         plot = plot_3, width = 12, height = PLOT_HEIGHT)
  
  # --- Plot 4: Combined Functional vs Pseudogene ---
  combined_gene_type_data <- source_data %>%
    filter(Gene_group != "TRNA") %>%
    mutate(
      Gene_Type_combined_simple = ifelse(Gene_Type == "Pseudogene", "Pseudogene", "Functional")
    )
  clip_info_combined <- panel_clip_info(combined_gene_type_data, y_scale_combined)
  combined_gene_type_data <- combined_gene_type_data %>% floor_for_plot(y_scale_combined)
  
  combined_gene_type_plot <- ggplot(combined_gene_type_data,
                                    aes(x = Gene_Type_combined_simple, y = enrichment_intergenic_plot, color = Gene_Type_combined_simple)) +
    geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
    stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +
    scale_color_manual(values = c("Functional" = "firebrick", "Pseudogene" = "cornflowerblue")) +
    log_y_scale(y_scale_combined) +
    coord_cartesian(ylim = y_scale_combined$ylim) +
    labs(
      title = paste("SNP Enrichment: Pseudogenes vs Functional Genes -", src),
      x = "Gene Type",
      y = "SNP Enrichment (intergenic, log10 scale)",
      caption = clip_caption(clip_info_combined, y_scale_combined)
    ) +
    single_source_theme +
    theme(axis.text.x = element_text(size = 24, angle = 0, hjust = 0.5))
  
  ggsave(filename = paste0("../../results/SNP_Enrichment_", src, "_Pseudogene_vs_Functional.pdf"),
         plot = combined_gene_type_plot, width = 12, height = PLOT_HEIGHT)
}