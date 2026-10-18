#!/usr/bin/env Rscript
# =============================================================================
# plot_motif_distributions.R
#
# Violin / "forest" style visualization of type1/type2/type3 motif counts and
# mean scores, split by gene_type (functional vs pseudogene). Motif count
# (top) and score (bottom) are stacked and share one x axis.
#
# Layout (matches the combined SNP / ENCODE / PhyloP figure):
#   - Promoter type (Type 1 / 2 / 3) on the X axis, value (count / score) on
#     the Y axis.
#   - Within each promoter type, Functional (F) is dodged left and
#     Pseudogene (P) is dodged right, with small "F" / "P" letters under each
#     violin and the promoter-type label beneath them.
#   - Violin outlines coloured by promoter type; light fill of the same hue.
#   - Only OUTLIER points (1.5*IQR rule) are drawn individually.
#   - Functional = solid dot (shape 16) + vibrant colour,
#     Pseudogene = hollow dot (shape 1) + pastel version of the SAME hue
#     (Type 1 = blue, Type 2 = red, Type 3 = yellow).
#   - No vertical gridlines.
#   - Motif counts and scores are side-by-side subplots (facets).
#
# Input : results/mcast_summary_grouped_annotated_with_snp.csv
# Output: PNG plots written to results/plots/
#
# Genes with no promoter of a given type (NA count/score) are treated as
# having a motif count of 0 and a score of 0 for that type.
# =============================================================================

## ---- 0. Setup -------------------------------------------------------------

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "readr", "stringr", "scales")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(scales)

input_csv  <- "../../results/mcast_summary_grouped_annotated_with_snp.csv"
output_dir <- "../../results/plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load data -----------------------------------------------------

df <- read_csv(input_csv, show_col_types = FALSE)

if (!"gene_type" %in% names(df)) {
  stop("Expected a 'gene_type' column (e.g. 'functional' / 'pseudogene') in the input CSV.")
}

## ---- 2. Clean / impute missing values ----------------------------------

df_clean <- df %>%
  mutate(
    type1_motif_counts = ifelse(is.na(type1_motif_counts), 0, type1_motif_counts),
    type2_motif_counts = ifelse(is.na(type2_motif_counts), 0, type2_motif_counts),
    type3_motif_counts = ifelse(is.na(type3_motif_counts), 0, type3_motif_counts),
    type1_score        = ifelse(is.na(type1_score), 0, type1_score),
    type2_score        = ifelse(is.na(type2_score), 0, type2_score),
    type3_score        = ifelse(is.na(type3_score), 0, type3_score)
  )

df_clean <- df_clean %>%
  mutate(gene_type = str_to_title(str_trim(gene_type)))

## ---- 3. Reshape to long format (counts + scores combined) -----------------

counts_long <- df_clean %>%
  select(gene_name, gene_type, type1_motif_counts, type2_motif_counts, type3_motif_counts) %>%
  pivot_longer(
    cols = c(type1_motif_counts, type2_motif_counts, type3_motif_counts),
    names_to = "promoter_type",
    values_to = "value"
  ) %>%
  mutate(
    promoter_type = recode(promoter_type,
                           type1_motif_counts = "Type 1",
                           type2_motif_counts = "Type 2",
                           type3_motif_counts = "Type 3"),
    metric = "Motif count"
  )

scores_long <- df_clean %>%
  select(gene_name, gene_type, type1_score, type2_score, type3_score) %>%
  pivot_longer(
    cols = c(type1_score, type2_score, type3_score),
    names_to = "promoter_type",
    values_to = "value"
  ) %>%
  mutate(
    promoter_type = recode(promoter_type,
                           type1_score = "Type 1",
                           type2_score = "Type 2",
                           type3_score = "Type 3"),
    metric = "Score"
  )

promoter_order <- c("Type 1", "Type 2", "Type 3")
n_types <- length(promoter_order)

combined_long <- bind_rows(counts_long, scores_long) %>%
  mutate(
    metric        = factor(metric, levels = c("Motif count", "Score")),
    gene_type     = factor(gene_type, levels = c("Functional", "Pseudogene")),
    promoter_type = factor(promoter_type, levels = promoter_order)
  )

## ---- 3b. Robust Z-scores vs. pooled pseudogene controls ------------------
# x* = pseudogene values pooled across all promoter types, computed
# separately for each metric (Motif count, Score), so every promoter type is
# scaled consistently against the same control distribution.
#   MAD(x*) >  1e-6:  Z = (x - median(x*)) / (1.4826 * MAD(x*))
#   MAD(x*) <= 1e-6:  Z = (x - median(x*)) / (1.2533 * MnAD(x*))
# (MnAD = mean absolute deviation about the median)

z_params <- combined_long %>%
  group_by(metric) %>%
  summarise(
    ctrl_median = median(value[gene_type == "Pseudogene"], na.rm = TRUE),
    ctrl_mad    = mad(value[gene_type == "Pseudogene"], constant = 1, na.rm = TRUE),
    ctrl_mnad   = mean(abs(value[gene_type == "Pseudogene"] -
                             median(value[gene_type == "Pseudogene"], na.rm = TRUE)),
                       na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    method       = ifelse(ctrl_mad > 1e-6, "MAD", "MnAD"),
    scale_factor = ifelse(ctrl_mad > 1e-6, 1.4826 * ctrl_mad, 1.2533 * ctrl_mnad)
  )

print(z_params)
write_csv(z_params, file.path(output_dir, "motif_zscore_parameters.csv"))

combined_long <- combined_long %>%
  left_join(z_params %>% select(metric, ctrl_median, scale_factor), by = "metric") %>%
  mutate(
    raw_value = value,
    value     = (raw_value - ctrl_median) / scale_factor   # robust Z used below
  ) %>%
  select(-ctrl_median, -scale_factor)

## ---- 4. Colour palette -----------------------------------------------

base_colors <- c(
  "Type 1" = "#1f77b4",  # blue
  "Type 2" = "#d62728",  # red
  "Type 3" = "#e8a800"   # yellow
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

combined_long <- combined_long %>%
  mutate(Type_combined = paste(promoter_type, gene_type, sep = "."))

## ---- 5. Layout: promoter type on X, F dodged left / P dodged right --------

dodge <- 0.15

combined_long <- combined_long %>%
  mutate(
    x_num    = as.numeric(promoter_type),
    x_offset = ifelse(gene_type == "Functional", -dodge, dodge),
    x_pos    = x_num + x_offset
  )

## ---- 6. Outlier detection --------------------------------------------
# Standard 1.5*IQR (Tukey) rule per promoter_type x gene_type x metric.

combined_long <- combined_long %>%
  group_by(metric, promoter_type, gene_type) %>%
  mutate(
    q1  = quantile(value, 0.25, na.rm = TRUE),
    q3  = quantile(value, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    is_outlier = value > (q3 + 1.5 * iqr) | value < (q1 - 1.5 * iqr)
  ) %>%
  ungroup()

outlier_points <- combined_long %>% filter(is_outlier)

median_data <- combined_long %>%
  group_by(metric, promoter_type, gene_type, Type_combined) %>%
  summarise(
    Median_value = median(value, na.rm = TRUE),
    x_pos = first(x_num) + ifelse(first(gene_type) == "Functional", -dodge, dodge),
    .groups = "drop"
  )

# Groups with <=2 points can't draw a meaningful violin -> plain points.
group_counts <- combined_long %>%
  count(metric, Type_combined, name = "n_points") %>%
  mutate(key = paste(metric, Type_combined))

thin_keys <- group_counts %>% filter(n_points <= 2) %>% pull(key)

combined_long <- combined_long %>% mutate(key = paste(metric, Type_combined))

violin_data <- combined_long %>% filter(!(key %in% thin_keys))
thin_data   <- combined_long %>% filter(key %in% thin_keys)

## ---- 6b. KS test: Functional vs Pseudogene, per promoter type -------------
format_p <- function(p) {
  paste0("p = ", formatC(p, digits = 2, format = "e"))
}

ks_results <- combined_long %>%
  group_by(metric, promoter_type) %>%
  summarise(
    ks = list(suppressWarnings(ks.test(
      value[gene_type == "Functional"],
      value[gene_type == "Pseudogene"]
    ))),
    .groups = "drop"
  ) %>%
  mutate(
    D_stat  = vapply(ks, function(k) unname(k$statistic), numeric(1)),
    p_value = vapply(ks, function(k) k$p.value, numeric(1))
  ) %>%
  select(metric, promoter_type, D_stat, p_value)

write_csv(ks_results, file.path(output_dir, "motif_dist_ks_test.csv"))

# One label per facet x promoter type, centred above that pair of violins.
ks_labels <- ks_results %>%
  mutate(
    x_pos = match(as.character(promoter_type), promoter_order),
    label = paste0("KS D = ", formatC(D_stat, digits = 3, format = "f"),
                   "\n", format_p(p_value))
  )

# F / P sub-labels under each violin (repeated in every facet).
fp_labels <- data.frame(
  x_pos = rep(seq_len(n_types), each = 2) + rep(c(-dodge, dodge), times = n_types),
  lab   = rep(c("F", "P"), times = n_types),
  # only the bottom panel (Score) gets the F/P letters; the x axis is shared
  metric = factor("Score", levels = levels(combined_long$metric))
)

## ---- 7. Combined plot (counts + scores as subplots) ----------------------
y_breaks <- c(-10, -1, 0, 1, 2, 3, 10, 100, 1000, 10000)

base_theme <- theme_minimal(base_size = 20) +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing.y = unit(1.2, "lines"),
    # invisible ticks, but their length pushes the type labels below the F/P letters
    axis.ticks.x = element_line(colour = "transparent"),
    axis.ticks.length.x = unit(16, "pt"),
    axis.text.x = element_text(face = "bold")
  )

forest_plot <- ggplot() +
  # reference lines: control median (0), Z = 2, Z = 3
  geom_hline(data = data.frame(yint = c(0, 2, 3),
                               lt   = c("solid", "dashed", "dotted")),
             aes(yintercept = yint, linetype = lt),
             color = "grey50", linewidth = 0.4, show.legend = FALSE) +
  scale_linetype_identity() +
  geom_violin(
    data = violin_data,
    aes(x = x_pos, y = value, group = Type_combined,
        fill = Type_combined, color = Type_combined),
    orientation = "x",
    scale = "width",
    width = 0.28,
    alpha = 0.35,
    trim = TRUE,
    linewidth = 0.6
  ) +
  geom_jitter(
    data = thin_data,
    aes(x = x_pos, y = value, color = Type_combined, shape = gene_type),
    width = 0.06, height = 0, size = 2.4, alpha = 0.8, stroke = 0.8
  ) +
  geom_jitter(
    data = outlier_points,
    aes(x = x_pos, y = value, color = Type_combined, shape = gene_type),
    width = 0.06, height = 0, size = 2.4, alpha = 0.8, stroke = 0.8
  ) +
  geom_point(
    data = median_data,
    aes(x = x_pos, y = Median_value),
    shape = 23, size = 2, fill = "white", color = "black"
  ) +
  # F / P letters just below the panel, one per violin
  geom_text(
    data = fp_labels,
    aes(x = x_pos, y = -Inf, label = lab),
    inherit.aes = FALSE, vjust = 1.4, size = 5,
    family = "Times New Roman", color = "grey25"
  ) +
  geom_text(
    data = ks_labels,
    aes(x = x_pos, y = Inf, label = label),
    vjust = 1.1, hjust = 0.5, size = 4, lineheight = 0.9,
    family = "Times New Roman", color = "black"
  ) +
  scale_shape_manual(values = c("Functional" = 16, "Pseudogene" = 1)) +  # solid vs hollow
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(
    breaks = seq_along(promoter_order),
    labels = promoter_order,
    expand = expansion(add = 0.6)
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = y_breaks,
    labels = y_breaks,
    expand = expansion(mult = c(0.02, 0.18))   # headroom for KS labels
  ) +
  coord_cartesian(clip = "off") +              # lets F/P letters draw below the panel
  # stacked panels, shared x axis; strip labels sit centred above each panel,
  # so "Score" falls in the middle between the two panels
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  labs(
    title = "Distribution of PolIII Promoter Motifs in Pseudogenes and Functional ncRNA",
    x = NULL,
    y = "Robust Z-score (pseudo-log10)"
  ) +
  base_theme

ggsave(
  filename = file.path(output_dir, "motif_distribution_forest_plot.png"),
  plot = forest_plot, width =12, height = 10, dpi = 300
)

## ---- 8. Export the outliers used in the plot ------------------------------

write_csv(
  outlier_points %>%
    select(gene_name, promoter_type, gene_type, metric,
           raw_value, Z_score = value),
  file.path(output_dir, "motif_distribution_outliers.csv")
)

## ---- 9. Export median stats ------------------------------------------------
write_csv(
  median_data %>%
    select(metric, promoter_type, gene_type, Median_value),
  file.path(output_dir, "motif_median_stats.csv")
)

## ---- 10. Done ---------------------------------------------------------------

cat("Plots written to:", normalizePath(output_dir), "\n")
cat(" - motif_distribution_forest_plot.png\n")
cat(" - motif_distribution_outliers.csv\n")
cat(" - motif_dist_ks_test.csv\n")
cat(" - motif_median_stats.csv\n")