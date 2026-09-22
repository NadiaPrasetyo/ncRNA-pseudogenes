#!/usr/bin/env Rscript
# =============================================================================
# plot_motif_distributions.R
#
# Forest-plot style visualization of type1/type2/type3 motif counts and
# mean scores, split by gene_type (functional vs pseudogene).
#
# Style notes (matches the ENCODE Z-score forest plot):
#   - Violin outlines are coloured by promoter type; fill is a light wash
#     of the same colour.
#   - Only OUTLIER points (1.5*IQR rule) are drawn individually, in colour.
#   - Functional = solid dot (shape 16), Pseudogene = hollow dot (shape 1).
#   - Functional = deep/vibrant colour, Pseudogene = pastel version of the
#     SAME hue (Type 1 = blue, Type 2 = red, Type 3 = yellow).
#   - Axes swapped vs. a standard boxplot: promoter type on the y-axis,
#     value (count / score) on the x-axis.
#   - No horizontal gridlines.
#   - Motif counts and motif scores are combined into ONE plot, shown as
#     side-by-side subplots (facets) rather than two separate files.
#
# Input : results/mcast_summary_grouped_annotated.csv
# Output: PNG plots written to results/plots/
#
# Genes with no promoter of a given type (contains_typeN_promoter == FALSE,
# or NA count/score) are treated as having a motif count of 0 and a score
# of 0 for that type, rather than being dropped.
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

input_csv  <- "../../results/mcast_summary_grouped_annotated.csv"
output_dir <- "../../results/plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load data -----------------------------------------------------

df <- read_csv(input_csv, show_col_types = FALSE)

# Basic sanity check on the column we need for grouping
if (!"gene_type" %in% names(df)) {
  stop("Expected a 'gene_type' column (e.g. 'functional' / 'pseudogene') in the input CSV.")
}

## ---- 2. Clean / impute missing values ----------------------------------
# Any gene lacking a given promoter type gets 0 for that type's
# motif count and score, instead of NA.

df_clean <- df %>%
  mutate(
    type1_motif_counts = ifelse(is.na(type1_motif_counts), 0, type1_motif_counts),
    type2_motif_counts = ifelse(is.na(type2_motif_counts), 0, type2_motif_counts),
    type3_motif_counts = ifelse(is.na(type3_motif_counts), 0, type3_motif_counts),
    type1_score        = ifelse(is.na(type1_score), 0, type1_score),
    type2_score        = ifelse(is.na(type2_score), 0, type2_score),
    type3_score        = ifelse(is.na(type3_score), 0, type3_score)
  )

# Normalize gene_type labels (trim whitespace / case), and title-case them
# ("Functional" / "Pseudogene") so they match the legend style used
# elsewhere.
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

# Single combined long data frame drives ONE faceted plot (counts + scores
# shown as subplots side by side) instead of two separate figures.
promoter_order <- c("Type 1", "Type 2", "Type 3")

combined_long <- bind_rows(counts_long, scores_long) %>%
  mutate(
    metric        = factor(metric, levels = c("Motif count", "Score")),
    gene_type     = factor(gene_type, levels = c("Functional", "Pseudogene")),
    promoter_type = factor(promoter_type, levels = promoter_order)
  )

## ---- 4. Colour palette -----------------------------------------------
# One vibrant base colour per promoter type (Type 1 = blue, Type 2 = red,
# Type 3 = yellow). Functional = the vibrant colour, Pseudogene = an
# auto-derived pastel ("washed out") version of the SAME hue.

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

## ---- 5. Forest-plot style layout ------------------------------------------
# y-axis = promoter type (categorical, like a forest plot), x-axis = value
# (motif count / score). Functional vs Pseudogene are dodged within the
# same row; colour + shape (solid vs hollow) both encode gene type.

group_levels_rev <- rev(promoter_order)  # so "Type 1" ends up at the top

combined_long <- combined_long %>%
  mutate(
    promoter_type_f = factor(promoter_type, levels = group_levels_rev),
    y_num    = as.numeric(promoter_type_f),
    y_offset = ifelse(gene_type == "Functional", 0.15, -0.15),
    y_pos    = y_num + y_offset
  )

## ---- 6. Outlier detection --------------------------------------------
# Only outlier points are drawn individually (in colour); everything else
# is summarized by the violin. Standard 1.5*IQR (Tukey) rule, computed
# separately per promoter_type x gene_type x metric.

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

# One median value per promoter_type x gene_type x metric, placed at the
# SAME fixed y_offset used for that gene_type's dodge, drawn as a diamond.
median_data <- combined_long %>%
  group_by(metric, promoter_type, gene_type, Type_combined) %>%
  summarise(
    Median_value = median(value, na.rm = TRUE),
    y_pos = first(y_num) + ifelse(first(gene_type) == "Functional", 0.15, -0.15),
    .groups = "drop"
  )

# Groups with too few points (<=2) can't draw a meaningful violin -> show
# them as plain points instead (same treatment as the outliers).
group_counts <- combined_long %>%
  count(metric, Type_combined, name = "n_points") %>%
  mutate(key = paste(metric, Type_combined))

thin_keys <- group_counts %>% filter(n_points <= 2) %>% pull(key)

combined_long <- combined_long %>% mutate(key = paste(metric, Type_combined))

violin_data <- combined_long %>% filter(!(key %in% thin_keys))
thin_data   <- combined_long %>% filter(key %in% thin_keys)

## ---- 7. Combined forest-style plot (counts + scores as subplots) ---------
x_breaks <- c(0, 25, 50, 100, 200, 400, 1000, 4000)

base_theme <- theme_minimal(base_size = 20) +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )

forest_plot <- ggplot() +
  geom_violin(
    data = violin_data,
    aes(x = value, y = y_pos, group = Type_combined,
        fill = Type_combined, color = Type_combined),
    orientation = "y",
    scale = "width",
    width = 0.28,
    alpha = 0.35,
    trim = TRUE,
    linewidth = 0.6
  ) +
  geom_jitter(
    data = thin_data,
    aes(x = value, y = y_pos, color = Type_combined, shape = gene_type),
    height = 0.06, width = 0, size = 2.4, alpha = 0.8, stroke = 0.8
  ) +
  geom_jitter(
    data = outlier_points,
    aes(x = value, y = y_pos, color = Type_combined, shape = gene_type),
    height = 0.06, width = 0, size = 2.4, alpha = 0.8, stroke = 0.8
  ) +
  geom_point(
    data = median_data,
    aes(x = Median_value, y = y_pos),
    shape = 23, size = 2, fill = "white", color = "black"
  ) +
  scale_shape_manual(values = c("Functional" = 16, "Pseudogene" = 1)) +  # solid vs hollow
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    breaks = x_breaks,
    labels = x_breaks
  ) +
  scale_y_continuous(
    breaks = seq_along(group_levels_rev),
    labels = group_levels_rev,
    expand = expansion(add = 0.6)
  ) +
  facet_wrap(~ metric, scales = "free_x") +
  labs(
    title = "Distribution of Motif Counts and Scores by Promoter Type",
    x = NULL,
    y = NULL
  ) +
  base_theme

ggsave(
  filename = file.path(output_dir, "motif_distribution_forest_plot.png"),
  plot = forest_plot, width = 12, height = 8, dpi = 300
)

## ---- 8. Export the outliers used in the plot ------------------------------

write_csv(
  outlier_points %>%
    select(gene_name, promoter_type, gene_type, metric, value),
  file.path(output_dir, "motif_distribution_outliers.csv")
)

## ---- 9. Summary bar plot (mean motif count, mean score) -------------------

summary_stats <- combined_long %>%
  group_by(gene_type, promoter_type, metric) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    n          = n(),
    .groups = "drop"
  ) %>%
  mutate(Type_combined = paste(promoter_type, gene_type, sep = "."))

write_csv(summary_stats, file.path(output_dir, "summary_stats.csv"))

p_summary <- ggplot(summary_stats, aes(x = promoter_type, y = mean_value,
                                       fill = Type_combined, color = Type_combined)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.75) +
  geom_errorbar(
    aes(ymin = pmax(mean_value - sd_value, 0), ymax = mean_value + sd_value),
    position = position_dodge(width = 0.8), width = 0.2
  ) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Mean Motif Count and Score by Promoter Type",
       x = "Promoter type", y = "Mean value (\u00B1 SD)") +
  base_theme +
  theme(legend.position = "none")

ggsave(
  filename = file.path(output_dir, "mean_summary_barplot.png"),
  plot = p_summary, width = 10, height = 5, dpi = 300
)

## ---- 10. Done ---------------------------------------------------------------

cat("Plots written to:", normalizePath(output_dir), "\n")
cat(" - motif_distribution_forest_plot.png\n")
cat(" - motif_distribution_outliers.csv\n")
cat(" - mean_summary_barplot.png\n")
cat(" - summary_stats.csv\n")