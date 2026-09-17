#!/usr/bin/env Rscript
# =============================================================================
# plot_motif_distributions.R
#
# Plots the distribution of type1/type2/type3 motif counts and mean scores,
# split by gene_type (functional vs pseudogene).
#
# Input : results/mcast_summary_grouped_annotated.csv
# Output: PNG plots written to results/plots/
#
# Genes with no promoter of a given type (contains_typeN_promoter == FALSE,
# or NA count/score) are treated as having a motif count of 0 and a score
# of 0 for that type, rather than being dropped.
# =============================================================================

## ---- 0. Setup -------------------------------------------------------------

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "readr", "stringr")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)

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

# Normalize gene_type labels (trim whitespace / case) just in case
df_clean <- df_clean %>%
  mutate(gene_type = str_trim(gene_type))

## ---- 3. Reshape to long format -------------------------------------------

counts_long <- df_clean %>%
  select(gene_name, gene_type, type1_motif_counts, type2_motif_counts, type3_motif_counts) %>%
  pivot_longer(
    cols = c(type1_motif_counts, type2_motif_counts, type3_motif_counts),
    names_to = "promoter_type",
    values_to = "motif_count"
  ) %>%
  mutate(promoter_type = recode(promoter_type,
                                type1_motif_counts = "Type 1",
                                type2_motif_counts = "Type 2",
                                type3_motif_counts = "Type 3"
  ))

scores_long <- df_clean %>%
  select(gene_name, gene_type, type1_score, type2_score, type3_score) %>%
  pivot_longer(
    cols = c(type1_score, type2_score, type3_score),
    names_to = "promoter_type",
    values_to = "score"
  ) %>%
  mutate(promoter_type = recode(promoter_type,
                                type1_score = "Type 1",
                                type2_score = "Type 2",
                                type3_score = "Type 3"
  ))

## ---- 4. Plotting helpers --------------------------------------------------

base_theme <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

fill_colors <- c("functional" = "#2E86AB", "pseudogene" = "#E76F51")

make_dist_plot <- function(data, y_var, y_label, title) {
  ggplot(data, aes(x = gene_type, y = .data[[y_var]], fill = gene_type)) +
    geom_violin(alpha = 0.5, trim = FALSE, scale = "width") +
    geom_jitter(width = 0.08, alpha = 0.1, size = 1) +
    facet_wrap(~ promoter_type, nrow = 1) +
    scale_fill_manual(values = fill_colors, na.value = "grey70") +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 1, base = 10),
    ) +
    labs(title = title, x = NULL, y = y_label, fill = "Gene type") +
    base_theme
}

## ---- 5. Motif count distribution plot -------------------------------------

p_counts <- make_dist_plot(
  counts_long,
  y_var = "motif_count",
  y_label = "Motif count",
  title = "Distribution of Motif Counts by Promoter Type"
)

ggsave(
  filename = file.path(output_dir, "motif_counts_distribution.png"),
  plot = p_counts, width = 10, height = 5, dpi = 300
)

## ---- 6. Score distribution plot -------------------------------------------

p_scores <- make_dist_plot(
  scores_long,
  y_var = "score",
  y_label = "Score",
  title = "Distribution of Motif Scores by Promoter Type"
)

ggsave(
  filename = file.path(output_dir, "motif_scores_distribution.png"),
  plot = p_scores, width = 10, height = 5, dpi = 300
)

## ---- 7. Combined summary bar plot (mean score, mean motif count) ---------

summary_stats <- bind_rows(
  counts_long %>%
    group_by(gene_type, promoter_type) %>%
    summarise(mean_value = mean(motif_count, na.rm = TRUE),
              sd_value = sd(motif_count, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    mutate(metric = "Mean motif count"),
  scores_long %>%
    group_by(gene_type, promoter_type) %>%
    summarise(mean_value = mean(score, na.rm = TRUE),
              sd_value = sd(score, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    mutate(metric = "Mean score")
)

write_csv(summary_stats, file.path(output_dir, "summary_stats.csv"))

p_summary <- ggplot(summary_stats, aes(x = promoter_type, y = mean_value, fill = gene_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = pmax(mean_value - sd_value, 0), ymax = mean_value + sd_value),
    position = position_dodge(width = 0.8), width = 0.2
  ) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = fill_colors, na.value = "grey70") +
  labs(title = "Mean Motif Count and Score by Promoter Type",
       x = "Promoter type", y = "Mean value (\u00B1 SD)", fill = "Gene type") +
  base_theme

ggsave(
  filename = file.path(output_dir, "mean_summary_barplot.png"),
  plot = p_summary, width = 10, height = 5, dpi = 300
)

## ---- 8. Done ---------------------------------------------------------------

cat("Plots written to:", normalizePath(output_dir), "\n")
cat(" - motif_counts_distribution.png\n")
cat(" - motif_scores_distribution.png\n")
cat(" - mean_summary_barplot.png\n")
cat(" - summary_stats.csv\n")