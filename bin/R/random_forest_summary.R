#!/usr/bin/env Rscript
# Random forest results: functional probability distribution + feature importance
# Usage: Rscript plot_rf_results.R input.csv [output_prefix]

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ragg)
  library(scales)
})

#!/usr/bin/env Rscript
# Random forest results: functional probability distribution + feature importance
# (with SD / SE error bars from leave-one-family-out models)

infile     <- "../../results/combined_gene_data_with_snp_promoters.csv"
imp_file   <- "../../results/rf_feature_importance.csv"
prefix     <- "../../results/rf"

ERROR_TYPE <- "sd"   # "se" (standard error) or "sd" (standard deviation)

FONT_FAMILY <- "serif"
# ---------------------------------------------------------------- data
df <- read_csv(infile, show_col_types = FALSE) %>%
  mutate(
    Gene_Type = factor(Gene_Type, levels = c("Functional", "Pseudogene")),
    label_txt = ifelse(label == 1, "Functional", "Pseudogene")
  )

# If your CSV uses a different pseudogene label in Gene_Type, fall back to `label`
if (all(is.na(df$Gene_Type)) || nlevels(droplevels(df$Gene_Type)) < 2) {
  df$Gene_Type <- factor(df$label_txt, levels = c("Functional", "Pseudogene"))
}

cols <- c(Functional = "#d62728", Pseudogene = "#1f77b4")

# Dataset counts (from the training log)
subtitle_txt <- paste0(
  "Full: 535 functional / 2919 pseudogenes  |  ",
  "Out-of-family ROC AUC Score: 0.93"
)

# ---------------------------------------------------------------- 1. score distribution
nb  <- 40
brk <- seq(0, 1, length.out = nb + 1)

bin_counts <- df %>%
  mutate(bin = factor(cut(functional_probability, brk,
                          include.lowest = TRUE, labels = FALSE),
                      levels = seq_len(nb))) %>%
  count(bin, Gene_Type, .drop = FALSE) %>%
  mutate(bin = as.integer(bin),
         mid = (brk[bin] + brk[bin + 1]) / 2)

# log2(count + 1) transform, applied AFTER stacking (via coord_transform), so the
# bars are true stacked counts drawn on a log2 axis. The +1 keeps empty bins
# at zero height instead of -Inf.
log2p1 <- trans_new("log2p1",
                    transform = function(x) log2(x + 1),
                    inverse   = function(x) 2^x - 1)

ymax   <- max(tapply(bin_counts$n, bin_counts$bin, sum))
y_brks <- 2^(0:ceiling(log2(ymax + 1))) - 1   # 0, 1, 3, 7, 15, ... (log2 spaced)

p_dist <- ggplot(bin_counts, aes(x = mid, y = n, fill = Gene_Type)) +
  geom_col(width = 1 / nb, position = "stack", alpha = 0.85,
           colour = "white", linewidth = 0.2) +
  geom_text(data = filter(bin_counts, n > 0),
            aes(label = n), position = position_stack(vjust = 0.5),
            size = 2.8, family = FONT_FAMILY, colour = "black",
            check_overlap = TRUE) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30") +
  scale_fill_manual(values = cols, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = y_brks, expand = expansion(mult = c(0, 0.02))) +
  coord_transform(y = log2p1) +
  labs(
    title    = "Random forest functional probability",
    subtitle = subtitle_txt,
    x        = "Functional probability (RF score)",
    y        = "Number of genes (log2 scale)"
  ) +
  theme_minimal(base_size = 18, base_family = FONT_FAMILY) +
  theme(legend.position = "top",
        plot.subtitle = element_text(size = 12, colour = "grey30"),
        panel.grid.minor = element_blank())

# Optional: counts version (log10 y) is useful because classes are imbalanced
p_dist_counts <- ggplot(df, aes(x = functional_probability, fill = Gene_Type)) +
  geom_histogram(bins = 40, alpha = 0.7, position = "identity",
                 colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30") +
  scale_fill_manual(values = cols, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_sqrt() +
  labs(title = "Functional probability (gene counts, sqrt scale)",
       x = "Functional probability (RF score)", y = "Number of genes") +
  theme_minimal(base_size = 18, base_family = FONT_FAMILY) +
  theme(legend.position = "top")

# ---------------------------------------------------------------- 2. feature importance
# ---------------------------------------------------------------- 2. feature importance
# Map raw feature names from the Python output to display labels. Any feature
# not listed here keeps its raw name.
feature_labels <- c(
  "PhyloP100_median"               = "Conservation",
  "enrichment_intergenic_dbsnp"    = "Variation",
  "ENCODE_max"                     = "Expression",
  "type1_score"                    = "Pol III Type 1 Score",
  "type1_motif_counts"             = "Pol III Type 1 Motif Count",
  "type2_score"                    = "Pol III Type 2 Score",
  "type2_motif_counts"             = "Pol III Type 2 Motif Count",
  "type3_score"                    = "Pol III Type 3 Score",
  "type3_motif_counts"             = "Pol III Type 3 Motif Count"
)

err_col   <- paste0("importance_", ERROR_TYPE)
err_label <- toupper(ERROR_TYPE)
n_models  <- NA_integer_

imp <- read_csv(imp_file, show_col_types = FALSE) %>%
  mutate(
    Feature    = ifelse(feature %in% names(feature_labels),
                        feature_labels[feature], feature),
    Importance = importance_mean,
    err        = .data[[err_col]],
    lo         = pmax(0, Importance - err),   # importance can't go below 0
    hi         = Importance + err
  ) %>%
  arrange(desc(Importance)) %>%
  mutate(Feature = factor(Feature, levels = rev(Feature)))  # largest at top

n_models <- imp$n_models[1]

p_imp <- ggplot(imp, aes(x = Importance, y = Feature)) +
  geom_col(fill = "orange", width = 0.7) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25,
                orientation = "y", linewidth = 0.6, colour = "black") +
  geom_text(aes(x = hi, label = sprintf("%#.3g", Importance)),
            family = FONT_FAMILY, hjust = -0.15, size = 6) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title    = "Random forest feature importance",
    subtitle = paste0("Mean \u00B1 ", err_label, " across ", n_models,
                      " leave-one-family-out models"),
    x        = "Importance (mean decrease in impurity)",
    y        = NULL
  ) +
  theme_classic(base_size = 18, base_family = FONT_FAMILY) +
  theme(plot.subtitle = element_text(size = 12, colour = "grey30"))

# ---------------------------------------------------------------- save
# save: use ragg for PNG, cairo_pdf for PDF
save_png <- function(name, plot, w, h) {
  ggsave(paste0(prefix, name, ".png"), plot, width = w, height = h,
         dpi = 300, device = ragg::agg_png)
}
save_pdf <- function(name, plot, w, h) {
  ggsave(paste0(prefix, name, ".pdf"), plot, width = w, height = h,
         device = cairo_pdf)
}

save_png("_score_distribution",        p_dist,        7, 4.5)
save_png("_score_distribution_counts", p_dist_counts, 7, 4.5)
save_png("_feature_importance",        p_imp,         7, 4.5)
save_pdf("_score_distribution",        p_dist,        7, 4.5)
save_pdf("_feature_importance",        p_imp,         7, 4.5)

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combined <- p_dist + p_imp + plot_annotation(tag_levels = "A")
  save_png("_combined", combined, 14, 5)
}

message("Done. Plots written with prefix '", prefix, "'.")