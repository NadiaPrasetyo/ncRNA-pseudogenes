# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Data-loading step
data <- read.csv("../../data/snp_enrichment_genes_pseudogenes.csv")

# ---------------------------------------------------------------------------
# CONFIG - adjust these if your data differs from these assumptions
# ---------------------------------------------------------------------------
# Metric to plot per source, used directly (no Z-score) since "enrichment"
# is already normalized against the flanking control region.
metric_col <- "enrichment"

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
    cols = matches("^(snp_count|snp_density|flank_count|flank_density|enrichment)_"),
    names_to = c(".value", "source"),
    names_pattern = paste0("(snp_count|snp_density|flank_count|flank_density|enrichment)_(",
                           paste(snp_sources, collapse = "|"), ")")
  )

# Clean the data: Remove rows with non-finite values for the chosen metric
clean_data <- long_data %>%
  filter(is.finite(.data[[metric_col]]))

# ---------------------------------------------------------------------------
# No Z-score normalization: "enrichment" is already normalized against the
# flanking control region, so it's plotted directly.
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
# Shared y-axis scale per plot type, computed ACROSS all sources, so the
# four "Major Spliceosomal" PDFs (one per source) share identical axis
# limits/ticks and are directly comparable, and likewise for each other
# plot type. Each helper returns a list with limits + breaks for use in
# scale_y_continuous().
# ---------------------------------------------------------------------------
make_y_scale <- function(values, n_breaks = 6) {
  values <- values[is.finite(values)]
  rng <- range(values, na.rm = TRUE)
  # Small padding so points at the extremes aren't drawn on the panel edge
  pad <- diff(rng) * 0.05
  limits <- c(rng[1] - pad, rng[2] + pad)
  list(limits = limits, breaks = pretty(rng, n = n_breaks))
}

y_scale_1 <- make_y_scale(normalized_data$enrichment[normalized_data$Gene_group %in% combined_groups])
y_scale_2 <- make_y_scale(normalized_data$enrichment[normalized_data$Gene_group %in% combined_groups_2])
y_scale_3 <- make_y_scale(normalized_data$enrichment[normalized_data$Gene_group %in% remaining_ncRNAs])
y_scale_combined <- make_y_scale(normalized_data$enrichment[normalized_data$Gene_group != "TRNA"])

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
  
  plot_1 <- ggplot(data_1, aes(x = Gene_Type_label, y = enrichment, color = Gene_Type_combined)) +
    geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
    stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +
    scale_color_manual(values = custom_colors) +
    scale_y_continuous(limits = y_scale_1$limits, breaks = y_scale_1$breaks) +
    labs(
      title = paste("SNP Enrichment of Major Spliceosomal RNAs -", src),
      x = "Gene Type",
      y = "SNP Enrichment"
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
  
  plot_2 <- ggplot(data_2, aes(x = Gene_Type_label, y = enrichment, color = Gene_Type_combined)) +
    geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
    stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +
    scale_color_manual(values = custom_colors) +
    scale_y_continuous(limits = y_scale_2$limits, breaks = y_scale_2$breaks) +
    labs(
      title = paste("SNP Enrichment of Minor Spliceosomal RNAs -", src),
      x = "Gene Type",
      y = "SNP Enrichment"
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
  
  plot_3 <- ggplot(data_3, aes(x = Gene_Type_label, y = enrichment, color = Gene_Type_combined)) +
    geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
    stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +
    scale_color_manual(values = custom_colors) +
    scale_y_continuous(limits = y_scale_3$limits, breaks = y_scale_3$breaks) +
    labs(
      title = paste("SNP Enrichment of Other sncRNAs -", src),
      x = "Gene Type",
      y = "SNP Enrichment"
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
                                    aes(x = Gene_Type_combined_simple, y = enrichment, color = Gene_Type_combined_simple)) +
    geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
    stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +
    scale_color_manual(values = c("Functional" = "firebrick", "Pseudogene" = "cornflowerblue")) +
    scale_y_continuous(limits = y_scale_combined$limits, breaks = y_scale_combined$breaks) +
    labs(
      title = paste("SNP Enrichment: Pseudogenes vs Functional Genes -", src),
      x = "Gene Type",
      y = "SNP Enrichment"
    ) +
    single_source_theme +
    theme(axis.text.x = element_text(size = 24, angle = 0, hjust = 0.5))
  
  ggsave(filename = paste0("../../results/SNP_Enrichment_", src, "_Pseudogene_vs_Functional.pdf"),
         plot = combined_gene_type_plot, width = 12, height = 7)
}