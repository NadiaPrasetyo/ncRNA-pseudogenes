# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)   # install.packages("patchwork") if needed

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
data <- read.csv("../../results/combined_gene_data_with_snp_promoters.csv")

snp_panel_source <- "dbsnp"      # "gnomad", "1000genomes", "pangenome" or "dbsnp"
snp_sources      <- c("1000genomes", "gnomad", "pangenome", "dbsnp")
snp_metric_col   <- "enrichment_intergenic"

# ---------------------------------------------------------------------------
# Colors (vibrant = Functional, pastel = Pseudogene)
# ---------------------------------------------------------------------------
base_colors <- c(
  "RNU1" = "#1f77b4", "RNU2" = "#d62728", "RNU4" = "#e8a800", "RNU5" = "#2ca02c",
  "RNU6" = "#9467bd", "RNU4ATAC" = "#8c564b", "RNU6ATAC" = "#e377c2",
  "RNU11" = "#555555", "RNU12" = "#8c8d00", "VTRNA" = "#17becf",
  "RNY" = "#e0507a", "TRNA" = "#7f5aa2", "RN7SL" = "#f4a300",
  "RNU7" = "#6a3fbf", "RN7SK" = "#3b8f6b", "Pooled" = "#000000"
)

make_pastel <- function(hex, amount = 0.6) {
  rgb_val <- grDevices::col2rgb(hex)
  pastel_val <- rgb_val + (255 - rgb_val) * amount
  grDevices::rgb(pastel_val[1, ], pastel_val[2, ], pastel_val[3, ], maxColorValue = 255)
}
pastel_colors <- setNames(make_pastel(base_colors), names(base_colors))

custom_colors <- c(
  setNames(base_colors,   paste0(names(base_colors),   ".Functional")),
  setNames(pastel_colors, paste0(names(pastel_colors), ".Pseudogene"))
)

outline_pastel <- setNames(make_pastel(base_colors, amount = 0.3), names(base_colors))
outline_colors <- c(
  setNames(base_colors,   paste0(names(base_colors),   ".Functional")),
  setNames(outline_pastel, paste0(names(outline_pastel), ".Pseudogene"))
)

# ---------------------------------------------------------------------------
# Shared preprocessing
# ---------------------------------------------------------------------------
data <- data %>%
  mutate(Gene_Type = ifelse(tolower(Gene_Type) == "pseudogene", "Pseudogene", "Functional"))

combined_groups   <- c("RNU1", "RNU2", "RNU4", "RNU5", "RNU6")
combined_groups_2 <- c("RNU4ATAC", "RNU6ATAC", "RNU11", "RNU12")
remaining_ncRNAs  <- setdiff(unique(data$Gene_group), c(combined_groups, combined_groups_2))
gene_group_order  <- c(combined_groups, combined_groups_2, remaining_ncRNAs, "Pooled")
n_groups          <- length(gene_group_order)

add_pooled_layout <- function(df) {
  bind_rows(df, df %>% mutate(Gene_group = "Pooled")) %>%
    filter(Gene_group %in% gene_group_order) %>%
    mutate(
      Gene_group = factor(Gene_group, levels = gene_group_order),
      Gene_Type  = factor(Gene_Type, levels = c("Functional", "Pseudogene")),
      Gene_Type_combined = paste(Gene_group, Gene_Type, sep = "."),
      x_num = as.numeric(Gene_group),
      x_pos = x_num + ifelse(Gene_Type == "Functional", -0.15, 0.15)
    )
}

make_medians <- function(df) {
  df %>%
    group_by(Gene_group, Gene_Type, Gene_Type_combined) %>%
    summarize(
      Median_val = median(val, na.rm = TRUE),
      x_pos = first(x_num) + ifelse(first(Gene_Type) == "Functional", -0.15, 0.15),
      .groups = "drop"
    )
}

thin_group_names <- function(df) {
  df %>% count(Gene_Type_combined, name = "n_points") %>%
    filter(n_points <= 2) %>% pull(Gene_Type_combined)
}

# ---------------------------------------------------------------------------
# FDR / p-value threshold
# Each dataset (ENCODE, PhyloP, SNP) is scored against its OWN pooled values
# (all genes, functional + pseudogene): robust median/MAD -> z -> p-value ->
# Benjamini-Hochberg. "upper" = one-sided (ENCODE, PhyloP); "two.sided" = either
# direction (SNP). Returns the p-value threshold and the least extreme passing z
# on each side (NA if none passes).
# ---------------------------------------------------------------------------
# Robust z-score parameters: median/MAD, falling back to mean absolute deviation
null_params <- function(x) {
  center <- median(x, na.rm = TRUE)
  mad_raw <- median(abs(x - center), na.rm = TRUE)
  s <- if (mad_raw > 1e-6) {
    1.4826 * mad_raw
  } else {
    1.2533 * mean(abs(x - center), na.rm = TRUE)
  }
  list(center = center, scale = s)
}

calc_fdr_threshold <- function(z, tail = c("upper", "two.sided"), alpha = 0.05) {
  tail <- match.arg(tail)
  p   <- if (tail == "upper") pnorm(z, lower.tail = FALSE) else 2 * pnorm(-abs(z))
  sig <- which(p.adjust(p, method = "BH") <= alpha)
  hi  <- sig[z[sig] > 0]
  lo  <- sig[z[sig] < 0]
  list(
    p_thr = if (length(sig)) max(p[sig]) else NA_real_,
    z_hi  = if (length(hi)) min(z[hi]) else NA_real_,
    z_lo  = if (length(lo)) max(z[lo]) else NA_real_
  )
}

# Write outlier CSV with the thresholds as "#" comment lines at the top
write_outliers <- function(df, path, thr, unit = "Z-score") {
  vals <- na.omit(thr$value_thr)
  writeLines(c(
    "# FDR (Benjamini-Hochberg) = 0.05",
    sprintf("# p-value threshold = %s", format(thr$p_thr, digits = 4)),
    sprintf("# %s threshold = %s", unit,
            if (length(vals)) paste(format(vals, digits = 4), collapse = ", ") else "none passed")
  ), path)
  write.table(df %>% select(-any_of("val")), path, sep = ",",
              row.names = FALSE, append = TRUE)
}

# ---------------------------------------------------------------------------
# Z-score panel data (ENCODE / PhyloP100): threshold from pooled data of that column
# ---------------------------------------------------------------------------
make_z_panel_data <- function(data, col, pooled_functional_cut,
                              log_transform = FALSE, floor_val = 0.001) {
  if (log_transform) {
    # NA and non-positive values (treated as 0) become a tiny pseudocount
    data[[col]][is.na(data[[col]]) | data[[col]] <= 0] <- floor_val
  }
  clean <- data %>% filter(is.finite(.data[[col]]))
  x     <- if (log_transform) log10(clean[[col]]) else clean[[col]]
  np    <- null_params(x)                    # pooled: all genes, F + P
  z     <- (x - np$center) / np$scale
  thr   <- calc_fdr_threshold(z, "upper")
  thr$value_thr <- thr$z_hi
  
  nd <- clean %>%
    mutate(Z_score = z, val = Z_score) %>%
    add_pooled_layout()
  
  thin <- thin_group_names(nd)
  
  outliers <- nd %>%
    filter(Gene_group != "Pooled",
           (Gene_Type == "Pseudogene" & Z_score >= thr$z_hi) |
             (Gene_Type == "Functional" & Z_score <= 0)) %>%
    mutate(Outlier_Reason = if_else(Gene_Type == "Pseudogene",
                                    "Pseudogene above FDR threshold",
                                    "Functional with Z <= 0")) %>%
    arrange(Gene_group, desc(Z_score))
  
  pooled_jitter <- nd %>%
    filter(Gene_group == "Pooled",
           (Gene_Type == "Pseudogene" & Z_score >= thr$z_hi) |
             (Gene_Type == "Functional" & Z_score > pooled_functional_cut))
  
  list(
    violin   = nd %>% filter(!(Gene_Type_combined %in% thin)),
    points   = bind_rows(nd %>% filter(Gene_Type_combined %in% thin), outliers, pooled_jitter),
    medians  = make_medians(nd),
    outliers = outliers,
    thr      = thr
  )
}

# ---------------------------------------------------------------------------
# SNP enrichment panel data: two-sided threshold from pooled SNP data
# ---------------------------------------------------------------------------
make_snp_panel_data <- function(data, src) {
  long_data <- data %>%
    pivot_longer(
      cols = matches("^(snp_density|enrichment_intergenic)_"),
      names_to = c(".value", "source"),
      names_pattern = paste0("(snp_density|enrichment_intergenic)_(",
                             paste(snp_sources, collapse = "|"), ")")
    ) %>%
    filter(is.finite(.data[[snp_metric_col]]), source == src) %>%
    mutate(val = .data[[snp_metric_col]])
  
  # pooled null (all genes) -> z -> two-sided FDR -> back to enrichment units
  np  <- null_params(long_data$val)
  thr <- calc_fdr_threshold((long_data$val - np$center) / np$scale, "two.sided")
  thr$value_thr <- np$center + c(thr$z_lo, thr$z_hi) * np$scale   # c(lower, upper)
  
  beyond <- function(v) {
    (!is.na(thr$value_thr[1]) & v <= thr$value_thr[1]) |
      (!is.na(thr$value_thr[2]) & v >= thr$value_thr[2])
  }
  
  nd <- add_pooled_layout(long_data) %>%
    group_by(Gene_Type_combined) %>%
    mutate(
      q1 = quantile(val, 0.25, na.rm = TRUE),
      q3 = quantile(val, 0.75, na.rm = TRUE),
      is_outlier = val < q1 - 1.5 * (q3 - q1) | val > q3 + 1.5 * (q3 - q1)
    ) %>%
    ungroup() %>%
    select(-q1, -q3) %>%
    mutate(flagged = Gene_Type == "Pseudogene" & beyond(val))
  
  thin <- thin_group_names(nd)
  
  list(
    violin   = nd %>% filter(!(Gene_Type_combined %in% thin) & !is_outlier),
    points   = nd %>% filter(Gene_Type_combined %in% thin | is_outlier | flagged),
    medians  = make_medians(nd),
    outliers = nd %>% filter(flagged, Gene_group != "Pooled") %>%
      arrange(Gene_group, val),
    thr      = thr
  )
}

#----------------------------------------------------------------------------
# HELPER: Y Axis superscript
#----------------------------------------------------------------------------
y_label_fmt <- function(x) {
  lab <- lapply(x, function(v) {
    if (is.na(v)) return(NA)
    if (v >= 1000) {
      e <- log10(v)
      if (abs(e - round(e)) < 1e-9) {
        bquote(10^.(round(e)))
      } else {
        e <- floor(e)
        bquote(.(v / 10^e) %*% 10^.(e))
      }
    } else {
      format(v, scientific = FALSE, trim = TRUE)
    }
  })
  as.expression(lab)
}

# ---------------------------------------------------------------------------
# Generic panel: gene groups on X, numeric value on Y
# ---------------------------------------------------------------------------
panel_theme <- theme_minimal() +
  theme(
    legend.position = "none",
    plot.title   = element_text(face = "bold", size = 15),
    axis.title.y = element_text(size = 14),
    axis.text.y  = element_text(size = 12),
    axis.text.y.right = element_text(size = 11, color = "grey30", hjust = 0, vjust = 1),
    axis.text.x  = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    text = element_text(family = "serif")
  )

fp_labels <- data.frame(
  x_pos = rep(seq_len(n_groups), each = 2) + rep(c(-0.15, 0.15), times = n_groups),
  lab   = rep(c("F", "P"), times = n_groups)
)

plot_panel <- function(pd, title, y_label, ref_lines, p_thr = NA,
                       y_trans = "identity", y_breaks = waiver(), y_view = NULL,
                       show_x = FALSE) {
  ref_lines <- ref_lines %>% filter(!is.na(yint))   # drop thresholds that don't exist
  # right-axis labels only for the FDR threshold line(s): value, FDR level, p-value
  thr_lines <- ref_lines %>% filter(lt == "dotted")
  thr_labs  <- paste0(format(signif(thr_lines$yint, 3), trim = TRUE),
                      "\nFDR = 0.05\np = ", formatC(p_thr, format = "e", digits = 2))
  p <- ggplot() +
    geom_hline(data = ref_lines, aes(yintercept = yint, linetype = lt),
               color = "grey30", linewidth = 0.4, show.legend = FALSE) +
    scale_linetype_identity() +
    geom_vline(xintercept = n_groups - 0.5, linetype = "dashed",
               color = "grey70", linewidth = 0.4) +
    geom_violin(
      data = pd$violin,
      aes(x = x_pos, y = val, group = Gene_Type_combined,
          fill = Gene_Type_combined, color = Gene_Type_combined),
      orientation = "x", scale = "width", width = 0.28,
      alpha = 0.55, trim = TRUE, linewidth = 0.7
    ) +
    geom_jitter(
      data = pd$points,
      aes(x = x_pos, y = val, color = Gene_Type_combined, shape = Gene_Type),
      width = 0.06, height = 0, size = 2.8, alpha = 0.9, stroke = 1.4
    ) +
    geom_point(
      data = pd$medians,
      aes(x = x_pos, y = Median_val),
      shape = 23, size = 2, fill = "white", color = "black"
    ) +
    geom_text(
      data = fp_labels,
      aes(x = x_pos, y = -Inf, label = lab),
      inherit.aes = FALSE, vjust = 1.2, size = 3.6,
      family = "serif", fontface = "plain", color = "grey25"
    ) +
    scale_shape_manual(values = c("Functional" = 16, "Pseudogene" = 1)) +
    scale_color_manual(values = outline_colors) +
    scale_fill_manual(values = custom_colors) +
    scale_x_continuous(
      breaks = seq_along(gene_group_order),
      labels = gene_group_order,
      expand = expansion(add = 0.6)
    ) +
    scale_y_continuous(
      trans = y_trans, breaks = y_breaks, labels = y_label_fmt,
      sec.axis = dup_axis(name = NULL, breaks = thr_lines$yint, labels = thr_labs)
    ) +
    coord_cartesian(ylim = y_view, clip = "off") +
    labs(title = title,
         x = if (show_x) "Gene family" else NULL,
         y = y_label) +
    panel_theme
  
  if (!show_x) {
    p <- p + theme(axis.text.x  = element_blank(),
                   axis.title.x = element_blank(),
                   plot.margin  = margin(t = 5.5, r = 5.5, b = 16, l = 5.5))
  } else {
    p <- p + theme(
      axis.title.x = element_text(size = 16, margin = margin(t = 8)),
      axis.text.x  = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
      axis.ticks.x = element_line(colour = "transparent"),
      axis.ticks.length.x = unit(16, "pt")
    )
  }
  p
}

# ---------------------------------------------------------------------------
# Build the three panels (each with its own pooled-data threshold)
# ---------------------------------------------------------------------------
snp_pd     <- make_snp_panel_data(data, snp_panel_source)
encode_pd  <- make_z_panel_data(data, "ENCODE_max",       pooled_functional_cut = 10, log_transform = TRUE)
phylop_pd  <- make_z_panel_data(data, "PhyloP100_median", pooled_functional_cut = 20)

p_snp <- plot_panel(
  snp_pd,
  title     = "SNP Enrichment",
  y_label   = "SNP enrichment",
  ref_lines = data.frame(yint = c(0, snp_pd$thr$value_thr),
                         lt   = c("solid", "dotted", "dotted")),
  p_thr = snp_pd$thr$p_thr,
  y_breaks  = seq(-10, 6, by = 2),
  y_view    = NULL,
  show_x    = TRUE
)

p_encode <- plot_panel(
  encode_pd,
  title   = "Expression",
  y_label =  expression(Z-score ~ (log[10] ~ ENCODE)),
  ref_lines = data.frame(yint = c(0, encode_pd$thr$z_hi), lt = c("solid", "dotted")),
  p_thr = encode_pd$thr$p_thr
)

p_phylop <- plot_panel(
  phylop_pd,
  title   = "Vertebrate Conservation",
  y_label = "Z-score",
  ref_lines = data.frame(yint = c(0, phylop_pd$thr$z_hi), lt = c("solid", "dotted")),
  p_thr = phylop_pd$thr$p_thr
)

# ---------------------------------------------------------------------------
# Join into one figure
# ---------------------------------------------------------------------------
combined_plot <- (p_phylop / p_encode / p_snp) +
  plot_layout(heights = c(1, 1, 1.25)) +
  plot_annotation(
    title = "Features of Pseudogenes and Functional ncRNA",
    tag_levels = "A"
  ) &
  theme(
    plot.tag   = element_text(face = "bold", size = 18, family = "serif"),
    plot.title = element_text(face = "bold", size = 22, family = "serif", hjust = 0.5)
  )

ggsave(
  filename = "../../results/combined_forest_plot_SNP_ENCODE_PhyloP100.pdf",
  plot = combined_plot,
  width = max(12, n_groups * 0.9),
  height = 16
)

# ---------------------------------------------------------------------------
# Outlier exports (pseudogenes past each dataset's own threshold; "Pooled" excluded)
# ---------------------------------------------------------------------------
write_outliers(snp_pd$outliers,    "../../results/SNP_Enrichment_outliers_dbsnp.csv",
               snp_pd$thr,    unit = "SNP enrichment (lower, upper)")
write_outliers(encode_pd$outliers, "../../results/ENCODE_Z_score_outliers.csv",
               encode_pd$thr)
write_outliers(phylop_pd$outliers, "../../results/phylop100_Z_score_outliers.csv",
               phylop_pd$thr)