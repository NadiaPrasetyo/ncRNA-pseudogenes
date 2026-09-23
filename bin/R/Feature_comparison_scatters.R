#!/usr/bin/env Rscript
# ==============================================================================
# gene_analysis_plots.R
#
# Single-input version. Reads ONE combined CSV:
#   results/combined_gene_data_with_snp_promoters.csv
#
# Uses only these columns from that file:
#   Gene, PhyloP100_median, ENCODE_max, enrichment_intergenic_dbsnp, Gene_Type
#
# Produces interactive (plotly) HTML-free plots (rendered via plotly print),
# colored by functional vs pseudogene, with gene name shown on hover:
#   - ENCODE max expression vs phyloP100 median conservation
#   - ENCODE max expression vs enrichment source
#   - enrichment source vs phyloP100 median conservation
#   - enrichment source: jittered comparison, Functional vs Pseudogene
#   - PCA (PC1 vs PC2) of the three metrics
#   - 3D scatter of all three metrics together
#   - multi-dimensional 2D plot (phyloP vs ENCODE, shape=category, color=enrichment)
#
# Usage:
#   Rscript gene_analysis_plots.R
# ==============================================================================

## ---- 0. CONFIG --------------------------------------------------------------

input_csv <- "../../results/combined_gene_data_with_snp_promoters.csv"

# Which enrichment column to use as "the" enrichment source for all plots below.
# Must match a column name present in input_csv.
enrichment_col   <- "enrichment_intergenic_dbsnp"

# Human-readable label used in plot titles/axes/hover text for the column above.
# Update this alongside enrichment_col so labels stay accurate.
enrichment_label <- "SNP Enrichment (dbSNP155)"

output_dir <- "../../results/interactive_plots"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## ---- 1. PACKAGES -------------------------------------------------------------

required_pkgs <- c("dplyr", "readr", "stringr", "tidyr",
                   "plotly", "rlang")

missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(plotly)
  library(rlang)
})

## ---- 2. LOAD COMBINED INPUT TABLE ---------------------------------------------
# One file now supplies enrichment, expression, and conservation metrics
# together, so there is no merging step -- just select the columns we need.

message("Reading combined input table: ", input_csv)

raw_df <- read_csv(input_csv, show_col_types = FALSE)

required_cols <- c("Gene", "PhyloP100_median", "ENCODE_max", enrichment_col, "Gene_Type")
missing_cols  <- setdiff(required_cols, names(raw_df))
if (length(missing_cols) > 0) {
  stop(
    "Missing required column(s) in ", input_csv, ": ", paste(missing_cols, collapse = ", "),
    ". Available columns: ", paste(names(raw_df), collapse = ", ")
  )
}

combined_df <- raw_df %>%
  transmute(
    gene_name           = as.character(Gene),
    Gene_Type           = Gene_Type,
    Max_TPM             = as.numeric(ENCODE_max),
    Median_Conservation = as.numeric(PhyloP100_median),
    !!enrichment_col    := as.numeric(.data[[enrichment_col]])
  ) %>%
  mutate(
    Max_TPM = replace_na(Max_TPM, 0)
  )

## ---- 3. CLASSIFY FUNCTIONAL vs PSEUDOGENE -------------------------------------
# Gene_Type is the definitive flag (values like "protein_coding", "pseudogene",
# etc.). Any value containing "pseudo" (case-insensitive) is labeled
# "Pseudogene"; everything else is labeled "Functional".

message("Classifying functional vs pseudogene using column: Gene_Type")
message("  Unique values in 'Gene_Type': ", paste(unique(combined_df$Gene_Type), collapse = " | "))

combined_df <- combined_df %>%
  mutate(
    category = if_else(
      str_detect(str_to_lower(Gene_Type), "pseudo"),
      "Pseudogene",
      "Functional"
    )
  )

message("Classification result -> Functional: ", sum(combined_df$category == "Functional"),
        ", Pseudogene: ", sum(combined_df$category == "Pseudogene"))

## ---- 4. CLEAN NUMERIC FIELDS ---------------------------------------------------
# Max_TPM is already filled with 0 for missing values (step 2).
# Only drop a row here if the configured enrichment column itself is
# missing/non-finite, since every gene should have it. Rows missing
# Median_Conservation are kept in plot_df (still shown in TPM vs enrichment
# plots) and are only dropped from the specific plots/analyses that require
# conservation (2D conservation plots, PCA, 3D), done per-plot below.

plot_df <- combined_df %>%
  filter(is.finite(.data[[enrichment_col]])) %>%
  mutate(log_Max_TPM = log10(Max_TPM + 1))

message("Genes available for plotting: ", nrow(plot_df))
message("  Functional: ", sum(plot_df$category == "Functional"))
message("  Pseudogene: ", sum(plot_df$category == "Pseudogene"))
message("Genes missing phyloP conservation data (excluded only from conservation/PCA/3D plots): ",
        sum(!is.finite(plot_df$Median_Conservation)))

# Save the cleaned table for reference
write_csv(combined_df, file.path(output_dir, "combined_gene_metrics.csv"))

## ---- 5. COLOR PALETTE ---------------------------------------------------------

category_colors <- c("Functional" = "#B22222", "Pseudogene" = "#6495ED")  # brick red / cornflower blue

## ---- 6. HELPER: 2D INTERACTIVE SCATTER ----------------------------------------

make_scatter <- function(df, x_var, y_var, x_label, y_label, title,
                         x_raw_var = NULL, x_raw_label = NULL) {
  hover_x_line <- if (!is.null(x_raw_var)) {
    # show both the raw value and the scaled/plotted value
    paste0("<br>", x_raw_label, ": ", round(df[[x_raw_var]], 3),
           "<br>", x_label, ": ", round(df[[x_var]], 3))
  } else {
    paste0("<br>", x_label, ": ", round(df[[x_var]], 3))
  }
  
  plot_ly(
    data   = df,
    x      = as.formula(paste0("~", x_var)),
    y      = as.formula(paste0("~", y_var)),
    color  = ~category,
    colors = category_colors,
    type   = "scatter",
    mode   = "markers",
    text   = ~paste0(
      "Gene: ", gene_name,
      "<br>Category: ", category,
      hover_x_line,
      "<br>", y_label, ": ", round(get(y_var), 3)
    ),
    hoverinfo = "text",
    marker = list(size = 8, opacity = 0.75, line = list(width = 0.5, color = "white"))
  ) %>%
    layout(
      title  = title,
      xaxis  = list(title = x_label),
      yaxis  = list(title = y_label),
      legend = list(title = list(text = "Gene category"))
    )
}

## ---- 7. 2D PLOTS ---------------------------------------------------------------

message("Building 2D interactive plots using enrichment source: ", enrichment_col,
        " (\"", enrichment_label, "\")")

plot_expr_vs_conservation <- make_scatter(
  plot_df %>% filter(is.finite(Max_TPM), is.finite(Median_Conservation)),
  "log_Max_TPM", "Median_Conservation",
  "log10(ENCODE Max Expression + 1)", "phyloP100 Median Conservation",
  "ENCODE Max Expression (log-scaled) vs phyloP100 Median Conservation",
  x_raw_var = "Max_TPM", x_raw_label = "ENCODE Max Expression"
)

plot_expr_vs_enrichment <- make_scatter(
  plot_df %>% filter(is.finite(Max_TPM), is.finite(.data[[enrichment_col]])),
  "log_Max_TPM", enrichment_col,
  "log10(ENCODE Max Expression + 1)", enrichment_label,
  paste0("ENCODE Max Expression (log-scaled) vs ", enrichment_label),
  x_raw_var = "Max_TPM", x_raw_label = "ENCODE Max Expression"
)

plot_enrichment_vs_conservation <- make_scatter(
  plot_df %>% filter(is.finite(.data[[enrichment_col]]), is.finite(Median_Conservation)),
  enrichment_col, "Median_Conservation",
  enrichment_label, "phyloP100 Median Conservation",
  paste0(enrichment_label, " vs phyloP100 Median Conservation")
)

## ---- 7b. JITTER PLOT: enrichment source, Functional vs Pseudogene -------------
# Plotly's R "scatter" trace has no native jitter mode, so we place each
# category at an integer x-position (1 = Functional, 2 = Pseudogene) and add
# a small random horizontal offset per point. This keeps individual genes
# distinguishable (and hoverable) instead of overplotting them into a single
# vertical line, while the x-axis is still labeled by category, not by number.
# Only enrichment_col is used here.

message("Building jitter plot for enrichment source: ", enrichment_col,
        " (\"", enrichment_label, "\") by gene category...")

set.seed(42)  # reproducible jitter offsets across runs

jitter_width <- 0.18

plot_enrichment_jitter_df <- plot_df %>%
  filter(is.finite(.data[[enrichment_col]])) %>%
  mutate(
    category = factor(category, levels = c("Functional", "Pseudogene")),
    x_center = as.numeric(category),
    x_jitter = x_center + runif(n(), -jitter_width, jitter_width)
  )

plot_enrichment_jitter <- plot_ly(
  data   = plot_enrichment_jitter_df,
  x      = ~x_jitter,
  y      = as.formula(paste0("~", enrichment_col)),
  color  = ~category,
  colors = category_colors,
  type   = "scatter",
  mode   = "markers",
  text   = ~paste0(
    "Gene: ", gene_name,
    "<br>Category: ", category,
    "<br>", enrichment_label, ": ", round(.data[[enrichment_col]], 3)
  ),
  hoverinfo = "text",
  marker = list(size = 7, opacity = 0.65, line = list(width = 0.5, color = "white"))
) %>%
  layout(
    title  = paste0(enrichment_label, ": Functional vs Pseudogene (jittered)"),
    xaxis  = list(
      title    = "Gene category",
      tickmode = "array",
      tickvals = c(1, 2),
      ticktext = c("Functional", "Pseudogene"),
      range    = c(1 - 3 * jitter_width, 2 + 3 * jitter_width)
    ),
    yaxis  = list(title = enrichment_label),
    legend = list(title = list(text = "Gene category"))
  )

## ---- 8. PCA -------------------------------------------------------------------

message("Running PCA on Max_TPM, ", enrichment_col, ", Median_Conservation...")

pca_input <- plot_df %>%
  select(gene_name, category, Max_TPM, all_of(enrichment_col), Median_Conservation) %>%
  filter(is.finite(Max_TPM), is.finite(.data[[enrichment_col]]), is.finite(Median_Conservation))

message("Genes included in PCA (require all 3 metrics): ", nrow(pca_input),
        " (Functional: ", sum(pca_input$category == "Functional"),
        ", Pseudogene: ", sum(pca_input$category == "Pseudogene"), ")")

pca_matrix <- pca_input %>%
  select(Max_TPM, all_of(enrichment_col), Median_Conservation) %>%
  as.matrix()

pca_result <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)

pca_scores <- as.data.frame(pca_result$x) %>%
  mutate(
    gene_name = pca_input$gene_name,
    category  = pca_input$category
  )

var_explained <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

plot_pca <- plot_ly(
  data   = pca_scores,
  x      = ~PC1,
  y      = ~PC2,
  color  = ~category,
  colors = category_colors,
  type   = "scatter",
  mode   = "markers",
  text   = ~paste0(
    "Gene: ", gene_name,
    "<br>Category: ", category,
    "<br>PC1: ", round(PC1, 3),
    "<br>PC2: ", round(PC2, 3)
  ),
  hoverinfo = "text",
  marker = list(size = 8, opacity = 0.75, line = list(width = 0.5, color = "white"))
) %>%
  layout(
    title  = paste0("PCA: Max Expression, ", enrichment_label, ", phyloP100 Median Conservation"),
    xaxis  = list(title = paste0("PC1 (", var_explained[1], "% variance)")),
    yaxis  = list(title = paste0("PC2 (", var_explained[2], "% variance)")),
    legend = list(title = list(text = "Gene category"))
  )

## ---- 9. 3D SCATTER --------------------------------------------------------------

message("Building 3D interactive plot...")

# ENCODE Max expression is typically heavily right-skewed (a few very highly
# expressed genes vs. many low/zero-expression genes), which compresses
# everything into one corner of a 3D plot on a raw linear scale. log10(x+1)
# keeps zeros well-defined (log1p(0) = 0) while spreading out the rest.
plot_3d_df <- plot_df %>%
  filter(is.finite(Max_TPM), is.finite(.data[[enrichment_col]]), is.finite(Median_Conservation))

plot_3d <- plot_ly(
  data   = plot_3d_df,
  x      = ~log_Max_TPM,
  y      = as.formula(paste0("~", enrichment_col)),
  z      = ~Median_Conservation,
  color  = ~category,
  colors = category_colors,
  type   = "scatter3d",
  mode   = "markers",
  text   = ~paste0(
    "Gene: ", gene_name,
    "<br>Category: ", category,
    "<br>Max Expression: ", round(Max_TPM, 3), " (log10+1: ", round(log_Max_TPM, 3), ")",
    "<br>", enrichment_label, ": ", round(.data[[enrichment_col]], 3),
    "<br>Median Conservation: ", round(Median_Conservation, 3)
  ),
  hoverinfo = "text",
  marker = list(size = 4, opacity = 0.75)
) %>%
  layout(
    title = paste0("3D: ENCODE Max Expression (log-scaled) vs ", enrichment_label, " vs phyloP100 Median Conservation"),
    scene = list(
      xaxis = list(title = "log10(ENCODE Max Expression + 1)"),
      yaxis = list(title = enrichment_label),
      zaxis = list(title = "phyloP100 Median Conservation")
    ),
    legend = list(title = list(text = "Gene category"))
  )

## ---- 10. MULTI-DIMENSIONAL PLOT: phyloP vs ENCODE, shape = category, color = enrichment -----
# Packs all four dimensions into a single 2D view:
#   x-axis  = phyloP100 Median Conservation
#   y-axis  = ENCODE Max Expression (log-scaled, same as other plots)
#   shape   = gene category (Functional = star, Pseudogene = circle)
#   color   = configured enrichment source (continuous scale)

message("Building multi-dimensional plot (phyloP vs ENCODE, shaped by category, colored by ",
        enrichment_label, ")...")

plot_multidim_df <- plot_df %>%
  filter(is.finite(Max_TPM), is.finite(Median_Conservation), is.finite(.data[[enrichment_col]]))

symbol_map <- c("Functional" = "star", "Pseudogene" = "circle")

# A red-white-blue diverging palette (classic ColorBrewer "RdBu", reversed so
# red = high enrichment, white = mid-range, blue = low enrichment). This gives
# much more visual spread than a sequential palette like viridis.
# To use a full rainbow instead, swap the line below for:
#   enrichment_palette <- grDevices::rainbow(256)
enrichment_palette <- grDevices::colorRampPalette(
  c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#EF8A62", "#B2182B")
)(256)

plot_multidim <- plot_ly(
  data    = plot_multidim_df,
  x       = ~Median_Conservation,
  y       = ~log_Max_TPM,
  symbol  = ~category,
  symbols = symbol_map,
  color   = as.formula(paste0("~", enrichment_col)),
  colors  = enrichment_palette,
  type    = "scatter",
  mode    = "markers",
  text    = ~paste0(
    "Gene: ", gene_name,
    "<br>Category: ", category,
    "<br>phyloP100 Median Conservation: ", round(Median_Conservation, 3),
    "<br>ENCODE Max Expression: ", round(Max_TPM, 3), " (log10+1: ", round(log_Max_TPM, 3), ")",
    "<br>", enrichment_label, ": ", round(.data[[enrichment_col]], 3)
  ),
  hoverinfo = "text",
  marker = list(size = 11, line = list(width = 1, color = "DarkSlateGrey"))
) %>%
  layout(
    title  = list(text = paste0(
      "phyloP100 Median Conservation vs ENCODE Max Expression<br>",
      "<sup>Shape: Functional = \u2605  |  Pseudogene = \u25CF   ---   Color: ", enrichment_label, "</sup>"
    )),
    xaxis  = list(title = "phyloP100 Median Conservation"),
    yaxis  = list(title = "log10(ENCODE Max Expression + 1)"),
    legend = list(title = list(text = "Gene category (shape)"))
  ) %>%
  colorbar(title = paste0(enrichment_label, "\n"))

## ---- 11. DISPLAY PLOTS (RStudio Viewer / browser via plotly) ------------------
# Each plot is a plotly widget object. Printing it renders it interactively
# in the RStudio "Viewer" pane (or your default browser if run outside RStudio).
# No HTML files are written to disk.

message("Rendering interactive plotly plots...")

print(plot_expr_vs_conservation)
print(plot_expr_vs_enrichment)
print(plot_enrichment_vs_conservation)
print(plot_enrichment_jitter)
print(plot_pca)
print(plot_3d)
print(plot_multidim)

# If you're running this via `Rscript` (not inside RStudio) and want the plots
# to pop open in your default browser one at a time, that happens automatically
# via httpuv/browser when each plot is printed -- no extra step needed.

message("Done. Seven interactive plotly objects have been rendered (enrichment source: ",
        enrichment_col, "):")
message("  plot_expr_vs_conservation, plot_expr_vs_enrichment, plot_enrichment_vs_conservation,")
message("  plot_enrichment_jitter, plot_pca, plot_3d, plot_multidim")
message("Each remains available as a variable in your R session for further tweaking,")
message("e.g. plot_3d %>% layout(title = 'My custom title')")