#!/usr/bin/env Rscript
# ==============================================================================
# gene_analysis_plots.R
#
# Combines:
#   1. snp_enrichment_genes_pseudogenes.csv   -> enrichment source (see CONFIG)
#   2. ENCODE-expr_summary/*.csv              -> Max_TPM / Max_FPKM (expression)
#   3. phyloP100_summary/*_summary_metrics.csv-> Median_Conservation / Max_Conservation
#
# Produces interactive (plotly) HTML plots, colored by functional vs pseudogene,
# with gene_name shown on hover:
#   - ENCODE max expression vs phyloP median conservation
#   - ENCODE max expression vs enrichment source
#   - enrichment source vs phyloP median conservation
#   - PCA (PC1 vs PC2) of the three metrics
#   - 3D scatter of all three metrics together
#
# Usage:
#   Rscript gene_analysis_plots.R
#
# Edit the CONFIG section below to point at your actual file/folder locations,
# and to choose which enrichment column to use (see `enrichment_col` below).
# ==============================================================================

## ---- 0. CONFIG --------------------------------------------------------------

snp_csv_path      <- "../../data/snp_enrichment_genes_pseudogenes.csv"
encode_dir        <- "../../data/ENCODE-expr_summary"      # contains *_expr.csv (Gene,Max_TPM,Max_FPKM)
phylop_dir        <- "../../data/phyloP100_summary"         # contains *_summary_metrics.csv

# Which enrichment column in snp_csv_path to use as "the" enrichment source
# for all plots below. Must match a column name present in snp_csv_path,
# e.g. "enrichment_dbsnp", "enrichment_gnomad", "enrichment_1000genomes", etc.
enrichment_col    <- "enrichment_dbsnp"

# Human-readable label used in plot titles/axes/hover text for the column above.
# Update this alongside enrichment_col so labels stay accurate.
enrichment_label  <- "dbSNP Enrichment"

output_dir <- "gene_analysis_output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## ---- 1. PACKAGES -------------------------------------------------------------

required_pkgs <- c("dplyr", "readr", "purrr", "stringr", "tidyr",
                   "plotly", "fs", "rlang")

missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(plotly)
  library(fs)
  library(rlang)
})

## ---- 2. LOAD SNP ENRICHMENT / GENE METADATA ----------------------------------

message("Reading SNP enrichment table: ", snp_csv_path)

snp_df <- read_csv(snp_csv_path, show_col_types = FALSE) %>%
  mutate(gene_name = as.character(gene_name))

if (!enrichment_col %in% names(snp_df)) {
  stop(
    "Configured enrichment_col '", enrichment_col, "' not found in ", snp_csv_path,
    ". Available columns: ", paste(names(snp_df), collapse = ", ")
  )
}

## ---- 3. LOAD & COMBINE ENCODE EXPRESSION SUMMARIES ---------------------------

message("Reading ENCODE expression summaries from: ", encode_dir)

encode_files <- dir_ls(encode_dir, glob = "*.csv")

if (length(encode_files) == 0) {
  stop("No CSV files found in encode_dir: ", encode_dir)
}

encode_df <- map_dfr(encode_files, function(f) {
  read_csv(f, show_col_types = FALSE) %>%
    mutate(across(everything(), as.character))
}) %>%
  mutate(
    Max_TPM  = as.numeric(Max_TPM),
    Max_FPKM = as.numeric(Max_FPKM)
  ) %>%
  rename(gene_name = Gene) %>%
  # if a gene appears in multiple files, keep the max values observed
  group_by(gene_name) %>%
  summarise(
    Max_TPM  = max(Max_TPM,  na.rm = TRUE),
    Max_FPKM = max(Max_FPKM, na.rm = TRUE),
    .groups = "drop"
  )

## ---- 4. LOAD & COMBINE phyloP100 CONSERVATION SUMMARIES ----------------------

message("Reading phyloP100 summaries from: ", phylop_dir)

phylop_files <- dir_ls(phylop_dir, glob = "*_summary_metrics.csv")

if (length(phylop_files) == 0) {
  stop("No '*_summary_metrics.csv' files found in phylop_dir: ", phylop_dir)
}

phylop_df <- map_dfr(phylop_files, function(f) {
  read_csv(f, show_col_types = FALSE) %>%
    mutate(across(everything(), as.character))
}) %>%
  mutate(
    Median_Conservation = as.numeric(Median_Conservation),
    Max_Conservation    = as.numeric(Max_Conservation)
  ) %>%
  rename(gene_name = Gene) %>%
  group_by(gene_name) %>%
  summarise(
    Median_Conservation = max(Median_Conservation, na.rm = TRUE),
    Max_Conservation    = max(Max_Conservation,    na.rm = TRUE),
    .groups = "drop"
  )

## ---- 5. MERGE ALL THREE SOURCES ----------------------------------------------
# IMPORTANT: use left_join, anchored on snp_df (which contains BOTH functional
# genes and pseudogenes). An inner_join here would silently drop any gene that
# doesn't have a matching row in the ENCODE or phyloP files -- which is common
# for pseudogenes (many aren't expressed / weren't summarized) -- and that is
# why earlier versions of this script only showed functional genes.
#
# Genes with no ENCODE record are treated as having zero expression (Max_TPM /
# Max_FPKM = 0) rather than being dropped. Genes with no phyloP record keep
# NA for conservation and are only excluded from the specific plots that need
# conservation values (handled per-plot below), not from the whole dataset.

message("Merging SNP enrichment + ENCODE expression + phyloP conservation...")

combined_df <- snp_df %>%
  left_join(encode_df, by = "gene_name") %>%
  left_join(phylop_df, by = "gene_name") %>%
  mutate(
    Max_TPM  = replace_na(Max_TPM,  0),
    Max_FPKM = replace_na(Max_FPKM, 0)
  )

## ---- 6. CLASSIFY FUNCTIONAL vs PSEUDOGENE -------------------------------------
# gene_type is the definitive flag (values like "protein_coding", "pseudogene",
# etc.). Any value containing "pseudo" (case-insensitive) is labeled
# "Pseudogene"; everything else is labeled "Functional".

if (!"gene_type" %in% names(combined_df)) {
  stop("Expected column 'gene_type' not found in the SNP enrichment table.")
}

message("Classifying functional vs pseudogene using column: gene_type")
message("  Unique values in 'gene_type': ", paste(unique(combined_df$gene_type), collapse = " | "))

combined_df <- combined_df %>%
  mutate(
    category = if_else(
      str_detect(str_to_lower(gene_type), "pseudo"),
      "Pseudogene",
      "Functional"
    )
  )

message("Classification result -> Functional: ", sum(combined_df$category == "Functional"),
        ", Pseudogene: ", sum(combined_df$category == "Pseudogene"))

## ---- 7. CLEAN NUMERIC FIELDS ---------------------------------------------------
# Max_TPM is already filled with 0 for genes with no ENCODE record (step 5).
# Only drop a row here if the configured enrichment column itself is
# missing/non-finite, since that comes straight from the SNP table and every
# gene should have it. Rows missing Median_Conservation are kept in plot_df
# (still shown in TPM vs enrichment plots) and are only dropped from the
# specific plots/analyses that require conservation (2D conservation plots,
# PCA, 3D), done per-plot below.

plot_df <- combined_df %>%
  filter(is.finite(.data[[enrichment_col]])) %>%
  mutate(log_Max_TPM = log10(Max_TPM + 1))

message("Genes available for plotting: ", nrow(plot_df))
message("  Functional: ", sum(plot_df$category == "Functional"))
message("  Pseudogene: ", sum(plot_df$category == "Pseudogene"))
message("Genes missing phyloP conservation data (excluded only from conservation/PCA/3D plots): ",
        sum(!is.finite(plot_df$Median_Conservation)))

# Save the merged table for reference
write_csv(combined_df, file.path(output_dir, "combined_gene_metrics.csv"))

## ---- 8. COLOR PALETTE ---------------------------------------------------------

category_colors <- c("Functional" = "#B22222", "Pseudogene" = "#6495ED")  # brick red / cornflower blue

## ---- 9. HELPER: 2D INTERACTIVE SCATTER ----------------------------------------

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

## ---- 10. 2D PLOTS ---------------------------------------------------------------

message("Building 2D interactive plots using enrichment source: ", enrichment_col,
        " (\"", enrichment_label, "\")")

plot_expr_vs_conservation <- make_scatter(
  plot_df %>% filter(is.finite(Max_TPM), is.finite(Median_Conservation)),
  "log_Max_TPM", "Median_Conservation",
  "log10(ENCODE Max TPM + 1)", "phyloP Median Conservation",
  "ENCODE Max Expression (log-scaled) vs phyloP Median Conservation",
  x_raw_var = "Max_TPM", x_raw_label = "ENCODE Max TPM"
)

plot_expr_vs_enrichment <- make_scatter(
  plot_df %>% filter(is.finite(Max_TPM), is.finite(.data[[enrichment_col]])),
  "log_Max_TPM", enrichment_col,
  "log10(ENCODE Max TPM + 1)", enrichment_label,
  paste0("ENCODE Max Expression (log-scaled) vs ", enrichment_label),
  x_raw_var = "Max_TPM", x_raw_label = "ENCODE Max TPM"
)

plot_enrichment_vs_conservation <- make_scatter(
  plot_df %>% filter(is.finite(.data[[enrichment_col]]), is.finite(Median_Conservation)),
  enrichment_col, "Median_Conservation",
  enrichment_label, "phyloP Median Conservation",
  paste0(enrichment_label, " vs phyloP Median Conservation")
)

## ---- 11. PCA -------------------------------------------------------------------

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
    title  = paste0("PCA: Max Expression, ", enrichment_label, ", phyloP Median Conservation"),
    xaxis  = list(title = paste0("PC1 (", var_explained[1], "% variance)")),
    yaxis  = list(title = paste0("PC2 (", var_explained[2], "% variance)")),
    legend = list(title = list(text = "Gene category"))
  )

## ---- 12. 3D SCATTER --------------------------------------------------------------

message("Building 3D interactive plot...")

# ENCODE Max_TPM is typically heavily right-skewed (a few very highly
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
    "<br>Max TPM: ", round(Max_TPM, 3), " (log10+1: ", round(log_Max_TPM, 3), ")",
    "<br>", enrichment_label, ": ", round(.data[[enrichment_col]], 3),
    "<br>Median Conservation: ", round(Median_Conservation, 3)
  ),
  hoverinfo = "text",
  marker = list(size = 4, opacity = 0.75)
) %>%
  layout(
    title = paste0("3D: ENCODE Max Expression (log-scaled) vs ", enrichment_label, " vs phyloP Median Conservation"),
    scene = list(
      xaxis = list(title = "log10(ENCODE Max TPM + 1)"),
      yaxis = list(title = enrichment_label),
      zaxis = list(title = "phyloP Median Conservation")
    ),
    legend = list(title = list(text = "Gene category"))
  )

## ---- 13. MULTI-DIMENSIONAL PLOT: phyloP vs ENCODE, shape = category, color = enrichment -----
# Packs all four dimensions into a single 2D view:
#   x-axis  = phyloP Median Conservation
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
    "<br>phyloP Median Conservation: ", round(Median_Conservation, 3),
    "<br>ENCODE Max TPM: ", round(Max_TPM, 3), " (log10+1: ", round(log_Max_TPM, 3), ")",
    "<br>", enrichment_label, ": ", round(.data[[enrichment_col]], 3)
  ),
  hoverinfo = "text",
  marker = list(size = 11, line = list(width = 1, color = "DarkSlateGrey"))
) %>%
  layout(
    title  = list(text = paste0(
      "phyloP Median Conservation vs ENCODE Max Expression<br>",
      "<sup>Shape: Functional = \u2605  |  Pseudogene = \u25CF   ---   Color: ", enrichment_label, "</sup>"
    )),
    xaxis  = list(title = "phyloP Median Conservation"),
    yaxis  = list(title = "log10(ENCODE Max TPM + 1)"),
    legend = list(title = list(text = "Gene category (shape)"))
  ) %>%
  colorbar(title = paste0(enrichment_label, "\n"))

## ---- 14. DISPLAY PLOTS (RStudio Viewer / browser via plotly) ------------------
# Each plot is a plotly widget object. Printing it renders it interactively
# in the RStudio "Viewer" pane (or your default browser if run outside RStudio).
# No HTML files are written to disk.

message("Rendering interactive plotly plots...")

print(plot_expr_vs_conservation)
print(plot_expr_vs_enrichment)
print(plot_enrichment_vs_conservation)
print(plot_pca)
print(plot_3d)
print(plot_multidim)

# If you're running this via `Rscript` (not inside RStudio) and want the plots
# to pop open in your default browser one at a time, that happens automatically
# via httpuv/browser when each plot is printed -- no extra step needed.

message("Done. Six interactive plotly objects have been rendered (enrichment source: ",
        enrichment_col, "):")
message("  plot_expr_vs_conservation, plot_expr_vs_enrichment, plot_enrichment_vs_conservation,")
message("  plot_pca, plot_3d, plot_multidim")
message("Each remains available as a variable in your R session for further tweaking,")
message("e.g. plot_3d %>% layout(title = 'My custom title')")