#!/usr/bin/env Rscript
# ============================================================================
# Multi-method dimension reduction of gene conservation/variation features
#
# Runs the dimension-reduction method(s) YOU choose on the same feature set
# and saves one plot per method to results/plots/, so they can be compared
# side-by-side for how well each separates pseudogenes from functional genes.
#
# Available methods (use these names on the command line):
#
#   pca       - PCA: unsupervised, linear, maximizes overall variance
#   lda       - LDA (supervised): linear, maximizes between-class separation.
#               Fit on 3 groups: Functional (all), Pseudogene with
#               functional_probability >= 0.5, Pseudogene with < 0.5. 3
#               classes -> 2 discriminant axes (LD1 vs LD2). Functional
#               genes with probability < 0.5 are not a class of their own
#               but are still plotted (orange outline) and labelled (blue).
#               (Also writes lda_pseudogene_plot_coordinates.csv, which the
#               interactive plotly script reads.)
#   tsne      - t-SNE: unsupervised, nonlinear, local neighborhoods
#   umap      - UMAP (unsupervised): nonlinear, neighbor-graph based
#   mds       - Classic MDS: unsupervised, linear (~ PCA on Euclidean dists)
#   lle       - LLE: unsupervised, nonlinear, preserves local linear weights
#   sammon    - Sammon mapping: nonlinear MDS variant that weights small
#               (nearby) distances much more heavily than large ones
#   umap_sup  - UMAP (supervised): blends in Gene_Type via `y` +
#               `target_weight`
#   all       - every method above (default when no argument is given)
#
# Aesthetics used in every plot:
#   1. Shape  -> Gene_Type (Functional = solid circle, Pseudogene = open diamond)
#                (LDA: 3-way LDA_Group instead: Functional = solid circle,
#                 Pseudogene >= 0.5 = hollow diamond, Pseudogene < 0.5 =
#                 solid diamond)
#      Functional genes with probability < 0.5 get a orange outline
#   2. Color  -> functional_probability (continuous, blue = low, red = high)
#   Fixed alpha = 0.8 on all points.
#   3. Text labels (with leader lines) for the genes listed in LABEL_GENES.
#   4. Font: Times New Roman.
#
# Input features: enrichment_intergenic_dbsnp, PhyloP100_median, ENCODE_max,
# and the type1/2/3 promoter score + motif-count columns. `label` and
# `functional_probability` are aesthetic-only and excluded as features.
#
# Requires: ggplot2, dplyr, readr, tibble, ggrepel
#   Per method: MASS (lda, sammon; ships with R), Rtsne (tsne),
#               uwot (umap, umap_sup), lle (lle)
#   Recommended: ragg (reliable Times New Roman in saved PNGs)
#   install.packages(c("ggrepel", "Rtsne", "uwot", "lle", "ragg"))
#   (if lle isn't on CRAN for you:
#    remotes::install_github("cran/lle"))
#
# Usage:
#   Rscript mds_pseudogene_plot_edited.R                # all methods
#   Rscript mds_pseudogene_plot_edited.R pca            # just PCA
#   Rscript mds_pseudogene_plot_edited.R pca lda umap   # any combination
#   (run interactively with no args -> a picker pops up)
#
#               Also writes lda_pseudogene_plot_interactive.html (plotly).
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
})

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
ALL_METHODS <- c("pca", "lda", "tsne", "umap", "mds", "lle", "sammon", "umap_sup")

FONT_FAMILY <- "Times New Roman"

# Text-label colours: Conflicted Functional genes get a blue label so they stand
# out; every other labelled gene stays black.
LABEL_COLOR_DEFAULT     <- "black"
LABEL_COLOR_FALSE_FUNC  <- "orange"

# Label offset: how far labels are pushed away from their dots (larger =
# further; units are lines of text).
LABEL_POINT_PADDING <- 0.3   # gap between a label and its dot (back to the original value)
LABEL_BOX_PADDING   <- 0.4   # gap between neighbouring labels (back to the original value)

LABEL_GENES <- c(
  "RNU6ATAC3P", "RN7SKP70", "RNU5E-6P", "RNU5A-8P", "RNU6-1334P", "RNU11-3P",
  "RNU6-1194P", "RNU11-6P", "RNU4ATAC15P", "RNU12-2P",
  "RNU1-28P", "RNU1-27P", "VTRNA2-1", "RN7SKP1", "VTRNA1-1",
  "RN7SL3", "VTRNA1-3", "RNU1-4", "RNU1-3"
)

# Genes to label on the TRNA-only plots (LABEL_GENES are non-TRNA, so they
# would never appear there). Leave empty for no labels.
TRNA_LABEL_GENES <- character(0)
# ---------------------------------------------------------------------------
# Choose which method(s) to run
# ---------------------------------------------------------------------------
args <- tolower(commandArgs(trailingOnly = TRUE))

if (length(args) == 0) {
  methods <- if (interactive()) {
    picked <- select.list(c(ALL_METHODS, "all"), multiple = TRUE,
                          title = "Dimension reduction method(s) to run:")
    if (length(picked) == 0) stop("No method selected.")
    picked
  } else {
    "all"
  }
} else {
  methods <- args
}
if ("all" %in% methods) methods <- ALL_METHODS

bad <- setdiff(methods, ALL_METHODS)
if (length(bad) > 0) {
  stop(sprintf("Unknown method(s): %s\nChoose from: %s, all",
               paste(bad, collapse = ", "), paste(ALL_METHODS, collapse = ", ")))
}
methods <- ALL_METHODS[ALL_METHODS %in% methods]  # canonical order, no dupes
run <- function(m) m %in% methods
message("Running method(s): ", paste(methods, collapse = ", "))

# ---------------------------------------------------------------------------
# Package checks (only for what was selected)
# ---------------------------------------------------------------------------
pkg_needed <- c("ggrepel")
if (run("lda") || run("sammon"))    pkg_needed <- c(pkg_needed, "MASS")
if (run("lda"))                     pkg_needed <- c(pkg_needed, "plotly", "htmlwidgets")
if (run("tsne"))                    pkg_needed <- c(pkg_needed, "Rtsne")
if (run("umap") || run("umap_sup")) pkg_needed <- c(pkg_needed, "uwot")
if (run("lle"))                     pkg_needed <- c(pkg_needed, "lle")

for (pkg in pkg_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "lle") {
      message('If lle is not on CRAN for you: remotes::install_github("cran/lle")')
    }
    stop(sprintf("This script requires the '%s' package. Install it with: install.packages(\"%s\")", pkg, pkg))
  }
}

# ragg gives reliable system-font (Times New Roman) rendering in PNGs
if (requireNamespace("ragg", quietly = TRUE)) {
  png_device <- ragg::agg_png
} else {
  warning("Package 'ragg' not installed; Times New Roman may not render in the saved PNGs. ",
          "install.packages('ragg') is recommended.")
  png_device <- "png"
}

# ---------------------------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------------------------
input_csv  <- "../../results/combined_gene_data_with_snp_promoters.csv"

if (!file.exists(input_csv)) {
  stop(sprintf("Input file not found: %s", input_csv))
}

output_dir <- "../../results/plots"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
df <- read_csv(input_csv, show_col_types = FALSE)

required_cols <- c(
  "Gene", "PhyloP100_median", "ENCODE_max", "enrichment_intergenic_dbsnp",
  "type1_score", "type1_motif_counts", "type2_score", "type2_motif_counts",
  "type3_score", "type3_motif_counts",
  "Gene_Type", "Gene_group", "label", "functional_probability"
)

missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing expected column(s): %s", paste(missing_cols, collapse = ", ")))
}

# fill in any NA in PhyloP100_median and ENCODE_max with 0
df <- df %>%
  mutate(
    PhyloP100_median = replace_na(PhyloP100_median, 0),
    ENCODE_max = replace_na(ENCODE_max, 0)
  )


# ---------------------------------------------------------------------------
# 2. Feature matrix shared by every method
# ---------------------------------------------------------------------------
feature_cols <- c(
  "PhyloP100_median", "ENCODE_max", "enrichment_intergenic_dbsnp",
  "type1_score", "type1_motif_counts",
  "type2_score", "type2_motif_counts",
  "type3_score", "type3_motif_counts"
)

df_complete <- df %>%
  filter(if_all(all_of(feature_cols), ~ !is.na(.))) %>%
  filter(!is.na(Gene_Type))

n_dropped <- nrow(df) - nrow(df_complete)
if (n_dropped > 0) {
  message(sprintf("Dropping %d row(s) with missing values in feature/Gene_Type columns.", n_dropped))
}

if (nrow(df_complete) < 5) {
  stop("Not enough complete rows to compute stable embeddings (need at least 5, ideally many more).")
}

df_complete <- df_complete %>%
  mutate(Gene_Type = trimws(Gene_Type)) %>%
  mutate(is_TRNA = toupper(Gene_group) == "TRNA") %>%
  mutate(Gene_Type = factor(Gene_Type, levels = c("Functional", "Pseudogene")))

if (any(is.na(df_complete$Gene_Type))) {
  warning("Some Gene_Type values did not match 'Functional' or 'Pseudogene' exactly; ",
          "check spelling/casing in the input CSV.")
}

feature_matrix <- df_complete %>%
  select(all_of(feature_cols)) %>%
  as.matrix()

# Z-scale each feature so no single column dominates distance/variance
# calculations just because of its scale.
feature_scaled <- scale(feature_matrix)

metadata <- df_complete %>%
  select(Gene, Gene_Type, Gene_group, is_TRNA,
         enrichment_intergenic_dbsnp, PhyloP100_median, ENCODE_max,
         functional_probability)

n_genes <- nrow(feature_scaled)

# Shared across several methods
dist_matrix     <- dist(feature_scaled, method = "euclidean")  # mds, sammon
n_neighbors_use <- min(15, n_genes - 1)                        # umap, umap_sup

# Report label genes that aren't in the data at all
missing_labels <- setdiff(LABEL_GENES, df_complete$Gene)
if (length(missing_labels) > 0) {
  message(sprintf("Label genes not found in the data (%d): %s",
                  length(missing_labels), paste(missing_labels, collapse = ", ")))
}

# ---------------------------------------------------------------------------
# 3. Shared plotting function
# ---------------------------------------------------------------------------
make_plot <- function(embed_df, x_col, y_col, title, subtitle, out_path,
                      shape_col = "Gene_Type",
                      shape_values = c("Functional" = 16, "Pseudogene" = 5),
                      shape_name = "Gene Type",
                      trna_only = FALSE,            # TRUE -> plot ONLY the TRNA dots
                      label_genes = LABEL_GENES,    # which genes to label
                      write_coords = TRUE) {        # FALSE for the TRNA plot so the
  # full-coordinates CSV isn't duplicated
  # Draw order (ggplot draws rows/layers in order, later ends up on top):
  #   1. everything else, low -> high functional_probability
  #   2. Conflicted Pseudogene (Pseudogene with probability >= 0.5), on top of 1
  #   3. Conflicted Functional (Functional gene with probability < 0.5), drawn
  #      LAST in its own layer as a SOLID circle (probability colour) with a
  #      orange outline (shape 21: border = black, fill = probability colour)
  #   (same definition as the LDA Match_Group, applied to every method)
  plot_df <- embed_df %>%
    filter(if (trna_only) coalesce(is_TRNA, FALSE) else !coalesce(is_TRNA, FALSE)) %>%
    mutate(
      .false_func = coalesce(Gene_Type == "Functional" & functional_probability <  0.5, FALSE),
      .false_pseu = coalesce(Gene_Type == "Pseudogene" & functional_probability >= 0.5, FALSE)
    )
  
  main_df <- plot_df %>%
    filter(!.false_func) %>%
    arrange(.false_pseu, functional_probability)
  
  ff_df <- plot_df %>%
    filter(.false_func) %>%
    arrange(functional_probability)
  
  label_df <- plot_df %>%
    filter(Gene %in% label_genes) %>%
    mutate(.label_col = ifelse(.false_func, LABEL_COLOR_FALSE_FUNC, LABEL_COLOR_DEFAULT),
           .label_face = ifelse(.false_func, "bold", "plain"))
  message(sprintf("  plotting %d point(s); labeling %d (%d Conflicted Functional in orange)",
                  nrow(plot_df), nrow(label_df), sum(label_df$.false_func)))
  
  p <- ggplot(main_df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(
      aes(shape = .data[[shape_col]], color = functional_probability),
      size = 3, stroke = 1.2, alpha = 0.8
    ) +
    # Conflicted Functional: last layer, solid circle + orange outline
    geom_point(
      data = ff_df,
      aes(fill = functional_probability),
      shape = 21, color = "orange",
      size = 3, stroke = 0.9, alpha = 0.8
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = Gene, colour = I(.label_col), fontface = .label_face),  # I(): literal colours, not the gradient scale
      family = FONT_FAMILY,      # geoms don't inherit the theme font
      size = 3.2,
      max.overlaps = Inf,        # label every requested gene, however crowded
      min.segment.length = 0,    # always draw the leader line
      segment.size = 0.3,
      segment.color = "grey35",
      box.padding = LABEL_BOX_PADDING,
      point.padding = LABEL_POINT_PADDING,
      seed = 42
    ) +
    scale_shape_manual(
      values = shape_values,  # solid circle vs open diamond by default
      name = shape_name
    ) +
    scale_color_gradient(
      name = "Functional\nProbability",
      low = "#1f77b4", high = "#d62728",
      limits = c(0, 1)
    ) +
    # Same blue -> red mapping for the fill of the outlined layer; hidden so
    # only one colorbar is shown.
    scale_fill_gradient(
      low = "#1f77b4", high = "#d62728",
      limits = c(0, 1),
      guide = "none"
    ) +
    guides(
      shape = guide_legend(override.aes = list(color = "black", alpha = 1))
    ) +
    labs(title = title, subtitle = subtitle, x = x_col, y = y_col) +
    theme_minimal(base_size = 13, base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.key = element_rect(fill = "white"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(out_path, plot = p, width = 9, height = 7, dpi = 300, device = png_device)
  message(sprintf("Saved plot to: %s", out_path))
  
  if (write_coords) {
    coords_out <- sub("\\.[^.]+$", "_coordinates.csv", out_path)
    write_csv(embed_df, coords_out)
    message(sprintf("Saved coordinates to: %s", coords_out))
  }
  
  invisible(p)
}

# ---------------------------------------------------------------------------
# 4. PCA (unsupervised, linear, maximizes overall variance)
# ---------------------------------------------------------------------------
if (run("pca")) {
  message("\n== PCA ==")
  pca_fit <- prcomp(feature_scaled, center = FALSE, scale. = FALSE)
  pca_var_pct <- round(100 * (pca_fit$sdev^2 / sum(pca_fit$sdev^2))[1:2], 1)
  
  pca_df <- bind_cols(metadata, as.data.frame(pca_fit$x[, 1:2])) %>%
    rename(PC1 = PC1, PC2 = PC2)
  
  make_plot(
    pca_df, "PC1", "PC2",
    title = "PCA of Gene Conservation / Variation Features",
    subtitle = sprintf("n = %d annotations | PC1: %.1f%% var | PC2: %.1f%% var", n_genes, pca_var_pct[1], pca_var_pct[2]),
    out_path = file.path(output_dir, "pca_pseudogene_plot.png")
  )
}

# ---------------------------------------------------------------------------
# 5. LDA (supervised, linear, maximizes between-class separation)
#    Fit on 3 groups (Conflicted Functional is far too rare -- 1 gene -- to be its
#    own class, so Functional genes are kept together as one class):
#      Functional       = every Functional gene (any functional_probability)
#      Conflicted Pseudogene = Pseudogene, functional_probability >= 0.5
#      Pseudogene  = Pseudogene, functional_probability <  0.5
#    (exactly 0.5 goes to the ">= 0.5" side.) 3 classes -> 2 discriminant
#    axes (LD1, LD2), both plotted.
#
#    Functional genes with functional_probability < 0.5 ("Conflicted Functional")
#    are NOT a separate class, but are still plotted and labelled: make_plot()
#    draws them last as solid circles with a orange outline, with blue labels.
#
#    The saved coordinates CSV keeps the old 4-way `Match_Group` column
#    (True/Conflicted Functional, Conflicted/Pseudogene) as a descriptive
#    annotation, so lda_pseudogene_plot_interactive.R keeps working. The
#    3-class label used for the fit and the plot shapes is `LDA_Group`.
# ---------------------------------------------------------------------------
if (run("lda")) {
  message("\n== LDA ==")
  lda_mask <- !is.na(df_complete$functional_probability)
  if (any(!lda_mask)) {
    message(sprintf(
      "LDA: dropping %d row(s) with missing functional_probability.",
      sum(!lda_mask)
    ))
  }
  
  feature_scaled_lda <- feature_scaled[lda_mask, , drop = FALSE]
  metadata_lda <- metadata[lda_mask, ]
  
  # 3-class label used for the fit and for the plot shapes
  lda_group <- case_when(
    metadata_lda$Gene_Type == "Functional" ~ "Functional",
    metadata_lda$Gene_Type == "Pseudogene" & metadata_lda$functional_probability >= 0.5 ~ "Conflicted Pseudogene",
    metadata_lda$Gene_Type == "Pseudogene" & metadata_lda$functional_probability <  0.5 ~ "Pseudogene"
  )
  lda_group <- factor(
    lda_group,
    levels = c("Functional", "Conflicted Pseudogene", "Pseudogene")
  )
  
  # Descriptive 4-way label (annotation only; NOT used for the fit)
  match_group <- case_when(
    metadata_lda$Gene_Type == "Pseudogene" & metadata_lda$functional_probability <  0.5 ~ "Pseudogene",
    metadata_lda$Gene_Type == "Pseudogene" & metadata_lda$functional_probability >= 0.5 ~ "Conflicted Pseudogene",
    metadata_lda$Gene_Type == "Functional" & metadata_lda$functional_probability <  0.5 ~ "Conflicted Functional",
    metadata_lda$Gene_Type == "Functional" & metadata_lda$functional_probability >= 0.5 ~ "True Functional"
  )
  match_group <- factor(
    match_group,
    levels = c("True Functional", "Conflicted Functional", "Conflicted Pseudogene", "Pseudogene")
  )
  
  group_counts <- table(lda_group)
  message("LDA group sizes: ",
          paste(names(group_counts), group_counts, sep = "=", collapse = ", "))
  if (any(group_counts < 2)) {
    warning("LDA: at least one group has fewer than 2 members ",
            "(counts: ", paste(names(group_counts), group_counts, sep = "=", collapse = ", "),
            "); LDA may be unstable or fail.")
  }
  
  lda_fit <- MASS::lda(x = feature_scaled_lda, grouping = lda_group)
  lda_scores <- predict(lda_fit)$x  # n x 2 matrix (LD1, LD2) for 3 classes
  
  # Proportion of between-class variance ("trace") captured by each axis
  lda_prop_trace <- round(100 * lda_fit$svd^2 / sum(lda_fit$svd^2), 1)
  
  lda_df <- bind_cols(
    metadata_lda,
    tibble(
      LD1 = as.numeric(lda_scores[, 1]),
      LD2 = as.numeric(lda_scores[, 2]),
      LDA_Group = lda_group,
      Match_Group = match_group
    )
  )
  
  lda_shape_values <- c(
    "Functional"            = 16,  # solid circle (prob < 0.5 gets an orange outline)
    "Conflicted Pseudogene" = 5,   # hollow diamond: Pseudogene, prob >= 0.5
    "Pseudogene"            = 18   # solid diamond:  Pseudogene, prob <  0.5
  )
  n_trna <- sum(lda_df$is_TRNA, na.rm = TRUE)
  
  # --- Static: non-TRNA (original plot) ---
  make_plot(
    lda_df, "LD1", "LD2",
    title = "LDA of Pseudogenes and Functional ncRNA",
    subtitle = sprintf(
      "n = %d annotations | 3-group LDA | LD1: %.1f%% | LD2: %.1f%% trace | orange = Functional with probability < 0.5",
      nrow(lda_df) - n_trna, lda_prop_trace[1], lda_prop_trace[2]
    ),
    out_path = file.path(output_dir, "lda_pseudogene_plot.png"),
    shape_col = "LDA_Group",
    shape_values = lda_shape_values,
    shape_name = "Classification"
  )
  
  # --- Static: TRNA only (same LDA fit / axes) ---
  if (n_trna > 0) {
    make_plot(
      lda_df, "LD1", "LD2",
      title = "LDA of tRNA Genes",
      subtitle = sprintf(
        "n = %d tRNA annotations | fit on all annotations | LD1: %.1f%% | LD2: %.1f%% trace | orange = Functional (p < 0.5)",
        n_trna, lda_prop_trace[1], lda_prop_trace[2]
      ),
      out_path = file.path(output_dir, "lda_trna_plot.png"),
      shape_col = "LDA_Group",
      shape_values = lda_shape_values,
      shape_name = "Classification",
      trna_only = TRUE,
      label_genes = TRNA_LABEL_GENES,
      write_coords = FALSE
    )
  } else {
    message("  no TRNA rows found; skipping TRNA-only LDA plot.")
  }
  
  # -------------------------------------------------------------------------
  # Interactive plotly versions (same colours / symbols as the static plots)
  # -------------------------------------------------------------------------
  prob_colorscale <- list(c(0, "#1f77b4"), c(1, "#d62728"))
  
  build_lda_interactive <- function(keep_trna, label_genes, title_main, html_out) {
    ly_df <- lda_df %>%
      filter(if (keep_trna) coalesce(is_TRNA, FALSE) else !coalesce(is_TRNA, FALSE)) %>%
      mutate(
        .false_func = coalesce(Gene_Type == "Functional" & functional_probability < 0.5, FALSE),
        .hover = sprintf(
          paste0("<b>%s</b><br>Group: %s<br>Match group: %s<br>",
                 "Functional probability: %.3f<br>LD1: %.3f<br>LD2: %.3f<br>",
                 "PhyloP100 median: %.3f<br>ENCODE max: %.3f<br>",
                 "dbSNP enrichment: %.3f"),
          Gene, LDA_Group, Match_Group, functional_probability, LD1, LD2,
          PhyloP100_median, ENCODE_max, enrichment_intergenic_dbsnp
        )
      )
    
    if (nrow(ly_df) == 0) {
      message("  no rows for ", basename(html_out), "; skipping.")
      return(invisible(NULL))
    }
    
    add_group_trace <- function(p, d, name, symbol, outline = NULL, show_bar = FALSE) {
      d <- d %>% arrange(functional_probability)
      marker <- list(
        symbol = symbol, size = 10, opacity = 0.8,
        color = d$functional_probability,
        colorscale = prob_colorscale, cmin = 0, cmax = 1,
        showscale = show_bar,
        line = if (is.null(outline)) list(width = 1.5, color = NULL)
        else list(width = 1.5, color = outline)
      )
      if (show_bar) {
        marker$colorbar <- list(title = list(text = "Functional<br>Probability"),
                                len = 0.6)
      }
      plotly::add_trace(
        p, data = d, x = ~LD1, y = ~LD2,
        type = "scatter", mode = "markers",
        name = name, text = ~.hover, hoverinfo = "text",
        marker = marker
      )
    }
    
    # Trace specs; the colorbar goes on the first NON-EMPTY trace (the TRNA
    # subset may have no plain "Functional" points).
    specs <- list(
      list(d = filter(ly_df, LDA_Group == "Functional" & !.false_func),
           name = "Functional", symbol = "circle", outline = NULL),
      list(d = filter(ly_df, LDA_Group == "Pseudogene"),
           name = "Pseudogene", symbol = "diamond", outline = NULL),
      list(d = filter(ly_df, LDA_Group == "Conflicted Pseudogene"),
           name = "Conflicted Pseudogene", symbol = "diamond-open", outline = NULL),
      list(d = filter(ly_df, .false_func),
           name = "Functional (prob < 0.5)", symbol = "circle", outline = "orange")
    )
    
    p_ly <- plotly::plot_ly()
    bar_done <- FALSE
    for (s in specs) {
      if (nrow(s$d) == 0) next
      p_ly <- add_group_trace(p_ly, s$d, s$name, s$symbol, s$outline,
                              show_bar = !bar_done)
      bar_done <- TRUE
    }
    
    lab_df <- ly_df %>% filter(Gene %in% label_genes)
    label_annots <- lapply(seq_len(nrow(lab_df)), function(i) {
      list(
        x = lab_df$LD1[i], y = lab_df$LD2[i],
        text = lab_df$Gene[i],
        showarrow = TRUE, arrowhead = 0, arrowwidth = 0.6, arrowcolor = "grey35",
        ax = 25, ay = -25,
        font = list(family = FONT_FAMILY, size = 11,
                    color = if (lab_df$.false_func[i]) LABEL_COLOR_FALSE_FUNC
                    else LABEL_COLOR_DEFAULT),
        bgcolor = "rgba(255,255,255,0.0)"
      )
    })
    
    p_ly <- p_ly %>%
      plotly::layout(
        title = list(
          text = sprintf(
            "<b>%s</b><br><sup>n = %d annotations | 3-group LDA | LD1: %.1f%% | LD2: %.1f%% trace | orange outline = Functional with probability &lt; 0.5</sup>",
            title_main, nrow(ly_df), lda_prop_trace[1], lda_prop_trace[2]),
          x = 0.01
        ),
        xaxis = list(title = "LD1", zeroline = FALSE),
        yaxis = list(title = "LD2", zeroline = FALSE),
        font = list(family = FONT_FAMILY, size = 13),
        legend = list(title = list(text = "Classification")),
        annotations = label_annots,
        plot_bgcolor = "white", paper_bgcolor = "white",
        margin = list(t = 90)
      ) %>%
      plotly::config(displaylogo = FALSE)
    
    tryCatch(
      htmlwidgets::saveWidget(p_ly, html_out, selfcontained = TRUE),
      error = function(e) {
        message("  selfcontained save failed (pandoc missing?): ", conditionMessage(e),
                "\n  saving with a companion _files folder instead.")
        htmlwidgets::saveWidget(p_ly, html_out, selfcontained = FALSE)
      }
    )
    message(sprintf("Saved interactive plot to: %s", html_out))
    invisible(p_ly)
  }
  
  message("  building interactive plotly figures")
  build_lda_interactive(
    keep_trna = FALSE, label_genes = LABEL_GENES,
    title_main = "LDA of Pseudogenes and Functional ncRNA",
    html_out = file.path(output_dir, "lda_pseudogene_plot_interactive.html")
  )
  build_lda_interactive(
    keep_trna = TRUE, label_genes = TRNA_LABEL_GENES,
    title_main = "LDA of tRNA Genes",
    html_out = file.path(output_dir, "lda_trna_plot_interactive.html")
  )
}


# ---------------------------------------------------------------------------
# 6. t-SNE (unsupervised, nonlinear, local neighborhood structure)
# ---------------------------------------------------------------------------
if (run("tsne")) {
  message("\n== t-SNE ==")
  perplexity_use <- max(2, min(30, floor((n_genes - 1) / 3)))
  set.seed(42)
  tsne_fit <- Rtsne::Rtsne(
    feature_scaled,
    dims = 2,
    perplexity = perplexity_use,
    check_duplicates = FALSE,
    verbose = FALSE
  )
  
  tsne_df <- bind_cols(metadata, as.data.frame(tsne_fit$Y))
  colnames(tsne_df)[(ncol(tsne_df) - 1):ncol(tsne_df)] <- c("tSNE1", "tSNE2")
  
  make_plot(
    tsne_df, "tSNE1", "tSNE2",
    title = "t-SNE of Gene Conservation / Variation Features",
    subtitle = sprintf("n = %d annotations | perplexity = %d | unsupervised", n_genes, perplexity_use),
    out_path = file.path(output_dir, "tsne_pseudogene_plot.png")
  )
}

# ---------------------------------------------------------------------------
# 7. UMAP (unsupervised, nonlinear neighbor-graph based)
# ---------------------------------------------------------------------------
if (run("umap")) {
  message("\n== UMAP ==")
  set.seed(42)
  umap_fit <- uwot::umap(
    feature_scaled,
    n_neighbors = n_neighbors_use,
    n_components = 2,
    metric = "euclidean",
    min_dist = 0.1,
    ret_model = FALSE
  )
  
  umap_df <- bind_cols(metadata, as.data.frame(umap_fit))
  colnames(umap_df)[(ncol(umap_df) - 1):ncol(umap_df)] <- c("UMAP1", "UMAP2")
  
  make_plot(
    umap_df, "UMAP1", "UMAP2",
    title = "UMAP of Gene Conservation / Variation Features",
    subtitle = sprintf("n = %d annotations | n_neighbors = %d | unsupervised", n_genes, n_neighbors_use),
    out_path = file.path(output_dir, "umap_pseudogene_plot.png")
  )
}

# ---------------------------------------------------------------------------
# 8. Classic MDS (unsupervised, linear; ~ PCA on Euclidean distances)
# ---------------------------------------------------------------------------
if (run("mds")) {
  message("\n== Classic MDS ==")
  mds_fit <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  mds_var_pct <- round(100 * mds_fit$eig[1:2] / sum(abs(mds_fit$eig)), 1)
  
  mds_df <- bind_cols(metadata, as.data.frame(mds_fit$points))
  colnames(mds_df)[(ncol(mds_df) - 1):ncol(mds_df)] <- c("MDS1", "MDS2")
  
  make_plot(
    mds_df, "MDS1", "MDS2",
    title = "Classic MDS of Gene Conservation / Variation Features",
    subtitle = sprintf("n = %d annotations | Axis 1: %.1f%% var | Axis 2: %.1f%% var", n_genes, mds_var_pct[1], mds_var_pct[2]),
    out_path = file.path(output_dir, "mds_pseudogene_plot.png")
  )
}

# ---------------------------------------------------------------------------
# 9. LLE (unsupervised, nonlinear manifold method)
#    Reconstructs each point as a weighted combination of its k nearest
#    neighbors, then finds a low-D embedding preserving those local weights.
# ---------------------------------------------------------------------------
if (run("lle")) {
  message("\n== LLE ==")
  lle_k_use <- min(15, n_genes - 1)
  set.seed(42)
  lle_fit <- lle::lle(
    X = feature_scaled,
    m = 2,
    k = lle_k_use,
    v = 0.9  # regularization for the reconstruction-weight solve
  )
  
  lle_df <- bind_cols(metadata, as.data.frame(lle_fit$Y))
  colnames(lle_df)[(ncol(lle_df) - 1):ncol(lle_df)] <- c("LLE1", "LLE2")
  
  make_plot(
    lle_df, "LLE1", "LLE2",
    title = "LLE of Gene Conservation / Variation Features",
    subtitle = sprintf("n = %d annotations | k = %d neighbors | unsupervised", n_genes, lle_k_use),
    out_path = file.path(output_dir, "lle_pseudogene_plot.png")
  )
}

# ---------------------------------------------------------------------------
# 10. Sammon mapping (unsupervised, nonlinear MDS variant)
#    Stress divides each squared error by the original distance, so nearby
#    points are weighted far more heavily. sammon() needs strictly positive
#    distances, so exact duplicate rows (if any) get a tiny jitter first.
# ---------------------------------------------------------------------------
if (run("sammon")) {
  message("\n== Sammon mapping ==")
  if (any(dist_matrix == 0)) {
    message("Sammon mapping: duplicate rows detected in feature_scaled; ",
            "adding a small jitter (1e-6 sd) to break exact ties.")
    set.seed(42)
    feature_scaled_sammon <- feature_scaled +
      matrix(rnorm(length(feature_scaled), sd = 1e-6), nrow = nrow(feature_scaled))
    dist_matrix_sammon <- dist(feature_scaled_sammon, method = "euclidean")
  } else {
    dist_matrix_sammon <- dist_matrix
  }
  
  sammon_init <- cmdscale(dist_matrix_sammon, k = 2)  # classic MDS start
  set.seed(42)
  sammon_fit <- MASS::sammon(dist_matrix_sammon, y = sammon_init, k = 2, trace = FALSE)
  
  sammon_df <- bind_cols(metadata, as.data.frame(sammon_fit$points))
  colnames(sammon_df)[(ncol(sammon_df) - 1):ncol(sammon_df)] <- c("Sammon1", "Sammon2")
  
  make_plot(
    sammon_df, "Sammon1", "Sammon2",
    title = "Sammon Mapping of Gene Conservation / Variation Features",
    subtitle = sprintf("n = %d annotations | final stress = %.4f | unsupervised", n_genes, sammon_fit$stress),
    out_path = file.path(output_dir, "sammon_pseudogene_plot.png")
  )
}

# ---------------------------------------------------------------------------
# 11. UMAP (supervised)
#    Same machinery as unsupervised UMAP, but Gene_Type is passed as `y`.
#    `target_weight`: 0 = ignore labels, 1 = collapse by label almost like a
#    classifier. 0.5 (uwot default) gives labels real influence without
#    discarding feature-space geometry.
# ---------------------------------------------------------------------------
if (run("umap_sup")) {
  message("\n== UMAP (supervised) ==")
  set.seed(42)
  umap_sup_fit <- uwot::umap(
    feature_scaled,
    y = df_complete$Gene_Type,
    n_neighbors = n_neighbors_use,
    n_components = 2,
    metric = "euclidean",
    target_metric = "categorical",
    target_weight = 0.5,
    min_dist = 0.1,
    ret_model = FALSE
  )
  
  umap_sup_df <- bind_cols(metadata, as.data.frame(umap_sup_fit))
  colnames(umap_sup_df)[(ncol(umap_sup_df) - 1):ncol(umap_sup_df)] <- c("UMAP_sup1", "UMAP_sup2")
  
  make_plot(
    umap_sup_df, "UMAP_sup1", "UMAP_sup2",
    title = "UMAP (Supervised) of Gene Conservation / Variation Features",
    subtitle = sprintf("n = %d annotations | n_neighbors = %d | target_weight = 0.5 | supervised on Gene_Type", n_genes, n_neighbors_use),
    out_path = file.path(output_dir, "umap_supervised_pseudogene_plot.png")
  )
}

message("\nDone. Plots saved to: ", output_dir, " (methods: ", paste(methods, collapse = ", "), ")")