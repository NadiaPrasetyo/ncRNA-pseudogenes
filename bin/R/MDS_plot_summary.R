#!/usr/bin/env Rscript
# ============================================================================
# Multi-method dimension reduction of gene conservation/variation features
#
# Runs several dimension-reduction methods on the same feature set and
# saves one plot per method to results/plots/, so they can be compared
# side-by-side for how well each separates pseudogenes from functional genes:
#
#   1. PCA               - unsupervised, linear, maximizes overall variance
#   2. LDA (supervised)  - supervised, linear, maximizes between-class
#                           separation. Grouping is a 4-way cross of
#                           Gene_Type x whether functional_probability
#                           agrees with that label (threshold 0.5): True/
#                           False Functional, True/False Pseudogene. 4
#                           classes -> up to 3 discriminant axes; LD1 vs
#                           LD2 plotted as a normal 2D scatter -- see note
#                           below
#   3. t-SNE              - unsupervised, nonlinear, emphasizes local
#                            neighborhoods
#   4. UMAP (unsupervised) - unsupervised, nonlinear neighbor-graph based
#   5. Classic MDS         - unsupervised, linear (mathematically ~ PCA on
#                            Euclidean distances)
#   6. LLE                 - unsupervised, nonlinear, preserves local
#                            linear reconstruction weights (manifold method)
#   7. Sammon mapping       - unsupervised, nonlinear MDS variant that
#                            weights small (nearby) distances much more
#                            heavily than large ones
#   8. UMAP (supervised)    - nonlinear neighbor-graph based, but blends in
#                            Gene_Type via `y` + `target_weight` so
#                            class label informs the embedding
#
# Aesthetics used in every plot (kept deliberately simple/legible):
#   1. Shape  -> Gene_Type (Functional = solid circle, Pseudogene = open diamond)
#   2. Color  -> functional_probability (continuous, blue = low, red = high)
#   Fixed alpha = 0.8 on all points (not mapped to data).
#
# enrichment_intergenic_dbsnp, PhyloP100_median, ENCODE_max, and the
# promoter-type (type1/2/3) score + motif-count columns are the input
# features used for every embedding below. They are not separately encoded
# as point aesthetics on top of that -- shape and color (below) come from
# Gene_Type and functional_probability instead.
#
# `label` and `functional_probability` are aesthetic-only (functional_
# probability drives point color; label isn't currently used) and are
# excluded as input features. Gene_Type is used as the class label for the
# two supervised methods (LDA, supervised UMAP); UMAP is also run in its
# unsupervised form for comparison.
#
# Requires: MASS (ships with base R), Rtsne, uwot, lle
#   install.packages(c("Rtsne", "uwot", "lle"))
#
# Usage:
#   Rscript mds_pseudogene_plot.R
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
})

for (pkg in c("MASS", "Rtsne", "uwot", "lle")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    print('Install lle from # install.packages("remotes") remotes::install_github("cran/lle")')
    stop(sprintf("This script requires the '%s' package. Install it with: install.packages(\"%s\")", pkg, pkg))
  }
}

# ---------------------------------------------------------------------------
# 0. Arguments
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

# ---------------------------------------------------------------------------
# 2. Feature matrix shared by every method
# ---------------------------------------------------------------------------
feature_cols <- c(
  "PhyloP100_median", "ENCODE_max", "enrichment_intergenic_dbsnp",
  "type1_score", "type1_motif_counts",
  "type2_score", "type2_motif_counts",
  "type3_score", "type3_motif_counts"
)
# Note: `label` and `functional_probability` are aesthetic-only (see header
# comment) and are excluded as input features -- functional_probability
# drives point color, label isn't currently used as an aesthetic at all.

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
  select(Gene, Gene_Type, Gene_group,
         enrichment_intergenic_dbsnp, PhyloP100_median, ENCODE_max,
         functional_probability)

n_genes <- nrow(feature_scaled)

# ---------------------------------------------------------------------------
# 3. Shared plotting function
# ---------------------------------------------------------------------------
make_plot <- function(embed_df, x_col, y_col, title, subtitle, out_path,
                      shape_col = "Gene_Type",
                      shape_values = c("Functional" = 16, "Pseudogene" = 5),
                      shape_name = "Gene Type") {
  # Sort so higher functional_probability points are drawn last (on top),
  # since those are the ones we want to highlight when points overlap.
  embed_df <- embed_df %>% arrange(functional_probability)
  
  
  p <- ggplot(embed_df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(
      aes(shape = .data[[shape_col]], color = functional_probability),
      size = 3, stroke = 1.2, alpha = 0.8
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
    guides(
      shape = guide_legend(override.aes = list(color = "black", alpha = 1))
    ) +
    labs(title = title, subtitle = subtitle, x = x_col, y = y_col) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.key = element_rect(fill = "white"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(out_path, plot = p, width = 9, height = 7, dpi = 300)
  message(sprintf("Saved plot to: %s", out_path))
  
  coords_out <- sub("\\.[^.]+$", "_coordinates.csv", out_path)
  write_csv(embed_df, coords_out)
  message(sprintf("Saved coordinates to: %s", coords_out))
  
  invisible(p)
}

# ---------------------------------------------------------------------------
# 4. Method 1: PCA (unsupervised, linear, maximizes overall variance)
# ---------------------------------------------------------------------------
pca_fit <- prcomp(feature_scaled, center = FALSE, scale. = FALSE)
pca_var_pct <- round(100 * (pca_fit$sdev^2 / sum(pca_fit$sdev^2))[1:2], 1)

pca_df <- bind_cols(metadata, as.data.frame(pca_fit$x[, 1:2])) %>%
  rename(PC1 = PC1, PC2 = PC2)

make_plot(
  pca_df, "PC1", "PC2",
  title = "PCA of Gene Conservation / Variation Features",
  subtitle = sprintf("n = %d genes | PC1: %.1f%% var | PC2: %.1f%% var", n_genes, pca_var_pct[1], pca_var_pct[2]),
  out_path = file.path(output_dir, "pca_pseudogene_plot.png")
)

# ---------------------------------------------------------------------------
# 5. Method 2: LDA (supervised, linear, maximizes between-class separation)
#    Grouping variable is now a 4-way cross of Gene_Type x whether
#    functional_probability agrees with that label (threshold 0.5), instead
#    of the binary Gene_Type used elsewhere in this script:
#      True Functional  = Functional gene, functional_probability >= 0.5
#      False Functional = Functional gene, functional_probability <  0.5
#      False Pseudogene = Pseudogene,      functional_probability >= 0.5
#      True Pseudogene  = Pseudogene,      functional_probability <  0.5
#    (functional_probability == 0.5 exactly is routed to the ">= 0.5" /
#    "true"-labeled side -- flag if you want the opposite convention.)
#    4 classes gives up to K-1 = 3 discriminant axes; we keep LD1 and LD2
#    and plot them as a normal 2D scatter (via make_plot) instead of the
#    1D strip plot the binary version needed. Shape now marks the 4-way
#    group (circle = Functional Gene_Type, diamond = Pseudogene Gene_Type,
#    filled = probability agrees with the label, open = it doesn't) and
#    color still shows the raw functional_probability gradient, since shape
#    only encodes which side of 0.5 a point falls on, not how far.
# ---------------------------------------------------------------------------
lda4_mask <- !is.na(df_complete$functional_probability)
if (any(!lda4_mask)) {
  message(sprintf(
    "4-group LDA: dropping %d row(s) with missing functional_probability.",
    sum(!lda4_mask)
  ))
}

feature_scaled_lda4 <- feature_scaled[lda4_mask, , drop = FALSE]
metadata_lda4 <- metadata[lda4_mask, ]

match_group <- case_when(
  metadata_lda4$Gene_Type == "Pseudogene" & metadata_lda4$functional_probability <  0.5 ~ "True Pseudogene",
  metadata_lda4$Gene_Type == "Pseudogene" & metadata_lda4$functional_probability >= 0.5 ~ "False Pseudogene",
  metadata_lda4$Gene_Type == "Functional" & metadata_lda4$functional_probability <  0.5 ~ "False Functional",
  metadata_lda4$Gene_Type == "Functional" & metadata_lda4$functional_probability >= 0.5 ~ "True Functional"
)
match_group <- factor(
  match_group,
  levels = c("True Functional", "False Functional", "False Pseudogene", "True Pseudogene")
)

group_counts <- table(match_group)
if (any(group_counts < 2)) {
  warning("4-group LDA: at least one group has fewer than 2 members ",
          "(counts: ", paste(names(group_counts), group_counts, sep = "=", collapse = ", "),
          "); LDA may be unstable or fail.")
}

lda_fit <- MASS::lda(x = feature_scaled_lda4, grouping = match_group)
lda_scores <- predict(lda_fit)$x  # n x up to 3 matrix (LD1..LD3, for 4 classes)

# Proportion of between-class variance ("trace") captured by each axis,
# the LDA analog of PCA's % variance explained.
lda_prop_trace <- round(100 * lda_fit$svd^2 / sum(lda_fit$svd^2), 1)

lda_df <- bind_cols(
  metadata_lda4,
  tibble(
    LD1 = as.numeric(lda_scores[, 1]),
    LD2 = as.numeric(lda_scores[, 2]),
    Match_Group = match_group
  )
)

make_plot(
  lda_df, "LD1", "LD2",
  title = "LDA of Gene Conservation / Variation Features",
  subtitle = sprintf(
    "n = %d genes | 4-group (Gene_Type x functional_probability >=/< 0.5) | LD1: %.1f%% | LD2: %.1f%% trace",
    nrow(lda_df), lda_prop_trace[1], lda_prop_trace[2]
  ),
  out_path = file.path(output_dir, "lda_pseudogene_plot.png"),
  shape_col = "Match_Group",
  shape_values = c(
    "True Functional"  = 16,  # filled circle: Functional gene, prob agrees
    "False Functional" = 1,   # open circle:   Functional gene, prob disagrees
    "False Pseudogene" = 5,   # open diamond:  Pseudogene, prob disagrees
    "True Pseudogene"  = 18   # filled diamond: Pseudogene, prob agrees
  ),
  shape_name = "Classification"
)

# ---------------------------------------------------------------------------
# 6. Method 3: t-SNE (unsupervised, nonlinear, local neighborhood structure)
# ---------------------------------------------------------------------------
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
  subtitle = sprintf("n = %d genes | perplexity = %d | unsupervised", n_genes, perplexity_use),
  out_path = file.path(output_dir, "tsne_pseudogene_plot.png")
)

# ---------------------------------------------------------------------------
# 7. Method 4: UMAP (unsupervised, nonlinear neighbor-graph based)
# ---------------------------------------------------------------------------
n_neighbors_use <- min(15, n_genes - 1)
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
  subtitle = sprintf("n = %d genes | n_neighbors = %d | unsupervised", n_genes, n_neighbors_use),
  out_path = file.path(output_dir, "umap_pseudogene_plot.png")
)

# ---------------------------------------------------------------------------
# 8. Method 5: Classic MDS (unsupervised, linear; ~ PCA on Euclidean distances)
# ---------------------------------------------------------------------------
dist_matrix <- dist(feature_scaled, method = "euclidean")
mds_fit <- cmdscale(dist_matrix, k = 2, eig = TRUE)
mds_var_pct <- round(100 * mds_fit$eig[1:2] / sum(abs(mds_fit$eig)), 1)

mds_df <- bind_cols(metadata, as.data.frame(mds_fit$points))
colnames(mds_df)[(ncol(mds_df) - 1):ncol(mds_df)] <- c("MDS1", "MDS2")

make_plot(
  mds_df, "MDS1", "MDS2",
  title = "Classic MDS of Gene Conservation / Variation Features",
  subtitle = sprintf("n = %d genes | Axis 1: %.1f%% var | Axis 2: %.1f%% var", n_genes, mds_var_pct[1], mds_var_pct[2]),
  out_path = file.path(output_dir, "mds_pseudogene_plot.png")
)

# ---------------------------------------------------------------------------
# 9. Method 6: LLE (unsupervised, nonlinear manifold method)
#    Locally Linear Embedding reconstructs each point as a weighted linear
#    combination of its k nearest neighbors in the original 9-D feature
#    space, then finds a low-D embedding that best preserves those local
#    reconstruction weights. Unlike PCA/MDS it only cares about local
#    neighborhoods, not global distances -- so a dense/homogeneous cluster
#    (e.g. pseudogenes, if their features are more uniform) can pull
#    together tightly while a more heterogeneous group stays spread out.
#    `k` (neighbors used for reconstruction) is chosen the same way as
#    n_neighbors_use for UMAP, capped at n_genes - 1.
# ---------------------------------------------------------------------------
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
  subtitle = sprintf("n = %d genes | k = %d neighbors | unsupervised", n_genes, lle_k_use),
  out_path = file.path(output_dir, "lle_pseudogene_plot.png")
)

# ---------------------------------------------------------------------------
# 10. Method 7: Sammon mapping (unsupervised, nonlinear MDS variant)
#    Like classic MDS, Sammon mapping tries to preserve pairwise distances,
#    but its stress function divides each squared error by the *original*
#    distance -- so small (nearby) distances are weighted far more heavily
#    than large ones. A tight, homogeneous cluster gets its internal
#    structure preserved with high fidelity; distances among more scattered
#    points are allowed more error. sammon() requires strictly positive
#    distances between distinct points, so exact duplicate feature rows (if
#    any) get a tiny jitter first purely to break ties -- it does not
#    materially change the embedding.
# ---------------------------------------------------------------------------
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

sammon_init <- cmdscale(dist_matrix_sammon, k = 2)  # classic MDS start, standard practice
set.seed(42)
sammon_fit <- MASS::sammon(dist_matrix_sammon, y = sammon_init, k = 2, trace = FALSE)

sammon_df <- bind_cols(metadata, as.data.frame(sammon_fit$points))
colnames(sammon_df)[(ncol(sammon_df) - 1):ncol(sammon_df)] <- c("Sammon1", "Sammon2")

make_plot(
  sammon_df, "Sammon1", "Sammon2",
  title = "Sammon Mapping of Gene Conservation / Variation Features",
  subtitle = sprintf("n = %d genes | final stress = %.4f | unsupervised", n_genes, sammon_fit$stress),
  out_path = file.path(output_dir, "sammon_pseudogene_plot.png")
)

# ---------------------------------------------------------------------------
# 11. Method 8: UMAP (supervised)
#    Same neighbor-graph machinery as the unsupervised UMAP above, but
#    Gene_Type is passed as `y` so class labels also inform which points get
#    pulled together. `target_weight` controls the balance between the
#    unsupervised (X) and label-driven (y) components: 0 = ignore labels
#    (identical to unsupervised UMAP), 1 = collapse points by label almost
#    like a classifier. 0.5 (the uwot default) gives labels real influence
#    without fully discarding the feature-space geometry -- adjust higher if
#    you want tighter, more label-driven clustering.
# ---------------------------------------------------------------------------
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
  subtitle = sprintf("n = %d genes | n_neighbors = %d | target_weight = 0.5 | supervised on Gene_Type", n_genes, n_neighbors_use),
  out_path = file.path(output_dir, "umap_supervised_pseudogene_plot.png")
)

message("All eight dimension-reduction plots saved to: ", output_dir)