#!/usr/bin/env Rscript
# ============================================================================
# Interactive (plotly) replica of the LDA 4-group plot
#
# Reuses the coordinates already written out by mds_pseudogene_plot.R
# (lda_pseudogene_plot_coordinates.csv) so this script doesn't need to
# recompute the LDA fit -- it just needs to have been run at least once
# first. Reproduces the same visual encoding as the static ggplot version
# as closely as plotly allows:
#
#   - Color   -> functional_probability, blue (low) -> red (high), 0-1 scale
#   - Symbol  -> Match_Group:
#                 True Functional  = filled circle
#                 False Functional = filled circle + thin BLACK outline
#                 False Pseudogene = open diamond
#                 True Pseudogene  = filled diamond
#   - Alpha   -> 0.6 on every point
#   - Draw order -> low -> high functional_probability, with
#                 False Functional / False Pseudogene drawn last (on top)
#
# Hover tooltip shows Gene, Gene_group, Gene_Type, Match_Group,
# functional_probability, and the LD1/LD2 coordinates -- use this to
# identify which points are which when you hover over the interactive plot.
#
# Requires: plotly, htmlwidgets, dplyr, readr
#   install.packages(c("plotly", "htmlwidgets"))
#
# Usage:
#   Rscript lda_pseudogene_plot_interactive.R
# (run mds_pseudogene_plot.R first if lda_pseudogene_plot_coordinates.csv
#  doesn't exist yet)
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(plotly)
  library(htmlwidgets)
})

# ---------------------------------------------------------------------------
# 0. Load the coordinates saved by mds_pseudogene_plot.R
# ---------------------------------------------------------------------------
coords_csv <- "../../results/plots/lda_pseudogene_plot_coordinates.csv"
out_html   <- "../../results/plots/lda_pseudogene_plot_interactive.html"

if (!file.exists(coords_csv)) {
  stop(sprintf(
    "Coordinates file not found: %s\nRun mds_pseudogene_plot.R first -- it writes this file as a side effect of the LDA plot.",
    coords_csv
  ))
}

lda_df <- read_csv(coords_csv, show_col_types = FALSE) %>%
  mutate(Match_Group = factor(
    Match_Group,
    levels = c("True Functional", "False Functional", "False Pseudogene", "True Pseudogene")
  ))

# ---------------------------------------------------------------------------
# 1. Draw order: low -> high functional_probability, with False Functional /
#    False Pseudogene always drawn last (on top), matching the static plot.
# ---------------------------------------------------------------------------
lda_df <- lda_df %>%
  mutate(.draw_last = Match_Group %in% c("False Functional", "False Pseudogene")) %>%
  arrange(.draw_last, functional_probability) %>%
  select(-.draw_last)

# ---------------------------------------------------------------------------
# 2. Per-point plotly symbol + color + outline (mirrors shape_values / the
#    black outline layer in the ggplot version).
#
#    IMPORTANT: colors are computed as literal hex strings in R (via
#    colorRamp), not left to plotly's `colorscale`. Passing `symbol = ~...`
#    to plot_ly() makes it silently split the data into one trace per
#    symbol level, and each split trace can end up normalizing a
#    `colorscale` mapping against its own local min/max instead of the
#    fixed 0-1 range -- that's what caused high-probability points to show
#    up blue. Computing literal per-point colors up front sidesteps that
#    entirely, and everything is drawn as a single trace so nothing gets
#    re-normalized.
# ---------------------------------------------------------------------------
symbol_map <- c(
  "True Functional"  = "circle",
  "False Functional" = "circle",
  "False Pseudogene" = "diamond-open",
  "True Pseudogene"  = "diamond"
)

prob_to_hex <- function(p, low = "#1f77b4", high = "#d62728") {
  p <- pmin(pmax(p, 0), 1)  # clamp to the same [0, 1] limits as the ggplot version
  ramp <- grDevices::colorRamp(c(low, high), space = "Lab")
  rgb_mat <- ramp(p)
  grDevices::rgb(rgb_mat[, 1], rgb_mat[, 2], rgb_mat[, 3], maxColorValue = 255)
}

lda_df <- lda_df %>%
  mutate(
    point_symbol = symbol_map[as.character(Match_Group)],
    point_color  = prob_to_hex(functional_probability),
    # Filled shapes (circle, diamond) show color as fill, no border needed.
    # The open diamond (False Pseudogene) has no fill at all, so it needs
    # its *own* border colored to match, or it's invisible. False Functional
    # gets a fixed black border on top of its normal colored fill.
    outline_color = case_when(
      Match_Group == "False Functional" ~ "black",
      Match_Group == "False Pseudogene" ~ point_color,
      TRUE ~ "rgba(0,0,0,0)"
    ),
    outline_width = case_when(
      Match_Group == "False Functional" ~ 2,
      Match_Group == "False Pseudogene" ~ 1.5,
      TRUE ~ 0
    ),
    hover_text = paste0(
      "Gene: ", Gene, "<br>",
      "Gene group: ", Gene_group, "<br>",
      "Gene type: ", Gene_Type, "<br>",
      "Classification: ", Match_Group, "<br>",
      "Functional probability: ", round(functional_probability, 3), "<br>",
      "LD1: ", round(LD1, 3), "<br>",
      "LD2: ", round(LD2, 3)
    )
  )

# ---------------------------------------------------------------------------
# 3. Build the interactive plot: one data trace (literal per-point colors/
#    symbols, so nothing gets split or re-normalized), plus a couple of
#    invisible helper traces so the color scale and the shape legend still
#    show up correctly.
# ---------------------------------------------------------------------------
p <- plot_ly() %>%
  add_markers(
    data = lda_df,
    x = ~LD1,
    y = ~LD2,
    text = ~hover_text,
    hoverinfo = "text",
    marker = list(
      size = 10,
      opacity = 0.6,
      symbol = ~point_symbol,
      color = ~point_color,
      line = list(
        color = ~outline_color,
        width = ~outline_width
      )
    ),
    showlegend = FALSE
  ) %>%
  # Invisible trace purely to draw the continuous colorbar at the true 0-1 scale.
  add_markers(
    data = lda_df,
    x = ~LD1, y = ~LD2,
    marker = list(
      size = 0.001,
      color = ~functional_probability,
      colorscale = list(c(0, "#1f77b4"), c(1, "#d62728")),
      cmin = 0, cmax = 1,
      showscale = TRUE,
      colorbar = list(title = "Functional\nProbability")
    ),
    hoverinfo = "skip",
    showlegend = FALSE
  ) %>%
  # Invisible dummy traces purely to build a shape/classification legend.
  add_trace(
    data = tibble(grp = names(symbol_map), sym = unname(symbol_map)),
    x = rep(NA, length(symbol_map)), y = rep(NA, length(symbol_map)),
    type = "scatter", mode = "markers",
    marker = list(symbol = ~sym, size = 10, color = "grey40",
                  line = list(color = "black", width = 1)),
    name = ~grp,
    split = ~grp,
    hoverinfo = "skip",
    showlegend = TRUE
  ) %>%
  layout(
    title = "LDA of Gene Conservation / Variation Features (interactive)",
    xaxis = list(title = "LD1"),
    yaxis = list(title = "LD2"),
    legend = list(title = list(text = "Classification"))
  )

saveWidget(p, out_html, selfcontained = TRUE)
message(sprintf("Saved interactive plot to: %s", out_html))
message("Open it in a browser and hover over points to see Gene / annotation info.")