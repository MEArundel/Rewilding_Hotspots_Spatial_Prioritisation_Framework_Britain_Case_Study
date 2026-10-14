############################################################
# Rewilding hotspot metric maps
#
# Three-panel figure showing:
#
#   A. Species richness
#   B. Faith's phylogenetic diversity
#   C. Mean viability headroom score
#
# Species richness has a discrete legend.
# Faith's PD and mean headroom have continuous legends.
############################################################


# ==========================================================
# 1. Packages
# ==========================================================

required_packages <- c(
  "terra",
  "dplyr"
)

packages_to_install <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

library(terra)
library(dplyr)


# ==========================================================
# 2. Paths
# ==========================================================

hotspot_dir <-
  "./results/Rewilding Hotspots/Barrier Polygon Summary/"

headroom_dir <- file.path(
  hotspot_dir,
  "Viability Headroom"
)

hotspot_geojson <- file.path(
  hotspot_dir,
  "Rewilding_Hotspots_FaithPD.geojson"
)

headroom_file <- file.path(
  headroom_dir,
  "Hotspot_Viability_Headroom_Summary.csv"
)

boundary_file <-
  "./data/Britain Boundary Polygon.geojson"

out_file <- file.path(
  headroom_dir,
  "Rewilding_Hotspots_SpeciesRichness_PD_Headroom_Multipanel.jpg"
)


# ==========================================================
# 3. Check files
# ==========================================================

if (!file.exists(hotspot_geojson)) {
  stop(
    "Could not find hotspot polygon file:\n",
    hotspot_geojson
  )
}

if (!file.exists(headroom_file)) {
  stop(
    "Could not find headroom summary file:\n",
    headroom_file,
    "\n\nRun the viability headroom analysis first."
  )
}

if (!file.exists(boundary_file)) {
  stop(
    "Could not find Britain boundary file:\n",
    boundary_file
  )
}


# ==========================================================
# 4. Read hotspot polygons
# ==========================================================

hotspot_polys <- terra::vect(
  hotspot_geojson
)


# ==========================================================
# 5. Read Britain boundary
# ==========================================================

boundary_vect <- terra::vect(
  boundary_file
)

boundary_vect <- terra::project(
  boundary_vect,
  terra::crs(hotspot_polys)
)

boundary_union <- terra::aggregate(
  boundary_vect
)


# ==========================================================
# 6. Read hotspot headroom scores
# ==========================================================

headroom_summary <- read.csv(
  headroom_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_headroom_fields <- c(
  "hotspot_id",
  "headroom_score_mean"
)

missing_headroom_fields <- setdiff(
  required_headroom_fields,
  names(headroom_summary)
)

if (length(missing_headroom_fields) > 0) {
  stop(
    paste0(
      "Headroom table is missing required field(s):\n",
      paste(
        missing_headroom_fields,
        collapse = "\n"
      )
    )
  )
}


# ==========================================================
# 7. Check hotspot polygon fields
# ==========================================================

required_hotspot_fields <- c(
  "barrier_poly_id",
  "species_richness",
  "median_faith_PD"
)

missing_hotspot_fields <- setdiff(
  required_hotspot_fields,
  names(hotspot_polys)
)

if (length(missing_hotspot_fields) > 0) {
  stop(
    paste0(
      "Hotspot polygons are missing required field(s):\n",
      paste(
        missing_hotspot_fields,
        collapse = "\n"
      )
    )
  )
}


# ==========================================================
# 8. Attach mean headroom score to hotspot polygons
# ==========================================================

hotspot_df <- as.data.frame(
  hotspot_polys
)

headroom_match <- match(
  hotspot_df$barrier_poly_id,
  headroom_summary$hotspot_id
)

hotspot_polys$headroom_score <-
  headroom_summary$headroom_score_mean[
    headroom_match
  ]

if (
  any(
    !is.finite(
      hotspot_polys$headroom_score
    )
  )
) {
  missing_ids <-
    hotspot_polys$barrier_poly_id[
      !is.finite(
        hotspot_polys$headroom_score
      )
    ]
  
  stop(
    paste0(
      "The following hotspots could not be matched to a ",
      "mean headroom score:\n",
      paste(
        missing_ids,
        collapse = ", "
      )
    )
  )
}


# ==========================================================
# 9. Continuous colour helper
# ==========================================================

make_metric_colours <- function(
    values,
    palette_colours,
    n_colours = 101
) {
  value_min <- min(
    values,
    na.rm = TRUE
  )
  
  value_max <- max(
    values,
    na.rm = TRUE
  )
  
  palette <- colorRampPalette(
    palette_colours
  )(
    n_colours
  )
  
  if (value_max > value_min) {
    scaled_values <- (
      values - value_min
    ) / (
      value_max - value_min
    )
  } else {
    scaled_values <- rep(
      1,
      length(values)
    )
  }
  
  colour_index <- round(
    scaled_values * (n_colours - 1)
  ) + 1
  
  colour_index <- pmax(
    1,
    pmin(n_colours, colour_index)
  )
  
  list(
    colours = palette[colour_index],
    palette = palette,
    minimum = value_min,
    maximum = value_max
  )
}


# ==========================================================
# 10. Species richness colours
# ==========================================================

richness_values <- sort(
  unique(
    hotspot_polys$species_richness
  )
)

richness_palette <- colorRampPalette(
  c(
    "#ffffcc",
    "#c2e699",
    "#78c679",
    "#31a354",
    "#006837"
  )
)(
  length(richness_values)
)

richness_colour_lookup <- setNames(
  richness_palette,
  richness_values
)

richness_polygon_colours <-
  richness_colour_lookup[
    as.character(
      hotspot_polys$species_richness
    )
  ]


# ==========================================================
# 11. Faith's PD colours
# ==========================================================

pd_colours <- make_metric_colours(
  hotspot_polys$median_faith_PD,
  c(
    "#fee8c8",
    "#fdbb84",
    "#e34a33",
    "#7f0000"
  )
)


# ==========================================================
# 12. Mean viability headroom colours
# ==========================================================

headroom_colours <- make_metric_colours(
  hotspot_polys$headroom_score,
  c(
    "#edf8fb",
    "#b3cde3",
    "#8c96c6",
    "#8856a7",
    "#810f7c"
  )
)


# ==========================================================
# 13. Helper to plot one map
# ==========================================================

plot_hotspot_panel <- function(
    hotspot_polys,
    boundary_union,
    polygon_colours,
    panel_title
) {
  map_extent <- terra::ext(
    boundary_union
  )
  
  dx <- map_extent$xmax - map_extent$xmin
  dy <- map_extent$ymax - map_extent$ymin
  
  # Small padding lets the maps use more of each panel.
  xpad <- 0.01 * dx
  ypad <- 0.01 * dy
  
  plot(
    boundary_union,
    col = NA,
    border = "black",
    lwd = 1.6,
    xlim = c(
      map_extent$xmin - xpad,
      map_extent$xmax + xpad
    ),
    ylim = c(
      map_extent$ymin - ypad,
      map_extent$ymax + ypad
    ),
    axes = FALSE,
    frame.plot = FALSE
  )
  
  plot(
    hotspot_polys,
    col = polygon_colours,
    border = "grey20",
    lwd = 0.9,
    add = TRUE
  )
  
  title(
    main = panel_title,
    cex.main = 1.85,
    font.main = 2,
    line = 0.5
  )
}


# ==========================================================
# 14. Helper for continuous legends
# ==========================================================

plot_continuous_legend <- function(
    metric_info,
    legend_title
) {
  legend_values <- seq(
    metric_info$minimum,
    metric_info$maximum,
    length.out = length(
      metric_info$palette
    )
  )
  
  image(
    x = 1,
    y = legend_values,
    z = matrix(
      legend_values,
      nrow = 1
    ),
    col = metric_info$palette,
    zlim = range(
      legend_values
    ),
    axes = FALSE,
    xlab = "",
    ylab = "",
    useRaster = TRUE
  )
  
  legend_breaks <- pretty(
    c(
      metric_info$minimum,
      metric_info$maximum
    ),
    n = 5
  )
  
  legend_breaks <- legend_breaks[
    legend_breaks >= metric_info$minimum &
      legend_breaks <= metric_info$maximum
  ]
  
  axis(
    side = 4,
    at = legend_breaks,
    labels = round(
      legend_breaks,
      1
    ),
    las = 1,
    cex.axis = 0.95,
    lwd = 0,
    lwd.ticks = 1
  )
  
  usr <- par("usr")
  
  segments(
    x0 = usr[2],
    y0 = metric_info$minimum,
    x1 = usr[2],
    y1 = metric_info$maximum,
    lwd = 1,
    col = "black",
    xpd = FALSE
  )
  
  mtext(
    legend_title,
    side = 4,
    line = 3,
    cex = 1.05
  )
  
  mtext(
    "Higher",
    side = 3,
    line = 1,
    cex = 0.85
  )
  
  mtext(
    "Lower",
    side = 1,
    line = 2,
    cex = 0.85
  )
}


# ==========================================================
# 15. Create three-panel figure
# ==========================================================

jpeg(
  filename = out_file,
  width = 7200,
  height = 3300,
  res = 300,
  quality = 100
)

layout(
  matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 1,
    byrow = TRUE
  ),
  widths = c(
    8, 1.3,
    8, 1.3,
    8, 1.3
  )
)

# Space for the overall title and the rightmost legend title.
par(
  oma = c(0, 0, 5, 6)
)


# ----------------------------------------------------------
# Panel A — Species richness
# ----------------------------------------------------------

par(
  mar = c(1, 0.5, 3.5, 0.5)
)

plot_hotspot_panel(
  hotspot_polys = hotspot_polys,
  boundary_union = boundary_union,
  polygon_colours = richness_polygon_colours,
  panel_title = "A. Species Richness"
)


# ----------------------------------------------------------
# Legend A — categorical species richness
# ----------------------------------------------------------

par(
  mar = c(4, 0, 3.5, 5)
)

plot.new()

legend(
  "center",
  legend = rev(richness_values),
  fill = rev(richness_palette),
  border = "grey20",
  bty = "n",
  cex = 1.55,
  y.intersp = 1.4
)

mtext(
  "Species Richness",
  side = 4,
  line = 3,
  cex = 1.4
)


# ----------------------------------------------------------
# Panel B — Faith's phylogenetic diversity
# ----------------------------------------------------------

par(
  mar = c(1, 0.5, 3.5, 0.5)
)

plot_hotspot_panel(
  hotspot_polys = hotspot_polys,
  boundary_union = boundary_union,
  polygon_colours = pd_colours$colours,
  panel_title = "B. Faith's Phylogenetic Diversity"
)


# ----------------------------------------------------------
# Legend B — Faith's PD
# ----------------------------------------------------------

par(
  mar = c(4, 0, 3.5, 5)
)

plot_continuous_legend(
  metric_info = pd_colours,
  legend_title = "Median Faith's PD"
)


# ----------------------------------------------------------
# Panel C — mean viability headroom
# ----------------------------------------------------------

par(
  mar = c(1, 0.5, 3.5, 0.5)
)

plot_hotspot_panel(
  hotspot_polys = hotspot_polys,
  boundary_union = boundary_union,
  polygon_colours = headroom_colours$colours,
  panel_title = "C. Viability Headroom"
)


# ----------------------------------------------------------
# Legend C — mean viability headroom
# ----------------------------------------------------------

par(
  mar = c(4, 0, 3.5, 5)
)

plot_continuous_legend(
  metric_info = headroom_colours,
  legend_title = "Mean Headroom Score"
)


# ----------------------------------------------------------
# Overall figure title
# ----------------------------------------------------------

mtext(
  "Candidate Rewilding Hotspots",
  side = 3,
  outer = TRUE,
  line = 1.5,
  cex = 2.8,
  font = 2
)

dev.off()


# ==========================================================
# 16. Completion message
# ==========================================================

message(
  "\nThree-panel hotspot metric figure complete.\n\n",
  "Saved to:\n",
  normalizePath(
    out_file,
    mustWork = FALSE
  )
)