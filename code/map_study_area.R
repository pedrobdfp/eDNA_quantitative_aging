# =============================================================================
# 06_map.R
# -----------------------------------------------------------------------------
# Study-area map: satellite view of the Hood Canal sampling site with an
# inset showing regional context (~50 km across).
#
# Marks the two sampling stations (Husbandry Area, NOAA Boat) and the
# inferred main eDNA source (captive dolphins), with annotated distances.
#
# Output: plots/map_study_area.tiff
#
# This script is self-contained; it does not depend on 00_setup.R's fits.
# =============================================================================

library(here)
library(sf)
library(terra)
library(maptiles)
library(ggplot2)
library(ggspatial)
library(grid)
library(cowplot)

dir.create(here("plots"), showWarnings = FALSE, recursive = TRUE)

# --- Study area polygon (lon, lat) -------------------------------------------

coords <- matrix(c(
  -122.75, 47.737635,
  -122.72, 47.734855,
  -122.75, 47.746697,
  -122.72, 47.743580,
  -122.75, 47.737635   # close polygon
), ncol = 2, byrow = TRUE)

area_sf <- st_sf(id = 1,
                 geometry = st_sfc(st_polygon(list(coords)), crs = 4326))

# --- Sampling points ---------------------------------------------------------

pts <- data.frame(
  id  = c("Husbandry (nearby)", "Boat (further)", "Main eDNA source"),
  lon = c(-122.729975, -122.743206, -122.729799),
  lat = c(  47.742236,   47.736441,   47.742805)
)
pts_sf <- st_as_sf(pts, coords = c("lon", "lat"), crs = 4326)

# --- Project everything to UTM zone 10N (metric CRS) ------------------------

crs_proj  <- 32610
area_proj <- st_transform(area_sf, crs_proj)
pts_proj  <- st_transform(pts_sf,  crs_proj)

# --- Satellite tiles for the main panel -------------------------------------

sat <- get_tiles(
  x        = area_proj,
  provider = "Esri.WorldImagery",
  crop     = TRUE,
  zoom     = 16,
  project  = TRUE,
  verbose  = TRUE
)

# --- Distances and connecting lines (for annotation) ------------------------

d12 <- as.numeric(st_distance(pts_proj[1, ], pts_proj[2, ]))  # Husb - Boat
d13 <- as.numeric(st_distance(pts_proj[1, ], pts_proj[3, ]))  # Husb - Source

line12 <- st_union(pts_proj[1, ], pts_proj[2, ]) |> st_cast("LINESTRING")
line13 <- st_union(pts_proj[1, ], pts_proj[3, ]) |> st_cast("LINESTRING")

midpoint <- function(line) st_line_sample(line, sample = 0.5) |> st_cast("POINT")
mid12 <- midpoint(line12)
mid13 <- midpoint(line13)

# --- Layout math for scale bar and north arrow -------------------------------

bb <- st_bbox(area_proj)
dx <- as.numeric(bb["xmax"] - bb["xmin"])
dy <- as.numeric(bb["ymax"] - bb["ymin"])

# Narrow the main map horizontally around its centre (80% of original width)
xc        <- (bb["xmin"] + bb["xmax"]) / 2
width_new <- dx * 0.8
xmin_new  <- as.numeric(xc - width_new / 2)
xmax_new  <- as.numeric(xc + width_new / 2)

# Scale bar (two-tone, bottom-right)
bar_len    <- 500
bar_height <- 0.015 * dy
x_right    <- as.numeric(xmax_new - 0.08 * width_new)
x_left     <- x_right - bar_len
x_mid      <- (x_left + x_right) / 2
y_bar      <- as.numeric(bb["ymin"] + 0.06 * dy)

# North arrow (triangle above the scale bar)
x_arrow_center <- x_mid
y_arrow_base   <- as.numeric(y_bar + 0.10 * dy)
y_arrow_tip    <- as.numeric(y_bar + 0.14 * dy)
na_half_width  <- 0.012 * width_new

north_arrow_head <- data.frame(
  x = c(x_arrow_center, x_arrow_center - na_half_width, x_arrow_center + na_half_width),
  y = c(y_arrow_tip,    y_arrow_base,                   y_arrow_base)
)
y_arrow_label <- as.numeric(y_arrow_tip + 0.025 * dy)

# Horizontal offset for distance-text labels (keep them off the yellow line)
offset_x <- 0.05 * width_new

# --- Inset panel: zoomed-out regional view -----------------------------------

x_center <- as.numeric(xc)
y_center <- as.numeric((bb["ymin"] + bb["ymax"]) / 2)
half_w   <- 25000   # half-width  in metres (-> ~50 km span)
half_h   <- 25000   # half-height in metres

inset_bbox_sf <- st_as_sfc(st_bbox(
  c(xmin = x_center - half_w, xmax = x_center + half_w,
    ymin = y_center - half_h, ymax = y_center + half_h),
  crs = crs_proj
))

sat_inset <- get_tiles(
  x        = inset_bbox_sf,
  provider = "Esri.WorldImagery",
  crop     = TRUE,
  zoom     = 11,
  project  = TRUE,
  verbose  = TRUE
)

main_bbox_sf <- st_as_sfc(st_bbox(area_proj))

# --- Main panel --------------------------------------------------------------

main_map <- ggplot() +
  layer_spatial(sat) +
  
  # Sampling points
  layer_spatial(pts_proj, aes(color = id), size = 5) +
  
  # Distance lines
  layer_spatial(line12, color = "yellow", linewidth = 1.5) +
  layer_spatial(line13, color = "yellow", linewidth = 1.5) +
  
  # Distance labels
  annotate("text",
           x = st_coordinates(mid12)[1] + offset_x, y = st_coordinates(mid12)[2],
           label = paste0(round(d12), " m"), color = "yellow", size = 6) +
  annotate("text",
           x = st_coordinates(mid13)[1] + offset_x, y = st_coordinates(mid13)[2],
           label = paste0(round(d13), " m"), color = "yellow", size = 6) +
  
  # Scale bar (two-tone)
  geom_rect(aes(xmin = x_left, xmax = x_mid,
                ymin = y_bar,  ymax = y_bar + bar_height),
            fill = "white", color = "black", linewidth = 0.5) +
  geom_rect(aes(xmin = x_mid, xmax = x_right,
                ymin = y_bar,  ymax = y_bar + bar_height),
            fill = "black", color = "black", linewidth = 0.5) +
  annotate("label",
           x = x_mid, y = y_bar + 2.5 * bar_height,
           label = paste0(bar_len, " m"),
           size = 5, fill = "white", alpha = 0.7,
           color = "black", fontface = "bold") +
  
  # North arrow
  geom_polygon(data = north_arrow_head, aes(x = x, y = y),
               fill = "white", color = "black", linewidth = 0.5) +
  annotate("label",
           x = x_arrow_center, y = y_arrow_label,
           label = "N", size = 5, fontface = "bold",
           fill = "white", alpha = 0.7, color = "black") +
  
  scale_color_manual(
    values = c("Husbandry (nearby)" = "#53CFDA",
               "Boat (further)"     = "#FF7994",
               "Main eDNA source"   = "#785eF1")
  ) +
  theme_minimal() +
  theme(
    panel.grid        = element_blank(),
    legend.position   = "top",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 14),
    legend.key.size   = unit(0.6, "cm"),
    axis.title        = element_blank(),
    axis.text         = element_text(size = 14, face = "bold",
                                     margin = margin(t = 2, r = 2)),
    axis.ticks.length = unit(0.15, "cm"),
    plot.margin       = margin(5, 10, 5, 10)
  ) +
  coord_sf(xlim = c(xmin_new, xmax_new),
           ylim = c(bb["ymin"], bb["ymax"]),
           expand = FALSE)

# --- Inset panel (small map with a box showing the main extent) -------------

inset_map <- ggplot() +
  layer_spatial(sat_inset) +
  layer_spatial(main_bbox_sf, fill = NA, color = "yellow", linewidth = 0.9) +
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin  = margin(0, 0, 0, 0)
  )

# --- Compose and save --------------------------------------------------------

combined_plot <- ggdraw(main_map) +
  draw_plot(inset_map, x = 0.05, y = 0.50, width = 0.4, height = 0.4)
combined_plot
ggsave(
  filename    = here("plots", "map_study_area.tiff"),
  plot        = combined_plot,
  width       = 12, height = 9, units = "in",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)

message("Saved: plots/map_study_area.tiff")
