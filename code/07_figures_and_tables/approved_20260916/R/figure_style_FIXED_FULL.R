suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(grid)
})

# Device-portable clean sans-serif family. Cairo maps this to a standard
# Helvetica/Arial-equivalent face without platform-specific font registration.
FIG_FAMILY <- "sans"

FIG_PALETTE <- list(
  neutral_dark = "#2C2C2C",
  neutral_mid = "#777777",
  neutral_light = "#D8DDE0",
  monthly = "#C7D0D6",
  centroid = "#985544",
  comparison = "#667F91",
  iwv = "#527DA3",
  zbar = "#C9854A",
  interaction = "#7F9992",
  mass = "#B56F43",
  height = "#6F91A6",
  residual = "#9A9A9A",
  diverging = c("#355F8D", "#6F8EAF", "#AABACB", "#D8DEE4", "#F7F6F2",
                "#E8D2C8", "#D3A08F", "#B86855", "#8E3B35"),
  separation = c("#F7F8F7", "#EFF3F2", "#E2ECE9", "#D1E2DE", "#BAD2CC", "#96B8B0")
)

# Low-saturation v5 palette. Kept separate so previous figure versions render
# exactly as before when their functions are called.
FIG_PALETTE_V5 <- list(
  neutral_dark = "#303234",
  neutral_mid = "#74797C",
  neutral_light = "#E6E8E7",
  neutral_wash = "#F5F5F2",
  iwv = "#587D9D",
  zbar = "#B77A43",
  interaction = "#81988D",
  positive = "#B66F4F",
  negative = "#6687A2",
  residual = "#94989A"
)

# Robinson's conventional width:height is approximately 1.97:1. The source
# projection is normalized to +/-1 on both axes, so this factor restores the
# natural projected aspect without adapting it to an arbitrary panel box.
ROBINSON_NATURAL_Y_FACTOR <- 1 / 1.97

# Single title constructor used by every v6 panel. The panel tag alone is bold;
# title text remains plain and shares the same baseline.
panel_title_v6 <- function(tag = NULL, title = NULL) {
  panel_title(tag, title)
}

natural_robinson_bundle <- function(data) {
  wide_map_bundle(data, y_factor = ROBINSON_NATURAL_Y_FACTOR)
}

panel_title <- function(tag = NULL, title = NULL) {
  if (is.null(title)) return(NULL)
  if (is.null(tag)) return(title)
  bquote(bold(.(tag)) ~ .(title))
}

theme_panel <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = FIG_FAMILY) +
    theme(
      text = element_text(family = FIG_FAMILY, colour = FIG_PALETTE$neutral_dark),
      plot.title = element_text(size = 9, face = "plain", hjust = 0,
                                margin = margin(b = 5)),
      axis.title = element_text(size = 7.8, face = "plain"),
      axis.text = element_text(size = 7.1, colour = "#444444"),
      axis.line = element_line(linewidth = 0.32, colour = "#555555"),
      axis.ticks = element_line(linewidth = 0.30, colour = "#555555"),
      axis.ticks.length = unit(1.4, "mm"),
      legend.title = element_text(size = 7.1, face = "plain"),
      legend.text = element_text(size = 6.7, face = "plain"),
      legend.key.height = unit(2.6, "mm"),
      legend.key.width = unit(9, "mm"),
      legend.spacing.x = unit(1.3, "mm"),
      plot.caption = element_text(size = 6.4, face = "plain", hjust = 0,
                                  colour = "#555555", margin = margin(t = 4)),
      plot.margin = margin(5, 7, 5, 6),
      panel.grid = element_blank()
    )
}

theme_map_panel <- function() {
  theme_void(base_family = FIG_FAMILY, base_size = 7.5) +
    theme(
      text = element_text(family = FIG_FAMILY, colour = FIG_PALETTE$neutral_dark),
      plot.title = element_text(size = 9, face = "plain", hjust = 0,
                                margin = margin(b = 4)),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 7.0, face = "plain", margin = margin(r = 3)),
      legend.text = element_text(size = 6.6, face = "plain"),
      legend.key.width = unit(11, "mm"),
      legend.key.height = unit(2.5, "mm"),
      legend.margin = margin(t = 1, b = 0),
      legend.box.spacing = unit(0.5, "mm"),
      plot.margin = margin(5, 7, 4, 6),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

# Robinson projection normalized to x in [-1, 1] at the equator and y in [-1, 1].
robinson_project <- function(lon, lat) {
  rob_x <- c(1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600,
             0.9427, 0.9216, 0.8962, 0.8679, 0.8350, 0.7986, 0.7597,
             0.7186, 0.6732, 0.6213, 0.5722, 0.5322)
  rob_y <- c(0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720,
             0.4340, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903,
             0.8435, 0.8936, 0.9394, 0.9761, 1.0000)
  alat <- pmin(abs(lat), 90)
  k <- pmin(floor(alat / 5) + 1L, 19L)
  k2 <- pmin(k + 1L, 19L)
  frac <- ifelse(k == 19L, 0, (alat - (k - 1L) * 5) / 5)
  xx <- rob_x[k] + frac * (rob_x[k2] - rob_x[k])
  yy <- rob_y[k] + frac * (rob_y[k2] - rob_y[k])
  data.frame(x = (lon / 180) * xx, y = sign(lat) * yy)
}

make_robinson_grid_geometry <- function(lon, lat) {
  stopifnot(all(diff(lon) > 0), all(diff(lat) > 0))
  lon_edge <- c(-180, (head(lon, -1) + tail(lon, -1)) / 2, 180)
  lat_edge <- c(-90, (head(lat, -1) + tail(lat, -1)) / 2, 90)
  cell <- expand.grid(i = seq_along(lon), j = seq_along(lat))
  x1 <- lon_edge[cell$i]; x2 <- lon_edge[cell$i + 1L]
  y1 <- lat_edge[cell$j]; y2 <- lat_edge[cell$j + 1L]
  plon <- as.vector(rbind(x1, x2, x2, x1, x1))
  plat <- as.vector(rbind(y1, y1, y2, y2, y1))
  pr <- robinson_project(plon, plat)
  data.frame(x = pr$x, y = pr$y,
             id = rep(seq_len(nrow(cell)), each = 5L),
             cell = rep(seq_len(nrow(cell)), each = 5L))
}

attach_map_field <- function(geometry, field, field_name = "value") {
  stopifnot(length(field) * 5L == nrow(geometry))
  out <- geometry
  out[[field_name]] <- rep(as.vector(field), each = 5L)
  out
}

load_robinson_coastline <- function(csv_file) {
  d <- read.csv(csv_file, stringsAsFactors = FALSE)
  n <- nrow(d)
  grp <- rep(NA_integer_, n)
  g <- 0L
  for (i in seq_len(n)) {
    valid <- is.finite(d$lon[i]) && is.finite(d$lat[i])
    if (!valid) next
    restart <- i == 1L || !is.finite(d$lon[i - 1L]) || !is.finite(d$lat[i - 1L]) ||
      abs(d$lon[i] - d$lon[i - 1L]) > 180
    if (restart) g <- g + 1L
    grp[i] <- g
  }
  ok <- is.finite(d$lon) & is.finite(d$lat) & !is.na(grp)
  pr <- robinson_project(d$lon[ok], d$lat[ok])
  data.frame(x = pr$x, y = pr$y, group = grp[ok])
}

make_robinson_graticule <- function() {
  parts <- list(); g <- 0L
  for (lo in seq(-150, 150, 30)) {
    g <- g + 1L
    la <- seq(-89.5, 89.5, length.out = 240)
    pr <- robinson_project(rep(lo, length(la)), la)
    parts[[length(parts) + 1L]] <- data.frame(x = pr$x, y = pr$y, group = g)
  }
  for (la in seq(-60, 60, 30)) {
    g <- g + 1L
    lo <- seq(-179.8, 179.8, length.out = 480)
    pr <- robinson_project(lo, rep(la, length(lo)))
    parts[[length(parts) + 1L]] <- data.frame(x = pr$x, y = pr$y, group = g)
  }
  do.call(rbind, parts)
}

make_robinson_border <- function() {
  n <- 500L
  lon <- c(seq(-180, 180, length.out = n), rep(180, n),
           seq(180, -180, length.out = n), rep(-180, n))
  lat <- c(rep(-90, n), seq(-90, 90, length.out = n),
           rep(90, n), seq(90, -90, length.out = n))
  pr <- robinson_project(lon, lat)
  data.frame(x = pr$x, y = pr$y, group = 1L)
}

map_panel <- function(map_df, fill_scale, coast, graticule, border,
                      tag = NULL, title = NULL, non_significant = NULL,
                      annotation = NULL) {
  p <- ggplot() +
    geom_polygon(data = map_df, aes(x, y, group = id, fill = value), colour = NA) +
    geom_path(data = graticule, aes(x, y, group = group),
              linewidth = 0.15, colour = "#B8B8B8", alpha = 0.65)
  if (!is.null(non_significant)) {
    wash <- map_df[rep(as.vector(non_significant), each = 5L), , drop = FALSE]
    p <- p + geom_polygon(data = wash, aes(x, y, group = id), inherit.aes = FALSE,
                          fill = "#F1F1F1", alpha = 0.48, colour = NA)
  }
  p <- p +
    geom_path(data = coast, aes(x, y, group = group),
              linewidth = 0.23, colour = "#4A4A4A") +
    geom_path(data = border, aes(x, y, group = group),
              linewidth = 0.26, colour = "#777777") +
    fill_scale +
    coord_equal(xlim = c(-1.015, 1.015), ylim = c(-1.02, 1.02), expand = FALSE) +
    labs(title = panel_title(tag, title), x = NULL, y = NULL) +
    theme_map_panel()
  if (!is.null(annotation)) {
    p <- p + annotate("text", x = -0.94, y = 0.84, label = annotation,
                      hjust = 0, vjust = 1, family = FIG_FAMILY,
                      size = 2.35, colour = "#3D3D3D")
  }
  p
}

wide_map_bundle <- function(data, y_factor = 0.52) {
  transform_y <- function(x) { x$y <- x$y * y_factor; x }
  list(
    z = transform_y(data$z_map_df),
    wgpe = transform_y(data$wgpe_map_df),
    coast = transform_y(data$coast),
    graticule = transform_y(data$graticule),
    border = transform_y(data$border)
  )
}

map_panel_wide <- function(map_df, fill_scale, coast, graticule, border,
                           tag = NULL, title = NULL, non_significant = NULL) {
  p <- ggplot() +
    geom_polygon(data = map_df, aes(x, y, group = id, fill = value), colour = NA) +
    geom_path(data = graticule, aes(x, y, group = group),
              linewidth = 0.14, colour = "#B8B8B8", alpha = 0.60)
  if (!is.null(non_significant)) {
    wash <- map_df[rep(as.vector(non_significant), each = 5L), , drop = FALSE]
    p <- p + geom_polygon(data = wash, aes(x, y, group = id), inherit.aes = FALSE,
                          fill = "#F2F2F2", alpha = 0.46, colour = NA)
  }
  p +
    geom_path(data = coast, aes(x, y, group = group),
              linewidth = 0.22, colour = "#4A4A4A") +
    geom_path(data = border, aes(x, y, group = group),
              linewidth = 0.25, colour = "#777777") +
    fill_scale +
    coord_equal(xlim = c(-1.015, 1.015), ylim = c(-0.535, 0.535), expand = FALSE) +
    labs(title = panel_title(tag, title), x = NULL, y = NULL) +
    theme_map_panel() +
    theme(plot.margin = margin(4, 7, 3, 6))
}

# Build an even-odd polygon that covers the plotting rectangle but leaves the
# Robinson globe as a transparent hole. Drawing this after all geographic
# layers provides strict globe-boundary clipping without changing the data.
make_robinson_outer_mask <- function(border, xlim = c(-1.02, 1.02),
                                     ylim = c(-0.54, 0.54)) {
  outer <- data.frame(
    x = c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]),
    y = c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1]),
    group = 1L, subgroup = 1L
  )
  inner <- border
  inner$group <- 1L
  inner$subgroup <- 2L
  rbind(outer[, c("x", "y", "group", "subgroup")],
        inner[, c("x", "y", "group", "subgroup")])
}

map_panel_wide_clipped <- function(map_df, fill_scale, coast, graticule, border,
                                   tag = NULL, title = NULL,
                                   non_significant = NULL) {
  xlim <- c(-1.015, 1.015)
  ypad <- 0.012
  ylim <- c(-max(abs(border$y)) - ypad, max(abs(border$y)) + ypad)
  mask <- make_robinson_outer_mask(border, xlim, ylim)
  p <- ggplot() +
    geom_polygon(data = map_df, aes(x, y, group = id, fill = value), colour = NA) +
    geom_path(data = graticule, aes(x, y, group = group),
              linewidth = 0.13, colour = "#B9BDBD", alpha = 0.55)
  if (!is.null(non_significant)) {
    wash <- map_df[rep(as.vector(non_significant), each = 5L), , drop = FALSE]
    p <- p + geom_polygon(data = wash, aes(x, y, group = id), inherit.aes = FALSE,
                          fill = "#F1F2F1", alpha = 0.47, colour = NA)
  }
  p +
    geom_path(data = coast, aes(x, y, group = group),
              linewidth = 0.20, colour = "#505456", lineend = "round") +
    geom_polygon(data = mask,
                 aes(x, y, group = group, subgroup = subgroup),
                 inherit.aes = FALSE, fill = "white", colour = NA,
                 rule = "evenodd") +
    geom_path(data = border, aes(x, y, group = group),
              linewidth = 0.30, colour = "#6A6E70", linejoin = "round") +
    fill_scale +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE, clip = "on") +
    labs(title = panel_title(tag, title), x = NULL, y = NULL) +
    theme_map_panel() +
    theme(
      plot.title = element_text(size = 8.8, face = "plain", margin = margin(b = 3)),
      plot.margin = margin(4, 6, 2, 6),
      legend.box.spacing = unit(0.2, "mm"),
      legend.margin = margin(t = 0, b = 0)
    )
}

theme_mechanism_panel <- function(base_size = 7.5) {
  theme_panel(base_size) +
    theme(
      panel.grid.major.y = element_line(linewidth = 0.20, colour = "#E7E7E7"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 7, 5, 7)
    )
}
