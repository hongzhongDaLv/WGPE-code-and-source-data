options(stringsAsFactors = FALSE, warn = 2)
.libPaths(.libPaths())

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(grid)
})

ROOT <- Sys.getenv('WGPE_PROJECT_ROOT', unset='.')
args <- commandArgs(trailingOnly = TRUE)
REF <- "absolute"
PALETTE_SET <- "reference2"
if (!REF %in% c("absolute", "surface_relative")) stop("REF must be absolute or surface_relative")
if (!PALETTE_SET %in% c("reference1", "reference2")) stop("Unknown PALETTE_SET")
V3 <- file.path(ROOT, "output_attribution_minimal", "FULL_RERUN_FINAL_absolute_relative_v3")
V4 <- file.path(ROOT, "output_attribution_minimal", "FULL_RERUN_TARGETED_REPAIR_v4")
V5 <- file.path(ROOT, "output_attribution_minimal", "FULL_RERUN_TARGETED_REPAIR_v5_VISUAL_QA")
ENSO_BASE <- file.path(ROOT, "WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard_RC2", "11_ENSO_ONI_STANDARD_UPDATE")
CLOSURE <- Sys.getenv("TOP100_CLOSURE_ROOT", unset=file.path(ROOT,"output_TOP100_revision","FINAL_FIGURE_LAYOUT_v5_panel_frame_fix"))
SOURCE_ROOT <- Sys.getenv("TOP100_SOURCE_ROOT", unset=file.path(ROOT,"output_TOP100_revision","FINAL_LOCKED_ARCHIVE_v2"))
COORD_ROOT <- file.path(ROOT, "output_TOP100_revision", "FINAL_COORDINATE_FIXED_CLOSURE")
ADAPTER_ROOT <- Sys.getenv("TOP100_ADAPTER_ROOT", unset=file.path(ROOT,"output_TOP100_revision","derived_data","figure_source_data"))
OUT <- file.path(CLOSURE, "figures", "main", "Figure2_TOP100_FINAL_v5")
NO_PANEL_LABELS <- identical(Sys.getenv("TOP100_NO_PANEL_LABELS", unset = "0"), "1")
FIG2A_ONLY <- identical(Sys.getenv("TOP100_FIG2A_ONLY", unset = "0"), "1")
if (FIG2A_ONLY) {
  OUT <- file.path(CLOSURE, "figures", "main",
                   "Figure2a_deseasonalized_only_TOP100")
}
if (NO_PANEL_LABELS) {
  OUT <- file.path(CLOSURE, "figures", "main_no_panel_labels", "Figure2_TOP100_FINAL_v5")
}
PANEL_DIR <- file.path(OUT, "Figure2_source_panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "R", "figure_style_FIXED_FULL.R"), encoding = "UTF-8")

FONT <- "sans"
pal_blue <- if (PALETTE_SET == "reference2") {
  c("#7F0D3A", "#C83D4D", "#E47B62", "#F6C5B3", "#F7F7F7",
    "#B9D5E5", "#68A8CF", "#2E73B8", "#174A7E")
} else {
  c("#794607", "#A66B1C", "#C89545", "#DFC685", "#F5EFE4",
    "#CCDEEF", "#9CC8DF", "#5EA6D1", "#0D56A0")
}
trial <- if (PALETTE_SET == "reference1") {
  list(blue = "#2E86C1", z = "#D2A23A", inter = "#7B1E7A")
} else {
  list(blue = "#2C7BB6", z = "#D66A55", inter = "#8E2A62")
}
COL <- list(
  blue = trial$blue, brown = trial$z, inter = trial$inter,
  actual = "#596B7A",
  none = "#F5F5F3", nonsig = "#ECECEA", ink = "#202020",
  axis = "#555555", grid = "#ECECEA"
)
PANEL_TITLE_SIZE <- 9.6
PANEL_SUBTITLE_SIZE <- PANEL_TITLE_SIZE
PANEL_TEXT_SIZE <- PANEL_TITLE_SIZE

inflate_all_text <- function(g, delta_pt = 2) {
  if (!is.null(g$gp) && !is.null(g$gp$fontsize)) g$gp$fontsize <- g$gp$fontsize + delta_pt
  if (!is.null(g$children) && length(g$children)) for (i in seq_along(g$children)) g$children[[i]] <- inflate_all_text(g$children[[i]], delta_pt)
  if (!is.null(g$grobs) && length(g$grobs)) for (i in seq_along(g$grobs)) g$grobs[[i]] <- inflate_all_text(g$grobs[[i]], delta_pt)
  g
}

save_png <- function(plot, filename, width_mm, height_mm, dpi = 500, dir = OUT) {
  plot <- ggdraw() + draw_plot(plot, x = .018, y = .018, width = .964, height = .964)
  gr <- grid::grid.force(cowplot::as_grob(plot))
  ggsave(file.path(dir, filename), gr, width = width_mm, height = height_mm,
         units = "mm", dpi = dpi, bg = "white", device = "png")
}

save_pdf <- function(plot, filename, width_mm, height_mm, dir = OUT) {
  plot <- ggdraw() + draw_plot(plot, x = .018, y = .018, width = .964, height = .964)
  gr <- grid::grid.force(cowplot::as_grob(plot))
  ggsave(file.path(dir, filename), gr, width = width_mm, height = height_mm,
         units = "mm", bg = "white", device = cairo_pdf)
}

ptitle <- function(tag, title) {
  if (NO_PANEL_LABELS) title else bquote(bold(.(tag)) ~ .(title))
}

theme_rec <- function(base = 6.2) {
  base <- base + 2
  theme_classic(base_size = base, base_family = FONT) +
    theme(
      text = element_text(family = FONT, colour = COL$ink),
      plot.title = element_text(size = PANEL_TITLE_SIZE, hjust = 0, margin = margin(b = 2), lineheight = .94),
      plot.subtitle = element_text(size = PANEL_SUBTITLE_SIZE, hjust = 0, margin = margin(b = 1.5), colour = COL$ink),
      axis.title = element_text(size = PANEL_TEXT_SIZE, colour = COL$ink, lineheight = .85),
      axis.text = element_text(size = PANEL_TEXT_SIZE, colour = COL$ink),
      axis.line = element_line(linewidth = .25, colour = COL$axis),
      axis.ticks = element_line(linewidth = .22, colour = COL$axis),
      axis.ticks.length = unit(.7, "mm"),
      panel.grid.major = element_line(colour = COL$grid, linewidth = .15),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = PANEL_TEXT_SIZE),
      legend.text = element_text(size = PANEL_TEXT_SIZE),
      legend.key.height = unit(1.9, "mm"),
      legend.key.width = unit(5.2, "mm"),
      plot.title.position = "plot",
      plot.margin = margin(2, 3, 2, 1)
    )
}

theme_map_rec <- function(base = 6.8) {
  base <- base + 2
  theme_void(base_family = FONT, base_size = base) +
    theme(
      text = element_text(family = FONT, colour = COL$ink),
      plot.title = element_text(size = PANEL_TITLE_SIZE, hjust = 0, margin = margin(b = 1), lineheight = .94),
      plot.subtitle = element_text(size = PANEL_SUBTITLE_SIZE, hjust = .5,
                                   margin = margin(t = 0, b = 1), colour = COL$ink),
      plot.margin = margin(1, 1, 0, 1),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

grid_spec <- function() list(lon = 0:359, lat = -90:90)

bundle_from_field <- function(field, q = NULL, y_factor = 0.52) {
  gs <- grid_spec()
  lon <- gs$lon; lat <- gs$lat
  stopifnot(nrow(field) == length(lon), ncol(field) == length(lat))
  if (!is.null(q) && !identical(dim(field), dim(q))) stop("field/q mismatch")
  lon2 <- ifelse(lon >= 180, lon - 360, lon)
  oi <- order(lon2); oj <- order(lat)
  fld <- field[oi, oj, drop = FALSE]
  q_c <- if (!is.null(q)) q[oi, oj, drop = FALSE] else NULL
  geom <- make_robinson_grid_geometry(lon2[oi], lat[oj])
  tf <- function(x) { x$y <- x$y * y_factor; x }
  list(
    lon = lon2[oi], lat = lat[oj], field = fld, q = q_c,
    map = tf(attach_map_field(geom, fld)),
    coast = tf(load_robinson_coastline(file.path(ROOT, "R", "coastlines_source_FIXED_FULL.csv"))),
    graticule = tf(make_robinson_graticule()),
    border = tf(make_robinson_border())
  )
}

bundle_from_lonlat_field <- function(field, lon, lat, y_factor = 0.52) {
  lon2 <- ifelse(lon >= 180, lon - 360, lon)
  oi <- order(lon2); oj <- order(lat)
  fld <- field[oi, oj, drop = FALSE]
  geom <- make_robinson_grid_geometry(lon2[oi], lat[oj])
  tf <- function(x) { x$y <- x$y * y_factor; x }
  list(
    lon = lon2[oi], lat = lat[oj], field = fld,
    map = tf(attach_map_field(geom, fld)),
    coast = tf(load_robinson_coastline(file.path(ROOT, "R", "coastlines_source_FIXED_FULL.csv"))),
    graticule = tf(make_robinson_graticule()),
    border = tf(make_robinson_border())
  )
}

outer_mask <- function(border, xlim, ylim) {
  outer <- data.frame(x = c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]),
                      y = c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1]),
                      group = 1L, subgroup = 1L)
  inner <- border; inner$group <- 1L; inner$subgroup <- 2L
  rbind(outer[, c("x", "y", "group", "subgroup")],
        inner[, c("x", "y", "group", "subgroup")])
}

lat_profile <- function(field, lon, lat, bin_width = 5) {
  br <- seq(-90, 90, by = bin_width)
  out <- lapply(seq_len(length(br) - 1L), function(i) {
    lo <- br[i]; hi <- br[i + 1L]
    idx <- which(if (i == length(br) - 1L) lat >= lo & lat <= hi else lat >= lo & lat < hi)
    vals <- field[, idx, drop = FALSE]
    w <- matrix(rep(cos(lat[idx] * pi / 180), each = length(lon)), nrow = length(lon), ncol = length(idx))
    data.frame(lat_mid = (lo + hi) / 2,
               value = weighted.mean(as.vector(vals), as.vector(w), na.rm = TRUE),
               n_grid = sum(is.finite(vals)))
  })
  do.call(rbind, out)
}

map_plot_cont <- function(bundle, lim, title_expr, fdr = FALSE, subtitle_text = NULL) {
  xlim <- c(-1.015, 1.015)
  ylim <- c(-max(abs(bundle$border$y)) - .012, max(abs(bundle$border$y)) + .012)
  mask <- outer_mask(bundle$border, xlim, ylim)
  mapdf <- bundle$map
  p <- ggplot()
  if (fdr) {
    valid <- is.finite(bundle$field)
    sig <- valid & is.finite(bundle$q) & bundle$q <= 0.05
    nonsig <- valid & !sig
    grey <- mapdf[rep(as.vector(nonsig), each = 5L), , drop = FALSE]
    sigdf <- mapdf
    sigdf$value[!rep(as.vector(sig), each = 5L)] <- NA_real_
    p <- p +
      geom_polygon(data = grey, aes(x, y, group = id), fill = COL$nonsig, colour = NA) +
      geom_polygon(data = sigdf, aes(x, y, group = id, fill = value), colour = NA)
  } else {
    p <- p + geom_polygon(data = mapdf, aes(x, y, group = id, fill = value), colour = NA)
  }
  p +
    geom_path(data = bundle$graticule, aes(x, y, group = group), linewidth = .12, colour = "#B9BDBD", alpha = .55) +
    geom_path(data = bundle$coast, aes(x, y, group = group), linewidth = .19, colour = "#505456", lineend = "round") +
    geom_polygon(data = mask, aes(x, y, group = group, subgroup = subgroup),
                 inherit.aes = FALSE, fill = "white", colour = NA, rule = "evenodd") +
    geom_path(data = bundle$border, aes(x, y, group = group), linewidth = .28, colour = "#6A6E70", linejoin = "round") +
    scale_fill_gradientn(colours = pal_blue, limits = c(-lim, lim), oob = squish, na.value = "transparent", guide = "none") +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE, clip = "on") +
    labs(title = title_expr, subtitle = subtitle_text, x = NULL, y = NULL) +
    theme_map_rec(6.6)
}

colorbar <- function(lim, title_expr, ticks, accuracy = .1, extend = FALSE) {
  d <- data.frame(x = seq(-lim, lim, length.out = 401), y = 0)
  td <- data.frame(x = ticks, label = label_number(accuracy = accuracy, big.mark = " ")(ticks))
  p <- ggplot(d, aes(x, y, fill = x)) +
    geom_tile(width = (2 * lim) / 400, height = .075) +
    annotate("rect", xmin = -lim, xmax = lim, ymin = -.037, ymax = .037,
             fill = NA, colour = "#111111", linewidth = .20) +
    geom_segment(data = td, aes(x = x, xend = x, y = -.040, yend = -.068),
                 inherit.aes = FALSE, linewidth = .18, colour = "#111111") +
    geom_text(data = td, aes(x = x, y = -.108, label = label),
              inherit.aes = FALSE, family = FONT, size = 3.35, colour = "#111111") +
    scale_fill_gradientn(colours = pal_blue, limits = c(-lim, lim), oob = squish, guide = "none") +
    coord_cartesian(xlim = c(-lim*ifelse(extend,1.055,1), lim*ifelse(extend,1.055,1)), ylim = c(-.13, .10), clip = "off") +
    labs(title = title_expr, x = NULL, y = NULL) +
    theme_void(base_family = FONT, base_size = 6) +
    theme(plot.title = element_text(size = PANEL_TEXT_SIZE, hjust = .5, margin = margin(b = .6), colour = "#111111"),
          plot.margin = margin(-3, 4, 0, 4))
  if (extend) p <- p +
    annotate("polygon", x=c(-lim*1.05,-lim,-lim), y=c(0,-.037,.037), fill=pal_blue[1], colour="#111111", linewidth=.20) +
    annotate("polygon", x=c(lim*1.05,lim,lim), y=c(0,-.037,.037), fill=pal_blue[length(pal_blue)], colour="#111111", linewidth=.20)
  p
}

side_profile_cont <- function(bundle, lim, xlabel, line_col = COL$blue, accuracy = .1, profile_lim = lim) {
  prof <- lat_profile(bundle$field, bundle$lon, bundle$lat, 5)
  SIDE_PROFILE_CALL <<- SIDE_PROFILE_CALL + 1L
  vn <- c("precip_r", "ENSO_WGPE")[SIDE_PROFILE_CALL]
  uu <- U_PROFILE[U_PROFILE$height_reference == REF & U_PROFILE$variable_name == vn, ]
  if (vn == "ENSO_WGPE") {
    prof$value <- prof$value / 1000
    uu$lower95 <- uu$lower95 / 1000
    uu$upper95 <- uu$upper95 / 1000
    profile_lim <- profile_lim / 1000
    xlabel <- expression(atop("W-GPE anomaly", "(10"^3~"J m"^{-2}*")"))
    accuracy <- 1
  }
  lim_vals <- c(prof$value, uu$lower95, uu$upper95)
  lim_vals <- lim_vals[is.finite(lim_vals)]
  if (length(lim_vals)) {
    if (min(lim_vals) >= 0) {
      xlim_use <- c(0, max(pretty(c(0, max(lim_vals) * 1.08), n = 3)))
      xbreaks <- c(0, xlim_use[2])
    } else if (max(lim_vals) <= 0) {
      xlim_use <- c(min(pretty(c(min(lim_vals) * 1.08, 0), n = 3)), 0)
      xbreaks <- c(xlim_use[1], 0)
    } else {
      profile_lim <- max(abs(pretty(lim_vals, n = 3)))
      xlim_use <- c(-profile_lim, profile_lim)
      xbreaks <- c(-profile_lim, profile_lim)
    }
  } else {
    xlim_use <- c(-profile_lim, profile_lim)
    xbreaks <- c(-profile_lim, profile_lim)
  }
  side_lab <- function(x) {
    out <- label_number(accuracy = accuracy, big.mark = " ")(x)
    out <- sub("^-0\\.", "-.", out)
    out <- sub("^0\\.", ".", out)
    out
  }
  ggplot(prof, aes(value, lat_mid)) +
    geom_ribbon(data = uu, aes(y = latitude_bin_mid, xmin = lower95, xmax = upper95),
                inherit.aes = FALSE, fill = alpha(line_col, .16), colour = NA) +
    geom_vline(xintercept = 0, colour = "#6D6D6D", linewidth = .19) +
    geom_path(colour = line_col, linewidth = .34, lineend = "round") +
    geom_point(colour = line_col, size = .30, stroke = 0) +
    scale_y_continuous(limits = c(-90, 90), breaks = c(-60, 0, 60),
                       expand = expansion(mult = c(.01, .01))) +
    scale_x_continuous(limits = xlim_use, breaks = xbreaks,
                       labels = side_lab,
                       expand = expansion(mult = c(0, 0), add = c(0, 0))) +
    coord_cartesian(xlim = xlim_use, clip = "on", expand = FALSE) +
    labs(x = xlabel, y = expression("Latitude ("*degree*")")) +
    theme_rec(5.0) +
    theme(panel.grid.major.y = element_line(colour = "#EEEEEA", linewidth = .14),
          panel.grid.major.x = element_blank(),
          axis.title.x = element_text(margin = margin(t = .45), lineheight = .78),
          axis.title.y = element_text(margin = margin(r = .4)),
          axis.text.x = element_text(margin = margin(t = .25)),
          axis.text.y = element_text(margin = margin(r = 0)),
          axis.ticks.length = unit(.48, "mm"),
          plot.margin = margin(.6, 10, .6, 0))
}

map_block_cont <- function(field, q, lim, title_expr, cbar_title_expr, cbar_ticks,
                           fdr = FALSE, accuracy = .1, profile_accuracy = accuracy,
                           profile_colour = COL$blue, profile_lim = lim,
                           subtitle_text = NULL, colorbar_extend = FALSE) {
  b <- bundle_from_field(field, q)
  mp <- map_plot_cont(b, lim, title_expr, fdr, subtitle_text)
  side <- plot_grid(NULL, side_profile_cont(b, lim, cbar_title_expr, profile_colour, profile_accuracy, profile_lim),
                    NULL, ncol = 1, rel_heights = c(.12, .76, .12))
  top <- plot_grid(mp, side, nrow = 1, rel_widths = c(.80, .20), align = "h", axis = "tb")
  cb <- plot_grid(NULL, colorbar(lim, cbar_title_expr, cbar_ticks, accuracy, colorbar_extend),
                  NULL, nrow = 1, rel_widths = c(.03, .74, .23))
  # Preserve the accepted map and colourbar grobs at exactly their existing
  # sizes. Only translate the complete colourbar grob upward slightly.
  ggdraw() +
    draw_plot(top, x = 0, y = .205, width = .965, height = .795) +
    draw_plot(cb,  x = 0, y = .035, width = .965, height = .150)
}

temporal_cols <- c("forward-only" = COL$blue, "reverse-only" = COL$brown,
                   "bidirectional" = COL$inter, "none" = COL$none)
TEMPORAL_LEGEND <- NULL

temporal_e_block <- function(fdr_grid) {
  td <- fdr_grid[fdr_grid$target == "precipitation" & fdr_grid$predictor == "W-GPE" &
                   fdr_grid$lag == 3 & fdr_grid$valid_TRUE_FALSE, , drop = FALSE]
  classes <- c("forward-only", "reverse-only", "bidirectional", "none")
  lon <- sort(unique(td$lon)); lat <- sort(unique(td$lat))
  code <- setNames(seq_along(classes), classes)
  mat <- matrix(NA_real_, nrow = length(lon), ncol = length(lat))
  mat[cbind(match(td$lon, lon), match(td$lat, lat))] <- code[as.character(td$class_q)]
  b <- bundle_from_lonlat_field(mat, lon, lat)
  xlim <- c(-1.015, 1.015); ylim <- c(-max(abs(b$border$y)) - .012, max(abs(b$border$y)) + .012)
  mask <- outer_mask(b$border, xlim, ylim)
  mapdf <- b$map
  mapdf$class <- factor(classes[round(mapdf$value)], levels = classes)
  mapdf <- mapdf[!is.na(mapdf$class), , drop = FALSE]
  mp <- ggplot() +
    geom_polygon(data = mapdf, aes(x, y, group = id, fill = class), colour = NA) +
    geom_path(data = b$graticule, aes(x, y, group = group), linewidth = .12, colour = "#B9BDBD", alpha = .55) +
    geom_path(data = b$coast, aes(x, y, group = group), linewidth = .19, colour = "#505456", lineend = "round") +
    geom_polygon(data = mask, aes(x, y, group = group, subgroup = subgroup),
                 inherit.aes = FALSE, fill = "white", colour = NA, rule = "evenodd") +
    geom_path(data = b$border, aes(x, y, group = group), linewidth = .28, colour = "#6A6E70", linejoin = "round") +
    scale_fill_manual(values = temporal_cols, name = NULL,
                      labels = c("forward-only", "reverse-only", "bidirectional", "none")) +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE, clip = "on") +
    labs(title = ptitle("e", "Lag-3 precipitation temporal precedence"),
         x = NULL, y = NULL) +
    theme_map_rec(6.4) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.text = element_text(size = 3.25),
          legend.key.width = unit(2.1, "mm"),
          legend.key.height = unit(1.45, "mm"),
          legend.spacing.x = unit(1.0, "mm"),
          legend.margin = margin(t = 1, b = 0),
          legend.box.margin = margin(t = 0, b = 0),
          plot.subtitle = element_text(size = 5.0, margin = margin(b = 1.0), colour = COL$ink))

  TEMPORAL_LEGEND <<- get_legend(mp + guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom", legend.spacing.x = unit(2.5, "mm"),
          legend.box.margin = margin(1, 4, 1, 4), legend.key.width = unit(5, "mm")))
  mp <- mp + theme(legend.position = "none")

  # Latitude profile: net forward-only area minus reverse-only area, per 5-degree bin.
  br <- seq(-90, 90, by = 5)
  prof <- lapply(seq_len(length(br) - 1), function(i) {
    lo <- br[i]; hi <- br[i + 1]
    idx <- td$lat >= lo & if (i == length(br) - 1) td$lat <= hi else td$lat < hi
    sub <- td[idx, , drop = FALSE]
    if (nrow(sub) == 0) return(data.frame(lat_mid = (lo + hi) / 2, net = NA_real_))
    w <- cos(sub$lat * pi / 180)
    f <- sum(w[sub$class_q == "forward-only"], na.rm = TRUE)
    r <- sum(w[sub$class_q == "reverse-only"], na.rm = TRUE)
    tot <- sum(w, na.rm = TRUE)
    data.frame(lat_mid = (lo + hi) / 2, net = 100 * (f - r) / tot)
  })
  prof <- do.call(rbind, prof)
  final_prof_file <- file.path(ROOT,"output_TOP100_revision","FINAL_SCIENTIFIC_ARCHIVE","figure_source_data","Figure2e_NET_ASYMMETRY_FINAL_SOURCE.csv")
  if (file.exists(final_prof_file)) {
    fp <- read.csv(final_prof_file,check.names=FALSE)
    prof <- data.frame(lat_mid=fp$latitude_bin_mid,net=fp$estimate_net_asymmetry)
    U_TEMP_PROF <- data.frame(height_reference=REF,latitude_bin_mid=fp$latitude_bin_mid,
                              lower95=fp$lower95_net_asymmetry,upper95=fp$upper95_net_asymmetry)
  }
  side_profile <- ggplot(prof, aes(net, lat_mid)) +
    geom_ribbon(data = U_TEMP_PROF[U_TEMP_PROF$height_reference == REF, ],
                aes(y = latitude_bin_mid, xmin = lower95, xmax = upper95),
                inherit.aes = FALSE, fill = alpha(COL$blue, .16), colour = NA) +
    geom_vline(xintercept = 0, colour = "#6D6D6D", linewidth = .19) +
    geom_path(colour = COL$blue, linewidth = .34, lineend = "round", na.rm = TRUE) +
    scale_y_continuous(limits = c(-90, 90), breaks = c(-60, 0, 60), expand = expansion(mult = c(.01, .01))) +
    scale_x_continuous(breaks = c(-50, 0, 50)) +
    coord_cartesian(xlim = c(-50, 50), clip = "on") +
    labs(x = "Net asym. (% area)", y = expression("Latitude ("*degree*")")) +
    theme_rec(5.0) +
    theme(panel.grid.major.y = element_line(colour = "#EEEEEA", linewidth = .14),
          panel.grid.major.x = element_blank(),
          axis.title.x = element_text(size = 3.95, margin = margin(t = .4), lineheight = .78),
          axis.title.y = element_text(size = 4.95, margin = margin(r = .35)),
          axis.text.x = element_text(size = 3.95, margin = margin(t = .2)),
          axis.text.y = element_text(size = 4.6, margin = margin(r = 0)),
          axis.ticks.length = unit(.45, "mm"),
          plot.margin = margin(.5, 2.8, .5, 0))

  # Area-fraction summary below profile.
  w_all <- cos(td$lat * pi / 180)
  sumdf <- data.frame(class = classes)
  sumdf$pct <- vapply(classes, function(cl) {
    100 * sum(w_all[td$class_q == cl], na.rm = TRUE) / sum(w_all, na.rm = TRUE)
  }, numeric(1))
  tu <- U_TEMP[U_TEMP$height_reference == REF & U_TEMP$target == "precipitation" & U_TEMP$lag == 3, ]
  sumdf$lower95 <- tu$lower95[match(sumdf$class, tu$variable)]
  sumdf$upper95 <- tu$upper95[match(sumdf$class, tu$variable)]
  sumdf$class <- factor(sumdf$class, levels = classes)
  small_bar <- ggplot(sumdf, aes(class, pct, fill = class)) +
    geom_col(width = .64, colour = "#555555", linewidth = .08) +
    geom_text(aes(y = 5, label = sprintf("%.1f", pct)), angle = 90,
              hjust = 0, vjust = .5, size = 1.35, colour = "#111111") +
    scale_fill_manual(values = temporal_cols, guide = "none") +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100), expand = expansion(mult = c(0, .03))) +
    labs(x = NULL, y = "% area") +
    theme_rec(4.8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 4.15, margin = margin(r = .2)),
          axis.text.y = element_text(size = 4.15),
          panel.grid.major.x = element_blank(),
          plot.margin = margin(0, 6, 2, 0))

  side <- plot_grid(NULL, side_profile, small_bar, ncol = 1, rel_heights = c(.06, .62, .32))
  plot_grid(mp, side, nrow = 1, rel_widths = c(5, 1), align = "h", axis = "tb")
}

guides_bg <- function(plot) {
  ggdraw(plot) +
    theme(plot.background = element_rect(fill = "white", colour = "#D0D0D0", linewidth = .4))
}

# Inputs: formal TOP100 values in the compatibility schema expected by the
# accepted archived layout.
if (REF != "absolute") stop("Candidate main Figure 2 is absolute-height only")
F2 <- readRDS(file.path(ADAPTER_ROOT, "Figure2_TOP100_adapter.rds"))
precip <- F2$precip
event <- F2$event
fdr_grid <- F2$fdr_grid
U_PROFILE <- F2$latitude_profile_uncertainty
U_EVENT <- F2$enso_uncertainty
U_EVENT_CENTERED <- F2$wet_uncertainty
U_TEMP <- F2$temporal_fraction_uncertainty
U_TEMP_PROF <- F2$temporal_profile_uncertainty
SIDE_PROFILE_CALL <- 0L

# a/b map blocks.
enso_global <- event$summary$anomaly[event$summary$event == "ENSO_difference" &
                                       event$summary$variable == "wgpe" &
                                       event$summary$region == "Global"][1]
pglob <- precip$summary[precip$summary$variable == "wgpe" &
                          precip$summary$lag_months_predictor_leads_precip == 0 &
                          precip$summary$region == "Global", ]
a <- map_block_cont(
  precip$precipitation$wgpe$lag0$r, precip$precipitation$wgpe$lag0$q, .8,
  ptitle("a", "W-GPE-precipitation association"),
  "Pearson r", c(-.8, -.4, 0, .4, .8), TRUE, .1, .1, COL$blue,
  subtitle_text = sprintf("Global mean r = %+.3f", pglob$mean_r_area[1])
)
if (FIG2A_ONLY) {
  save_png(a, "Figure2a_deseasonalized_only_TOP100.png", 92, 62, dpi = 600)
  writeLines(c(
    "# Figure 2a deseasonalized-only plotting QA",
    "",
    "- Main field: calendar-month anomalies only; no linear detrending.",
    "- Longitude alignment: one-to-one 0-359 to -179-180 before correlation.",
    paste0("- Displayed global mean r: ", sprintf("%+.6f", pglob$mean_r_area[1]), "."),
    "- Detrended result is sensitivity only.",
    "- Only the standalone Figure 2a PNG was generated; the full Figure 2 was not redrawn."
  ), file.path(OUT, "Figure2a_deseasonalized_only_TOP100_plot_QA.md"), useBytes = TRUE)
  quit(save = "no", status = 0L)
}
b <- map_block_cont(
  event$event_maps$wgpe$ENSO_difference, NULL, 50000,
  ptitle("b", "ENSO W-GPE anomaly"),
  expression("W-GPE anomaly (J m"^{-2}*")"), c(-50000, -25000, 0, 25000, 50000),
  FALSE, 25000, 90000, COL$blue, 180000,
  subtitle_text = bquote("Global mean = "*.(sprintf("%+.0f", enso_global))*" J m"^{-2}),
  colorbar_extend = TRUE
)

# c wet/dry total + channel decomposition.  Read the eight locked result_ids
# directly; the plotted 10^4 scale is a display transform only.
wet_source_file <- file.path(SOURCE_ROOT,"figure_source_data","Fig2c_source.csv")
ws <- read.csv(wet_source_file,check.names=FALSE)
metric_label <- c(actual="Total",IWV_channel="IWV",z_channel="<z>",interaction_channel="Interaction")
wetdry_long <- data.frame(event=ws$event,event_label=ws$event,
                          component=unname(metric_label[ws$metric]),
                          value=ws$estimate/10000,lower95=ws$lower95/10000,
                          upper95=ws$upper95/10000,result_id=ws$result_id)
wetdry_long$event_label <- factor(wetdry_long$event_label,levels=c("Dry","Wet"))
wetdry_long$component <- factor(wetdry_long$component, levels = c("Total", "IWV", "<z>", "Interaction"))
wetdry_long$ci_matches_estimate <- is.finite(wetdry_long$lower95) & is.finite(wetdry_long$upper95) &
  wetdry_long$value >= wetdry_long$lower95 & wetdry_long$value <= wetdry_long$upper95
wetdry_long$lower95_plot <- ifelse(wetdry_long$ci_matches_estimate, wetdry_long$lower95, NA_real_)
wetdry_long$upper95_plot <- ifelse(wetdry_long$ci_matches_estimate, wetdry_long$upper95, NA_real_)
ev_order <- c("Dry", "Wet")
base_y_w <- setNames(seq_along(ev_order), ev_order)
offset_w <- c("Total" = .30, "IWV" = .10, "<z>" = -.10, "Interaction" = -.30)
height_w <- c("Total" = .145, "IWV" = .145, "<z>" = .145, "Interaction" = .145)
wetdry_long$y <- base_y_w[as.character(wetdry_long$event_label)] + offset_w[as.character(wetdry_long$component)]
wetdry_long$h <- height_w[as.character(wetdry_long$component)]
wetdry_long$xmin <- pmin(0, wetdry_long$value)
wetdry_long$xmax <- pmax(0, wetdry_long$value)
xlim_w <- ceiling(max(abs(wetdry_long$value), na.rm = TRUE) * 1.14)

c <- ggplot(wetdry_long) +
  geom_vline(xintercept = 0, colour = "#606060", linewidth = .28) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y - h / 2, ymax = y + h / 2, fill = component),
            colour = "#222222", linewidth = .08) +
  geom_errorbar(aes(xmin = lower95_plot, xmax = upper95_plot, y = y), orientation = "y",
                width = .07, linewidth = .22, colour = "#111111", na.rm = TRUE) +
  scale_fill_manual(values = c("Total" = COL$actual, "IWV" = COL$blue, "<z>" = COL$brown, "Interaction" = COL$inter),
                    breaks = c("Total", "IWV", "<z>", "Interaction"),
                    labels = c("Total", "IWV", "<z>", "Int."),
                    name = NULL) +
  scale_y_continuous(breaks = base_y_w[ev_order], labels = ev_order,
                     limits = c(.45, length(ev_order) + .55), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-xlim_w, xlim_w),
                     breaks = pretty(c(-xlim_w, xlim_w), n = 5),
                     labels = label_number(accuracy=.1)) +
  labs(title = ptitle("c", "Wet/dry decomposition"),
       x = expression("Anomaly / contribution (10"^4~"J m"^{-2}*")"), y = NULL) +
  theme_rec(5.65) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.justification = "center", legend.box.spacing = unit(.1, "mm"),
        legend.margin = margin(3, 0, 0, 0),
        legend.box.margin = margin(t = 3, b = 0),
        legend.text = element_text(size = PANEL_TEXT_SIZE), legend.key.width = unit(4.6, "mm"),
        legend.key.height = unit(1.9, "mm"),
        plot.title = element_text(size = PANEL_TITLE_SIZE, hjust = 0, margin = margin(b = 1.0)),
        plot.subtitle = element_text(size = PANEL_TEXT_SIZE, hjust = 0, margin = margin(b = .8)),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = PANEL_TEXT_SIZE, colour = COL$ink),
        axis.text.x = element_text(size = PANEL_TEXT_SIZE),
        axis.title.x = element_text(size = PANEL_TEXT_SIZE, margin = margin(t = 3)),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(2.0, 4, 7.5, 4))

# d ENSO standard-climate-band signed stacked channel bars, following the existing candidate panel.
enso_dec <- F2$enso_decomposition
band_order <- c("Global", "N high latitude", "N temperate", "Tropics", "S temperate", "S high latitude")
enso_dec$region <- factor(enso_dec$region, levels = band_order)
enso_dec <- enso_dec[order(enso_dec$region), , drop = FALSE]
enso_long <- rbind(
  data.frame(region = enso_dec$region, component = "IWV", value = enso_dec$IWV_channel_J_m2 / 1000),
  data.frame(region = enso_dec$region, component = "<z>", value = enso_dec$z_channel_J_m2 / 1000),
  data.frame(region = enso_dec$region, component = "Interaction", value = enso_dec$interaction_channel_J_m2 / 1000)
)
enso_long$component <- factor(enso_long$component, levels = c("IWV", "<z>", "Interaction"))
enso_long$region <- factor(as.character(enso_long$region), levels = rev(band_order))
enso_actual_ci <- U_EVENT[U_EVENT$height_reference == REF & U_EVENT$event_type == "ENSO_difference" &
                            U_EVENT$metric == "actual", c("region","estimate","lower95","upper95")]
enso_actual_ci$region <- factor(enso_actual_ci$region, levels = rev(band_order))
xlim_enso <- c(-8, 36)
d <- ggplot(enso_long, aes(x = value, y = region, fill = component)) +
  geom_vline(xintercept = 0, colour = "#606060", linewidth = .25) +
  geom_col(width = .50, colour = "#222222", linewidth = .055) +
  scale_fill_manual(values = c("IWV" = COL$blue, "<z>" = COL$brown, "Interaction" = COL$inter),
                    breaks = c("IWV", "<z>", "Interaction"), name = NULL) +
  scale_x_continuous(breaks = c(-5, 0, 10, 20, 30),
                     labels = label_number(big.mark = " ")) +
  coord_cartesian(xlim = xlim_enso, clip = "on") +
  labs(title = ptitle("d", "ENSO decomposition"),
       x = expression("Channel anomaly (10"^3~"J m"^{-2}*")"), y = NULL) +
  theme_rec(5.45) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.justification = "center", legend.box.spacing = unit(.1, "mm"),
        legend.margin = margin(3, 0, 0, 0),
        legend.box.margin = margin(t = 3, b = 0),
        legend.text = element_text(size = PANEL_TEXT_SIZE),
        legend.key.width = unit(4.6, "mm"),
        legend.key.height = unit(1.9, "mm"),
        plot.title = element_text(size = PANEL_TITLE_SIZE, hjust = 0, margin = margin(b = 1.0)),
        plot.subtitle = element_text(size = PANEL_TEXT_SIZE, hjust = 0, margin = margin(b = .8)),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = PANEL_TEXT_SIZE, colour = COL$ink),
        axis.text.x = element_text(size = PANEL_TEXT_SIZE),
        axis.title.x = element_text(size = PANEL_TEXT_SIZE, margin = margin(t = 3)),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(2.0, 4, 7.5, 4))

# e Incremental information beyond IWV.  This panel reads the corrected
# canonical true-increment source:
# M0=IWV, M1=IWV+<z>, M2=IWV+<z>+IWV:<z>, M3=IWV+W-GPE.
# The old single-predictor M3 (Y ~ W-GPE) is not used in the main figure.
TRUE_INC_ROOT <- file.path(ROOT, "CANONICAL_WGPE_TOP100_v2", "04_RESULTS_CANONICAL", "FIG2E_TRUE_INCREMENT")
true_inc_file <- file.path(TRUE_INC_ROOT, "FIG2E_TRUE_INCREMENT_RESULTS.csv")
if (!file.exists(true_inc_file)) stop("Missing corrected Fig2e true-increment source: ", true_inc_file)
e_model <- read.csv(true_inc_file, check.names = FALSE)
e_model <- e_model[e_model$region == "Global" & e_model$metric == "blocked-CV R2", , drop = FALSE]
e_model$source_file <- true_inc_file
e_model$model <- factor(e_model$model, levels = paste0("M", 0:3))
e_model$target <- factor(e_model$target, levels = c("Precipitation", "Continuous Nino3.4"))
e_model$model_definition <- c("IWV", "IWV + <z>", "IWV + <z> + IWV:<z>", "IWV + W-GPE")[match(as.character(e_model$model), paste0("M", 0:3))]
dir.create(file.path(CLOSURE, "figure_source_data"), recursive = TRUE, showWarnings = FALSE)
write.csv(e_model, file.path(CLOSURE, "figure_source_data", "Fig2e_four_model_source.csv"), row.names = FALSE)

e <- ggplot(e_model, aes(model, estimate, colour = target, group = target)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = .24) +
  geom_line(position = position_dodge(width = .22), linewidth = .42) +
  geom_point(position = position_dodge(width = .22), size = 1.25) +
  scale_colour_manual(values = c("Precipitation" = COL$brown, "Continuous Nino3.4" = COL$blue),
                      name = NULL) +
  scale_x_discrete(labels = c("M0", "M1", "M2", "M3")) +
  scale_y_continuous(labels = label_number(accuracy = .01),
                     breaks = c(0, .1, .2, .3),
                     limits = c(0, .305),
                     expand = expansion(mult = c(0, .03))) +
  labs(title = ptitle("e", "Incremental information beyond IWV"),
       x = NULL, y = expression("Blocked-CV "*R^2)) +
  theme_rec(6.0) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.justification = "center", legend.key.width = unit(6.5, "mm"),
        legend.key.height = unit(2.0, "mm"), legend.spacing.x = unit(4.0, "mm"),
        legend.margin = margin(3, 0, 0, 0), legend.box.margin = margin(t = 3, b = 0),
        axis.text.x = element_text(size = PANEL_TEXT_SIZE),
        axis.text.y = element_text(size = PANEL_TEXT_SIZE, margin = margin(r = 2)),
        axis.title.y = element_text(size = PANEL_TEXT_SIZE, margin = margin(r = 5)),
        plot.margin = margin(2, 6, 8, 9),
        panel.grid.major.x = element_blank())
e_individual <- e

# Main layout: row 1 maps; bottom row keeps c/d/e as peer panels with the
# same row frame height.  This avoids the earlier c+d wrapper making panel e
# read as a larger, separate figure.
top_row <- plot_grid(a, b, nrow = 1, rel_widths = c(1, 1), align = "h", axis = "tb")
bottom_row <- plot_grid(c, d, e, nrow = 1, rel_widths = c(1, 1, 1), align = "h", axis = "tb")
main <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(.50, .50), align = "v", axis = "lr")

save_png(main, "Figure2_TOP100_FINAL_v5.png", 225, 160, dpi = 600)
save_pdf(main, "Figure2_TOP100_FINAL_v5.pdf", 225, 160)
save_png(a, "Figure2_TOP100_panel_a_precip.png", 92, 55, dpi = 600)
save_png(b, "Figure2_TOP100_panel_b_ENSO.png", 92, 55, dpi = 600)

guided <- plot_grid(
  plot_grid(guides_bg(a), guides_bg(b), nrow = 1, rel_widths = c(1, 1), align = "h", axis = "tb"),
  plot_grid(guides_bg(c), guides_bg(d), guides_bg(e), nrow = 1, rel_widths = c(1, 1, 1), align = "h", axis = "tb"),
  ncol = 1, rel_heights = c(.50, .50), align = "v", axis = "lr"
)
save_png(guided, "Figure2_TOP100_v5_with_guides.png", 225, 160, dpi = 600)

save_png(a, "a.png", 92, 55, dir = PANEL_DIR)
save_png(b, "b.png", 92, 55, dir = PANEL_DIR)
save_png(c, "c.png", 78, 50, dir = PANEL_DIR)
save_png(d, "d.png", 78, 50, dir = PANEL_DIR)
save_png(e_individual, "e.png", 78, 50, dir = PANEL_DIR)
save_pdf(a, "a.pdf", 92, 55, dir = PANEL_DIR)
save_pdf(b, "b.pdf", 92, 55, dir = PANEL_DIR)
save_pdf(c, "c.pdf", 78, 50, dir = PANEL_DIR)
save_pdf(d, "d.pdf", 78, 50, dir = PANEL_DIR)
save_pdf(e_individual, "e.pdf", 78, 50, dir = PANEL_DIR)

source_values <- rbind(
  data.frame(panel = "a", group = "Global", component = "area_weighted_global_mean_r",
             value = pglob$mean_r_area[1], unit = "Pearson r",
             source = "output_TOP100_revision/results/precip_wetdry/Fig2a_deseasoned_only_TOP100/Fig2a_WGPE_precip_correlation_summary_TOP100.csv"),
  data.frame(panel = "b", group = "Global", component = "ENSO_global_WGPE_anomaly",
             value = enso_global, unit = "J m^-2",
             source = "output_TOP100_revision/results/ENSO/ONI_standard_TOP100/03_absolute_results/ENSO_ONIstandard_absolute_global_summary.csv"),
  data.frame(panel = "c", group = as.character(wetdry_long$event_label), component = as.character(wetdry_long$component),
             value = wetdry_long$value*10000, unit = "J m^-2",
             source = "FINAL_SCIENTIFIC_ARCHIVE/figure_source_data/Figure2c_WET_DRY_FINAL_SOURCE.csv"),
  data.frame(panel = "d", group = as.character(enso_long$region), component = as.character(enso_long$component),
             value = enso_long$value, unit = "J m^-2",
             source = "output_TOP100_revision/results/ENSO/ONI_standard_TOP100/03_absolute_results/ENSO_ONIstandard_absolute_channel_decomposition.csv")
)
write.csv(source_values, file.path(OUT, "Fig3_wetdry_cdmerge_plus_enso_v1_source_values.csv"), row.names = FALSE)
write.csv(source_values, file.path(OUT, "Figure2_TOP100_source_values.csv"), row.names = FALSE)

qa <- c(
  "# Figure 2 TOP100 FINAL v5 true-increment QA",
  "",
  "Scope: TOP100 Figure 2 plotting and source-data audit. Core calculations were not rerun during plotting.",
  "",
  "1. Core scientific calculations rerun during plotting: NO; all panels read formal TOP100 outputs.",
  "2. Panel e corrected scientifically: YES; the plotting frame is retained, but the source is the corrected true-increment M0-M3 table.",
  "3. Wet/dry total anomaly plus wet/dry channel decomposition merged into panel c: YES; grouped horizontal bars, not stacked.",
  "4. ENSO standard-climate-band signed stacked bar added as panel d: YES.",
  "4b. Panels c and d arranged side-by-side: YES.",
  "5. Panel a uses calendar-month deseasonalized TOP100 W-GPE and precipitation anomalies without linear detrending; IWV is not controlled. The detrended result is retained as sensitivity only.",
  "6. PDF output generated: YES, as requested by the v3 layout task.",
  "7. Old project figures overwritten: NO; outputs are in `output_TOP100_revision/figures/main/Figure2_TOP100_FINAL`.",
  "8. DOCX/SI/workbook modified: NO.",
  "",
  "Panel e four-model definitions:",
  "- M0 = IWV.",
  "- M1 = IWV + <z>.",
  "- M2 = IWV + <z> + IWV:<z>.",
  "- M3 = IWV + W-GPE.",
  "- Old single-predictor W-GPE-only M3 is excluded from the main figure and retained only as a legacy/SI benchmark.",
  "- Metric = area-weighted mean blocked-CV R2 with 95% bootstrap confidence intervals.",
  "- Source file: CANONICAL_WGPE_TOP100_v2/04_RESULTS_CANONICAL/FIG2E_TRUE_INCREMENT/FIG2E_TRUE_INCREMENT_RESULTS.csv.",
  "",
  "Panel c checks:",
  "- Rows: Wet events and Dry events.",
  "- Bars per row: Total, IWV, <z>, Interaction.",
  "- Total bar height matches the other channel bars, and c bar height is increased to align better with panel d.",
  "- Values are signed J m^-2 anomalies.",
  "- Stacked bars: NO.",
  "- Legend is placed below the x-axis title with increased spacing.",
  "- In the merged figure, the c/d block height is compressed to 0.8 of the previous local row height and centred.",
  "",
  "Panel d checks:",
  "- Data source: formal TOP100 ONI-standard exact decomposition under `output_TOP100_revision/results/ENSO/ONI_standard_TOP100/03_absolute_results/`.",
  "- Bands: Global, S high latitude, S temperate, Tropics, N temperate, N high latitude.",
  "- Bars per band: IWV, <z>, Interaction; signed stacked bars following the existing candidate panel.",
  "- Values are signed J m^-2 anomalies.",
  "- x-axis range adjusted to -5,000 to +32,000 J m^-2 to reduce excess negative whitespace while preserving negative bars.",
  "- Legend is placed below the x-axis title with increased spacing.",
  "",
  "Panel e checks:",
  "- Uses the original temporal-precedence block from Fig3_layout_fix_v3.",
  "- Remains labelled e."
)
writeLines(qa, file.path(OUT, "Fig3_wetdry_cdmerge_plus_enso_v1_QA.md"), useBytes = TRUE)
writeLines(qa, file.path(OUT, "Figure2_TOP100_FINAL_QA.md"), useBytes = TRUE)

cat("DONE Figure2 TOP100 layout v5 panel-frame-fix\n")
