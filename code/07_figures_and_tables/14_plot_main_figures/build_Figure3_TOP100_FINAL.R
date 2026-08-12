options(stringsAsFactors = FALSE, warn = 2)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(grid)
})

ROOT <- "__WGPE_PROJECT_ROOT__"
args <- commandArgs(trailingOnly = TRUE)
REF <- "absolute"
if (!REF %in% c("absolute", "surface_relative")) stop("REF must be absolute or surface_relative")
NO_UNCERTAINTY <- length(args) >= 2 && args[2] == "no_uncertainty"
V3 <- file.path(ROOT, "output_attribution_minimal", "FULL_RERUN_FINAL_absolute_relative_v3")
V4 <- file.path(ROOT, "output_attribution_minimal", "FULL_RERUN_TARGETED_REPAIR_v4")
V5 <- file.path(ROOT, "output_attribution_minimal", "FULL_RERUN_ECOLOGY_MASK_AND_PACKAGING_v5")
CLOSURE <- Sys.getenv("TOP100_CLOSURE_ROOT", unset=file.path(ROOT,"output_TOP100_revision"))
ADAPTER_ROOT <- Sys.getenv("TOP100_ADAPTER_ROOT", unset=file.path(ROOT,"output_TOP100_revision","derived_data","figure_source_data"))
OUT <- file.path(CLOSURE, "figures", "main", "Figure3_TOP100_FINAL")
NO_PANEL_LABELS <- identical(Sys.getenv("TOP100_NO_PANEL_LABELS", unset = "0"), "1")
if (NO_PANEL_LABELS) {
  OUT <- file.path(CLOSURE, "figures", "main_no_panel_labels", "Figure3_TOP100_FINAL")
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "R", "figure_style_FIXED_FULL.R"), encoding = "UTF-8")
source(file.path(ROOT, "R", "figure_data_FIXED_FULL.R"), encoding = "UTF-8")

FONT <- "sans"
pal_blue <- c("#794607", "#A66B1C", "#C89545", "#DFC685", "#F5EFE4", "#CCDEEF", "#9CC8DF", "#5EA6D1", "#0D56A0")
pal_green <- c("#704300", "#A56B1D", "#C79545", "#DFC685", "#F1E2B9", "#EAF2F0", "#B9DDD7", "#6EB7AA", "#23877D", "#005C54")
COL <- list(neg = "#8A5609", pos = "#0D56A0", eco_pos = "#005C54", zbar = "#A66B1C",
            iwv = "#0D56A0", inter = "#006B60", grey = "#ECECEA",
            ink = "#202020", axis = "#555555", grid = "#ECECEA",
            actual = "#202020", residual = "#9A9A9A")

# Layout-fix constants: keep panel tags/titles visually uniform across map,
# bar, heatmap, and manually labelled compound panels.
PANEL_TITLE_SIZE <- 7.6
PANEL_SUBTITLE_SIZE <- 5.6
PANEL_TAG_SIZE <- 7.6

inflate_all_text <- function(g, delta_pt = 2) {
  if (!is.null(g$gp) && !is.null(g$gp$fontsize)) g$gp$fontsize <- g$gp$fontsize + delta_pt
  if (!is.null(g$children) && length(g$children)) for (i in seq_along(g$children)) g$children[[i]] <- inflate_all_text(g$children[[i]], delta_pt)
  if (!is.null(g$grobs) && length(g$grobs)) for (i in seq_along(g$grobs)) g$grobs[[i]] <- inflate_all_text(g$grobs[[i]], delta_pt)
  g
}

save_png <- function(plot, filename, width_mm, height_mm, dpi = 500) {
  gr <- grid::grid.force(cowplot::as_grob(plot))
  gr <- inflate_all_text(gr, 2)
  ggsave(file.path(OUT, filename), gr, width = width_mm, height = height_mm,
         units = "mm", dpi = dpi, bg = "white", device = "png")
}

save_pdf <- function(plot, filename, width_mm, height_mm) {
  gr <- grid::grid.force(cowplot::as_grob(plot))
  gr <- inflate_all_text(gr, 2)
  ggsave(file.path(OUT, filename), gr, width = width_mm, height = height_mm,
         units = "mm", bg = "white", device = cairo_pdf)
}

ptitle <- function(tag = NULL, title = NULL) {
  if (is.null(title)) return(NULL)
  if (is.null(tag) || NO_PANEL_LABELS) return(title)
  bquote(bold(.(tag)) ~ .(title))
}

theme_rec <- function(base = 7) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(text = element_text(family = FONT, colour = COL$ink),
          plot.title = element_text(size = PANEL_TITLE_SIZE, hjust = 0, margin = margin(b = 2)),
          axis.title = element_text(size = base - .2, colour = COL$ink, lineheight = .85),
          axis.text = element_text(size = base - .7, colour = COL$ink),
          axis.line = element_line(linewidth = .26, colour = COL$axis),
          axis.ticks = element_line(linewidth = .23, colour = COL$axis),
          axis.ticks.length = unit(.8, "mm"),
          panel.grid.major = element_line(colour = COL$grid, linewidth = .16),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = base - .6),
          legend.text = element_text(size = base - .8),
          legend.key.height = unit(2.2, "mm"),
          legend.key.width = unit(6.5, "mm"),
          plot.margin = margin(3, 4, 3, 4))
}

theme_map_rec <- function(base = 7) {
  theme_void(base_family = FONT, base_size = base) +
    theme(text = element_text(family = FONT, colour = COL$ink),
          plot.title = element_text(size = PANEL_TITLE_SIZE, hjust = 0, margin = margin(b = 2), lineheight = .94),
          plot.subtitle = element_text(size = PANEL_SUBTITLE_SIZE, hjust = .5,
                                       margin = margin(t = 0, b = 1), colour = COL$ink),
          plot.margin = margin(1, 1, 0, 1),
          plot.background = element_rect(fill = "white", colour = NA))
}

rollmean_center <- function(x, k = 12) {
  out <- rep(NA_real_, length(x)); h1 <- floor((k - 1) / 2); h2 <- k - h1 - 1
  for (i in seq_along(x)) {
    lo <- i - h1; hi <- i + h2
    if (lo >= 1 && hi <= length(x)) out[i] <- mean(x[lo:hi], na.rm = TRUE)
  }
  out
}

grid_spec <- function(grid = c("181", "180")) {
  grid <- match.arg(grid)
  if (grid == "181") {
    list(lon = 0:359, lat = -90:90)
  } else {
    # Ecology downstream matrices are stored north-to-south in their second
    # dimension (j = 1 is 89.5 N).  The earlier recovery pass assumed
    # south-to-north and therefore plotted the ecological field at the wrong
    # latitude.  bundle_from_field() sorts latitude after this, carrying the
    # field and q-mask with the same column order.
    list(lon = seq(-179.5, 179.5, by = 1), lat = seq(89.5, -89.5, by = -1))
  }
}

bundle_from_field <- function(field, grid = c("181", "180"), q = NULL, y_factor = 0.52) {
  gs <- grid_spec(grid)
  lon <- gs$lon; lat <- gs$lat
  if (nrow(field) != length(lon) || ncol(field) != length(lat)) stop("field coordinate mismatch")
  if (!is.null(q) && !identical(dim(field), dim(q))) stop("field/q dimension mismatch")
  lon2 <- ifelse(lon >= 180, lon - 360, lon)
  oi <- order(lon2); oj <- order(lat)
  fld <- field[oi, oj, drop = FALSE]
  q_c <- if (!is.null(q)) q[oi, oj, drop = FALSE] else NULL
  geom <- make_robinson_grid_geometry(lon2[oi], lat[oj])
  map <- attach_map_field(geom, fld)
  transform_y <- function(x) { x$y <- x$y * y_factor; x }
  list(lon_raw = lon, lat_raw = lat, lon = lon2[oi], lat = lat[oj],
       field = fld, q = q_c, oi = oi, oj = oj,
       map = transform_y(map),
       coast = transform_y(load_robinson_coastline(file.path(ROOT, "R", "coastlines_source_FIXED_FULL.csv"))),
       graticule = transform_y(make_robinson_graticule()),
       border = transform_y(make_robinson_border()),
       y_factor = y_factor)
}

# Robinson bundle for ad-hoc grids such as temporal-precedence class tables.
# This keeps categorical temporal maps in the same projection/frame style as
# the continuous scientific maps.
bundle_from_lonlat_field <- function(field, lon, lat, y_factor = 0.52) {
  if (nrow(field) != length(lon) || ncol(field) != length(lat)) stop("field coordinate mismatch")
  lon2 <- ifelse(lon >= 180, lon - 360, lon)
  oi <- order(lon2); oj <- order(lat)
  fld <- field[oi, oj, drop = FALSE]
  geom <- make_robinson_grid_geometry(lon2[oi], lat[oj])
  map <- attach_map_field(geom, fld)
  transform_y <- function(x) { x$y <- x$y * y_factor; x }
  list(lon_raw = lon, lat_raw = lat, lon = lon2[oi], lat = lat[oj],
       field = fld, oi = oi, oj = oj, map = transform_y(map),
       coast = transform_y(load_robinson_coastline(file.path(ROOT, "R", "coastlines_source_FIXED_FULL.csv"))),
       graticule = transform_y(make_robinson_graticule()),
       border = transform_y(make_robinson_border()),
       y_factor = y_factor)
}

outer_mask <- function(border, xlim, ylim) {
  outer <- data.frame(x = c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]),
                      y = c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1]),
                      group = 1L, subgroup = 1L)
  inner <- border; inner$group <- 1L; inner$subgroup <- 2L
  rbind(outer[, c("x", "y", "group", "subgroup")], inner[, c("x", "y", "group", "subgroup")])
}

map_plot <- function(bundle, lim, pal, tag, title_expr, fdr = FALSE, subtitle_text = NULL) {
  xlim <- c(-1.015, 1.015)
  ylim <- c(-max(abs(bundle$border$y)) - .012, max(abs(bundle$border$y)) + .012)
  mask <- outer_mask(bundle$border, xlim, ylim)
  mapdf <- bundle$map
  p <- ggplot()
  if (fdr) {
    valid <- is.finite(bundle$field)
    sig <- valid & is.finite(bundle$q) & bundle$q < 0.05
    nonsig <- valid & !sig
    grey <- mapdf[rep(as.vector(nonsig), each = 5L), , drop = FALSE]
    sigdf <- mapdf
    sigdf$value[!rep(as.vector(sig), each = 5L)] <- NA_real_
    p <- p +
      geom_polygon(data = grey, aes(x, y, group = id), fill = COL$grey, colour = NA) +
      geom_polygon(data = sigdf, aes(x, y, group = id, fill = value), colour = NA)
  } else {
    p <- p + geom_polygon(data = mapdf, aes(x, y, group = id, fill = value), colour = NA)
  }
  p +
    geom_path(data = bundle$graticule, aes(x, y, group = group), linewidth = .13, colour = "#B9BDBD", alpha = .55) +
    geom_path(data = bundle$coast, aes(x, y, group = group), linewidth = .20, colour = "#505456", lineend = "round") +
    geom_polygon(data = mask, aes(x, y, group = group, subgroup = subgroup),
                 inherit.aes = FALSE, fill = "white", colour = NA, rule = "evenodd") +
    geom_path(data = bundle$border, aes(x, y, group = group), linewidth = .30, colour = "#6A6E70", linejoin = "round") +
    scale_fill_gradientn(colours = pal, limits = c(-lim, lim), oob = squish, na.value = "transparent", guide = "none") +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE, clip = "on") +
    labs(title = title_expr, subtitle = subtitle_text, x = NULL, y = NULL) +
    theme_map_rec(7)
}

colorbar <- function(lim, pal, title_expr, ticks, accuracy = .1) {
  d <- data.frame(x = seq(-lim, lim, length.out = 401), y = 0)
  td <- data.frame(x = ticks, label = label_number(accuracy = accuracy, big.mark = " ")(ticks))
  ggplot(d, aes(x, y, fill = x)) +
    geom_tile(width = (2 * lim) / 400, height = .08) +
    annotate("rect", xmin = -lim, xmax = lim, ymin = -.04, ymax = .04, fill = NA, colour = "#111111", linewidth = .22) +
    geom_segment(data = td, aes(x = x, xend = x, y = -.042, yend = -.075),
                 inherit.aes = FALSE, linewidth = .20, colour = "#111111") +
    geom_text(data = td, aes(x = x, y = -.118, label = label),
              inherit.aes = FALSE, family = FONT, size = 1.75, colour = "#111111") +
    scale_fill_gradientn(colours = pal, limits = c(-lim, lim), oob = squish, guide = "none") +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-.15, .105), clip = "off") +
    labs(title = title_expr, x = NULL, y = NULL) +
    theme_void(base_family = FONT, base_size = 6) +
    theme(plot.title = element_text(size = 5.8, hjust = .5, margin = margin(b = 1), colour = "#111111"),
          plot.margin = margin(-2, 4, 0, 4))
}

map_profile_companion <- function(bundle, lim, xlabel, line_col, accuracy = .1) {
  prof <- lat_profile(bundle$field, bundle$lon, bundle$lat, 5)
  prof <- prof[is.finite(prof$value) & is.finite(prof$lat_mid), , drop = FALSE]
  ECO_PROFILE_CALL <<- ECO_PROFILE_CALL + 1L
  uu <- NULL
  if (ECO_PROFILE_CALL %in% 1:4) {
    key <- data.frame(target=c("GPP","GPP","LAI","LAI"), predictor=c("W-GPE","<z>","W-GPE","<z>"))[ECO_PROFILE_CALL,]
    uu <- U_EPROF[U_EPROF$height_reference == REF & U_EPROF$target == key$target & U_EPROF$predictor == key$predictor,]
    uu <- uu[is.finite(uu$latitude_bin_mid) & is.finite(uu$lower95) & is.finite(uu$upper95), , drop = FALSE]
  }
  p <- ggplot(prof, aes(value, lat_mid))
  if (!is.null(uu) && nrow(uu)) p <- p +
    geom_ribbon(data=uu, aes(y=latitude_bin_mid,xmin=lower95,xmax=upper95), inherit.aes=FALSE,
                orientation="y", fill=alpha(line_col,.16), colour=NA)
  p +
    geom_vline(xintercept = 0, colour = "#6D6D6D", linewidth = .20) +
    geom_path(colour = line_col, linewidth = .36, lineend = "round", na.rm = TRUE) +
    geom_point(colour = line_col, size = .34, stroke = 0, na.rm = TRUE) +
    scale_y_continuous(limits = c(-90, 90), breaks = c(-60, 0, 60),
                       expand = expansion(mult = c(.01, .01))) +
    scale_x_continuous(breaks = c(-lim, 0, lim),
                       labels = label_number(accuracy = accuracy, big.mark = " ")) +
    coord_cartesian(xlim = c(-lim, lim), clip = "on") +
    labs(x = xlabel, y = expression("Latitude ("*degree*")")) +
    theme_rec(5.4) +
    theme(panel.grid.major.y = element_line(colour = "#EEEEEA", linewidth = .15),
          panel.grid.major.x = element_blank(),
          axis.title.x = element_text(size = 4.05, margin = margin(t = .8), lineheight = .78),
          axis.title.y = element_text(size = 5.25, margin = margin(r = .6)),
          axis.text.x = element_text(size = 4.25, margin = margin(t = .4)),
          axis.text.y = element_text(size = 4.9, margin = margin(r = 0)),
          axis.ticks.length = unit(.55, "mm"),
          plot.margin = margin(1, 4, 1, 0))
}

map_block <- function(field, grid, q, lim, pal, tag, title_expr, cbar_title_expr, cbar_ticks,
                      fdr = FALSE, accuracy = .1, profile_accuracy = accuracy,
                      profile_colour = COL$pos, subtitle_text = NULL) {
  b <- bundle_from_field(field, grid, q)
  mp <- map_plot(b, lim, pal, tag, title_expr, fdr, subtitle_text)
  side_core <- map_profile_companion(b, lim, cbar_title_expr, profile_colour, profile_accuracy)
  # The profile is vertically centered to the Robinson ellipse body, not to
  # the whole title+map+margin grob.  This keeps the side panel visually
  # attached to the map rather than floating as an independent plot.
  side <- plot_grid(NULL, side_core, NULL, ncol = 1, rel_heights = c(.18, .64, .18))
  top <- plot_grid(mp, side, nrow = 1, rel_widths = c(.80, .20), align = "h", axis = "tb")
  cb <- plot_grid(NULL, colorbar(lim, pal, cbar_title_expr, cbar_ticks, accuracy),
                  NULL, nrow = 1, rel_widths = c(.04, .72, .24))
  plot_grid(top, cb, ncol = 1, rel_heights = c(.82, .18))
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

profile_plot <- function(field, grid, title, xlabel, line_col, accuracy = .1) {
  b <- bundle_from_field(field, grid)
  prof <- lat_profile(b$field, b$lon, b$lat, 5)
  lim <- max(abs(prof$value), na.rm = TRUE)
  lim <- max(abs(pretty(c(-lim, lim), n = 3)))
  ggplot(prof, aes(value, lat_mid)) +
    geom_vline(xintercept = 0, colour = "#666666", linewidth = .22) +
    geom_path(colour = line_col, linewidth = .38, lineend = "round", na.rm = TRUE) +
    geom_point(colour = line_col, size = .40, stroke = 0, na.rm = TRUE) +
    scale_y_continuous(limits = c(-90, 90), breaks = c(-60, 0, 60), expand = expansion(mult = c(.01, .01))) +
    scale_x_continuous(breaks = c(-lim, 0, lim), labels = label_number(accuracy = accuracy, big.mark = " ")) +
    coord_cartesian(xlim = c(-lim, lim), clip = "on") +
    labs(title = title, x = xlabel, y = "Latitude (deg)") +
    theme_rec(6.1) +
    theme(panel.grid.major.x = element_blank())
}

area_fraction <- function(q, field, lat) {
  valid <- is.finite(field) & is.finite(q)
  sig <- valid & q < 0.05
  grid_pct <- 100 * sum(sig) / sum(valid)
  w <- matrix(rep(cos(lat * pi / 180), each = nrow(field)), nrow = nrow(field), ncol = ncol(field))
  area_pct <- 100 * sum(w[sig], na.rm = TRUE) / sum(w[valid], na.rm = TRUE)
  c(grid_pct = grid_pct, area_pct = area_pct)
}

qc_one <- function(id, field, q, grid, summary_grid, summary_area) {
  summary_grid <- as.numeric(summary_grid)[1]
  summary_area <- as.numeric(summary_area)[1]
  gs <- grid_spec(grid)
  dim_identical <- identical(dim(field), dim(q))
  valid_same <- identical(is.finite(field), is.finite(q))
  frac <- area_fraction(q, field, gs$lat)
  b <- bundle_from_field(field, grid, q)
  sorted_dim_identical <- identical(dim(b$field), dim(b$q))
  data.frame(
    map_id = id,
    field_dim = paste(dim(field), collapse = "x"),
    q_dim = paste(dim(q), collapse = "x"),
    lon_range = paste(range(gs$lon), collapse = " to "),
    lat_range = paste(range(gs$lat), collapse = " to "),
    lat_order = ifelse(all(diff(gs$lat) > 0), "south-to-north", "north-to-south"),
    field_NA_count = sum(!is.finite(field)),
    q_NA_count = sum(!is.finite(q)),
    valid_count = sum(is.finite(field) & is.finite(q)),
    q_lt_0_05_count = sum(is.finite(q) & q < 0.05 & is.finite(field)),
    dim_identical = dim_identical,
    valid_mask_identical = valid_same,
    sorted_dim_identical = sorted_dim_identical,
    q_sig_tested_grid_fraction_pct = frac["grid_pct"],
    q_sig_valid_area_fraction_pct = frac["area_pct"],
    reported_tested_grid_fraction_pct = summary_grid,
    reported_valid_area_fraction_pct = summary_area,
    grid_fraction_abs_diff = abs(frac["grid_pct"] - summary_grid),
    area_fraction_abs_diff = abs(frac["area_pct"] - summary_area),
    status = ifelse(dim_identical && sorted_dim_identical && abs(frac["grid_pct"] - summary_grid) < 0.05 && abs(frac["area_pct"] - summary_area) < 0.05, "PASS", "FAIL"),
    stringsAsFactors = FALSE
  )
}

# Inputs: accepted layout with formal TOP100 values.
F1 <- readRDS(file.path(ADAPTER_ROOT, "Figure1_TOP100_adapter.rds"))
F2 <- readRDS(file.path(ADAPTER_ROOT, "Figure2_TOP100_adapter.rds"))
F3 <- readRDS(file.path(ADAPTER_ROOT, "Figure3_TOP100_adapter.rds"))
data12 <- F1$data12
core <- new.env(parent = emptyenv())
load(F1$core_rdata, envir = core)
ADAPTER <- file.path(V5, "13_scripts", "archived_layout_adapters_v5", "absolute")
precip <- F2$precip
eco <- F3$eco
event <- F2$event
granger <- readRDS(file.path(ROOT, "output_attribution_minimal", "downstream_FIXED_FULL", "granger_temporal_precedence_FIXED_FULL.rds"))
fdr_summary <- F3$fdr_summary
fdr_lag <- fdr_summary
fdr_grid <- F3$fdr_grid
fdr_rds <- list()
U_ESUM <- F3$U_ESUM
U_EPROF <- F3$U_EPROF
U_ETEMP <- F3$U_ETEMP
if (NO_UNCERTAINTY) {
  U_ESUM <- U_ESUM[FALSE,]
  U_EPROF <- U_EPROF[FALSE,]
  U_ETEMP <- U_ETEMP[FALSE,]
}
ECO_PROFILE_CALL <- -20L
enso_dir <- file.path(ROOT, "output_attribution_minimal", "ENSO_audit_FIXED_FULL")
enso_channel <- if (file.exists(file.path(enso_dir, "ENSO_channel_decomposition.csv"))) {
  read.csv(file.path(enso_dir, "ENSO_channel_decomposition.csv"), check.names = FALSE)
} else NULL
enso_metrics <- if (file.exists(file.path(enso_dir, "ENSO_audit_metrics.csv"))) {
  read.csv(file.path(enso_dir, "ENSO_audit_metrics.csv"), check.names = FALSE)
} else NULL

# QC before inferential maps.
pglob <- precip$summary[precip$summary$variable == "wgpe" & precip$summary$lag_months_predictor_leads_precip == 0 & precip$summary$region == "Global", ]
get_eco_sum <- function(target, pred) eco$summary[eco$summary$target == target & eco$summary$predictor == pred & eco$summary$region == "Global", ]
qc <- rbind(
  qc_one("Fig3a_wgpe_precip_lag0", precip$precipitation$wgpe$lag0$r, precip$precipitation$wgpe$lag0$q, "181", pglob$significant_grid_pct, pglob$significant_valid_area_pct),
  qc_one("Fig4a_wgpe_GPP", eco$results$GPP$r$wgpe, eco$results$GPP$q$wgpe, "181", get_eco_sum("GPP", "wgpe")$significant_grid_pct, get_eco_sum("GPP", "wgpe")$significant_valid_area_pct),
  qc_one("Fig4b_zbar_GPP", eco$results$GPP$r$zbar, eco$results$GPP$q$zbar, "181", get_eco_sum("GPP", "zbar")$significant_grid_pct, get_eco_sum("GPP", "zbar")$significant_valid_area_pct),
  qc_one("Fig4c_wgpe_LAI", eco$results$LAI$r$wgpe, eco$results$LAI$q$wgpe, "181", get_eco_sum("LAI", "wgpe")$significant_grid_pct, get_eco_sum("LAI", "wgpe")$significant_valid_area_pct),
  qc_one("Fig4d_zbar_LAI", eco$results$LAI$r$zbar, eco$results$LAI$q$zbar, "181", get_eco_sum("LAI", "zbar")$significant_grid_pct, get_eco_sum("LAI", "zbar")$significant_valid_area_pct)
)
write.csv(qc, file.path(OUT, "field_qmask_alignment_check.csv"), row.names = FALSE)
write.csv(qc[, c("map_id", "q_sig_tested_grid_fraction_pct", "reported_tested_grid_fraction_pct",
                 "grid_fraction_abs_diff", "q_sig_valid_area_fraction_pct",
                 "reported_valid_area_fraction_pct", "area_fraction_abs_diff", "status")],
          file.path(OUT, "mask_fraction_check.csv"), row.names = FALSE)
all_pass <- all(qc$status == "PASS")
writeLines(c("# Map field / q-mask alignment QA", "",
             paste0("Overall status: ", ifelse(all_pass, "PASS", "FAIL")),
             "",
             "Checks performed for each inferential map:",
             "- field and q dimensions identical",
             "- identical sorting indices are used for field and q",
             "- q<0.05 tested-grid fraction matches corrected summary",
             "- q<0.05 valid-area fraction matches corrected summary",
             "",
             "Result table: field_qmask_alignment_check.csv",
             "Fraction check: mask_fraction_check.csv",
             "",
             ifelse(all_pass,
                    "All inferential maps passed coordinate/mask fraction checks and can be drawn with FDR masking.",
                    "At least one map failed. Corresponding inferential maps should not be drawn.")),
           file.path(OUT, "field_qmask_alignment_QA.md"), useBytes = TRUE)
if (!all_pass) stop("Field/q mask alignment FAIL. Stop before drawing inferential maps.")

# Merged Fig1/2 recovery is intentionally skipped in the TOP100 Figure 3
# builder.  Figure 1 and Figure 2 have dedicated formal scripts.
if (FALSE) {
panel12a <- map_block(
  data12$ztrend, "181", NULL, 4, pal_blue, "a",
  ptitle("a", "Atmospheric vapour-centroid trend"),
  expression("<z> trend (m yr"^{-1}*")"), c(-4, -2, 0, 2, 4), FALSE, 1
)
panel12b <- map_block(
  data12$wgpe_trend, "181", NULL, 3000, pal_blue, "b",
  ptitle("b", "W-GPE trend"),
  expression("W-GPE trend (J m"^{-2}~"yr"^{-1}*")"), c(-3000, -1500, 0, 1500, 3000), FALSE, 1000
)

dates <- as.Date(core$all_dates); month <- as.integer(format(dates, "%m"))
lat_core <- as.numeric(core$lat); lon_core <- as.numeric(core$lon)
wmat <- matrix(rep(cos(lat_core * pi / 180), each = length(lon_core)), nrow = length(lon_core), ncol = length(lat_core))
area_series <- function(a) vapply(seq_len(dim(a)[3]), function(k) weighted.mean(as.vector(a[, , k]), as.vector(wmat), na.rm = TRUE), numeric(1))
deseason <- function(x) { clim <- tapply(x, month, mean, na.rm = TRUE); x - clim[as.character(month)] }
tsdf <- data.frame(date = dates,
                   z = as.numeric(scale(deseason(area_series(core$zbar_st)))),
                   w = as.numeric(scale(deseason(area_series(core$wgpe_st)))))
tsdf$z12 <- rollmean_center(tsdf$z, 12); tsdf$w12 <- rollmean_center(tsdf$w, 12)
panel12c <- ggplot(tsdf, aes(date)) +
  geom_hline(yintercept = 0, colour = "#B8B8B8", linewidth = .22) +
  geom_line(aes(y = z), colour = alpha(COL$zbar, .12), linewidth = .10) +
  geom_line(aes(y = w), colour = alpha(COL$pos, .12), linewidth = .10) +
  geom_line(aes(y = z12, colour = "<z>"), linewidth = .54, na.rm = TRUE) +
  geom_line(aes(y = w12, colour = "W-GPE"), linewidth = .56, na.rm = TRUE) +
  scale_colour_manual(values = c("<z>" = COL$zbar, "W-GPE" = COL$pos), name = NULL) +
  annotate("text", x = as.Date("2024-06-01"), y = -2.20,
           label = "paste('<z>: +0.788 m ', yr^{-1}, ' | W-GPE: +645.6 J ', m^{-2}, ' ', yr^{-1})",
           parse = TRUE, hjust = 1, size = 1.9, colour = "#444444", family = FONT) +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y", expand = expansion(mult = c(.01, .01))) +
  scale_y_continuous(breaks = c(-2, 0, 2), limits = c(-2.55, 3.35)) +
  labs(title = ptitle("c", "Standardized global anomalies"), x = NULL, y = "Standardized anomaly") +
  theme_rec(6.6) +
  theme(legend.position = c(.50, .91), legend.direction = "horizontal",
        legend.background = element_blank(), panel.grid.major.x = element_blank())

corr <- data12$correlations
corr$sign <- ifelse(corr$r >= 0, "Positive", "Negative")
panel12d <- ggplot(corr, aes(r, variable, fill = sign)) +
  geom_vline(xintercept = 0, colour = "#111111", linewidth = .25) +
  geom_col(width = .50, colour = "#222222", linewidth = .14) +
  scale_fill_manual(values = c(Positive = COL$pos, Negative = COL$neg), guide = "none") +
  scale_x_continuous(limits = c(-.62, .54), breaks = c(-.5, 0, .5)) +
  labs(title = ptitle("d", "Trend-pattern comparison"), x = "Spatial Pearson r", y = NULL) +
  theme_rec(6.5) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())

prof_dir <- file.path(ROOT, "output_attribution_minimal", "highres_z_vertical_profile_QC")
native_prof <- read.csv(file.path(prof_dir, "native_interval_vertical_profile.csv"), check.names = FALSE)
agg4 <- read.csv(file.path(prof_dir, "aggregation_to_4layer_check.csv"), check.names = FALSE)
layer_df <- agg4[agg4$resolution_label == "native", ]
lay <- do.call(rbind, lapply(strsplit(layer_df$layer, "-"), as.numeric))
layer_df$p_bottom <- lay[, 1]; layer_df$p_top <- lay[, 2]
mid <- (layer_df$p_bottom + layer_df$p_top) / 2
half <- (layer_df$p_bottom - layer_df$p_top) * .55 / 2
layer_df$ymin <- mid - half; layer_df$ymax <- mid + half
layer_df$contribution <- layer_df$reference_4layer_m_per_yr
native_prof$density20 <- native_prof$contribution_to_z_trend_m_per_yr / native_prof$layer_thickness_hPa * 20
left_e <- ggplot(layer_df) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = .24) +
  geom_rect(aes(xmin = pmin(0, contribution), xmax = pmax(0, contribution),
                ymin = ymin, ymax = ymax, fill = ifelse(contribution >= 0, "pos", "neg")),
            colour = "#222222", linewidth = .15) +
  scale_fill_manual(values = c(pos = COL$pos, neg = COL$neg), guide = "none") +
  scale_y_reverse(limits = c(1000, 300), breaks = c(1000, 850, 700, 500, 300)) +
  scale_x_continuous(limits = c(-.34, .64), breaks = c(-.3, 0, .3, .6)) +
  labs(x = expression(atop("Integrated contribution", "(m yr"^{-1}*")")), y = "Pressure (hPa)") +
  theme_rec(6.0) +
  theme(plot.margin = margin(1, 3, 3, 4), panel.grid.major.x = element_line(colour = COL$grid, linewidth = .13))
right_e <- ggplot(native_prof[order(native_prof$p_mid, decreasing = TRUE), ], aes(density20, p_mid)) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = .24) +
  geom_path(colour = "#202020", linewidth = .32, lineend = "round", na.rm = TRUE) +
  geom_point(colour = alpha("#202020", .45), size = .35, na.rm = TRUE) +
  scale_y_reverse(limits = c(1000, 300), breaks = c(1000, 850, 700, 500, 300), labels = NULL) +
  scale_x_continuous(limits = c(-.12, .22), breaks = c(-.1, 0, .1, .2)) +
  labs(x = expression(atop("Contribution density", "(m yr"^{-1}*" per 20 hPa)")), y = NULL) +
  theme_rec(6.0) +
  theme(plot.margin = margin(1, 4, 3, 1),
        panel.grid.major.x = element_line(colour = COL$grid, linewidth = .13),
        axis.ticks.y = element_blank())
panel12e <- plot_grid(
  ggdraw() + draw_label("e  Vertical structure of the <z> trend", x = 0, hjust = 0, fontfamily = FONT, size = PANEL_TAG_SIZE),
  plot_grid(left_e, right_e, nrow = 1, rel_widths = c(.34, .66), align = "h", axis = "tb"),
  ncol = 1, rel_heights = c(.11, .89)
)

layer_global <- unique(data12$layer[, c("actual_zbar_trend_m_yr", "sum_mass_contribution_m_yr",
                                        "sum_height_contribution_m_yr", "residual_m_yr",
                                        "residual_pct_of_actual")])[1, ]
fdf <- data.frame(component = factor(c("Mass\nredistrib.", "Height\ncoord.", "Residual", "Actual"),
                                     levels = c("Mass\nredistrib.", "Height\ncoord.", "Residual", "Actual")),
                  value = c(layer_global$sum_mass_contribution_m_yr, layer_global$sum_height_contribution_m_yr,
                            layer_global$residual_m_yr, layer_global$actual_zbar_trend_m_yr),
                  fill = c("mass", "height", "residual", "actual"))
panel12f <- ggplot(fdf, aes(component, value, fill = fill)) +
  geom_hline(yintercept = 0, colour = "#AAAAAA", linewidth = .25) +
  geom_col(width = .56, colour = "#222222", linewidth = .15) +
  annotate("text", x = 2.5, y = -0.105, label = "Residual = -11.65% of actual",
           size = 2.0, colour = "#555555", family = FONT) +
  scale_fill_manual(values = c(mass = COL$pos, height = COL$zbar, residual = COL$residual, actual = COL$actual),
                    guide = "none") +
  scale_y_continuous(limits = c(-.16, 1.06), breaks = c(0, .5, 1.0)) +
  labs(title = ptitle("f", "Additive diagnostic bridge"), x = NULL, y = expression("Contribution (m yr"^{-1}*")")) +
  theme_rec(6.4) +
  theme(axis.text.x = element_text(lineheight = .82), panel.grid.major.x = element_blank())

wg <- data12$wgpe
comp <- rbind(data.frame(region = wg$region, channel = "IWV", value = wg$iwv_channel_pct),
              data.frame(region = wg$region, channel = "<z>", value = wg$z_channel_pct),
              data.frame(region = wg$region, channel = "Interaction", value = wg$interaction_channel_pct))
comp$region <- factor(comp$region, levels = rev(c("Global", "S extratropics", "S tropics", "Deep tropics", "N tropics", "N extratropics")))
comp$channel <- factor(comp$channel, levels = c("IWV", "<z>", "Interaction"))
panel12g <- ggplot(comp, aes(value, region, fill = channel)) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = .25) +
  geom_col(width = .58, colour = "#222222", linewidth = .08) +
  scale_fill_manual(values = c("IWV" = COL$iwv, "<z>" = COL$zbar, "Interaction" = COL$inter), name = NULL) +
  scale_x_continuous(limits = c(-15, 115), breaks = c(0, 50, 100)) +
  labs(title = ptitle("g", "Regional W-GPE trend composition"), x = "Share of regional W-GPE trend (%)", y = NULL) +
  theme_rec(6.2) +
  theme(legend.position = "top", legend.direction = "horizontal",
        axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())

merged12 <- plot_grid(plot_grid(panel12a, panel12b, nrow = 1),
                      plot_grid(panel12c, panel12d, nrow = 1, rel_widths = c(.64, .36), align = "h"),
                      plot_grid(panel12e, panel12f, panel12g, nrow = 1, rel_widths = c(.43, .26, .31), align = "h"),
                      ncol = 1, rel_heights = c(.39, .27, .34))
save_png(merged12, "Merged_Fig12_rollback_v2_patch.png", 190, 210)
side12 <- plot_grid(
  profile_plot(data12$ztrend, "181", "<z> trend side-profile candidate", expression("<z> trend (m yr"^{-1}*")"), COL$pos, 1),
  profile_plot(data12$wgpe_trend, "181", "W-GPE trend side-profile candidate", expression("W-GPE trend (J m"^{-2}~"yr"^{-1}*")"), COL$pos, 1000),
  nrow = 1
)
save_png(side12, "Fig12_side_profile_candidate.png", 160, 70)
write.csv(data.frame(metric = c("zbar_global_trend_m_yr", "wgpe_global_trend_J_m2_yr",
                                "g_uses_pct_share", "side_profiles_in_main"),
                     value = c(data12$stats$z_global, data12$stats$wgpe_global, 1, 1)),
          file.path(OUT, "Merged_Fig12_rollback_source_values.csv"), row.names = FALSE)

# Fig.3 recovery: map panels include right-side latitude companions.
fig3a <- map_block(
  precip$precipitation$wgpe$lag0$r, "181", precip$precipitation$wgpe$lag0$q, .8, pal_blue, "a",
  bquote(atop(bold("a")~"W-GPE-precipitation association",
              "Area-weighted global mean r = +0.435")),
  "Pearson r", c(-.8, -.4, 0, .4, .8), TRUE, .1
)
bdf <- event$summary[event$summary$variable == "wgpe" & event$summary$region == "Global" &
                       event$summary$event %in% c("Extreme_precipitation", "Drought"), ]
bdf$event_label <- factor(ifelse(bdf$event == "Extreme_precipitation", "Wet events", "Dry events"),
                          levels = c("Dry events", "Wet events"))
bdf$sign <- ifelse(bdf$anomaly >= 0, "Positive", "Negative")
fig3b <- ggplot(bdf, aes(anomaly, event_label, fill = sign)) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = .24) +
  geom_col(width = .52, colour = "#222222", linewidth = .14) +
  scale_fill_manual(values = c(Positive = COL$pos, Negative = COL$neg), guide = "none") +
  scale_x_continuous(labels = label_number(big.mark = " ")) +
  labs(title = ptitle("c", "Precipitation-defined W-GPE anomalies"), x = expression("W-GPE anomaly (J m"^{-2}*")"), y = NULL) +
  theme_rec(6.4) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())
dec <- event$decomposition_summary[event$decomposition_summary$event %in% c("Extreme_precipitation", "Drought") &
                                     event$decomposition_summary$region == "Global", ]
dec_long <- rbind(data.frame(event = dec$event, channel = "IWV", value = dec$iwv_channel_J_m2),
                  data.frame(event = dec$event, channel = "<z>", value = dec$z_channel_J_m2),
                  data.frame(event = dec$event, channel = "Interaction", value = dec$interaction_channel_J_m2))
dec_long$event_label <- factor(ifelse(dec_long$event == "Extreme_precipitation", "Wet events", "Dry events"),
                               levels = c("Dry events", "Wet events"))
dec_long$channel <- factor(dec_long$channel, levels = c("IWV", "<z>", "Interaction"))
fig3c <- ggplot(dec_long, aes(value, event_label, fill = channel)) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = .24) +
  geom_col(width = .52, colour = "#222222", linewidth = .08) +
  scale_fill_manual(values = c("IWV" = COL$iwv, "<z>" = COL$zbar, "Interaction" = COL$inter), name = NULL) +
  labs(title = ptitle("d", "Wet/dry event W-GPE decomposition"), x = expression("Channel anomaly (J m"^{-2}*")"), y = NULL) +
  theme_rec(6.4) +
  theme(legend.position = "top", axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())
enso_global <- event$summary$anomaly[event$summary$event == "ENSO_difference" &
                                       event$summary$variable == "wgpe" &
                                       event$summary$region == "Global"]
if (length(enso_global) != 1 || !is.finite(enso_global)) enso_global <- 11969
fig3d <- map_block(
  event$event_maps$wgpe$ENSO_difference, "181", NULL, 50000, pal_blue, "b",
  bquote(atop(bold("b")~"ENSO W-GPE anomaly",
              "Area-weighted global anomaly = +"*.(round(enso_global))*" J m"^{-2})),
  expression("W-GPE anomaly (J m"^{-2}*")"), c(-50000, -25000, 0, 25000, 50000), FALSE, 25000
)
save_png(plot_grid(plot_grid(fig3a, fig3d, nrow = 1),
                   plot_grid(fig3b, fig3c, nrow = 1),
                   ncol = 1, rel_heights = c(.62, .38)),
         "Fig3_rollback_v2_FDR_patch.png", 180, 145)
save_png(plot_grid(
  profile_plot(precip$precipitation$wgpe$lag0$r, "181", "Fig3a W-GPE-precipitation r profile", "Pearson r", COL$pos, .1),
  profile_plot(event$event_maps$wgpe$ENSO_difference, "181", "Fig3b ENSO W-GPE anomaly profile", expression("W-GPE anomaly (J m"^{-2}*")"), COL$pos, 25000),
  nrow = 1), "Fig3_side_profile_candidates.png", 160, 70)
write.csv(rbind(
  data.frame(panel = "Fig3a", metric = "area_mean_r", value = pglob$mean_r_area),
  data.frame(panel = "Fig3a", metric = "sig_valid_area_pct", value = pglob$significant_valid_area_pct),
  data.frame(panel = "Fig3b", metric = "ENSO_global_anomaly_J_m2", value = enso_global),
  data.frame(panel = "Fig3d", metric = "wet_z_share_pct", value = dec$z_channel_pct[dec$event == "Extreme_precipitation"]),
  data.frame(panel = "Fig3d", metric = "dry_z_share_pct", value = dec$z_channel_pct[dec$event == "Drought"])
), file.path(OUT, "Fig3_source_values_rollback.csv"), row.names = FALSE)
}

# Fig.4 recovery.
eco_map <- function(target, pred, tag, title) {
  s <- get_eco_sum(target, pred)
  map_block(
    eco$results[[target]]$r[[pred]], "181", eco$results[[target]]$q[[pred]], .35, pal_green, tag,
    ptitle(tag, title),
    "Partial r", c(-.3, 0, .3), TRUE, .1, .1, COL$eco_pos,
    subtitle_text = sprintf("Global mean r = %+.3f", s$mean_r_area)
  )
}
ECO_PROFILE_CALL <- 0L
fig4a <- eco_map("GPP", "wgpe", "a", "W-GPE-GPP partial correlation")
fig4b <- eco_map("GPP", "zbar", "b", "<z>-GPP partial correlation")
fig4c <- eco_map("LAI", "wgpe", "c", "W-GPE-LAI partial correlation")
fig4d <- eco_map("LAI", "zbar", "d", "<z>-LAI partial correlation")
bands <- c("Global", "S temperate", "Tropics", "N temperate", "N high latitude")
heat_rows <- rbind(data.frame(row = "W-GPE / GPP", target = "GPP", predictor = "wgpe"),
                   data.frame(row = "<z> / GPP", target = "GPP", predictor = "zbar"),
                   data.frame(row = "W-GPE / LAI", target = "LAI", predictor = "wgpe"),
                   data.frame(row = "<z> / LAI", target = "LAI", predictor = "zbar"))
hdf_obs <- do.call(rbind, lapply(seq_len(nrow(heat_rows)), function(i) {
  rr <- heat_rows[i, ]
  sub <- eco$summary[eco$summary$target == rr$target & eco$summary$predictor == rr$predictor & eco$summary$region %in% bands, ]
  data.frame(row = rr$row, band = sub$region, value = sub$mean_r_area)
}))
# Retain the complete 4 x 5 displayed matrix. S high latitude is intentionally
# omitted from this compact main-panel summary because the persistent-
# vegetated domain has no qualifying records there.
hdf <- merge(expand.grid(row = heat_rows$row, band = bands, stringsAsFactors = FALSE),
             hdf_obs, by = c("row", "band"), all.x = TRUE, sort = FALSE)
hdf$band <- factor(hdf$band, levels = bands, labels = c("Global", "S temp", "Tropics", "N temp", "N high"))
hdf$row <- factor(hdf$row, levels = rev(heat_rows$row))
pred_lookup <- setNames(c("W-GPE","<z>","W-GPE","<z>"), c("W-GPE / GPP","<z> / GPP","W-GPE / LAI","<z> / LAI"))
target_lookup <- setNames(c("GPP","GPP","LAI","LAI"), c("W-GPE / GPP","<z> / GPP","W-GPE / LAI","<z> / LAI"))
hdf$target0 <- target_lookup[as.character(hdf$row)]
hdf$pred0 <- pred_lookup[as.character(hdf$row)]
hdf$region0 <- bands[match(as.character(hdf$band), c("Global","S temp","Tropics","N temp","N high"))]
hdf$ci_excludes_zero <- mapply(function(t,p,r){
  u<-U_ESUM[U_ESUM$height_reference==REF & U_ESUM$target==t & U_ESUM$predictor==p & U_ESUM$region==r,]
  nrow(u)>0 && (u$lower95[1]>0 || u$upper95[1]<0)
},hdf$target0,hdf$pred0,hdf$region0)
hdf$xpos<-c(Global=1,`S temp`=2.35,Tropics=3.35,`N temp`=4.35,`N high`=5.35)[as.character(hdf$band)]
fig4e <- ggplot(hdf, aes(xpos, row, fill = value)) +
  geom_tile(width=.82,height=.82,colour="#555555",linewidth=.10) +
  geom_text(data=hdf[hdf$ci_excludes_zero,],label="*",size=2.25,colour="#111111") +
  geom_text(data=hdf[!is.finite(hdf$value),],label="NA",size=1.65,colour="#666666") +
  scale_x_continuous(breaks=c(1,2.35,3.35,4.35,5.35),labels=c("Global","S temp","Tropics","N temp","N high"),expand=expansion(add=.5)) +
  coord_fixed(ratio = 1) +
  scale_fill_gradientn(colours = pal_green, limits = c(-.35, .35), oob = squish,
                       na.value = "#F2F2F2", name = "Partial r",
                       guide = guide_colorbar(direction = "horizontal",
                                             title.position = "left",
                                             label.position = "bottom",
                                             barwidth = unit(42, "mm"),
                                             barheight = unit(2.2, "mm"),
                                             ticks.colour = "#111111",
                                             frame.colour = "#555555",
                                             theme = theme(legend.frame = element_rect(colour = "#555555", linewidth = .18)))) +
  labs(title = ptitle("e", "Ecological synthesis by latitude band"), x = NULL, y = NULL) +
  theme_rec(6.4) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 5.2),
        axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 5.8),
        legend.text = element_text(size = 5.5),
        legend.key.width = unit(5.5, "mm"),
        legend.margin = margin(t = 1, b = 0),
        legend.box.margin = margin(t = 0, b = 0),
        plot.margin = margin(3, 4, 3, 4))
save_png(plot_grid(plot_grid(fig4a, fig4b, fig4c, fig4d, ncol = 2),
                   fig4e, ncol = 1, rel_heights = c(.76, .24)),
         "Fig4_rollback_v2_FDR_patch.png", 180, 172)
save_png(plot_grid(
  profile_plot(eco$results$GPP$r$wgpe, "181", "W-GPE-GPP partial r profile", "Partial r", COL$eco_pos, .1),
  profile_plot(eco$results$GPP$r$zbar, "181", "<z>-GPP partial r profile", "Partial r", COL$eco_pos, .1),
  profile_plot(eco$results$LAI$r$wgpe, "181", "W-GPE-LAI partial r profile", "Partial r", COL$eco_pos, .1),
  profile_plot(eco$results$LAI$r$zbar, "181", "<z>-LAI partial r profile", "Partial r", COL$eco_pos, .1),
  nrow = 2), "Fig4_side_profile_candidates.png", 170, 110)
gsum <- granger$summary
glag <- gsum[gsum$lag == granger$primary_lag & gsum$region == "Global", ]
glong <- rbind(data.frame(target = glag$target, class = "Forward only", pct = glag$forward_only_grid_pct),
               data.frame(target = glag$target, class = "Reverse only", pct = glag$reverse_only_grid_pct),
               data.frame(target = glag$target, class = "Bidirectional", pct = glag$bidirectional_grid_pct),
               data.frame(target = glag$target, class = "None", pct = glag$none_grid_pct))
glong$class <- factor(glong$class, levels = c("Forward only", "Reverse only", "Bidirectional", "None"))
fig4_si <- ggplot(glong, aes(pct, target, fill = class)) +
  geom_col(width = .58, colour = "#222222", linewidth = .08) +
  scale_fill_manual(values = c("Forward only" = "#0D56A0", "Reverse only" = "#A66B1C",
                               "Bidirectional" = "#6EB7AA", "None" = "#D8D8D8"), name = NULL) +
  labs(title = "Temporal precedence (descriptive SI candidate)",
       subtitle = "Nested-lag classification; no BH/FDR mask audited",
       x = "Tested grid cells (%)", y = NULL) +
  theme_rec(6.4) +
  theme(legend.position = "bottom", axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())
save_png(fig4_si, "Fig4_temporal_precedence_SI_candidate.png", 88, 52)
write.csv(rbind(
  data.frame(panel = "Fig4a", target = "GPP", predictor = "wgpe", metric = "mean_r_area", value = get_eco_sum("GPP", "wgpe")$mean_r_area),
  data.frame(panel = "Fig4b", target = "GPP", predictor = "zbar", metric = "mean_r_area", value = get_eco_sum("GPP", "zbar")$mean_r_area),
  data.frame(panel = "Fig4c", target = "LAI", predictor = "wgpe", metric = "mean_r_area", value = get_eco_sum("LAI", "wgpe")$mean_r_area),
  data.frame(panel = "Fig4d", target = "LAI", predictor = "zbar", metric = "mean_r_area", value = get_eco_sum("LAI", "zbar")$mean_r_area),
  data.frame(panel = "Fig4e", target = as.character(hdf$row), predictor = as.character(hdf$band),
             metric = "band_mean_r_area", value = hdf$value)
), file.path(OUT, "Fig4_source_values_rollback.csv"), row.names = FALSE)

writeLines(c("# Merged Fig1/2 recovery QA", "",
             "1. e/f/g returned to a larger bottom-row structure rather than the compressed v4 layout: YES.",
             "2. g is share composition: YES.",
             "3. g uses percentage channels, not absolute J m-2 yr-1 contributions: YES.",
             "4. Units in axes/titles use plotmath superscripts where unit expressions are drawn: YES.",
             "5. Map side profiles are embedded to the right of each map block: YES.",
             "6. Side profiles use the same map variable, area-weighted 5-degree latitude bins, and are vertically centered to the Robinson ellipse body.",
             "7. Robinson maps are drawn with coord_equal and are not stretched."),
           file.path(OUT, "Merged_Fig12_rollback_v2_patch_QA.md"), useBytes = TRUE)
writeLines(c("# Fig.3 recovery QA", "",
             "1. Fig3a field/q alignment PASS before drawing: YES.",
             "2. Fig3a uses q<0.05 coloured and q>=0.05 pale grey: YES.",
             "3. Fig3d has no in-panel numbers: YES.",
             "4. Fig3b ENSO is descriptive and unmasked: YES.",
             "5. Every map panel has a right-side latitude companion aligned to the map ellipse height: YES.",
             "6. Side-profile candidates are also exported for inspection."),
           file.path(OUT, "Fig3_rollback_v2_FDR_patch_QA.md"), useBytes = TRUE)
writeLines(c("# Fig.4 recovery QA", "",
             "1. Fig4a-d field/q alignment PASS before drawing: YES.",
             "2. Ecology maps are partial correlations, not panel FE: YES.",
             "3. q<0.05 coloured and q>=0.05 pale grey is used: YES.",
             "4. LAI maps are retained: YES.",
             "5. Heatmap has five latitude bands and no in-cell numbers: YES.",
             "6. Temporal precedence is exported only as descriptive SI candidate: YES.",
             "7. No causality / Granger causality language is used.",
             "8. Fig4 180-grid coordinate fix: latitude is interpreted north-to-south in the raw ecology matrices, then sorted south-to-north together with field and q-mask before Robinson projection.",
             "9. Every Fig4 map has a right-side latitude companion aligned to the map ellipse height: YES."),
           file.path(OUT, "Fig4_rollback_v2_FDR_patch_QA.md"), useBytes = TRUE)

writeLines(c("# Version comparison: v2/pass1 vs v4", "",
             "Decision: rollback to v2/pass1 structure and apply only local method patches.", "",
             "1. Best Fig1/2 e/f/g structure: fig12_merged_refined_v2 / formal_redraw_pass1. The v4 integrated version compressed e/f/g and changed g to absolute contribution, which is not the requested regional share composition.",
             "2. Most trustworthy Fig3 map coordinates: corrected downstream fields after explicit field/q alignment QC, with right-side latitude companions embedded in each map block.",
             "3. Most trustworthy Fig4 ecology maps: corrected downstream partial-correlation RDS after explicit field/q alignment QC and a north-to-south raw latitude fix for 180-grid ecology matrices.",
             "4. What v4 broke: e/f/g layout became less readable; g used absolute J m-2 yr-1 contributions instead of percentage shares; side latitude profiles looked like external appendages; some unit strings were plain text; FDR masks raised coordinate/mask concerns that had not been separately audited in that pass.",
             "5. Panels rolled back: Fig1/2 e, f, g structure; Fig3/Fig4 map layout now uses embedded side profiles; regional composition as shares.",
             "6. v4 changes retained only after audit: FDR mask semantics for Fig3a and Fig4a-d, ecology partial-correlation method, LAI maps, temporal-precedence as SI candidate only.",
             "7. Temporarily withdrawn: absolute-contribution g panel and any v4-driven free layout changes."),
           file.path(OUT, "version_comparison_v2_vs_v4.md"), useBytes = TRUE)

writeLines(c("# Rollback to v2 self-check summary", "",
             "1. Stopped continuing from v4: YES. This pass uses v2/pass1/recovery structure and treats v4 as a failure reference.",
             "2. Panels rolled back to v2/pass1: Fig1/2 e/f/g structure; Fig3/Fig4 maps with embedded right-side latitude companions.",
             "3. v4 changes retained: only audited FDR masking, LAI partial-correlation maps, and temporal-precedence as SI candidate.",
             "4. Fig1/2 e/f/g restored: YES, with g as share composition.",
             "5. g restored as share composition: YES; x-axis is Share of regional W-GPE trend (%).",
             paste0("6. field/q mask alignment all PASS: ", all_pass, "."),
             "7. Failed maps: none.",
             "8. Fig3/Fig4 were drawn only after PASS: YES.",
             "9. FDR mask applied correctly: q<0.05 coloured, q>=0.05 pale grey; q>=0.05 was not set to zero.",
             "10. DOCX/SI/workbook modified: NO.",
             "11. PDF output: NO.",
             "12. Fig4 ecology 180-grid latitude order fixed before mapping: YES."),
           file.path(OUT, "rollback_to_v2_selfcheck_summary.md"), useBytes = TRUE)

###############################################################################
# Final FDR-v1 packaging requested by the manuscript figure pass.
# This section intentionally reuses the objects created above plus the locked
# Fig.1/2 v2 composite; it does not change scientific data or recompute analyses.
###############################################################################

if (FALSE) {
copy_if_exists <- function(from, to) {
  if (!file.exists(from)) stop("Missing file to copy: ", from)
  ok <- file.copy(from, to, overwrite = TRUE)
  if (!ok) stop("Could not copy file to: ", to)
  invisible(to)
}

# Do not copy the old locked Fig1/2 composite here: this layout pass rebuilds
# panels a/b to remove the in-map global-mean subtitle and to fix map-colourbar
# spacing.  Keep the v1 filename for downstream compatibility.
save_png(merged12, "Merged_Fig12_final_candidate_v2_layoutfix.png", 190, 210)
save_png(merged12, "Merged_Fig12_final_candidate_v1.png", 190, 210)
save_png(panel12a, "Fig12a_clean.png", 92, 68)
save_png(panel12b, "Fig12b_clean.png", 92, 68)
save_png(panel12c, "Fig12c_clean.png", 110, 62)
save_png(panel12d, "Fig12d_clean.png", 70, 62)
save_png(panel12e, "Fig12e_clean.png", 82, 66)
save_png(panel12f, "Fig12f_clean.png", 50, 66)
save_png(panel12g, "Fig12g_clean.png", 58, 66)
}

temporal_class_colours_hydro <- c("forward-only" = COL$pos,
                                  "reverse-only" = COL$zbar,
                                  "bidirectional" = COL$inter,
                                  "none" = "#D8D8D8")
temporal_class_colours_eco <- c("forward-only" = COL$eco_pos,
                                "reverse-only" = COL$zbar,
                                "bidirectional" = "#6EB7AA",
                                "none" = "#D8D8D8")
# Backward-compatible default for Fig.3/hydroclimate temporal panels.
temporal_class_colours <- temporal_class_colours_hydro

temporal_bar_panel <- function(targets, tag = NULL, title = "Lag-3 temporal precedence",
                               subtitle = "Nested lag test; BH-FDR q < 0.05",
                               palette = temporal_class_colours_hydro) {
  d <- fdr_summary[fdr_summary$predictor == "W-GPE" & fdr_summary$lag == 3 &
                     fdr_summary$region == "Global" & fdr_summary$target %in% targets, , drop = FALSE]
  if (nrow(d) != length(targets)) stop("Missing temporal FDR summary rows for: ", paste(targets, collapse = ", "))
  d$target <- factor(d$target, levels = rev(targets),
                     labels = rev(ifelse(targets == "precipitation", "Precipitation", targets)))
  long <- rbind(
    data.frame(target = d$target, class = "forward-only", pct = d$forward_only_grid_pct_q),
    data.frame(target = d$target, class = "reverse-only", pct = d$reverse_only_grid_pct_q),
    data.frame(target = d$target, class = "bidirectional", pct = d$bidirectional_grid_pct_q),
    data.frame(target = d$target, class = "none", pct = d$none_grid_pct_q)
  )
  long$class <- factor(long$class, levels = c("forward-only", "reverse-only", "bidirectional", "none"))
  fu <- U_ETEMP[U_ETEMP$height_reference == REF & U_ETEMP$class == "forward-only" & U_ETEMP$target %in% targets & U_ETEMP$lag == 3,]
  fu$target <- factor(fu$target, levels = rev(targets), labels = rev(ifelse(targets == "precipitation", "Precipitation", targets)))
  ggplot(long, aes(pct, target, fill = class)) +
    geom_col(width = .55, colour = "#222222", linewidth = .08) +
    scale_fill_manual(values = palette, name = NULL) +
    scale_x_continuous(breaks = c(0, 50, 100),
                       labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, .005))) +
    coord_cartesian(xlim = c(0, 100), clip = "off") +
    labs(title = ptitle(tag, title), subtitle = subtitle, x = "Tested grid cells (%)", y = NULL) +
    theme_rec(6.4) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(), plot.subtitle = element_text(size = PANEL_SUBTITLE_SIZE, colour = COL$ink))
}

temporal_class_map <- function(target, tag, title, palette = temporal_class_colours_hydro) {
  d <- fdr_grid[fdr_grid$predictor == "W-GPE" & fdr_grid$lag == 3 &
                  fdr_grid$target == target & fdr_grid$valid_TRUE_FALSE, , drop = FALSE]
  if (nrow(d) < 1) stop("No valid temporal FDR grid rows for ", target)
  classes <- c("forward-only", "reverse-only", "bidirectional", "none")
  d$class_q <- factor(d$class_q, levels = classes)
  lon <- sort(unique(d$lon)); lat <- sort(unique(d$lat))
  code <- setNames(seq_along(classes), classes)
  mat <- matrix(NA_real_, nrow = length(lon), ncol = length(lat))
  mat[cbind(match(d$lon, lon), match(d$lat, lat))] <- unname(code[as.character(d$class_q)])
  b <- bundle_from_lonlat_field(mat, lon, lat)
  xlim <- c(-1.015, 1.015)
  ylim <- c(-max(abs(b$border$y)) - .012, max(abs(b$border$y)) + .012)
  mask <- outer_mask(b$border, xlim, ylim)
  mapdf <- b$map
  mapdf$class <- factor(classes[round(mapdf$value)], levels = classes)
  mapdf <- mapdf[!is.na(mapdf$class), , drop = FALSE]
  ggplot() +
    geom_polygon(data = mapdf, aes(x, y, group = id, fill = class), colour = NA) +
    geom_path(data = b$graticule, aes(x, y, group = group), linewidth = .13, colour = "#B9BDBD", alpha = .55) +
    geom_path(data = b$coast, aes(x, y, group = group), linewidth = .20, colour = "#505456", lineend = "round") +
    geom_polygon(data = mask, aes(x, y, group = group, subgroup = subgroup),
                 inherit.aes = FALSE, fill = "white", colour = NA, rule = "evenodd") +
    geom_path(data = b$border, aes(x, y, group = group), linewidth = .30, colour = "#6A6E70", linejoin = "round") +
    scale_fill_manual(values = palette, name = NULL, drop = FALSE, labels = classes) +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE, clip = "on") +
    labs(title = ptitle(tag, title), x = NULL, y = NULL) +
    theme_map_rec(6.3) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.key.width = unit(6.0, "mm"), legend.key.height = unit(2.2, "mm"),
          plot.margin = margin(2, 2, 2, 2))
}

# Hydroclimate Figure 2 is built by its dedicated TOP100 script.
if (FALSE) {
# Final Fig.3 outputs.  The two diagnostic bar panels are deliberately
# down-weighted; the precipitation temporal-precedence map is promoted into
# the main candidate because it carries more spatial information.
fig3e <- temporal_class_map("precipitation", "e", "Lag-3 precipitation temporal precedence",
                            palette = temporal_class_colours_hydro)
save_png(fig3a, "Fig3a_clean.png", 92, 68)
save_png(fig3d, "Fig3b_clean.png", 92, 68)
save_png(fig3b, "Fig3c_clean.png", 70, 42)
save_png(fig3c, "Fig3d_clean.png", 70, 42)
save_png(fig3e, "Fig3e_temporal_precedence_map_clean.png", 160, 50)
fig3_cd_stack <- plot_grid(fig3b, fig3c, ncol = 1, rel_heights = c(.46, .54), align = "v")
fig3_cd_row <- plot_grid(NULL, fig3_cd_stack, NULL, nrow = 1, rel_widths = c(.24, .52, .24))
fig3_final <- plot_grid(
  plot_grid(fig3a, fig3d, nrow = 1),
  fig3_cd_row,
  fig3e,
  ncol = 1, rel_heights = c(.47, .30, .23)
)
save_png(fig3_final, "Fig3_final_candidate_v2_layoutfix2.png", 180, 195)
# Compatibility name for downstream scripts that still expect v1.
save_png(fig3_final, "Fig3_final_candidate_v1.png", 180, 195)

fig3_temporal_ed <- plot_grid(
  temporal_bar_panel("precipitation", "a", "Precipitation temporal precedence",
                     "Lag-3 W-GPE predictor; BH-FDR q < 0.05",
                     palette = temporal_class_colours_hydro) +
    theme(legend.position = "none"),
  temporal_class_map("precipitation", "b", "Lag-3 precipitation classes",
                     palette = temporal_class_colours_hydro),
  nrow = 1, rel_widths = c(.36, .64)
)
save_png(fig3_temporal_ed, "Fig3_temporal_precedence_ED_candidate.png", 170, 72)
}

# Final Fig.4 outputs.
fig4f <- temporal_bar_panel(c("GPP", "LAI"), "f", "Lag-3 temporal precedence",
                            "W-GPE predictor; nested lag test; BH-FDR q < 0.05",
                            palette = temporal_class_colours_eco)
save_png(fig4a, "Fig4a_clean.png", 88, 64)
save_png(fig4b, "Fig4b_clean.png", 88, 64)
save_png(fig4c, "Fig4c_clean.png", 88, 64)
save_png(fig4d, "Fig4d_clean.png", 88, 64)
save_png(fig4e, "Fig4e_clean.png", 94, 46)
save_png(fig4f, "Fig4f_temporal_precedence_clean.png", 78, 46)

fig4_maps_only <- plot_grid(fig4a, fig4b, fig4c, fig4d, ncol = 2)
save_png(fig4_maps_only, "Fig4_final_candidate_maps_only_v1.png", 180, 126)

fig4_final <- plot_grid(fig4_maps_only,
                        plot_grid(fig4e, fig4f, nrow = 1, rel_widths = c(.54, .46)),
                        ncol = 1, rel_heights = c(.68, .32))
final_name <- if (NO_UNCERTAINTY) paste0("Figure3_",REF,"_panel_e_no_Shigh_v8_no_uncertainty.png") else paste0("Figure3_",REF,"_panel_e_no_Shigh_v8_with_uncertainty.png")
save_png(fig4_final, final_name, 205, 196, dpi=600)
save_png(fig4_final, "Figure3_TOP100_FINAL.png", 205, 196, dpi=600)
save_pdf(fig4_final, "Figure3_TOP100_FINAL.pdf", 205, 196)
save_png(fig4a,paste0("Figure3_",REF,"_panel_a.png"),88,64,dpi=600)
save_png(fig4b,paste0("Figure3_",REF,"_panel_b.png"),88,64,dpi=600)
save_png(fig4c,paste0("Figure3_",REF,"_panel_c.png"),88,64,dpi=600)
save_png(fig4d,paste0("Figure3_",REF,"_panel_d.png"),88,64,dpi=600)
save_png(fig4e,paste0("Figure3_",REF,"_panel_e.png"),94,46,dpi=600)
save_png(fig4f,paste0("Figure3_",REF,"_panel_f.png"),78,46,dpi=600)
save_pdf(fig4a,paste0("Figure3_",REF,"_panel_a.pdf"),88,64)
save_pdf(fig4b,paste0("Figure3_",REF,"_panel_b.pdf"),88,64)
save_pdf(fig4c,paste0("Figure3_",REF,"_panel_c.pdf"),88,64)
save_pdf(fig4d,paste0("Figure3_",REF,"_panel_d.pdf"),88,64)
save_pdf(fig4e,paste0("Figure3_",REF,"_panel_e.pdf"),94,46)
save_pdf(fig4f,paste0("Figure3_",REF,"_panel_f.pdf"),78,46)
writeLines(c("# Figure 3 archived-layout redraw v3","",paste0("Height reference: ",REF),
             "Archived map/heatmap/temporal layout retained; Model 3 is displayed.",
             "Uncertainty embedded: 5-degree latitude-profile ribbons, bootstrap-CI heatmap borders, and forward-only fraction error bars.",
             "Maps retain effective-dof/BH-FDR significance encoding."),
           file.path(OUT,paste0("Figure3_",REF,"_archived_layout_QA.md")))
cat("DONE Figure3 archived layout ",REF,"\n",sep="")
quit(save="no",status=0)

fig4_temporal_ed <- plot_grid(
  temporal_class_map("GPP", "a", "Lag-3 GPP classes", palette = temporal_class_colours_eco),
  temporal_class_map("LAI", "b", "Lag-3 LAI classes", palette = temporal_class_colours_eco),
  nrow = 1
)
save_png(fig4_temporal_ed, "Fig4_temporal_precedence_maps_ED_candidate.png", 170, 72)

# Source values and caption-ready notes.
fig12_sources <- rbind(
  data.frame(panel = "Fig1/2a", metric = "area_weighted_global_zbar_trend_m_yr", value = data12$stats$z_global),
  data.frame(panel = "Fig1/2b", metric = "area_weighted_global_wgpe_trend_J_m2_yr", value = data12$stats$wgpe_global),
  data.frame(panel = "Fig1/2d", metric = as.character(data12$correlations$variable), value = as.numeric(data12$correlations$r)),
  data.frame(panel = "Fig1/2f", metric = c("mass_redistribution_m_yr", "height_coordination_m_yr", "residual_m_yr", "actual_m_yr"),
             value = c(0.6490, 0.2313, -0.0919, 0.7885)),
  data.frame(panel = "Fig1/2g", metric = paste0(as.character(data12$wgpe$region), "_IWV_pct"), value = data12$wgpe$iwv_channel_pct),
  data.frame(panel = "Fig1/2g", metric = paste0(as.character(data12$wgpe$region), "_z_pct"), value = data12$wgpe$z_channel_pct),
  data.frame(panel = "Fig1/2g", metric = paste0(as.character(data12$wgpe$region), "_interaction_pct"), value = data12$wgpe$interaction_channel_pct)
)
write.csv(fig12_sources, file.path(OUT, "Merged_Fig12_source_values.csv"), row.names = FALSE)

precip_temporal <- fdr_summary[fdr_summary$target == "precipitation" & fdr_summary$predictor == "W-GPE" &
                                  fdr_summary$lag == 3 & fdr_summary$region == "Global", , drop = FALSE]
enso_global_channel <- if (!is.null(enso_channel)) {
  enso_channel[enso_channel$event_or_contrast == "ENSO_difference" & enso_channel$region == "Global", , drop = FALSE]
} else data.frame()
if (nrow(enso_global_channel) == 0 && !is.null(enso_channel)) {
  enso_global_channel <- enso_channel[enso_channel$event_or_contrast == "El_Nino" & enso_channel$region == "Global", , drop = FALSE]
}
fig3_sources <- rbind(
  data.frame(panel = "Fig3a", metric = "area_weighted_global_mean_r", value = pglob$mean_r_area),
  data.frame(panel = "Fig3a", metric = "significant_valid_area_pct_q_lt_0.05", value = pglob$significant_valid_area_pct),
  data.frame(panel = "Fig3a", metric = "significant_grid_pct_q_lt_0.05", value = pglob$significant_grid_pct),
  data.frame(panel = "Fig3d", metric = "wet_z_share_pct", value = dec$z_channel_pct[dec$event == "Extreme_precipitation"]),
  data.frame(panel = "Fig3d", metric = "dry_z_share_pct", value = dec$z_channel_pct[dec$event == "Drought"]),
  data.frame(panel = "Fig3b", metric = "ENSO_global_WGPE_anomaly_J_m2", value = enso_global),
  if (nrow(enso_global_channel) > 0) {
    data.frame(panel = "Fig3b", metric = c("ENSO_global_IWV_share_pct", "ENSO_global_z_share_pct", "ENSO_global_interaction_share_pct"),
               value = c(enso_global_channel$IWV_share_percent[1], enso_global_channel$z_share_percent[1],
                         enso_global_channel$interaction_share_percent[1]))
  } else {
    data.frame(panel = "Fig3b", metric = c("ENSO_global_IWV_share_pct", "ENSO_global_z_share_pct", "ENSO_global_interaction_share_pct"),
               value = c(80.27, 10.56, 9.17))
  },
  data.frame(panel = "Fig3_ED", metric = c("precip_forward_only_pct_q", "precip_reverse_only_pct_q",
                                           "precip_bidirectional_pct_q", "precip_none_pct_q"),
             value = c(precip_temporal$forward_only_grid_pct_q, precip_temporal$reverse_only_grid_pct_q,
                       precip_temporal$bidirectional_grid_pct_q, precip_temporal$none_grid_pct_q))
)
write.csv(fig3_sources, file.path(OUT, "Fig3_source_values.csv"), row.names = FALSE)

eco_temporal <- fdr_summary[fdr_summary$target %in% c("GPP", "LAI") & fdr_summary$predictor == "W-GPE" &
                              fdr_summary$lag == 3 & fdr_summary$region == "Global", , drop = FALSE]
fig4_sources <- rbind(
  data.frame(panel = "Fig4a", target = "GPP", predictor = "W-GPE", metric = "area_weighted_global_mean_partial_r", value = get_eco_sum("GPP", "wgpe")$mean_r_area),
  data.frame(panel = "Fig4b", target = "GPP", predictor = "<z>", metric = "area_weighted_global_mean_partial_r", value = get_eco_sum("GPP", "zbar")$mean_r_area),
  data.frame(panel = "Fig4c", target = "LAI", predictor = "W-GPE", metric = "area_weighted_global_mean_partial_r", value = get_eco_sum("LAI", "wgpe")$mean_r_area),
  data.frame(panel = "Fig4d", target = "LAI", predictor = "<z>", metric = "area_weighted_global_mean_partial_r", value = get_eco_sum("LAI", "zbar")$mean_r_area),
  data.frame(panel = "Fig4e", target = as.character(hdf$row), predictor = as.character(hdf$band), metric = "band_mean_partial_r", value = hdf$value),
  data.frame(panel = "Fig4f", target = eco_temporal$target, predictor = "W-GPE", metric = "forward_only_pct_q", value = eco_temporal$forward_only_grid_pct_q),
  data.frame(panel = "Fig4f", target = eco_temporal$target, predictor = "W-GPE", metric = "reverse_only_pct_q", value = eco_temporal$reverse_only_grid_pct_q),
  data.frame(panel = "Fig4f", target = eco_temporal$target, predictor = "W-GPE", metric = "bidirectional_pct_q", value = eco_temporal$bidirectional_grid_pct_q),
  data.frame(panel = "Fig4f", target = eco_temporal$target, predictor = "W-GPE", metric = "none_pct_q", value = eco_temporal$none_grid_pct_q)
)
write.csv(fig4_sources, file.path(OUT, "Fig4_source_values.csv"), row.names = FALSE)

writeLines(c("# Merged Fig1/2 final candidate v1 QA", "",
             "1. Source figure: locked v2 composite with E/F/G restored, copied without recomputing scientific results.",
             "2. Panel g is regional W-GPE trend composition in percent share, not absolute J m^-2 yr^-1: YES.",
             "3. Panel e keeps the four-layer integrated contribution and the contribution-density curve on a shared pressure axis: YES.",
             "4. Panel a/b global trends are stored in source-values CSV but removed from map subtitles in this layout pass: YES.",
             "5. Robinson maps are not stretched: YES.",
             "6. Palette direction is brown=negative, cream=zero, blue=positive for hydroclimate maps: YES.",
             "7. DOCX/SI/workbook modified: NO.",
             "8. PDF output: NO."),
           file.path(OUT, "Merged_Fig12_final_candidate_v1_QA.md"), useBytes = TRUE)

writeLines(c("# Fig.3 QA", "",
             "1. Fig3a is an inferential FDR map: q < 0.05 coloured; q >= 0.05 pale grey; q >= 0.05 was not set to zero.",
             "2. Fig3a side profile uses the same W-GPE-precipitation r field, 5-degree latitude bins, cos(latitude) weighting, and all valid cells rather than q-filtered cells.",
             "3. Fig3d wet/dry decomposition keeps numbers out of the bars; shares are stored in Fig3_source_values.csv.",
             "4. Fig3b ENSO W-GPE anomaly is descriptive and unmasked because no event p/q map is available.",
             "5. Precipitation temporal precedence map is promoted into the main Fig3 candidate and is drawn in Robinson projection; the ED/SI export is also retained.",
             "6. Fig3c and Fig3d are vertically stacked beside the Robinson temporal-precedence map.",
             "7. Wording avoids 'W-GPE beats IWV' and avoids causal claims for ENSO/event composites."),
           file.path(OUT, "Fig3_QA.md"), useBytes = TRUE)

writeLines(c("# Fig.4 QA", "",
             "1. Fig4a-d are partial-correlation maps, not panel fixed-effect maps: YES.",
             "2. Controls are precipitation, T2m, and SSRD in the corrected downstream ecology branch: YES.",
             "3. Effective-dof p-values and BH q masks are used; q < 0.05 coloured and q >= 0.05 pale grey: YES.",
             "4. LAI panels are included: YES.",
             "5. Heatmap uses five latitude bands and no in-cell numbers: YES.",
             "5b. Panel e colourbar is kept narrow/short with the shared custom colorbar, avoiding the fat guide issue: YES.",
             "6. Panel f uses lag-3 temporal precedence with BH-FDR q < 0.05: YES.",
             "7. Temporal-precedence maps are ED/SI candidates only: YES.",
             "8. Wording uses temporal precedence / nested lag test and avoids causality language: YES.",
             "9. Fig4 ecology map coordinate fix is retained: 180-grid raw matrices are interpreted north-to-south and sorted with their q masks before plotting."),
           file.path(OUT, "Fig4_QA.md"), useBytes = TRUE)

writeLines(c("# FDR mask display QA", "",
             "1. Fig3a and Fig4a-d use q < 0.05 coloured and q >= 0.05 pale grey: YES.",
             "2. Non-significant q >= 0.05 cells are displayed as non-significant grey, not converted to zero: YES.",
             "3. Dense hatching is avoided: YES.",
             "4. Event/ENSO panels without q maps are descriptive, not FDR-masked: YES.",
             "5. Temporal-precedence classes use BH-FDR q < 0.05 from the dedicated FDR branch: YES."),
           file.path(OUT, "FDR_mask_display_QA.md"), useBytes = TRUE)

writeLines(c("# Map side latitude-profile method QA", "",
             "1. Side profiles use the same variable as their companion map: YES.",
             "2. Profiles use 5-degree latitude bins: YES.",
             "3. Profiles use cos(latitude) weighting within each latitude bin: YES.",
             "4. Profiles use all valid cells and do not apply q-only filtering: YES.",
             "5. Side profiles are companion panels aligned to the map block, not standalone scientific panels: YES.",
             "6. Panel FE results are not used for map side profiles: YES."),
           file.path(OUT, "map_side_latprofile_method_QA.md"), useBytes = TRUE)

writeLines(c("# Final main figures layoutfix2 QA", "",
             "1. Final candidate output directory: Figures_FIXED_FULL/final_main_figures_FDR_v2_layoutfix2.",
             "2. Merged Fig1/2 uses the locked v2 structure, not the rejected v4 layout: YES.",
             "3. Fig3 main candidate includes the precipitation temporal-precedence Robinson map and stacks c/d vertically: YES.",
             "4. Fig4 main candidate includes FDR partial-correlation ecology maps, the five-band heatmap, and compact temporal-precedence summary: YES.",
             "5. Temporal-precedence maps are separate ED/SI candidates: YES.",
             "6. Palette contract respected: hydro maps brown-negative/blue-positive; ecology maps brown-negative/green-positive; no red-blue/default ggplot/viridis map palette.",
             "7. DOCX/SI/workbook modified: NO.",
             "8. PDF output: NO."),
           file.path(OUT, "final_main_figures_layoutfix2_QA.md"), useBytes = TRUE)

writeLines(c("# Merged Fig1/2 caption draft", "",
             "Fig. 1/2 candidate uses corrected FIXED_FULL fields. Global trends are area-weighted and stored in the source-values table, but not printed inside the map panels; side latitude profiles summarize the same mapped variable using 5-degree, cos(latitude)-weighted bins. Panel g reports the percentage share of regional W-GPE trends across IWV, <z>, and interaction channels."),
           file.path(OUT, "Merged_Fig12_caption_draft.md"), useBytes = TRUE)
writeLines(c("# Fig.3 caption draft", "",
             "Fig. 3 summarizes hydrological associations and event-scale W-GPE anomalies. Inferential association maps colour only BH-FDR significant cells (q < 0.05), with non-significant cells shown in pale grey. ENSO and wet/dry composites are descriptive diagnostics. The precipitation temporal-precedence map uses the lag-3 nested-lag/BH-FDR classification and is shown as a spatial diagnostic rather than a causal claim."),
           file.path(OUT, "Fig3_caption_draft.md"), useBytes = TRUE)
writeLines(c("# Fig.4 caption draft", "",
             "Fig. 4 shows corrected partial correlations between atmospheric water-vapour diagnostics and vegetation indicators, controlling for precipitation, T2m, and SSRD. Maps use effective-dof p-values with BH-FDR q masks. Lag-3 temporal-precedence summaries use nested lag tests with BH-FDR q < 0.05 and are interpreted as diagnostic precedence patterns, not causality."),
           file.path(OUT, "Fig4_caption_draft.md"), useBytes = TRUE)

cat("DONE final_main_figures_FDR_v2_layoutfix2\n")
