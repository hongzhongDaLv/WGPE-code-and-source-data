options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))
if (!requireNamespace("ncdf4", quietly = TRUE)) stop("ncdf4 is required")
library(ncdf4)

ROOT <- "__WGPE_PROJECT_ROOT__"
CANONICAL <- file.path(ROOT, "CANONICAL_WGPE_TOP100_v2")
OUT <- file.path(CANONICAL, "08_REPRODUCTION_OUTPUT", "Figure1_sensitivity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PLEV_DIR <- file.path(ROOT, "era5_plev_ascii")
SP_FILE <- file.path(ROOT, "output_reviewer_P0_resolution", "ERA5_surface_pressure_aligned_0_359.nc")
CORE_FILE <- file.path(ROOT, "output_TOP100_revision", "derived_data", "ERA5_TOP100", "ERA5_TOP100_monthly_core_metrics.nc")
A_GRID_FILE <- file.path(CANONICAL, "FINAL_LARGE_DATA", "exact_layer_grid.csv.gz")
G <- 9.80665

LAYERS <- data.frame(layer = c("1000-850", "850-700", "700-500", "500-300", "300-200", "200-100"),
                     p_bottom = c(1000, 850, 700, 500, 300, 200),
                     p_top = c(850, 700, 500, 300, 200, 100))

files <- sort(list.files(PLEV_DIR, pattern = "\\.nc$", full.names = TRUE))
if (!length(files)) stop("No pressure-level files")
first <- nc_open(files[1])
lon <- as.numeric(ncvar_get(first, "longitude")); lat <- as.numeric(ncvar_get(first, "latitude"))
plev_all <- as.numeric(ncvar_get(first, "pressure_level")); nc_close(first)
use <- which(plev_all <= 1000 & plev_all >= 100); plev <- plev_all[use]
if (any(diff(plev) >= 0)) stop("Pressure levels must descend")
nlon <- length(lon); nlat <- length(lat); nlay <- nrow(LAYERS)
expected_dates <- seq(as.Date("1979-01-01"), as.Date("2024-12-01"), by = "month")
nt <- length(expected_dates)
tdec <- as.numeric(format(expected_dates, "%Y")) + (as.numeric(format(expected_dates, "%m")) - .5) / 12

get_dates <- function(nc, name = "valid_time") {
  as.Date(as.POSIXct(ncvar_get(nc, name), origin = "1970-01-01", tz = "UTC"))
}
file_meta <- do.call(rbind, lapply(files, function(p) {
  n <- nc_open(p); on.exit(nc_close(n), add = TRUE); dd <- get_dates(n)
  data.frame(path = p, first = min(dd), last = max(dd), n_time = length(dd))
}))
file_meta <- file_meta[order(file_meta$first), ]

sp_nc <- nc_open(SP_FILE)
core_nc <- nc_open(CORE_FILE)
on.exit({try(nc_close(sp_nc), silent = TRUE); try(nc_close(core_nc), silent = TRUE)}, add = TRUE)
stopifnot(all.equal(as.numeric(ncvar_get(core_nc, "longitude")), lon, tolerance = 0),
          all.equal(as.numeric(ncvar_get(core_nc, "latitude")), lat, tolerance = 0))
core_dates <- as.Date(ncvar_get(core_nc, "time"), origin = "1970-01-01")
if (!identical(format(core_dates, "%Y-%m"), format(expected_dates, "%Y-%m"))) {
  stop("Canonical core months do not match 1979-2024 monthly period: ",
       min(core_dates), " to ", max(core_dates))
}

find_layer <- function(p0, p1) {
  idx <- which(p0 <= LAYERS$p_bottom & p1 >= LAYERS$p_top)
  if (length(idx) != 1L) stop("Cannot assign pressure interval")
  idx
}
integrate_block <- function(nc, ls, le, sp_start) {
  nb <- le - ls + 1L
  q <- ncvar_get(nc, "q", start = c(1, 1, min(use), ls), count = c(nlon, nlat, length(use), nb))
  zh <- ncvar_get(nc, "z", start = c(1, 1, min(use), ls), count = c(nlon, nlat, length(use), nb)) / G
  q[!is.finite(q) | q < 0] <- NA_real_; zh[!is.finite(zh) | zh < 0] <- NA_real_
  sp <- ncvar_get(sp_nc, "sp", start = c(1, 1, sp_start), count = c(nlon, nlat, nb))
  mass <- array(0, c(nlon, nlat, nb, nlay)); moment <- array(0, c(nlon, nlat, nb, nlay))
  for (k in seq_len(length(plev) - 1L)) {
    p0 <- plev[k]; p1 <- plev[k + 1L]; li <- find_layer(p0, p1)
    q0 <- q[, , k, ]; q1 <- q[, , k + 1L, ]; z0 <- zh[, , k, ]; z1 <- zh[, , k + 1L, ]
    ok <- is.finite(q0) & is.finite(q1) & is.finite(z0) & is.finite(z1) & is.finite(sp)
    full <- ok & sp >= p0 * 100; partial <- ok & sp < p0 * 100 & sp > p1 * 100
    active <- full | partial; start_p <- ifelse(full, p0, sp / 100)
    alpha <- pmax(0, pmin(1, (p0 - start_p) / (p0 - p1)))
    qs <- q0 + (q1 - q0) * alpha; zs <- z0 + (z1 - z0) * alpha
    dp <- pmax(start_p - p1, 0) * 100
    mass[, , , li] <- mass[, , , li] + ifelse(active, .5 * (qs + q1) * dp / G, 0)
    moment[, , , li] <- moment[, , , li] + ifelse(active, .5 * (qs * zs + q1 * z1) * dp / G, 0)
  }
  list(mass = mass, moment = moment)
}

sum_mass <- sum_moment <- sum_w <- sum_Z <- array(0, c(nlon, nlat, nlay))
count_mass <- count_w <- count_Z <- array(0, c(nlon, nlat, nlay))
recon_global <- canon_global <- rep(NA_real_, nt)
sum_abs_diff <- sum_sq_diff <- sum_diff <- 0; n_diff <- 0; max_abs_diff <- 0
sy_r <- sty_r <- st_r <- st2_r <- n_r <- matrix(0, nlon, nlat)
sy_c <- sty_c <- matrix(0, nlon, nlat)
weights <- matrix(rep(pmax(0, cos(lat * pi / 180)), each = nlon), nrow = nlon, ncol = nlat)
wmean <- function(x) { ok <- is.finite(x) & weights > 0; if (!any(ok)) NA_real_ else sum(x[ok] * weights[ok]) / sum(weights[ok]) }

message("Single pass: references plus canonical bridge")
for (fp in file_meta$path) {
  nc <- nc_open(fp); dd <- get_dates(nc)
  for (ls in seq(1, length(dd), by = 12)) {
    le <- min(length(dd), ls + 11L); start_idx <- match(dd[ls], expected_dates)
    blk <- integrate_block(nc, ls, le, start_idx)
    canonical <- ncvar_get(core_nc, "zbar", start = c(1, 1, start_idx), count = c(nlon, nlat, le - ls + 1L))
    for (b in seq_len(dim(blk$mass)[3])) {
      idx <- start_idx + b - 1L; m <- blk$mass[, , b, ]; mm <- blk$moment[, , b, ]
      ok_m <- is.finite(m) & m > 0 & is.finite(mm)
      sum_mass <- sum_mass + ifelse(ok_m, m, 0); sum_moment <- sum_moment + ifelse(ok_m, mm, 0)
      count_mass <- count_mass + ok_m
      total_m <- apply(m, c(1, 2), sum, na.rm = TRUE); total_ok <- is.finite(total_m) & total_m > 0
      wi <- sweep(m, c(1, 2), total_m, "/"); wi[!is.finite(wi)] <- 0
      total_ok3 <- array(rep(total_ok, nlay), c(nlon, nlat, nlay))
      sum_w <- sum_w + ifelse(total_ok3, wi, 0); count_w <- count_w + total_ok3
      Zi <- mm / m; ok_Z <- is.finite(Zi) & ok_m
      sum_Z <- sum_Z + ifelse(ok_Z, Zi, 0); count_Z <- count_Z + ok_Z
      recon <- apply(mm, c(1, 2), sum, na.rm = TRUE) / total_m
      recon[!total_ok] <- NA_real_; canon <- canonical[, , b]
      common <- is.finite(recon) & is.finite(canon)
      recon_global[idx] <- wmean(ifelse(common, recon, NA_real_))
      canon_global[idx] <- wmean(ifelse(common, canon, NA_real_))
      dif <- recon[common] - canon[common]
      sum_abs_diff <- sum_abs_diff + sum(abs(dif)); sum_sq_diff <- sum_sq_diff + sum(dif^2)
      sum_diff <- sum_diff + sum(dif); n_diff <- n_diff + length(dif)
      if (length(dif)) max_abs_diff <- max(max_abs_diff, max(abs(dif)))
      tt <- tdec[idx]
      sy_r <- sy_r + ifelse(common, recon, 0); sty_r <- sty_r + ifelse(common, recon * tt, 0)
      sy_c <- sy_c + ifelse(common, canon, 0); sty_c <- sty_c + ifelse(common, canon * tt, 0)
      st_r <- st_r + common * tt; st2_r <- st2_r + common * tt^2; n_r <- n_r + common
    }
    message(format(dd[ls], "%Y-%m"), " to ", format(dd[le], "%Y-%m"))
  }
  nc_close(nc)
}
nc_close(sp_nc); nc_close(core_nc)

mean_mass <- sum_mass / count_mass; mean_moment <- sum_moment / count_mass
mean_mass[!is.finite(mean_mass)] <- NA; mean_moment[!is.finite(mean_moment)] <- NA
mean_total_mass <- apply(mean_mass, c(1, 2), sum, na.rm = TRUE)
w0A <- sweep(mean_mass, c(1, 2), mean_total_mass, "/"); Z0A <- mean_moment / mean_mass
w0B <- sum_w / count_w; Z0B <- sum_Z / count_Z
w0A[!is.finite(w0A)] <- NA; Z0A[!is.finite(Z0A)] <- NA
w0B[!is.finite(w0B)] <- NA; Z0B[!is.finite(Z0B)] <- NA

den <- n_r * st2_r - st_r^2
trend_recon <- (n_r * sty_r - st_r * sy_r) / den
trend_canon <- (n_r * sty_c - st_r * sy_c) / den
trend_recon[n_r < 24 | !is.finite(trend_recon)] <- NA
trend_canon[n_r < 24 | !is.finite(trend_canon)] <- NA

deseason <- function(x) {
  mo <- as.integer(format(expected_dates, "%m")); clim <- tapply(x, mo, mean, na.rm = TRUE)
  x - clim[as.character(mo)]
}
global_trend <- function(x) coef(lm(deseason(x) ~ tdec))[2]
global_monthly <- data.frame(date = expected_dates, canonical_zbar_m = canon_global,
                             reconstructed_zbar_m = recon_global,
                             difference_m = recon_global - canon_global)
write.csv(global_monthly, file.path(OUT, "FIG1_CANONICAL_BRIDGE_MONTHLY.csv"), row.names = FALSE)

grid_bridge <- expand.grid(lon = lon, lat = lat)
grid_bridge$reconstructed_trend_m_yr <- as.vector(trend_recon)
grid_bridge$canonical_trend_m_yr <- as.vector(trend_canon)
grid_bridge$trend_difference_m_yr <- grid_bridge$reconstructed_trend_m_yr - grid_bridge$canonical_trend_m_yr
write.csv(grid_bridge, file.path(OUT, "FIG1_CANONICAL_BRIDGE_GRID_TRENDS.csv"), row.names = FALSE)

common_map <- is.finite(trend_recon) & is.finite(trend_canon)
bridge_summary <- data.frame(
  metric = c("grid_month_mean_absolute_difference_m", "grid_month_max_absolute_difference_m",
             "grid_month_RMSE_m", "grid_month_mean_bias_m", "global_monthly_correlation",
             "trend_map_correlation", "global_reconstructed_trend_m_yr",
             "global_canonical_trend_m_yr", "global_trend_difference_m_yr",
             "grid_trend_RMSE_m_yr", "grid_trend_max_absolute_difference_m_yr"),
  value = c(sum_abs_diff / n_diff, max_abs_diff, sqrt(sum_sq_diff / n_diff), sum_diff / n_diff,
            cor(recon_global, canon_global, use = "complete.obs"),
            cor(as.vector(trend_recon[common_map]), as.vector(trend_canon[common_map])),
            global_trend(recon_global), global_trend(canon_global),
            global_trend(recon_global) - global_trend(canon_global),
            sqrt(mean((trend_recon[common_map] - trend_canon[common_map])^2)),
            max(abs(trend_recon[common_map] - trend_canon[common_map]))))
write.csv(bridge_summary, file.path(OUT, "FIG1_CANONICAL_BRIDGE_SUMMARY.csv"), row.names = FALSE)

# Algebraically transform the existing exact A trends to Reference B.  The
# exact layer total trend is invariant; only the allocation among the three
# centered product terms changes with the reference state.
if (requireNamespace("data.table", quietly = TRUE)) {
  gd <- data.table::fread(A_GRID_FILE, data.table = FALSE)
} else gd <- read.csv(A_GRID_FILE, check.names = FALSE)
ii <- match(gd$lon, lon); jj <- match(gd$lat, lat); kk <- match(gd$layer, LAYERS$layer)
idx3 <- cbind(ii, jj, kk)
gd$w0A <- w0A[idx3]; gd$Z0A <- Z0A[idx3]; gd$w0B <- w0B[idx3]; gd$Z0B <- Z0B[idx3]
gd$trend_w <- gd$mass_redistribution_m_yr / gd$Z0A
gd$trend_Z <- gd$height_coordinate_m_yr / gd$w0A
gd$mass_B <- gd$Z0B * gd$trend_w
gd$height_B <- gd$w0B * gd$trend_Z
gd$interaction_B <- gd$exact_layer_total_m_yr - gd$mass_B - gd$height_B
gd$area_weight <- pmax(0, cos(gd$lat * pi / 180))

wm_df <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0; sum(x[ok] * w[ok]) / sum(w[ok]) }
summarize_ref <- function(ref = c("A", "B"), groups, labels) {
  ref <- match.arg(ref); out <- list()
  for (h in seq_along(groups)) {
    ss <- gd[gd$layer %in% groups[[h]], ]
    if (ref == "A") {
      mass <- ss$mass_redistribution_m_yr; height <- ss$height_coordinate_m_yr
      inter <- ss$higher_order_interaction_m_yr
    } else { mass <- ss$mass_B; height <- ss$height_B; inter <- ss$interaction_B }
    # Aggregate multiple source layers within each grid cell before weighting.
    ss$mass_component <- mass
    ss$height_component <- height
    ss$interaction_component <- inter
    ss$exact_component <- ss$exact_layer_total_m_yr
    tmp <- aggregate(cbind(mass_component, height_component, interaction_component, exact_component) ~ cell_id + lat,
                     data = ss,
                     FUN = function(z) if (all(!is.finite(z))) NA_real_ else sum(z, na.rm = TRUE))
    ww <- pmax(0, cos(tmp$lat * pi / 180))
    out[[h]] <- data.frame(reference = paste0("Reference_", ref), layer = labels[h],
                           mass_redistribution_m_yr = wm_df(tmp$mass_component, ww),
                           within_layer_mean_height_m_yr = wm_df(tmp$height_component, ww),
                           higher_order_interaction_m_yr = wm_df(tmp$interaction_component, ww),
                           exact_layer_total_m_yr = wm_df(tmp$exact_component, ww))
  }
  do.call(rbind, out)
}
groups5 <- list("1000-850", "850-700", "700-500", "500-300", c("300-200", "200-100"))
labels5 <- c("1000-850", "850-700", "700-500", "500-300", "300-100")
sens <- rbind(summarize_ref("A", groups5, labels5), summarize_ref("B", groups5, labels5))
actual <- bridge_summary$value[bridge_summary$metric == "global_reconstructed_trend_m_yr"]
global_rows <- do.call(rbind, lapply(split(sens, sens$reference), function(z) {
  data.frame(reference = z$reference[1], layer = "GLOBAL",
             mass_redistribution_m_yr = sum(z$mass_redistribution_m_yr),
             within_layer_mean_height_m_yr = sum(z$within_layer_mean_height_m_yr),
             higher_order_interaction_m_yr = sum(z$higher_order_interaction_m_yr),
             exact_layer_total_m_yr = sum(z$exact_layer_total_m_yr))
}))
sens <- rbind(sens, global_rows)
sens$mass_share_pct <- 100 * sens$mass_redistribution_m_yr / sens$exact_layer_total_m_yr
sens$mean_height_share_pct <- 100 * sens$within_layer_mean_height_m_yr / sens$exact_layer_total_m_yr
sens$interaction_share_pct <- 100 * sens$higher_order_interaction_m_yr / sens$exact_layer_total_m_yr
sens$canonical_actual_trend_m_yr <- actual
write.csv(sens, file.path(OUT, "FIG1_REFERENCE_STATE_SENSITIVITY.csv"), row.names = FALSE)

a <- sens[sens$reference == "Reference_A" & sens$layer != "GLOBAL", ]
b <- sens[sens$reference == "Reference_B" & sens$layer != "GLOBAL", ]
exact_change <- max(abs(a$exact_layer_total_m_yr - b$exact_layer_total_m_yr))
ga <- sens[sens$reference == "Reference_A" & sens$layer == "GLOBAL", ]
gb <- sens[sens$reference == "Reference_B" & sens$layer == "GLOBAL", ]
share_change <- max(abs(c(ga$mass_share_pct - gb$mass_share_pct,
                          ga$mean_height_share_pct - gb$mean_height_share_pct,
                          ga$interaction_share_pct - gb$interaction_share_pct)))

writeLines(c(
  "# Figure 1 canonical bridge QA", "",
  paste0("- Monthly global correlation: ", signif(bridge_summary$value[bridge_summary$metric == "global_monthly_correlation"], 12)),
  paste0("- Grid-month mean absolute difference: ", signif(sum_abs_diff / n_diff, 8), " m."),
  paste0("- Grid-month RMSE: ", signif(sqrt(sum_sq_diff / n_diff), 8), " m."),
  paste0("- Grid-month maximum absolute difference: ", signif(max_abs_diff, 8), " m."),
  paste0("- Trend-map correlation: ", signif(bridge_summary$value[bridge_summary$metric == "trend_map_correlation"], 12)),
  paste0("- Global reconstructed trend: ", signif(global_trend(recon_global), 12), " m yr-1."),
  paste0("- Global canonical trend: ", signif(global_trend(canon_global), 12), " m yr-1."),
  paste0("- Global trend difference: ", signif(global_trend(recon_global) - global_trend(canon_global), 12), " m yr-1."),
  if (abs(global_trend(recon_global) - global_trend(canon_global)) < 1e-6) "- PASS: canonical and reconstructed global trends agree to numerical precision (<1e-6 m yr-1)." else "- REVIEW: canonical and reconstructed trends differ by more than 1e-6 m yr-1."
), file.path(OUT, "FIG1_CANONICAL_BRIDGE_QA.md"))

writeLines(c(
  "# Figure 1 reference-state sensitivity", "",
  "Reference A uses mean layer mass divided by mean total mass and mean layer moment divided by mean layer mass.",
  "Reference B uses the temporal means of instantaneous layer mass fraction and instantaneous layer mean height.",
  "Exact layer totals are algebraically invariant; the mass/mean-height/interaction allocation can depend on the reference state.", "",
  paste0("- Maximum absolute exact-layer change: ", signif(exact_change, 8), " m yr-1."),
  paste0("- Maximum global mechanism-share change: ", signif(share_change, 6), " percentage points."),
  if (exact_change <= 0.05) "- Exact layer-total stability threshold: PASS." else "- Exact layer-total stability threshold: FAIL.",
  if (share_change <= 10) "- Mechanism-share stability threshold: PASS." else "- Mechanism-share stability threshold: SENSITIVE.",
  "The term formerly called height-coordinate change is more accurately labelled within-layer mean-height change because Z_i is moisture-weighted within each pressure layer."
), file.path(OUT, "FIG1_REFERENCE_STATE_SENSITIVITY.md"))

cat("DONE reference-state and canonical bridge analysis\n")
