canonicalize_grid <- function(lon, lat, field) {
  lon2 <- ifelse(lon >= 180, lon - 360, lon)
  oi <- order(lon2)
  oj <- order(lat)
  list(lon = lon2[oi], lat = lat[oj], field = field[oi, oj, drop = FALSE])
}

area_mean_field <- function(field, lat) {
  w <- rep(cos(lat * pi / 180), each = nrow(field))
  weighted.mean(as.vector(field), w, na.rm = TRUE)
}

profile_separation_map <- function(zbar, iwv, dates, block_size = 750L) {
  d <- dim(zbar)
  stopifnot(identical(dim(iwv), d), d[3] == length(dates))
  ncell <- d[1] * d[2]
  nt <- d[3]
  dim(zbar) <- c(ncell, nt)
  dim(iwv) <- c(ncell, nt)
  month <- as.integer(format(dates, "%m"))
  out <- rep(NA_real_, ncell)
  for (start in seq.int(1L, ncell, by = block_size)) {
    ii <- start:min(ncell, start + block_size - 1L)
    z <- zbar[ii, , drop = FALSE]
    q <- iwv[ii, , drop = FALSE]
    for (m in 1:12) {
      jj <- which(month == m)
      z[, jj] <- z[, jj, drop = FALSE] - rowMeans(z[, jj, drop = FALSE])
      q[, jj] <- q[, jj, drop = FALSE] - rowMeans(q[, jj, drop = FALSE])
    }
    z <- z - rowMeans(z)
    q <- q - rowMeans(q)
    den <- sqrt(rowSums(z ^ 2) * rowSums(q ^ 2))
    r <- rowSums(z * q) / den
    out[ii] <- 100 * (1 - r ^ 2)
  }
  matrix(out, nrow = d[1], ncol = d[2])
}

load_figure_data_FIXED_FULL <- function(root) {
  fig_dir <- file.path(root, "Figures_FIXED_FULL")
  out_dir <- file.path(root, "output_attribution_minimal")
  core_file <- file.path(out_dir, "global_wgpe_1deg_surface_truncated_300hPa_FIXED_FULL.RData")
  mask_file <- file.path(fig_dir, "Fig1_zbar_trend_effective_n_BH_mask_FIXED_FULL.rds")
  sig_file <- file.path(out_dir, "WGPE_effective_n_BH_grid_FIXED_FULL.rds")
  time_file <- file.path(fig_dir, "Fig1b_ERA5_monthly_series_FIXED_FULL.csv")
  layer_file <- file.path(out_dir, "WGPE_ztrend_layer_contribution_FIXED_FULL.csv")
  wgpe_file <- file.path(out_dir, "WGPE_trend_decomposition_FIXED_FULL.csv")
  coast_file <- file.path(root, "R", "coastlines_source_FIXED_FULL.csv")
  needed <- c(core_file, mask_file, sig_file, time_file, layer_file, wgpe_file, coast_file)
  if (any(!file.exists(needed))) stop("Missing FIXED_FULL input(s): ", paste(needed[!file.exists(needed)], collapse = "; "))

  mask <- readRDS(mask_file)
  sig <- readRDS(sig_file)
  time <- read.csv(time_file, stringsAsFactors = FALSE)
  time$date <- as.Date(time$date)
  layer <- read.csv(layer_file, stringsAsFactors = FALSE, check.names = FALSE)
  wgpe <- read.csv(wgpe_file, stringsAsFactors = FALSE, check.names = FALSE)

  e <- new.env(parent = emptyenv())
  load(core_file, envir = e)
  required <- c("lon", "lat", "all_dates", "zbar_st", "iwv_st")
  if (!all(required %in% ls(e))) stop("Corrected core is missing required objects")
  sep_map <- profile_separation_map(e$zbar_st, e$iwv_st, as.Date(e$all_dates))

  lon <- as.numeric(mask$lon)
  lat <- as.numeric(mask$lat)
  if (!identical(as.numeric(e$lon), lon) || !identical(as.numeric(e$lat), lat)) stop("Core/mask coordinate mismatch")
  if (!identical(as.numeric(sig$lon), lon) || !identical(as.numeric(sig$lat), lat)) stop("Significance/core coordinate mismatch")

  lon2 <- ifelse(lon >= 180, lon - 360, lon)
  oi <- order(lon2)
  oj <- order(lat)
  lon_c <- lon2[oi]
  lat_c <- lat[oj]
  canon <- function(x) x[oi, oj, drop = FALSE]

  ztrend <- canon(mask$trend_m_yr)
  nonsig <- canon(mask$non_significant)
  wgpe_trend <- canon(sig$variables$W_GPE$trend)
  separation <- canon(sep_map)
  rm(e, sep_map)
  gc(verbose = FALSE)

  geometry <- make_robinson_grid_geometry(lon_c, lat_c)
  coast <- load_robinson_coastline(coast_file)
  graticule <- make_robinson_graticule()
  border <- make_robinson_border()

  z_map_df <- attach_map_field(geometry, ztrend)
  sep_map_df <- attach_map_field(geometry, separation)
  wgpe_map_df <- attach_map_field(geometry, wgpe_trend)

  weights <- rep(cos(lat_c * pi / 180), each = length(lon_c))
  profile_area <- weighted.mean(as.vector(separation), weights, na.rm = TRUE)
  profile_grid <- mean(separation, na.rm = TRUE)
  z_global <- weighted.mean(as.vector(ztrend), weights, na.rm = TRUE)
  wgpe_global <- weighted.mean(as.vector(wgpe_trend), weights, na.rm = TRUE)
  sig_pos_area <- mask$fractions$positive_effective_BH_area_pct[1]
  sig_neg_area <- mask$fractions$negative_effective_BH_area_pct[1]
  sig_pos_grid <- mask$fractions$positive_effective_BH_grid_pct[1]
  sig_neg_grid <- mask$fractions$negative_effective_BH_grid_pct[1]

  time_fit <- lm(zbar_anom_m ~ decimal_year, data = time)
  time_trend <- unname(coef(time_fit)[2])

  layer_ds <- layer[layer$region == "Global" & layer$method_raw_or_deseasonalized == "deseasonalized", , drop = FALSE]
  layer_ds$layer <- factor(layer_ds$layer, levels = c("1000-850", "850-700", "700-500", "500-300"))
  layer_ds <- layer_ds[order(layer_ds$layer), , drop = FALSE]
  wgpe_ds <- wgpe[wgpe$method_raw_or_deseasonalized == "deseasonalized", , drop = FALSE]
  wgpe_ds$region <- factor(wgpe_ds$region,
    levels = c("Global", "S extratropics", "S tropics", "Deep tropics", "N tropics", "N extratropics"))
  wgpe_ds <- wgpe_ds[order(wgpe_ds$region), , drop = FALSE]

  if (abs(z_global - 0.788114) > 0.001) stop("Fig. 1 <z> global-trend QC failed")
  if (abs(time_trend - z_global) > 1e-5) stop("Fig. 1 time-series/map trend mismatch")
  if (abs(profile_area - 64.96) > 0.08 || abs(profile_grid - 74.89) > 0.08) stop("Profile-separation QC failed")
  if (abs(layer_ds$sum_mass_contribution_m_yr[1] - 0.649037) > 0.001) stop("Layer mass-summary QC failed")
  if (abs(layer_ds$sum_height_contribution_m_yr[1] - 0.231348) > 0.001) stop("Layer height-summary QC failed")
  if (abs(layer_ds$residual_pct_of_actual[1] + 11.65) > 0.03) stop("Layer residual QC failed")
  global_row <- wgpe_ds[wgpe_ds$region == "Global", ]
  deep_row <- wgpe_ds[wgpe_ds$region == "Deep tropics", ]
  if (abs(global_row$iwv_channel_pct - 59.9745) > 0.02 ||
      abs(global_row$z_channel_pct - 37.2458) > 0.02 ||
      deep_row$z_channel_pct <= deep_row$iwv_channel_pct) stop("W-GPE decomposition QC failed")

  list(
    lon = lon_c, lat = lat_c,
    geometry = geometry, coast = coast, graticule = graticule, border = border,
    ztrend = ztrend, nonsig = nonsig, separation = separation, wgpe_trend = wgpe_trend,
    z_map_df = z_map_df, sep_map_df = sep_map_df, wgpe_map_df = wgpe_map_df,
    time = time, layer = layer_ds, wgpe = wgpe_ds,
    stats = list(
      z_global = z_global, wgpe_global = wgpe_global,
      sig_pos_area = sig_pos_area, sig_neg_area = sig_neg_area,
      sig_pos_grid = sig_pos_grid, sig_neg_grid = sig_neg_grid,
      n_grid = length(ztrend), profile_area = profile_area, profile_grid = profile_grid,
      time_trend = time_trend,
      layer_actual = layer_ds$actual_zbar_trend_m_yr[1],
      layer_mass = layer_ds$sum_mass_contribution_m_yr[1],
      layer_height = layer_ds$sum_height_contribution_m_yr[1],
      layer_residual = layer_ds$residual_m_yr[1],
      layer_residual_pct = layer_ds$residual_pct_of_actual[1],
      wgpe_global_pct = c(IWV = global_row$iwv_channel_pct,
                          zbar = global_row$z_channel_pct,
                          interaction = global_row$interaction_channel_pct),
      wgpe_deep_pct = c(IWV = deep_row$iwv_channel_pct,
                        zbar = deep_row$z_channel_pct,
                        interaction = deep_row$interaction_channel_pct)
    )
  )
}

load_figure_data_v3_FIXED_FULL <- function(root) {
  fig_dir <- file.path(root, "Figures_FIXED_FULL")
  out_dir <- file.path(root, "output_attribution_minimal")
  mask_file <- file.path(fig_dir, "Fig1_zbar_trend_effective_n_BH_mask_FIXED_FULL.rds")
  sig_file <- file.path(out_dir, "WGPE_effective_n_BH_grid_FIXED_FULL.rds")
  time_file <- file.path(fig_dir, "Fig1b_ERA5_monthly_series_FIXED_FULL.csv")
  layer_file <- file.path(out_dir, "WGPE_ztrend_layer_contribution_FIXED_FULL.csv")
  wgpe_file <- file.path(out_dir, "WGPE_trend_decomposition_FIXED_FULL.csv")
  corr_file <- file.path(fig_dir, "corrected_trend_correlations_1979_2024.csv")
  recheck_file <- file.path(fig_dir, "corrected_core_recheck_table.csv")
  coast_file <- file.path(root, "R", "coastlines_source_FIXED_FULL.csv")
  needed <- c(mask_file, sig_file, time_file, layer_file, wgpe_file,
              corr_file, recheck_file, coast_file)
  if (any(!file.exists(needed))) stop("Missing audited FIXED_FULL input(s): ",
                                      paste(needed[!file.exists(needed)], collapse = "; "))

  mask <- readRDS(mask_file)
  sig <- readRDS(sig_file)
  time <- read.csv(time_file, stringsAsFactors = FALSE)
  time$date <- as.Date(time$date)
  layer <- read.csv(layer_file, stringsAsFactors = FALSE, check.names = FALSE)
  wgpe <- read.csv(wgpe_file, stringsAsFactors = FALSE, check.names = FALSE)
  corr <- read.csv(corr_file, stringsAsFactors = FALSE, check.names = FALSE)
  recheck <- read.csv(recheck_file, stringsAsFactors = FALSE, check.names = FALSE)

  lon <- as.numeric(mask$lon)
  lat <- as.numeric(mask$lat)
  if (!identical(as.numeric(sig$lon), lon) || !identical(as.numeric(sig$lat), lat)) {
    stop("Significance/core coordinate mismatch")
  }
  lon2 <- ifelse(lon >= 180, lon - 360, lon)
  oi <- order(lon2); oj <- order(lat)
  lon_c <- lon2[oi]; lat_c <- lat[oj]
  canon <- function(x) x[oi, oj, drop = FALSE]
  ztrend <- canon(mask$trend_m_yr)
  nonsig <- canon(mask$non_significant)
  wgpe_trend <- canon(sig$variables$W_GPE$trend)

  geometry <- make_robinson_grid_geometry(lon_c, lat_c)
  coast <- load_robinson_coastline(coast_file)
  graticule <- make_robinson_graticule()
  border <- make_robinson_border()
  z_map_df <- attach_map_field(geometry, ztrend)
  wgpe_map_df <- attach_map_field(geometry, wgpe_trend)

  weights <- rep(cos(lat_c * pi / 180), each = length(lon_c))
  z_global <- weighted.mean(as.vector(ztrend), weights, na.rm = TRUE)
  wgpe_global <- weighted.mean(as.vector(wgpe_trend), weights, na.rm = TRUE)
  time_trend <- unname(coef(lm(zbar_anom_m ~ decimal_year, data = time))[2])

  get_metric <- function(name) {
    z <- recheck$value[recheck$metric == name]
    if (length(z) != 1L || !is.finite(z)) stop("Missing recheck metric: ", name)
    z
  }
  profile_area <- get_metric("one_minus_R2_zbar_IWV_area_weighted_mean")
  profile_grid <- get_metric("one_minus_R2_zbar_IWV_unweighted_gridcell_mean")

  wanted <- c(
    IWV = "<z> trend vs IWV trend",
    T2m = "<z> trend vs T2m trend (K)",
    Precipitation = "<z> trend vs total precipitation trend (mm day^-1)",
    SSRD = "<z> trend vs ssrd trend (W m^-2)"
  )
  corr_rows <- lapply(names(wanted), function(label) {
    row <- corr[corr$comparison == wanted[[label]], , drop = FALSE]
    if (nrow(row) != 1L) stop("Missing trend-correlation row: ", wanted[[label]])
    data.frame(variable = label, r = row$spatial_r_unweighted, stringsAsFactors = FALSE)
  })
  correlations <- do.call(rbind, corr_rows)
  correlations$variable <- factor(correlations$variable,
                                  levels = c("SSRD", "Precipitation", "T2m", "IWV"))

  layer_ds <- layer[layer$region == "Global" &
                      layer$method_raw_or_deseasonalized == "deseasonalized", , drop = FALSE]
  layer_ds$layer <- factor(layer_ds$layer,
                          levels = c("1000-850", "850-700", "700-500", "500-300"))
  layer_ds <- layer_ds[order(layer_ds$layer), , drop = FALSE]
  wgpe_ds <- wgpe[wgpe$method_raw_or_deseasonalized == "deseasonalized", , drop = FALSE]
  wgpe_ds$region <- factor(wgpe_ds$region,
    levels = c("Global", "S extratropics", "S tropics", "Deep tropics", "N tropics", "N extratropics"))
  wgpe_ds <- wgpe_ds[order(wgpe_ds$region), , drop = FALSE]
  global_row <- wgpe_ds[wgpe_ds$region == "Global", ]
  deep_row <- wgpe_ds[wgpe_ds$region == "Deep tropics", ]

  expected_r <- c(IWV = 0.393883671914035, T2m = -0.522518723409532,
                  Precipitation = 0.248376991677403, SSRD = -0.146119707334996)
  observed_r <- setNames(correlations$r, as.character(correlations$variable))
  if (max(abs(observed_r[names(expected_r)] - expected_r)) > 1e-12) stop("Correlation-strip QC failed")
  if (abs(z_global - 0.788114) > 0.001 || abs(time_trend - z_global) > 1e-5) stop("Fig. 1 v3 trend QC failed")
  if (abs(profile_area - 64.96393) > 0.001 || abs(profile_grid - 74.88834) > 0.001) stop("Stored profile-statistic QC failed")
  if (abs(layer_ds$residual_pct_of_actual[1] + 11.65) > 0.03) stop("Layer residual QC failed")
  if (deep_row$z_channel_pct <= deep_row$iwv_channel_pct) stop("Deep-tropical contrast QC failed")

  list(
    lon = lon_c, lat = lat_c, geometry = geometry,
    coast = coast, graticule = graticule, border = border,
    ztrend = ztrend, nonsig = nonsig, wgpe_trend = wgpe_trend,
    z_map_df = z_map_df, wgpe_map_df = wgpe_map_df,
    time = time, correlations = correlations, layer = layer_ds, wgpe = wgpe_ds,
    stats = list(
      z_global = z_global, wgpe_global = wgpe_global,
      sig_pos_area = mask$fractions$positive_effective_BH_area_pct[1],
      sig_neg_area = mask$fractions$negative_effective_BH_area_pct[1],
      sig_pos_grid = mask$fractions$positive_effective_BH_grid_pct[1],
      sig_neg_grid = mask$fractions$negative_effective_BH_grid_pct[1],
      n_grid = length(ztrend), profile_area = profile_area, profile_grid = profile_grid,
      time_trend = time_trend,
      layer_actual = layer_ds$actual_zbar_trend_m_yr[1],
      layer_mass = layer_ds$sum_mass_contribution_m_yr[1],
      layer_height = layer_ds$sum_height_contribution_m_yr[1],
      layer_residual = layer_ds$residual_m_yr[1],
      layer_residual_pct = layer_ds$residual_pct_of_actual[1],
      wgpe_global_pct = c(IWV = global_row$iwv_channel_pct,
                          zbar = global_row$z_channel_pct,
                          interaction = global_row$interaction_channel_pct),
      wgpe_deep_pct = c(IWV = deep_row$iwv_channel_pct,
                        zbar = deep_row$z_channel_pct,
                        interaction = deep_row$interaction_channel_pct)
    )
  )
}
