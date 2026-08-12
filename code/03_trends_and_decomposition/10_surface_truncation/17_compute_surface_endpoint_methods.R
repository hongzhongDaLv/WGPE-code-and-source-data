options(stringsAsFactors = FALSE)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))
suppressPackageStartupMessages(library(ncdf4))

ROOT <- "__WGPE_PROJECT_ROOT__"
OUT <- file.path(ROOT, "output_reviewer_P0_resolution")
PLEV_DIR <- file.path(ROOT, "era5_plev_ascii")
G <- 9.80665
EPS <- 0.622
P_TOP <- 100

decode_nc_time <- function(nc, nm = NULL) {
  if (is.null(nm)) nm <- intersect(c("valid_time", "time"), names(nc$dim))[1]
  d <- nc$dim[[nm]]
  vals <- d$vals
  units <- d$units
  origin <- sub("^[A-Za-z]+ since ", "", units)
  mult <- if (grepl("hours since", units)) 3600 else if (grepl("days since", units)) 86400 else 1
  as.POSIXct(origin, tz = "UTC") + vals * mult
}

nearest_index <- function(src, target, tol = 1e-6) {
  idx <- vapply(target, function(x) which.min(abs(src - x)), integer(1))
  if (max(abs(src[idx] - target)) > tol) stop("Target grid is not an exact subset of source grid")
  idx
}

q_from_dewpoint <- function(td_k, sp_pa) {
  # Bolton-style saturation vapour pressure over liquid water.  The endpoint
  # is a sensitivity construction, not a replacement for ERA5 model levels.
  tc <- td_k - 273.15
  e <- 611.2 * exp(17.67 * tc / (tc + 243.5))
  pmax(0, pmin(0.1, EPS * e / (sp_pa - (1 - EPS) * e)))
}

integrate_surface_method <- function(q, zh, sp_hpa, qsurf, zsurf, plev, method) {
  nx <- dim(q)[1]; ny <- dim(q)[2]
  mass <- matrix(0, nx, ny)
  moment <- matrix(0, nx, ny)
  used <- matrix(FALSE, nx, ny)
  pbot <- pmin(sp_hpa, 1000)
  for (k in seq_len(length(plev) - 1L)) {
    pb <- plev[k]
    pt <- plev[k + 1L]
    if (pt < P_TOP) next
    overlap_bottom <- pmin(pb, pbot)
    active <- is.finite(pbot) & overlap_bottom > pt
    if (!any(active)) next
    dp <- (overlap_bottom - pt) * 100
    q_top <- q[, , k + 1L]
    z_top <- zh[, , k + 1L]
    full <- active & pbot > pb + 1e-8
    exact <- active & abs(pbot - pb) <= 1e-8
    crossing <- active & pbot < pb - 1e-8
    q_bottom <- q[, , k]
    z_bottom <- zh[, , k]
    if (method == "P1_surface_anchored") {
      q_bottom[exact | crossing] <- qsurf[exact | crossing]
      z_bottom[exact | crossing] <- zsurf[exact | crossing]
    } else if (method == "P2_one_sided") {
      # For a partial crossing interval, extend the first real above-surface
      # pressure-level endpoint downward without reading the below-ground one.
      q_bottom[crossing] <- q_top[crossing]
      z_bottom[crossing] <- z_top[crossing]
      # At an exact pressure-level boundary use that real level as endpoint.
    } else stop("Unknown method")
    valid <- active & is.finite(q_bottom) & is.finite(q_top) & is.finite(z_bottom) & is.finite(z_top)
    if (!any(valid)) next
    q_bottom[q_bottom < 0] <- NA_real_
    q_top[q_top < 0] <- NA_real_
    z_bottom[z_bottom < 0] <- NA_real_
    z_top[z_top < 0] <- NA_real_
    valid <- valid & is.finite(q_bottom) & is.finite(q_top) & is.finite(z_bottom) & is.finite(z_top)
    dm <- 0.5 * (q_bottom + q_top) * dp / G
    dM <- 0.5 * (q_bottom * z_bottom + q_top * z_top) * dp / G
    mass[valid] <- mass[valid] + dm[valid]
    moment[valid] <- moment[valid] + dM[valid]
    used[valid] <- TRUE
  }
  mass[!used | mass <= 0] <- NA_real_
  moment[!used | mass <= 0] <- NA_real_
  zbar <- moment / mass
  list(iwv = mass, zbar = zbar, wgpe = G * moment)
}

sp_nc <- nc_open(file.path(OUT, "ERA5_surface_pressure_aligned_0_359.nc"))
lon <- sp_nc$dim$longitude$vals
lat <- sp_nc$dim$latitude$vals
dates <- decode_nc_time(sp_nc, "time")
nt <- length(dates)

d2_nc <- nc_open(file.path(ROOT, "2m dewpoint temperature.nc"))
d2_lon <- d2_nc$dim$longitude$vals
d2_lat <- d2_nc$dim$latitude$vals
d2_dates <- decode_nc_time(d2_nc, "valid_time")
ix_d2 <- nearest_index(d2_lon, lon)
iy_d2 <- nearest_index(d2_lat, lat)
d2_lookup <- match(format(dates, "%Y-%m"), format(d2_dates, "%Y-%m"))
if (anyNA(d2_lookup)) stop("Dewpoint does not cover all analysis months")

zg_nc <- nc_open(file.path(ROOT, "geopotential.nc"))
zg_lon <- zg_nc$dim$longitude$vals
zg_lat <- zg_nc$dim$latitude$vals
ix_zg <- nearest_index(zg_lon, lon)
iy_zg <- nearest_index(zg_lat, lat)
zg0 <- ncvar_get(zg_nc, "z", start = c(1, 1, 1), count = c(-1, -1, 1))
zsurf <- zg0[ix_zg, iy_zg] / G
zg_range <- 0
if (zg_nc$dim$valid_time$len > 1) {
  for (ti in 2:zg_nc$dim$valid_time$len) {
    zgi <- ncvar_get(zg_nc, "z", start = c(1, 1, ti), count = c(-1, -1, 1))
    zg_range <- max(zg_range, max(abs(zgi[ix_zg, iy_zg] / G - zsurf), na.rm = TRUE))
  }
}
nc_close(zg_nc)
if (zg_range > 1e-4) warning("Surface geopotential is not static across its 12 records")

plev_files <- list.files(PLEV_DIR, pattern = "\\.nc$", full.names = TRUE)
plev_ncs <- lapply(plev_files, nc_open)
plev_dates <- lapply(plev_ncs, decode_nc_time)
file_month <- do.call(rbind, lapply(seq_along(plev_ncs), function(fi) {
  data.frame(fi = fi, ti = seq_along(plev_dates[[fi]]), ym = format(plev_dates[[fi]], "%Y-%m"))
}))
file_month <- file_month[match(format(dates, "%Y-%m"), file_month$ym), ]
if (anyNA(file_month$fi)) stop("Pressure-level files do not cover all months")

plev_all <- plev_ncs[[1]]$dim$pressure_level$vals
keep <- which(plev_all <= 1000 & plev_all >= P_TOP)
plev <- plev_all[keep]
if (!all(diff(plev) < 0)) stop("Pressure levels must be descending")

lon_dim <- ncdim_def("longitude", "degrees_east", lon)
lat_dim <- ncdim_def("latitude", "degrees_north", lat)
time_days <- as.numeric(difftime(dates, as.POSIXct("1970-01-01", tz = "UTC"), units = "days"))
time_dim <- ncdim_def("time", "days since 1970-01-01", time_days, unlim = FALSE, calendar = "standard")
defs <- list()
for (method in c("P1", "P2")) for (variable in c("iwv", "zbar", "wgpe")) {
  unit <- c(iwv = "kg m-2", zbar = "m", wgpe = "J m-2")[[variable]]
  defs[[length(defs) + 1L]] <- ncvar_def(paste0(method, "_", variable), unit, list(lon_dim, lat_dim, time_dim), -9999, prec = "float", compression = NA)
}
out_file <- file.path(OUT, "ERA5_TOP100_SURFACE_ENDPOINT_METHODS_monthly.nc")
out_nc <- nc_create(out_file, defs, force_v4 = TRUE)
ncatt_put(out_nc, 0, "analysis_period", "1979-01 to 2024-12")
ncatt_put(out_nc, 0, "pressure_top_hPa", P_TOP)
ncatt_put(out_nc, 0, "P1_definition", "surface endpoint from surface geopotential and q diagnosed from 2 m dewpoint plus surface pressure; no below-ground q/z used")
ncatt_put(out_nc, 0, "P2_definition", "crossing partial interval uses first above-surface pressure-level q/z at both endpoints; no below-ground q/z used")
ncatt_put(out_nc, 0, "surface_q_formula", "e=611.2*exp(17.67*Td_C/(Td_C+243.5)); q=0.622e/(sp-0.378e)")

progress <- file.path(OUT, "SURFACE_ENDPOINT_BUILD_PROGRESS_TOP100.txt")
for (tt in seq_len(nt)) {
  fi <- file_month$fi[tt]; ti <- file_month$ti[tt]
  nc <- plev_ncs[[fi]]
  q <- ncvar_get(nc, "q", start = c(1, 1, min(keep), ti), count = c(-1, -1, length(keep), 1))
  geop <- ncvar_get(nc, "z", start = c(1, 1, min(keep), ti), count = c(-1, -1, length(keep), 1))
  zh <- geop / G
  sp <- ncvar_get(sp_nc, "sp", start = c(1, 1, tt), count = c(-1, -1, 1))
  td_full <- ncvar_get(d2_nc, "d2m", start = c(1, 1, d2_lookup[tt]), count = c(-1, -1, 1))
  td <- td_full[ix_d2, iy_d2]
  qsurf <- q_from_dewpoint(td, sp)
  p1 <- integrate_surface_method(q, zh, sp / 100, qsurf, zsurf, plev, "P1_surface_anchored")
  p2 <- integrate_surface_method(q, zh, sp / 100, qsurf, zsurf, plev, "P2_one_sided")
  for (v in names(p1)) ncvar_put(out_nc, paste0("P1_", v), p1[[v]], start = c(1, 1, tt), count = c(-1, -1, 1))
  for (v in names(p2)) ncvar_put(out_nc, paste0("P2_", v), p2[[v]], start = c(1, 1, tt), count = c(-1, -1, 1))
  if (tt %% 12 == 0 || tt == nt) {
    writeLines(sprintf("completed %d/%d (%s)", tt, nt, format(dates[tt], "%Y-%m")), progress)
    message(sprintf("completed %d/%d", tt, nt))
  }
}

nc_close(out_nc)
nc_close(sp_nc)
nc_close(d2_nc)
invisible(lapply(plev_ncs, nc_close))

manifest <- data.frame(
  output = normalizePath(out_file, winslash = "/", mustWork = TRUE),
  months = nt,
  grid = paste(length(lon), length(lat), sep = "x"),
  top_hPa = P_TOP,
  surface_geopotential_temporal_max_abs_difference_m = zg_range,
  methods = "P1_surface_anchored;P2_one_sided",
  below_ground_qz_used = FALSE
)
write.csv(manifest, file.path(OUT, "SURFACE_ENDPOINT_MONTHLY_BUILD_MANIFEST_TOP100.csv"), row.names = FALSE)
