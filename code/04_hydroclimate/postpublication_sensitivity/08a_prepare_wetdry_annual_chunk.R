options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !grepl("^[0-9]+$", args[1])) {
  stop("Usage: Rscript 08a_prepare_wetdry_annual_chunk.R CHUNK_ID")
}
chunk_id <- as.integer(args[1])

ROOT <- "__WGPE_PROJECT_ROOT__"
CANON <- file.path(ROOT, "CANONICAL_WGPE_TOP100_v2")
OUT <- file.path(
  CANON, "08_REPRODUCTION_OUTPUT",
  "wetdry_spatiotemporal_uncertainty", "annual_chunks"
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

CORE_FILE <- file.path(
  ROOT, "output_TOP100_revision", "derived_data", "ERA5_TOP100",
  "global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData"
)
SINGLE_FILE <- file.path(
  ROOT, "output_attribution_minimal",
  "ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds"
)
COORD_FILE <- file.path(
  ROOT, "scripts_TOP100_revision", "18_coordinate_harmonization",
  "coordinate_harmonization.R"
)
source(COORD_FILE)

G <- 9.80665
CHUNK_SIZE <- 10000L
BATCH <- 1000L

e <- new.env(parent = emptyenv())
load(CORE_FILE, envir = e)
lon <- as.numeric(e$lon)
lat <- as.numeric(e$lat)
dates <- as.Date(e$all_dates)
nt <- length(dates)
ncell <- length(lon) * length(lat)

chunk_start <- (chunk_id - 1L) * CHUNK_SIZE + 1L
chunk_end <- min(ncell, chunk_id * CHUNK_SIZE)
if (chunk_start > ncell) stop("Chunk ID exceeds the grid.")
ids_chunk <- chunk_start:chunk_end

Imat <- matrix(
  as.numeric(e$iwv_st),
  nrow = ncell,
  ncol = nt
)[ids_chunk, , drop = FALSE]
Zmat <- matrix(
  as.numeric(e$zbar_st),
  nrow = ncell,
  ncol = nt
)[ids_chunk, , drop = FALSE]
rm(e)
gc()

s <- readRDS(SINGLE_FILE)
p_aligned <- reorder_exact_grid(
  s$tp_native_m, s$lon, s$lat, lon, lat
)
Pmat <- matrix(
  as.numeric(p_aligned),
  nrow = ncell,
  ncol = nt
)[ids_chunk, , drop = FALSE]
rm(s, p_aligned)
gc()

month <- as.integer(format(dates, "%m"))
year <- as.integer(format(dates, "%Y"))
years <- sort(unique(year))
ny <- length(years)
year_id <- match(year, years)
td <- as.numeric(dates)
tc <- td - mean(td)
n_chunk <- length(ids_chunk)

adjust_and_anomaly <- function(m) {
  an <- m
  for (mm in 1:12) {
    ii <- which(month == mm)
    an[, ii] <- sweep(
      m[, ii, drop = FALSE],
      1,
      rowMeans(m[, ii, drop = FALSE], na.rm = TRUE),
      "-"
    )
  }
  slope <- as.numeric(an %*% tc / sum(tc^2))
  adjusted <- m - outer(slope, tc)
  climatology <- matrix(NA_real_, nrow(m), 12)
  for (mm in 1:12) {
    climatology[, mm] <- rowMeans(
      adjusted[, month == mm, drop = FALSE],
      na.rm = TRUE
    )
  }
  list(
    adjusted = adjusted,
    climatology = climatology,
    anomaly = adjusted - climatology[, month]
  )
}

annual_count <- function(mask) {
  out <- matrix(0L, nrow(mask), ny)
  for (yy in seq_len(ny)) {
    out[, yy] <- rowSums(mask[, year_id == yy, drop = FALSE], na.rm = TRUE)
  }
  out
}

annual_sum <- function(x, mask) {
  x[!mask] <- 0
  x[!is.finite(x)] <- 0
  out <- matrix(0, nrow(x), ny)
  for (yy in seq_len(ny)) {
    out[, yy] <- rowSums(x[, year_id == yy, drop = FALSE])
  }
  out
}

new_store <- function(integer = FALSE) {
  if (integer) matrix(0L, n_chunk, ny) else matrix(0, n_chunk, ny)
}

wet_count <- new_store(TRUE)
dry_count <- new_store(TRUE)
wet_actual <- new_store()
wet_iwv <- new_store()
wet_z <- new_store()
wet_interaction <- new_store()
dry_actual <- new_store()
dry_iwv <- new_store()
dry_z <- new_store()
dry_interaction <- new_store()

for (local_start in seq.int(1L, n_chunk, by = BATCH)) {
  local_end <- min(n_chunk, local_start + BATCH - 1L)
  ii <- local_start:local_end
  message(
    "Chunk ", chunk_id, ": cells ",
    ids_chunk[local_start], "-", ids_chunk[local_end]
  )

  I <- Imat[ii, , drop = FALSE]
  Z <- Zmat[ii, , drop = FALSE]
  P <- Pmat[ii, , drop = FALSE]
  ai <- adjust_and_anomaly(I)
  az <- adjust_and_anomaly(Z)
  ap <- adjust_and_anomaly(P)

  covariance_climatology <- matrix(NA_real_, nrow(I), 12)
  product_anomaly <- ai$anomaly * az$anomaly
  for (mm in 1:12) {
    covariance_climatology[, mm] <- rowMeans(
      product_anomaly[, month == mm, drop = FALSE],
      na.rm = TRUE
    )
  }

  actual <- G * (
    ai$adjusted * az$adjusted -
      ai$climatology[, month] * az$climatology[, month] -
      covariance_climatology[, month]
  )
  iwv_channel <- G * ai$anomaly * az$climatology[, month]
  z_channel <- G * az$anomaly * ai$climatology[, month]
  interaction_channel <- G * (
    product_anomaly - covariance_climatology[, month]
  )

  q90 <- apply(ap$anomaly, 1, quantile, 0.9, na.rm = TRUE)
  q10 <- apply(ap$anomaly, 1, quantile, 0.1, na.rm = TRUE)
  wet_mask <- ap$anomaly > q90
  dry_mask <- ap$anomaly < q10

  wet_count[ii, ] <- annual_count(wet_mask)
  dry_count[ii, ] <- annual_count(dry_mask)
  wet_actual[ii, ] <- annual_sum(actual, wet_mask)
  wet_iwv[ii, ] <- annual_sum(iwv_channel, wet_mask)
  wet_z[ii, ] <- annual_sum(z_channel, wet_mask)
  wet_interaction[ii, ] <- annual_sum(interaction_channel, wet_mask)
  dry_actual[ii, ] <- annual_sum(actual, dry_mask)
  dry_iwv[ii, ] <- annual_sum(iwv_channel, dry_mask)
  dry_z[ii, ] <- annual_sum(z_channel, dry_mask)
  dry_interaction[ii, ] <- annual_sum(interaction_channel, dry_mask)

  rm(
    I, Z, P, ai, az, ap, covariance_climatology, product_anomaly,
    actual, iwv_channel, z_channel, interaction_channel,
    wet_mask, dry_mask
  )
  gc()
}

cell_all <- expand.grid(lon = lon, lat = lat)
cell <- cell_all[ids_chunk, , drop = FALSE]
cell$cell_id <- ids_chunk
cell$area_weight <- pmax(0, cos(cell$lat * pi / 180))
cell$raw_block <- paste(
  floor(((cell$lon + 180) %% 360) / 10),
  floor((cell$lat + 90) / 10),
  sep = "_"
)

chunk <- list(
  chunk_id = chunk_id,
  cell = cell,
  years = years,
  dates = dates,
  wet_count = wet_count,
  wet_actual = wet_actual,
  wet_iwv = wet_iwv,
  wet_z = wet_z,
  wet_interaction = wet_interaction,
  dry_count = dry_count,
  dry_actual = dry_actual,
  dry_iwv = dry_iwv,
  dry_z = dry_z,
  dry_interaction = dry_interaction
)
out_file <- file.path(
  OUT, sprintf("WETDRY_ANNUAL_CHUNK_%02d.rds", chunk_id)
)
saveRDS(chunk, out_file, compress = FALSE)
message("WETDRY_ANNUAL_CHUNK_COMPLETE=", out_file)
