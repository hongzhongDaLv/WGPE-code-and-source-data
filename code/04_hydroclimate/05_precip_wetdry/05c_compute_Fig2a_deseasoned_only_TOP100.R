# Purpose: Recompute the Figure 2a W-GPE--precipitation association after
#          removing the calendar-month climatology only (no linear detrending).
#          A correctly aligned detrended result is retained as sensitivity.
# Inputs:  formal ERA5 TOP100 absolute-height core and ERA5 precipitation.
# Outputs: output_TOP100_revision/results/precip_wetdry/
#          Fig2a_deseasoned_only_TOP100/*
# Important: the TOP100 core uses 0--359 degree longitude whereas the single-
#            level file uses -179--180 degrees.  The TOP100 field is explicitly
#            reordered before any grid-cell pairing.

options(stringsAsFactors = FALSE, warn = 1)

ROOT <- "__WGPE_PROJECT_ROOT__"
source(file.path(ROOT, "scripts_TOP100_revision", "18_coordinate_harmonization", "coordinate_harmonization.R"))
CORE <- file.path(ROOT, "output_TOP100_revision", "derived_data", "ERA5_TOP100",
                  "global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData")
SINGLE <- file.path(ROOT, "output_attribution_minimal",
                    "ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds")
OUT <- Sys.getenv("TOP100_CLOSURE_FIG2A_OUT", unset = file.path(ROOT, "output_TOP100_revision", "results", "precip_wetdry",
                 "Fig2a_deseasoned_only_TOP100"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(CORE), file.exists(SINGLE))

message("Loading TOP100 W-GPE and precipitation")
e <- new.env(parent = emptyenv())
load(CORE, envir = e)
lon_core <- as.numeric(e$lon)
lat_core <- as.numeric(e$lat)
dates <- as.Date(e$all_dates)
wgpe <- e$wgpe_st
rm(e); gc()

s <- readRDS(SINGLE)
lon_ref <- as.numeric(s$lon)
lat_ref <- as.numeric(s$lat)
dates_ref <- as.Date(s$dates)
tp <- s$tp_native_m
rm(s); gc()

stopifnot(identical(as.character(dates), as.character(dates_ref)))
stopifnot(identical(lat_core, lat_ref))
mapping <- exact_grid_mapping(lon_core, lat_core, lon_ref, lat_ref)
idx_lon <- mapping$lon_index

nlon <- length(lon_ref)
nlat <- length(lat_ref)
nt <- length(dates)
ncell <- nlon * nlat
stopifnot(identical(dim(wgpe), c(length(lon_core), nlat, nt)))
stopifnot(identical(dim(tp), c(nlon, nlat, nt)))

# Matrix row order is longitude-fastest, then latitude.  row_core maps each
# -179--180 reference-grid cell to its physically identical 0--359 TOP100 cell.
row_core <- rep(idx_lon, times = nlat) +
  (rep(seq_len(nlat), each = nlon) - 1L) * length(lon_core)
dim(wgpe) <- c(length(lon_core) * nlat, nt)
dim(tp) <- c(ncell, nt)

month <- as.integer(format(dates, "%m"))
time_centered <- seq_len(nt) - mean(seq_len(nt))
time_ss <- sum(time_centered^2)

deseason_chunk <- function(x) {
  out <- x
  for (m in 1:12) {
    j <- month == m
    out[, j] <- sweep(x[, j, drop = FALSE], 1,
                      rowMeans(x[, j, drop = FALSE], na.rm = TRUE), "-")
  }
  out
}

detrend_chunk <- function(x) {
  out <- x
  complete <- rowSums(!is.finite(x)) == 0L
  if (any(complete)) {
    slope <- as.numeric(x[complete, , drop = FALSE] %*% time_centered) / time_ss
    out[complete, ] <- x[complete, , drop = FALSE] - slope %o% time_centered
  }
  if (any(!complete)) {
    for (i in which(!complete)) {
      ok <- is.finite(x[i, ])
      out[i, ] <- NA_real_
      if (sum(ok) >= 36L) out[i, ok] <- residuals(lm(x[i, ok] ~ time_centered[ok]))
    }
  }
  out
}

row_cor_ar1 <- function(x, y) {
  nr <- nrow(x)
  r <- p <- neff <- rep(NA_real_, nr)
  for (i in seq_len(nr)) {
    ok <- is.finite(x[i, ]) & is.finite(y[i, ])
    if (sum(ok) < 36L) next
    xx <- x[i, ok]; yy <- y[i, ok]
    rr <- suppressWarnings(cor(xx, yy))
    if (!is.finite(rr) || abs(rr) >= 1) next
    nn <- length(xx)
    a1 <- suppressWarnings(cor(xx[-nn], xx[-1]))
    b1 <- suppressWarnings(cor(yy[-nn], yy[-1]))
    if (!is.finite(a1)) a1 <- 0
    if (!is.finite(b1)) b1 <- 0
    a1 <- max(-0.99, min(0.99, a1))
    b1 <- max(-0.99, min(0.99, b1))
    ne <- max(4, min(nn, nn * (1 - a1 * b1) / (1 + a1 * b1)))
    stat <- rr * sqrt((ne - 2) / max(1e-12, 1 - rr^2))
    r[i] <- rr
    p[i] <- 2 * pt(-abs(stat), df = ne - 2)
    neff[i] <- ne
  }
  list(r = r, p = p, neff = neff)
}

methods <- c("deseasonalized_only", "deseasonalized_and_detrended")
res <- lapply(methods, function(z) {
  list(r = rep(NA_real_, ncell), p = rep(NA_real_, ncell),
       neff = rep(NA_real_, ncell))
})
names(res) <- methods

chunk_size <- 1000L
starts <- seq.int(1L, ncell, by = chunk_size)
for (z in seq_along(starts)) {
  i1 <- starts[z]
  ii <- i1:min(ncell, i1 + chunk_size - 1L)
  message("Chunk ", z, "/", length(starts), " cells ", min(ii), "--", max(ii))
  x <- deseason_chunk(wgpe[row_core[ii], , drop = FALSE])
  y <- deseason_chunk(tp[ii, , drop = FALSE])

  a <- row_cor_ar1(x, y)
  res$deseasonalized_only$r[ii] <- a$r
  res$deseasonalized_only$p[ii] <- a$p
  res$deseasonalized_only$neff[ii] <- a$neff

  xd <- detrend_chunk(x)
  yd <- detrend_chunk(y)
  b <- row_cor_ar1(xd, yd)
  res$deseasonalized_and_detrended$r[ii] <- b$r
  res$deseasonalized_and_detrended$p[ii] <- b$p
  res$deseasonalized_and_detrended$neff[ii] <- b$neff
  rm(x, y, xd, yd, a, b); gc(FALSE)
}

cell <- expand.grid(lon = lon_ref, lat = lat_ref)
cell$weight <- cos(cell$lat * pi / 180)
band_levels <- c("Global", "S high latitude", "S temperate", "Tropics",
                 "N temperate", "N high latitude")
band_mask <- function(region) {
  if (region == "Global") rep(TRUE, ncell)
  else if (region == "S high latitude") cell$lat < -66.5
  else if (region == "S temperate") cell$lat >= -66.5 & cell$lat < -23.5
  else if (region == "Tropics") cell$lat >= -23.5 & cell$lat <= 23.5
  else if (region == "N temperate") cell$lat > 23.5 & cell$lat <= 66.5
  else cell$lat > 66.5
}

grid_rows <- list()
summary_rows <- list()
k <- 0L
for (method in methods) {
  d <- res[[method]]
  q <- rep(NA_real_, ncell)
  okp <- is.finite(d$p)
  q[okp] <- p.adjust(d$p[okp], method = "BH")
  grid_rows[[method]] <- data.frame(
    lon = cell$lon, lat = cell$lat, r = d$r, p = d$p, q = q,
    neff = d$neff, height_reference = "absolute", method = method,
    period = "1979-2024", top_hPa = 100,
    longitude_alignment = "TOP100 0-359 reordered one-to-one to -179-180"
  )
  for (region in band_levels) {
    use <- band_mask(region) & is.finite(d$r)
    sig <- use & is.finite(q) & q <= 0.05
    w <- cell$weight
    k <- k + 1L
    summary_rows[[k]] <- data.frame(
      height_reference = "absolute", method = method, region = region,
      mean_r_area = weighted.mean(d$r[use], w[use]),
      positive_significant_area_pct = 100 * sum(w[sig & d$r > 0]) / sum(w[use]),
      negative_significant_area_pct = 100 * sum(w[sig & d$r < 0]) / sum(w[use]),
      total_significant_area_pct = 100 * sum(w[sig]) / sum(w[use]),
      n_grid = sum(use), period = "1979-2024", top_hPa = 100,
      area_weighting = "cos(latitude)",
      significance = "AR1 effective sample size; BH-FDR q <= 0.05"
    )
  }
}

grid_out <- do.call(rbind, grid_rows)
summary_out <- do.call(rbind, summary_rows)
write.csv(grid_out, file.path(OUT, "Fig2a_WGPE_precip_correlation_grid_TOP100.csv"), row.names = FALSE)
write.csv(summary_out, file.path(OUT, "Fig2a_WGPE_precip_correlation_summary_TOP100.csv"), row.names = FALSE)

gm <- summary_out[summary_out$region == "Global", ]
qa <- c(
  "# Fig. 2a deseasonalized-only TOP100 QA",
  "",
  "- Main method: remove the 1979-2024 calendar-month climatology only; no linear detrending.",
  "- Sensitivity method: remove calendar-month climatology and then linear trend.",
  "- TOP100 longitude was reordered one-to-one from 0-359 to the precipitation -179-180 grid before pairing.",
  paste0("- Longitude mapping unique: ", length(unique(idx_lon)), "/", length(idx_lon), "."),
  paste0("- Latitude coordinates identical: ", identical(lat_core, lat_ref), "."),
  paste0("- Dates identical: ", identical(as.character(dates), as.character(dates_ref)), "."),
  paste0("- Global mean r, deseasonalized only: ", sprintf("%+.6f", gm$mean_r_area[gm$method == "deseasonalized_only"]), "."),
  paste0("- Global mean r, detrended sensitivity: ", sprintf("%+.6f", gm$mean_r_area[gm$method == "deseasonalized_and_detrended"]), "."),
  "- Existing detrended/misaligned files were not overwritten.",
  "- No DOCX, SI, workbook, or figure was modified."
)
writeLines(qa, file.path(OUT, "Fig2a_deseasonalized_only_TOP100_QA.md"), useBytes = TRUE)
message("Wrote Fig. 2a corrected source data to: ", OUT)
