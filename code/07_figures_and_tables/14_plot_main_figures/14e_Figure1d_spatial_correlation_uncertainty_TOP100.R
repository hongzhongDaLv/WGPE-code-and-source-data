options(stringsAsFactors = FALSE, warn = 1)
ROOT <- "__WGPE_PROJECT_ROOT__"
source(file.path(ROOT, "scripts_TOP100_revision", "18_coordinate_harmonization", "coordinate_harmonization.R"))
OUT <- Sys.getenv("TOP100_CLOSURE_FIG1D_OUT", unset = file.path(ROOT, "output_TOP100_revision", "derived_data", "figure_source_data"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
B <- 1000L

e <- new.env(parent = emptyenv())
load(file.path(ROOT, "output_TOP100_revision", "derived_data", "ERA5_TOP100",
               "global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData"), envir = e)
lon <- as.numeric(e$lon)
lat <- as.numeric(e$lat)
dates <- as.Date(e$all_dates)
z <- e$zbar_st
iwv <- e$iwv_st
rm(e)
gc()

single <- readRDS(file.path(ROOT, "output_attribution_minimal",
                            "ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds"))

# The TOP100 pressure-level core is stored on 0--359 degrees, whereas the
# single-level ERA5 bundle is stored on -179--180 degrees.  Pairing the raw
# array positions shifts the single-level fields by about 180 degrees.  Reorder
# every external field to the TOP100 coordinates before any grid-cell statistic.
align_single <- function(a) {
  reorder_exact_grid(a, single$lon, single$lat, lon, lat)
}

decimal_year <- as.integer(format(dates, "%Y")) + (as.integer(format(dates, "%m")) - .5) / 12
month <- as.integer(format(dates, "%m"))
t0 <- decimal_year - mean(decimal_year)

trend_complete <- function(arr) {
  m <- matrix(as.numeric(arr), nrow = length(lon) * length(lat))
  for (mm in 1:12) {
    jj <- which(month == mm)
    m[, jj] <- m[, jj, drop = FALSE] - rowMeans(m[, jj, drop = FALSE], na.rm = TRUE)
  }
  ok <- rowSums(is.finite(m)) == ncol(m)
  out <- rep(NA_real_, nrow(m))
  out[ok] <- as.numeric(m[ok, , drop = FALSE] %*% t0) / sum(t0^2)
  matrix(out, nrow = length(lon), ncol = length(lat))
}

ztrend <- trend_complete(z)
rm(z)
gc()
predictors <- list(
  IWV = trend_complete(iwv),
  T2m = trend_complete(align_single(single$t2m_K)),
  Precipitation = trend_complete(align_single(single$tp_native_m)),
  SSRD = trend_complete(align_single(single$ssrd_native_J_m2))
)
rm(iwv, single)
gc()

cell <- expand.grid(lon = ifelse(lon > 180, lon - 360, lon), lat = lat)
block <- interaction(floor((cell$lon + 180) / 10), floor((cell$lat + 90) / 10), drop = TRUE)

weighted_cor <- function(x, y, ww) {
  sw <- sum(ww)
  mx <- sum(ww * x) / sw
  my <- sum(ww * y) / sw
  vx <- sum(ww * (x - mx)^2) / sw
  vy <- sum(ww * (y - my)^2) / sw
  cv <- sum(ww * (x - mx) * (y - my)) / sw
  cv / sqrt(vx * vy)
}

bootstrap_one <- function(y, variable, seed) {
  x <- as.vector(ztrend)
  y <- as.vector(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  bb <- droplevels(block[ok])
  bi <- as.integer(bb)
  nb <- nlevels(bb)
  estimate <- cor(x, y)
  set.seed(seed)
  sims <- numeric(B)
  for (b in seq_len(B)) {
    cnt <- tabulate(sample.int(nb, nb, replace = TRUE), nb)
    sims[b] <- weighted_cor(x, y, cnt[bi])
  }
  data.frame(variable = variable, estimate = estimate,
             lower95 = unname(quantile(sims, .025, na.rm = TRUE)),
             upper95 = unname(quantile(sims, .975, na.rm = TRUE)),
             n_grid = length(x), n_blocks = nb, n_boot = B,
             uncertainty_method = "10 degree x 10 degree spatial-block bootstrap")
}

out <- do.call(rbind, lapply(seq_along(predictors), function(k) {
  bootstrap_one(predictors[[k]], names(predictors)[k], 20260714L + k)
}))
locked <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "core",
                             "Trend_pattern_correlations_FINAL_v3.csv"), check.names = FALSE)
locked <- locked[locked$height_reference == "absolute", c("variable", "pearson_r")]
out$locked_estimate <- locked$pearson_r[match(out$variable, locked$variable)]
out$abs_difference_from_locked <- abs(out$estimate - out$locked_estimate)
out$longitude_alignment <- ifelse(out$variable == "IWV", "same TOP100 grid",
                                  "single-level -179--180 reordered one-to-one to TOP100 0--359")
write.csv(out, file.path(OUT, "Figure1d_spatial_block_bootstrap_LONGITUDE_FIXED_TOP100.csv"), row.names = FALSE)
write.csv(out[, c("variable", "estimate")],
          file.path(OUT, "Trend_pattern_correlations_LONGITUDE_FIXED_TOP100.csv"), row.names = FALSE)
print(out[, c("variable", "estimate", "locked_estimate", "lower95", "upper95")])
