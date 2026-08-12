options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))

ROOT <- "__WGPE_PROJECT_ROOT__"
CANONICAL <- file.path(ROOT, "CANONICAL_WGPE_TOP100_v2")
SRC <- file.path(CANONICAL, "FINAL_LARGE_DATA", "exact_layer_grid.csv.gz")
OUTDIR <- file.path(CANONICAL, "08_REPRODUCTION_OUTPUT", "Figure1_spatial_block_CI")
OUT <- file.path(OUTDIR, "FIG1_EXACT_SPATIAL_BLOCK_CI.csv")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(SRC)) stop("Missing exact grid-cell source: ", SRC)
if (requireNamespace("data.table", quietly = TRUE)) {
  d <- data.table::fread(SRC, data.table = FALSE, showProgress = TRUE)
} else {
  d <- read.csv(SRC, check.names = FALSE)
}

d$lon180 <- ifelse(d$lon >= 180, d$lon - 360, d$lon)
d$block_lon <- pmin(35L, pmax(0L, floor((d$lon180 + 180) / 10)))
d$block_lat <- pmin(17L, pmax(0L, floor((d$lat + 90) / 10)))
d$block_id <- sprintf("%02d_%02d", d$block_lon, d$block_lat)
d$area_weight <- pmax(0, cos(d$lat * pi / 180))

weighted_block_boot <- function(x, w, block, B = 2000L, seed = 20260722L) {
  ok <- is.finite(x) & is.finite(w) & w > 0 & !is.na(block)
  x <- x[ok]; w <- w[ok]; block <- block[ok]
  if (!length(x)) return(c(estimate = NA, lower95 = NA, upper95 = NA,
                           n_grid = 0, n_blocks = 0))
  num <- tapply(x * w, block, sum, na.rm = TRUE)
  den <- tapply(w, block, sum, na.rm = TRUE)
  ids <- intersect(names(num), names(den))
  num <- num[ids]; den <- den[ids]
  estimate <- sum(num) / sum(den)
  set.seed(seed)
  boot <- replicate(B, {
    ii <- sample.int(length(ids), length(ids), replace = TRUE)
    sum(num[ii]) / sum(den[ii])
  })
  q <- unname(quantile(boot, c(.025, .975), na.rm = TRUE, names = FALSE))
  c(estimate = estimate, lower95 = q[1], upper95 = q[2],
    n_grid = length(x), n_blocks = length(ids))
}

layer_defs <- list(
  total_contrib_1000_850 = "1000-850",
  total_contrib_850_700  = "850-700",
  total_contrib_700_500  = "700-500",
  total_contrib_500_300  = "500-300",
  total_contrib_300_100  = c("300-200", "200-100")
)

rows <- list()
for (metric in names(layer_defs)) {
  ss <- d[d$layer %in% layer_defs[[metric]], ]
  if (length(layer_defs[[metric]]) == 2L) {
    ss <- aggregate(cbind(exact_layer_total_m_yr, area_weight) ~ cell_id + lon180 + lat + block_id,
                    data = ss,
                    FUN = function(z) if (all(!is.finite(z))) NA_real_ else sum(z, na.rm = TRUE))
    ss$area_weight <- pmax(0, cos(ss$lat * pi / 180))
  }
  z <- weighted_block_boot(ss$exact_layer_total_m_yr, ss$area_weight, ss$block_id)
  rows[[length(rows) + 1L]] <- data.frame(metric = metric, t(z), stringsAsFactors = FALSE)
}

cells <- aggregate(cbind(mass_redistribution_m_yr,
                         height_coordinate_m_yr,
                         higher_order_interaction_m_yr) ~ cell_id + lon180 + lat + block_id,
                   data = d,
                   FUN = function(z) if (all(!is.finite(z))) NA_real_ else sum(z, na.rm = TRUE))
actual <- aggregate(actual_zbar_trend_m_yr ~ cell_id + lon180 + lat + block_id,
                    data = d, FUN = function(z) z[which(is.finite(z))[1]])
cells <- merge(cells, actual, by = c("cell_id", "lon180", "lat", "block_id"), all = TRUE)
cells$area_weight <- pmax(0, cos(cells$lat * pi / 180))

global_metrics <- c(
  sum_mass = "mass_redistribution_m_yr",
  sum_height = "height_coordinate_m_yr",
  residual = "higher_order_interaction_m_yr",
  actual = "actual_zbar_trend_m_yr"
)
for (metric in names(global_metrics)) {
  v <- cells[[global_metrics[[metric]]]]
  z <- weighted_block_boot(v, cells$area_weight, cells$block_id)
  rows[[length(rows) + 1L]] <- data.frame(metric = metric, t(z), stringsAsFactors = FALSE)
}

out <- do.call(rbind, rows)
out$estimate <- as.numeric(out$estimate)
out$lower95 <- as.numeric(out$lower95)
out$upper95 <- as.numeric(out$upper95)
out$n_grid <- as.integer(round(as.numeric(out$n_grid)))
out$n_blocks <- as.integer(round(as.numeric(out$n_blocks)))
out$method <- "paired 10x10-degree spatial-block bootstrap of exact grid-cell trends; 2000 replicates"
out$source_file <- SRC
write.csv(out, OUT, row.names = FALSE)
cat("Wrote ", OUT, "\n", sep = "")
