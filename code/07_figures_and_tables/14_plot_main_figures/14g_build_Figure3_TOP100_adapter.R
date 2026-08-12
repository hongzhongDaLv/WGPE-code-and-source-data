options(stringsAsFactors = FALSE, warn = 1)
ROOT <- "__WGPE_PROJECT_ROOT__"
OUT <- file.path(ROOT, "output_TOP100_revision", "derived_data", "figure_source_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
B <- 1000L

grid <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "ecology",
                           "Ecology_partial_correlations_grid_FINAL_v3.csv"), check.names = FALSE)
grid <- grid[grid$height_reference == "absolute" & grid$model == "Model3" &
               grid$predictor %in% c("W-GPE", "<z>"), ]
summary0 <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "ecology",
                               "Ecology_partial_correlations_summary_FINAL_v3.csv"), check.names = FALSE)
summary0 <- summary0[summary0$height_reference == "absolute" & summary0$model == "Model3" &
                       summary0$predictor %in% c("W-GPE", "<z>"), ]

lon <- seq(-179.5, 179.5, 1)
lat_raw <- seq(89.5, -89.5, -1)
to_matrix <- function(d, value) {
  m <- matrix(NA_real_, length(lon), length(lat_raw))
  ii <- match(d$lon, lon); jj <- match(d$lat, lat_raw)
  ok <- is.finite(ii) & is.finite(jj)
  m[cbind(ii[ok], jj[ok])] <- d[[value]][ok]
  m
}
pred_key <- c("W-GPE" = "wgpe", "<z>" = "zbar")
eco <- list(results = list(), summary = NULL)
for (target in c("GPP", "LAI")) {
  eco$results[[target]] <- list(r = list(), p = list(), q = list())
  for (predictor in names(pred_key)) {
    d <- grid[grid$target == target & grid$predictor == predictor, ]
    key <- pred_key[[predictor]]
    eco$results[[target]]$r[[key]] <- to_matrix(d, "partial_r")
    eco$results[[target]]$p[[key]] <- to_matrix(d, "p_value")
    eco$results[[target]]$q[[key]] <- to_matrix(d, "q_value")
  }
}

summary <- summary0
summary$predictor <- unname(pred_key[summary$predictor])
summary$mean_r_area <- summary$mean_partial_r_area
summary$significant_valid_area_pct <- summary$positive_significant_area_pct + summary$negative_significant_area_pct
summary$significant_grid_pct <- NA_real_
for (i in seq_len(nrow(summary))) {
  d <- grid[grid$target == summary$target[i] &
              pred_key[grid$predictor] == summary$predictor[i], ]
  summary$significant_grid_pct[i] <- 100 * sum(is.finite(d$q_value) & d$q_value < .05) /
    sum(is.finite(d$q_value))
}
eco$summary <- summary

# Formal map-side latitude-profile uncertainty.
profile_unc <- list()
u <- 0L
br <- seq(-90, 90, 5)
for (target in c("GPP", "LAI")) for (predictor in c("W-GPE", "<z>")) {
  d0 <- grid[grid$target == target & grid$predictor == predictor, ]
  for (i in seq_len(length(br) - 1)) {
    lo <- br[i]; hi <- br[i + 1]
    d <- d0[d0$lat >= lo & if (i == length(br) - 1) d0$lat <= hi else d0$lat < hi, ]
    if (!nrow(d)) next
    lb <- factor(floor((d$lon + 180) / 10)); li <- as.integer(lb); nl <- nlevels(lb)
    aw <- cos(d$lat * pi / 180)
    est <- weighted.mean(d$partial_r, aw, na.rm = TRUE)
    set.seed(20260714L + u + i)
    sims <- replicate(B, {
      cnt <- tabulate(sample.int(nl, nl, replace = TRUE), nl)
      weighted.mean(d$partial_r, aw * cnt[li], na.rm = TRUE)
    })
    u <- u + 1L
    profile_unc[[u]] <- data.frame(height_reference = "absolute", target = target,
                                   predictor = predictor, latitude_bin_start = lo,
                                   latitude_bin_end = hi, latitude_bin_mid = (lo + hi) / 2,
                                   estimate = est, lower95 = quantile(sims, .025, na.rm = TRUE),
                                   upper95 = quantile(sims, .975, na.rm = TRUE),
                                   uncertainty_method = "longitude 10-degree block bootstrap; 1000 resamples")
  }
}
U_EPROF <- do.call(rbind, profile_unc)

# Region-mean block-bootstrap CIs used to mark heatmap cells.
region_rule <- function(lat, region) {
  switch(region, Global = rep(TRUE, length(lat)),
         `S high latitude` = lat < -66.5,
         `S temperate` = lat >= -66.5 & lat < -23.5,
         Tropics = lat >= -23.5 & lat <= 23.5,
         `N temperate` = lat > 23.5 & lat <= 66.5,
         `N high latitude` = lat > 66.5)
}
regions <- c("Global", "S high latitude", "S temperate", "Tropics", "N temperate", "N high latitude")
esum <- list(); z <- 0L
for (target in c("GPP", "LAI")) for (predictor in c("W-GPE", "<z>")) for (region in regions) {
  d <- grid[grid$target == target & grid$predictor == predictor, ]
  d <- d[region_rule(d$lat, region), ]
  if (!nrow(d)) next
  block <- factor(interaction(floor((d$lon + 180) / 10), floor((d$lat + 90) / 10), drop = TRUE))
  bi <- as.integer(block); nb <- nlevels(block); aw <- cos(d$lat * pi / 180)
  est <- weighted.mean(d$partial_r, aw, na.rm = TRUE)
  set.seed(20261000L + z)
  sims <- replicate(B, {
    cnt <- tabulate(sample.int(nb, nb, replace = TRUE), nb)
    weighted.mean(d$partial_r, aw * cnt[bi], na.rm = TRUE)
  })
  z <- z + 1L
  esum[[z]] <- data.frame(height_reference = "absolute", target = target,
                          predictor = predictor, region = region, estimate = est,
                          lower95 = quantile(sims, .025, na.rm = TRUE),
                          upper95 = quantile(sims, .975, na.rm = TRUE), n_grid = nrow(d),
                          uncertainty_method = "spatial 10-degree block bootstrap; 1000 resamples")
}
U_ESUM <- do.call(rbind, esum)

# Conditional Model-3 temporal precedence, formal TOP100 branch.
tp <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "temporal_precedence",
                         "TOP100_conditional", "conditional_temporal_precedence_absolute_TOP100.csv"),
               check.names = FALSE)
tp$valid_TRUE_FALSE <- is.finite(tp$q_forward) & is.finite(tp$q_reverse)
tp$class_q <- tp$class
fdr_grid <- tp[, c("lon", "lat", "target", "predictor", "lag", "valid_TRUE_FALSE", "class_q")]
ts <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "temporal_precedence",
                         "TOP100_conditional", "conditional_temporal_precedence_summary_TOP100.csv"),
               check.names = FALSE)
ts <- ts[ts$height_reference == "absolute", ]
wide <- reshape(ts[, c("target", "predictor", "lag", "class", "area_fraction_pct")],
                idvar = c("target", "predictor", "lag"), timevar = "class", direction = "wide")
names(wide) <- sub("area_fraction_pct\\.", "", names(wide))
wide$region <- "Global"
wide$forward_only_grid_pct_q <- wide$`forward-only`
wide$reverse_only_grid_pct_q <- wide$`reverse-only`
wide$bidirectional_grid_pct_q <- wide$bidirectional
wide$none_grid_pct_q <- wide$none
fdr_summary <- wide

U_ETEMP <- data.frame(height_reference = character(), target = character(), predictor = character(),
                      lag = integer(), region = character(), class = character(), estimate = numeric(),
                      lower95 = numeric(), upper95 = numeric())

write.csv(U_EPROF, file.path(OUT, "Figure3_ecology_profile_uncertainty_TOP100.csv"), row.names = FALSE)
write.csv(U_ESUM, file.path(OUT, "Figure3_ecology_summary_uncertainty_TOP100.csv"), row.names = FALSE)
saveRDS(list(eco = eco, fdr_grid = fdr_grid, fdr_summary = fdr_summary,
             U_EPROF = U_EPROF, U_ESUM = U_ESUM, U_ETEMP = U_ETEMP),
        file.path(OUT, "Figure3_TOP100_adapter.rds"))
cat("Saved Figure3_TOP100_adapter.rds\n")
