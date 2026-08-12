options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))
suppressPackageStartupMessages(library(ncdf4))

ROOT <- "__WGPE_PROJECT_ROOT__"
OUT <- file.path(ROOT, "output_TOP100_revision", "derived_data", "figure_source_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
lon <- -180:179
lat <- -90:90

grid_matrix <- function(d, value) {
  m <- matrix(NA_real_, length(lon), length(lat))
  dl <- ifelse(d$lon >= 180, d$lon - 360, d$lon)
  ii <- match(dl, lon)
  jj <- match(d$lat, lat)
  ok <- is.finite(ii) & is.finite(jj)
  m[cbind(ii[ok], jj[ok])] <- d[[value]][ok]
  m
}

fig2a_dir <- file.path(ROOT, "output_TOP100_revision", "results", "precip_wetdry",
                       "Fig2a_deseasoned_only_TOP100")
pg <- read.csv(file.path(fig2a_dir, "Fig2a_WGPE_precip_correlation_grid_TOP100.csv"),
               check.names = FALSE)
ps <- read.csv(file.path(fig2a_dir, "Fig2a_WGPE_precip_correlation_summary_TOP100.csv"),
               check.names = FALSE)
pg <- pg[pg$method == "deseasonalized_only", , drop = FALSE]
ps <- ps[ps$method == "deseasonalized_only", , drop = FALSE]
ps <- ps[ps$height_reference == "absolute", , drop = FALSE]
ps$variable <- "wgpe"
ps$lag_months_predictor_leads_precip <- 0
ps$significant_valid_area_pct <- ps$total_significant_area_pct
ps$significant_grid_pct <- NA_real_
for (i in seq_len(nrow(ps))) {
  in_region <- switch(ps$region[i], Global = rep(TRUE, nrow(pg)),
                      `S high latitude` = pg$lat < -66.5,
                      `S temperate` = pg$lat >= -66.5 & pg$lat < -23.5,
                      Tropics = pg$lat >= -23.5 & pg$lat <= 23.5,
                      `N temperate` = pg$lat > 23.5 & pg$lat <= 66.5,
                      `N high latitude` = pg$lat > 66.5)
  ps$significant_grid_pct[i] <- 100 * sum(is.finite(pg$q[in_region]) & pg$q[in_region] < .05) /
    sum(is.finite(pg$q[in_region]))
}
precip <- list(
  precipitation = list(wgpe = list(lag0 = list(
    r = grid_matrix(pg, "r"), p = grid_matrix(pg, "p"),
    q = grid_matrix(pg, "q"), neff = grid_matrix(pg, "neff")
  ))),
  summary = ps
)

enso_nc_path <- file.path(ROOT, "output_TOP100_revision", "results", "ENSO",
                          "ONI_standard_TOP100", "03_absolute_results",
                          "ENSO_ONIstandard_absolute_maps.nc")
n <- nc_open(enso_nc_path)
lon_nc <- n$dim$lon$vals
lat_nc <- n$dim$lat$vals
read_map <- function(variable) {
  a <- ncvar_get(n, variable)
  lon2 <- ifelse(lon_nc >= 180, lon_nc - 360, lon_nc)
  a[order(lon2), order(lat_nc), drop = FALSE]
}
enso_map <- read_map("ENSO_difference_actual")
elnino_map <- read_map("ElNino_actual")
lanina_map <- read_map("LaNina_actual")
nc_close(n)

wet <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "precip_wetdry",
                          "exact_TOP100", "WET_DRY_EXACT_SUMMARY_TOP100.csv"), check.names = FALSE)
ens <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "ENSO",
                          "ONI_standard_TOP100", "03_absolute_results",
                          "ENSO_ONIstandard_absolute_channel_decomposition.csv"), check.names = FALSE)
ens_global <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "ENSO",
                                 "ONI_standard_TOP100", "03_absolute_results",
                                 "ENSO_ONIstandard_absolute_global_summary.csv"), check.names = FALSE)

wet_event <- wet[wet$metric == "actual", c("event", "region", "estimate")]
wet_event$event <- ifelse(wet_event$event == "Wet", "Extreme_precipitation", "Drought")
names(wet_event)[3] <- "anomaly"
wet_event$variable <- "wgpe"
enso_event <- data.frame(event = "ENSO_difference", region = ens_global$region[ens_global$metric == "actual"],
                         anomaly = ens_global$estimate[ens_global$metric == "actual"], variable = "wgpe")
event_summary <- rbind(wet_event[, c("event", "region", "anomaly", "variable")], enso_event)

wet_dec <- wet[wet$metric %in% c("IWV_channel", "z_channel", "interaction_channel"), ]
wet_wide <- reshape(wet_dec[, c("event", "region", "metric", "estimate")],
                    idvar = c("event", "region"), timevar = "metric", direction = "wide")
names(wet_wide) <- sub("estimate\\.", "", names(wet_wide))
wet_wide$event <- ifelse(wet_wide$event == "Wet", "Extreme_precipitation", "Drought")
names(wet_wide)[names(wet_wide) == "IWV_channel"] <- "iwv_channel_J_m2"
names(wet_wide)[names(wet_wide) == "z_channel"] <- "z_channel_J_m2"
names(wet_wide)[names(wet_wide) == "interaction_channel"] <- "interaction_channel_J_m2"

event <- list(
  event_maps = list(wgpe = list(ENSO_difference = enso_map,
                                El_Nino = elnino_map, La_Nina = lanina_map)),
  summary = event_summary,
  decomposition_summary = wet_wide
)

tp_grid <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "temporal_precedence",
                              "Precip_temporal_precedence_grid_FINAL_v3.csv"), check.names = FALSE)
tp_summary <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "temporal_precedence",
                                 "Precip_temporal_precedence_FINAL_v3.csv"), check.names = FALSE)
enso_dec <- ens[ens$metric %in% c("IWV_channel", "z_channel", "interaction_channel"),
                c("region", "metric", "estimate", "lower95", "upper95")]
enso_wide <- reshape(enso_dec[, c("region", "metric", "estimate")], idvar = "region",
                     timevar = "metric", direction = "wide")
names(enso_wide) <- sub("estimate\\.", "", names(enso_wide))
names(enso_wide)[names(enso_wide) == "IWV_channel"] <- "IWV_channel_J_m2"
names(enso_wide)[names(enso_wide) == "z_channel"] <- "z_channel_J_m2"
names(enso_wide)[names(enso_wide) == "interaction_channel"] <- "interaction_channel_J_m2"

tp_long <- do.call(rbind, lapply(1:4, function(k) {
  data.frame(lon = tp_grid$lon, lat = tp_grid$lat, target = "precipitation",
             predictor = "W-GPE", lag = k, valid_TRUE_FALSE = !is.na(tp_grid[[paste0("lag", k, "_class")]]),
             class_q = tp_grid[[paste0("lag", k, "_class")]])
}))

# Formal TOP100 lag-3 uncertainty for the map-side summaries.
td <- tp_long[tp_long$lag == 3 & tp_long$valid_TRUE_FALSE, ]
classes <- c("forward-only", "reverse-only", "bidirectional", "none")
block <- interaction(floor((td$lon + 180) / 10), floor((td$lat + 90) / 10), drop = TRUE)
bi <- as.integer(block)
nb <- nlevels(block)
aw <- cos(td$lat * pi / 180)
set.seed(20260714)
B <- 1000L
boot_fraction <- matrix(NA_real_, B, length(classes))
for (b in 1:B) {
  cnt <- tabulate(sample.int(nb, nb, replace = TRUE), nb)
  ww <- aw * cnt[bi]
  boot_fraction[b, ] <- vapply(classes, function(cl) 100 * sum(ww[td$class_q == cl]) / sum(ww), numeric(1))
}
temp_unc <- data.frame(
  variable = classes,
  estimate = vapply(classes, function(cl) 100 * sum(aw[td$class_q == cl]) / sum(aw), numeric(1)),
  lower95 = apply(boot_fraction, 2, quantile, .025, na.rm = TRUE),
  upper95 = apply(boot_fraction, 2, quantile, .975, na.rm = TRUE),
  height_reference = "absolute", target = "precipitation", lag = 3, region = "Global",
  uncertainty_method = "spatial 10-degree block bootstrap of cos(latitude)-weighted class fractions; 1000 resamples"
)

br <- seq(-90, 90, by = 5)
temp_prof <- do.call(rbind, lapply(seq_len(length(br) - 1), function(i) {
  lo <- br[i]; hi <- br[i + 1]
  sub <- td[td$lat >= lo & if (i == length(br) - 1) td$lat <= hi else td$lat < hi, ]
  if (!nrow(sub)) return(NULL)
  lon_block <- factor(floor((sub$lon + 180) / 10))
  li <- as.integer(lon_block); nl <- nlevels(lon_block); sw <- cos(sub$lat * pi / 180)
  net_fun <- function(ww) 100 * (sum(ww[sub$class_q == "forward-only"]) -
                                  sum(ww[sub$class_q == "reverse-only"])) / sum(ww)
  set.seed(20260714L + i)
  sims <- replicate(B, {
    cnt <- tabulate(sample.int(nl, nl, replace = TRUE), nl)
    net_fun(sw * cnt[li])
  })
  data.frame(variable = "net", estimate = net_fun(sw), lower95 = quantile(sims, .025),
             upper95 = quantile(sims, .975), height_reference = "absolute",
             target = "precipitation", lag = 3, latitude_bin_start = lo,
             latitude_bin_end = hi, latitude_bin_mid = (lo + hi) / 2,
             uncertainty_method = "longitude 10-degree block bootstrap; 1000 resamples")
}))

profile_boot <- function(field, variable_name) {
  do.call(rbind, lapply(seq_len(length(br) - 1), function(i) {
    lo <- br[i]; hi <- br[i + 1]
    jj <- which(lat >= lo & if (i == length(br) - 1) lat <= hi else lat < hi)
    sub <- field[, jj, drop = FALSE]
    lat_sub <- lat[jj]
    ww <- matrix(rep(cos(lat_sub * pi / 180), each = length(lon)), nrow = length(lon))
    lon_block <- factor(floor((lon + 180) / 10)); li <- as.integer(lon_block); nl <- nlevels(lon_block)
    estimate <- weighted.mean(as.vector(sub), as.vector(ww), na.rm = TRUE)
    set.seed(20260800L + i)
    sims <- replicate(B, {
      cnt <- tabulate(sample.int(nl, nl, replace = TRUE), nl)
      bw <- ww * cnt[li]
      weighted.mean(as.vector(sub), as.vector(bw), na.rm = TRUE)
    })
    data.frame(variable = "value", estimate = estimate, lower95 = quantile(sims, .025, na.rm = TRUE),
               upper95 = quantile(sims, .975, na.rm = TRUE), height_reference = "absolute",
               variable_name = variable_name, latitude_bin_start = lo, latitude_bin_end = hi,
               latitude_bin_mid = (lo + hi) / 2,
               uncertainty_method = "longitude 10-degree block bootstrap; 1000 resamples")
  }))
}
precip_prof_unc <- profile_boot(precip$precipitation$wgpe$lag0$r, "precip_r")
enso_prof_raw <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "ENSO",
                                    "ONI_standard_TOP100", "03_absolute_results",
                                    "ENSO_ONIstandard_absolute_lat5_profile.csv"), check.names = FALSE)
enso_prof_raw <- enso_prof_raw[enso_prof_raw$event_type == "ENSO_difference" &
                                enso_prof_raw$metric == "actual", ]
enso_prof_unc <- data.frame(variable = "value", estimate = enso_prof_raw$estimate,
                            lower95 = enso_prof_raw$lower95, upper95 = enso_prof_raw$upper95,
                            height_reference = "absolute", variable_name = "ENSO_WGPE",
                            latitude_bin_start = enso_prof_raw$lat_mid - 2.5,
                            latitude_bin_end = enso_prof_raw$lat_mid + 2.5,
                            latitude_bin_mid = enso_prof_raw$lat_mid,
                            uncertainty_method = "ONI event-block bootstrap")
profile_unc <- rbind(precip_prof_unc, enso_prof_unc)

write.csv(enso_wide, file.path(OUT, "Figure2_ENSO_standard_band_decomposition_TOP100.csv"), row.names = FALSE)
write.csv(temp_unc, file.path(OUT, "Figure2_temporal_fraction_uncertainty_TOP100.csv"), row.names = FALSE)
write.csv(temp_prof, file.path(OUT, "Figure2_temporal_profile_uncertainty_TOP100.csv"), row.names = FALSE)
write.csv(profile_unc, file.path(OUT, "Figure2_latitude_profile_uncertainty_TOP100.csv"), row.names = FALSE)
saveRDS(list(precip = precip, event = event, fdr_grid = tp_long, temporal_summary = tp_summary,
             wet_uncertainty = wet, enso_uncertainty = ens,
             enso_decomposition = enso_wide, temporal_fraction_uncertainty = temp_unc,
             temporal_profile_uncertainty = temp_prof, latitude_profile_uncertainty = profile_unc,
             enso_nc = enso_nc_path),
        file.path(OUT, "Figure2_TOP100_adapter.rds"))
cat("Saved Figure2_TOP100_adapter.rds\n")
cat("Precip global mean r:", ps$mean_r_area[ps$region == "Global"], "\n")
cat("ENSO global anomaly:", ens_global$estimate[ens_global$metric == "actual"], "\n")
