options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))
suppressPackageStartupMessages(library(ncdf4))

ROOT <- "__WGPE_PROJECT_ROOT__"
OUT <- file.path(ROOT, "output_TOP100_revision", "derived_data", "figure_source_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "R", "figure_style_FIXED_FULL.R"), encoding = "UTF-8")
source(file.path(ROOT, "R", "figure_data_FIXED_FULL.R"), encoding = "UTF-8")

read_nc_slice <- function(nc, variable, selector = list()) {
  v <- nc$var[[variable]]
  if (is.null(v)) stop("Missing NetCDF variable: ", variable)
  dn <- vapply(v$dim, function(d) d$name, character(1))
  a <- ncvar_get(nc, variable)
  for (nm in names(selector)) {
    k <- match(nm, dn)
    if (is.na(k)) stop("Selector dimension not present: ", nm)
    vals <- v$dim[[k]]$vals
    # Some files use an integer top index as the dimension coordinate and
    # store the physical pressure boundary in a companion top_hPa variable.
    if (nm == "top" && "top_hPa" %in% names(nc$var)) vals <- ncvar_get(nc, "top_hPa")
    ii <- match(selector[[nm]], vals)
    if (is.na(ii)) stop("Selector value not present for ", nm, ": ", selector[[nm]])
    subs <- rep(list(quote(expr = )), length(dim(a)))
    subs[[k]] <- ii
    a <- do.call(`[`, c(list(a), subs, list(drop = FALSE)))
    a <- drop(a)
    dn <- dn[-k]
  }
  if (!identical(dn, c("longitude", "latitude"))) {
    perm <- match(c("longitude", "latitude"), dn)
    if (anyNA(perm)) stop("Cannot orient ", variable, "; dimensions are ", paste(dn, collapse = ", "))
    a <- aperm(a, perm)
  }
  a
}

canonicalize_lonlat <- function(a, lon, lat) {
  lon2 <- ifelse(lon >= 180, lon - 360, lon)
  oi <- order(lon2)
  oj <- order(lat)
  list(field = a[oi, oj, drop = FALSE], lon = lon2[oi], lat = lat[oj])
}

# Start from the locked plotting geometry only; every scientific value below is
# replaced by the formal TOP100 result.
data12 <- load_figure_data_v3_FIXED_FULL(ROOT)

trend_nc_path <- file.path(ROOT, "output_TOP100_revision", "derived_data", "trend_fields",
                           "TOP_BOUNDARY_GRIDCELL_FIELDS.nc")
n <- nc_open(trend_nc_path)
lon_nc <- n$dim$longitude$vals
lat_nc <- n$dim$latitude$vals
z100 <- canonicalize_lonlat(read_nc_slice(n, "zbar_trend", list(top = 100)), lon_nc, lat_nc)
w100 <- canonicalize_lonlat(read_nc_slice(n, "wgpe_trend", list(top = 100)), lon_nc, lat_nc)
i100 <- canonicalize_lonlat(read_nc_slice(n, "iwv_trend", list(top = 100)), lon_nc, lat_nc)
nc_close(n)

stopifnot(identical(as.numeric(z100$lon), as.numeric(data12$lon)))
stopifnot(identical(as.numeric(z100$lat), as.numeric(data12$lat)))
data12$ztrend <- z100$field
data12$wgpe_trend <- w100$field
data12$iwv_trend <- i100$field

trend <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "core",
                            "Trend_results_FINAL_v3.csv"), check.names = FALSE)
trend <- trend[trend$height_reference == "absolute", , drop = FALSE]
global <- trend[trend$region == "Global", , drop = FALSE]
data12$stats$z_global <- global$z_trend_m_yr[1]
data12$stats$wgpe_global <- global$wgpe_trend_J_m2_yr[1]
data12$stats$iwv_global <- global$iwv_trend_kg_m2_yr[1]

decomp <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "core",
                             "WGPE_decomposition_FINAL_v3.csv"), check.names = FALSE)
decomp <- decomp[decomp$height_reference == "absolute", , drop = FALSE]
data12$wgpe <- transform(
  decomp,
  iwv_channel_J_m2_yr = IWV,
  z_channel_J_m2_yr = z,
  interaction_channel_J_m2_yr = interaction,
  iwv_channel_pct = 100 * IWV / actual,
  z_channel_pct = 100 * z / actual,
  interaction_channel_pct = 100 * interaction / actual
)

layers <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "top_boundary",
                             "TOP_BOUNDARY_LAYER_CONTRIBUTIONS.csv"), check.names = FALSE)
layers_global <- layers[layers$region == "Global", , drop = FALSE]
layer_summary <- data.frame(
  actual_zbar_trend_m_yr = layers_global$actual_TOP100_zbar_trend_m_yr[1],
  sum_mass_contribution_m_yr = sum(layers_global$mass_contribution_m_yr, na.rm = TRUE),
  sum_height_contribution_m_yr = sum(layers_global$height_contribution_m_yr, na.rm = TRUE),
  residual_m_yr = layers_global$actual_TOP100_zbar_trend_m_yr[1] -
    sum(layers_global$total_first_order_contribution_m_yr, na.rm = TRUE)
)
layer_summary$residual_pct_of_actual <- 100 * layer_summary$residual_m_yr /
  layer_summary$actual_zbar_trend_m_yr
data12$layer <- layer_summary
data12$layers6 <- layers_global
upper <- layers_global[layers_global$layer %in% c("300-200", "200-100"), , drop = FALSE]
upper_agg <- upper[1, , drop = FALSE]
upper_agg$layer <- "300-100"
upper_agg$p_bottom_hPa <- 300
upper_agg$p_top_hPa <- 100
for (nm in c("mass_contribution_m_yr", "height_contribution_m_yr",
             "total_first_order_contribution_m_yr")) {
  upper_agg[[nm]] <- sum(upper[[nm]], na.rm = TRUE)
}
data12$layers5 <- rbind(
  layers_global[!layers_global$layer %in% c("300-200", "200-100"), , drop = FALSE],
  upper_agg
)
data12$correlations <- read.csv(file.path(ROOT, "output_TOP100_revision", "results", "core",
                                        "Trend_pattern_correlations_FINAL_v3.csv"), check.names = FALSE)
data12$correlations <- data12$correlations[data12$correlations$height_reference == "absolute",
                                           c("variable", "pearson_r"), drop = FALSE]
names(data12$correlations)[2] <- "r"

U_CORE <- read.csv(file.path(ROOT, "output_TOP100_revision", "results",
                             "baseline_FULL_RERUN_TOP100", "05_uncertainty",
                             "Figure1_core_uncertainty_FINAL_v3.csv"), check.names = FALSE)
U_LAT <- read.csv(file.path(ROOT, "output_TOP100_revision", "results",
                            "baseline_FULL_RERUN_TOP100", "05_uncertainty",
                            "latitude_profile_uncertainty_FINAL_v3.csv"), check.names = FALSE)

# Spatial 10-degree block bootstrap for the six-layer vertical diagnostic.
layer_nc_path <- file.path(ROOT, "output_TOP100_revision", "derived_data", "decomposition_fields",
                           "TOP100_LAYER_GRIDCELL_FIELDS.nc")
ln <- nc_open(layer_nc_path)
layer_names <- as.character(ncvar_get(ln, "layer_name"))
if (length(layer_names) != 6L || any(!nzchar(layer_names))) {
  layer_names <- c("1000-850", "850-700", "700-500", "500-300", "300-200", "200-100")
}
read_layer_stack <- function(variable) {
  v <- ln$var[[variable]]
  dn <- vapply(v$dim, function(d) d$name, character(1))
  a <- ncvar_get(ln, variable)
  perm <- match(c("longitude", "latitude", "layer"), dn)
  if (anyNA(perm)) stop("Unexpected layer field dimensions for ", variable)
  a <- aperm(a, perm)
  lon2 <- ifelse(ln$dim$longitude$vals >= 180, ln$dim$longitude$vals - 360,
                 ln$dim$longitude$vals)
  a[order(lon2), order(ln$dim$latitude$vals), , drop = FALSE]
}
mass_stack <- read_layer_stack("mass_redistribution_contribution")
height_stack <- read_layer_stack("height_coordinate_contribution")
total_stack <- read_layer_stack("total_first_order_contribution")
nc_close(ln)

lon <- z100$lon
lat <- z100$lat
w <- outer(rep(1, length(lon)), cos(lat * pi / 180))
block_id <- outer(floor((lon + 180) / 10), floor((lat + 90) / 10),
                  function(x, y) paste(x, y, sep = "_"))
block_levels <- unique(as.vector(block_id))

block_stats <- function(field) {
  valid <- is.finite(field)
  num <- tapply(as.vector(field * w), as.vector(block_id), sum, na.rm = TRUE)
  den <- tapply(as.vector(w * valid), as.vector(block_id), sum, na.rm = TRUE)
  data.frame(block = names(num), num = as.numeric(num), den = as.numeric(den))
}

bootstrap_mean <- function(field, B = 1000L, seed = 20260714L) {
  bs <- block_stats(field)
  bs <- bs[is.finite(bs$den) & bs$den > 0, , drop = FALSE]
  estimate <- sum(bs$num) / sum(bs$den)
  set.seed(seed)
  sims <- replicate(B, {
    ii <- sample.int(nrow(bs), nrow(bs), replace = TRUE)
    sum(bs$num[ii]) / sum(bs$den[ii])
  })
  c(estimate = estimate, lower95 = unname(quantile(sims, .025, na.rm = TRUE)),
    upper95 = unname(quantile(sims, .975, na.rm = TRUE)))
}

metric_fields <- list(
  actual = z100$field,
  sum_mass = apply(mass_stack, c(1, 2), sum, na.rm = TRUE),
  sum_height = apply(height_stack, c(1, 2), sum, na.rm = TRUE),
  residual = z100$field - apply(total_stack, c(1, 2), sum, na.rm = TRUE)
)
for (k in seq_along(layer_names)) {
  metric_fields[[paste0("total_contrib_", gsub("-", "_", layer_names[k]))]] <- total_stack[, , k]
}
metric_fields[["total_contrib_300_100"]] <- total_stack[, , 5] + total_stack[, , 6]
U_VB <- do.call(rbind, lapply(seq_along(metric_fields), function(k) {
  ci <- bootstrap_mean(metric_fields[[k]], seed = 20260714L + k)
  data.frame(variable = names(metric_fields)[k], estimate = ci["estimate"],
             lower95 = ci["lower95"], upper95 = ci["upper95"],
             metric = names(metric_fields)[k], height_reference = "absolute",
             region = "Global",
             uncertainty_method = "spatial 10-degree block bootstrap; 1000 resamples")
}))
row.names(U_VB) <- NULL

write.csv(U_VB, file.path(OUT, "Figure1_vertical_bridge_uncertainty_TOP100.csv"), row.names = FALSE)
write.csv(layers_global, file.path(OUT, "Figure1_six_layer_contributions_TOP100.csv"), row.names = FALSE)
write.csv(decomp, file.path(OUT, "Figure1_regional_WGPE_decomposition_TOP100.csv"), row.names = FALSE)
saveRDS(list(data12 = data12, U_CORE = U_CORE, U_LAT = U_LAT, U_VB = U_VB,
             trend_nc = trend_nc_path, layer_nc = layer_nc_path,
             core_rdata = file.path(ROOT, "output_TOP100_revision", "derived_data", "ERA5_TOP100",
                                    "global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData")),
        file.path(OUT, "Figure1_TOP100_adapter.rds"))

cat("Saved Figure1_TOP100_adapter.rds\n")
cat("z trend:", data12$stats$z_global, "\n")
cat("W-GPE trend:", data12$stats$wgpe_global, "\n")
cat("six-layer sum:", sum(layers_global$total_first_order_contribution_m_yr), "\n")
cat("vertical residual:", layer_summary$residual_m_yr, "\n")
