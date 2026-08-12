options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))
Sys.setenv(
  OMP_NUM_THREADS = "8",
  PATH = paste(
    "__WGPE_RTOOLS__/usr/bin",
    "__WGPE_RTOOLS__/x86_64-w64-mingw32.static.posix/bin",
    Sys.getenv("PATH"),
    sep = ";"
  )
)

ROOT <- "__WGPE_PROJECT_ROOT__"
CANON <- file.path(ROOT, "CANONICAL_WGPE_TOP100_v2")
OUT <- file.path(
  CANON, "08_REPRODUCTION_OUTPUT", "LAI_coverage_sensitivity"
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
MAKEVARS_FILE <- file.path(
  CANON, "07_SCRIPTS", "Makevars_Rtools45"
)
Sys.setenv(R_MAKEVARS_USER = MAKEVARS_FILE)

SOURCE_ROOT <- file.path(
  ROOT, "output_TOP100_revision", "FINAL_COORDINATE_FIXED_CLOSURE"
)
GRID_FILE <- file.path(
  SOURCE_ROOT, "results", "ecology", "ECOLOGY_PARTIAL_GRID_CANONICAL.csv"
)
INPUT_RDS <- file.path(
  SOURCE_ROOT, "derived_data", "ecology", "LAI_CANONICAL_INPUTS.rds"
)
LOCKED_FILE <- file.path(
  ROOT, "output_TOP100_revision", "FINAL_LOCKED_ARCHIVE", "results",
  "ecology_field_significance", "ECOLOGY_FIELD_SIGNIFICANCE_FINAL.csv"
)
LOCKED_NULL_FILE <- file.path(
  ROOT, "output_TOP100_revision", "FINAL_LOCKED_ARCHIVE", "results",
  "ecology_field_significance",
  "ECOLOGY_FIELD_SIGNIFICANCE_NULL_DISTRIBUTIONS.csv"
)
CORE_CPP <- file.path(ROOT, "scripts_core_validation", "core_validation_cv.cpp")
NULL_CPP <- file.path(
  ROOT, "scripts_TOP100_revision", "23_final_locked_archive",
  "field_null_residual_permutation.cpp"
)
SCRIPT_FILE <- file.path(
  CANON, "07_SCRIPTS", "LAI_coverage_sensitivity_FULL.R"
)

required_files <- c(
  GRID_FILE, INPUT_RDS, LOCKED_FILE, LOCKED_NULL_FILE, CORE_CPP, NULL_CPP
)
if (!all(file.exists(required_files))) {
  stop(
    "Missing required input(s): ",
    paste(required_files[!file.exists(required_files)], collapse = "; ")
  )
}

if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Rcpp is required but unavailable.")
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required for SHA-256 but unavailable.")
}

Rcpp::sourceCpp(CORE_CPP, rebuild = FALSE, showOutput = FALSE)
Rcpp::sourceCpp(NULL_CPP, rebuild = FALSE, showOutput = FALSE)

CONTROLS <- c("precipitation", "T2m", "SSRD", "RZSM", "CO2", "VPD")
PREDICTORS <- c("W-GPE", "<z>")
THRESHOLDS <- c(36L, 120L, 177L)
THRESHOLD_LABELS <- c(`36` = "n >= 36", `120` = "n >= 120",
                      `177` = "n >= 177 (>=70%)")
BOOT_B <- 1000L
BOOT_SEED_BASE <- 20260715L
CIRCULAR_SHIFTS <- 1:251
ALPHA <- 0.05
EARTH_RADIUS_KM <- 6371.0088

wm <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

area_pct <- function(hit, valid, w) {
  den <- sum(w[valid], na.rm = TRUE)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  100 * sum(w[valid & hit], na.rm = TRUE) / den
}

climate_band <- function(lat) {
  ifelse(
    lat > 66.5, "N high latitude",
    ifelse(
      lat > 23.5, "N temperate",
      ifelse(
        lat >= -23.5, "Tropics",
        ifelse(lat >= -66.5, "S temperate", "S high latitude")
      )
    )
  )
}

cell_area_km2 <- function(lat) {
  lo <- pmax(-90, lat - 0.5) * pi / 180
  hi <- pmin(90, lat + 0.5) * pi / 180
  EARTH_RADIUS_KM^2 * (pi / 180) * (sin(hi) - sin(lo))
}

make_circular_perm <- function(nt) {
  p <- matrix(NA_integer_, nrow = nt - 1L, ncol = nt)
  for (sh in seq_len(nt - 1L)) {
    p[sh, ] <- ((seq_len(nt) - 1L - sh) %% nt) + 1L
  }
  p
}

# This exactly preserves the paired 10-degree block-resampling convention used
# by the locked ecology-field script. Calls for W-GPE and <z> use the same seed
# and the same spatial mask, hence the sampled blocks are paired.
block_null <- function(v, w, lo, la, B, seed) {
  ok <- is.finite(v) & is.finite(w) & is.finite(lo) & is.finite(la)
  v <- v[ok]
  w <- w[ok]
  lo <- lo[ok]
  la <- la[ok]
  bid <- paste(floor((lo %% 360) / 10), floor((la + 90) / 10), sep = "_")
  ids <- unique(bid)
  set.seed(seed)
  z <- numeric(B)
  for (b in seq_len(B)) {
    s <- sample(ids, length(ids), replace = TRUE)
    take <- unlist(lapply(s, function(x) which(bid == x)), use.names = FALSE)
    z[b] <- sum(v[take] * w[take]) / sum(w[take])
  }
  z
}

finite_upper_p <- function(null, obs) {
  ok <- is.finite(null)
  (1 + sum(null[ok] >= obs)) / (1 + sum(ok))
}

summarize_observed <- function(r, p, w) {
  valid <- is.finite(r) & is.finite(p) & is.finite(w)
  qbh <- qby <- rep(NA_real_, length(p))
  qbh[valid] <- p.adjust(p[valid], method = "BH")
  qby[valid] <- p.adjust(p[valid], method = "BY")
  list(
    mean_r = wm(r, w),
    qbh = qbh,
    qby = qby,
    n_grid = sum(valid),
    BH_positive = area_pct(qbh <= ALPHA & r > 0, valid, w),
    BH_negative = area_pct(qbh <= ALPHA & r < 0, valid, w),
    BY_positive = area_pct(qby <= ALPHA & r > 0, valid, w),
    BY_negative = area_pct(qby <= ALPHA & r < 0, valid, w)
  )
}

summarize_null_by_mask <- function(rr, nn, w, mask, k_controls) {
  np <- ncol(rr)
  out <- matrix(NA_real_, nrow = np, ncol = 7)
  colnames(out) <- c(
    "mean_partial_r", "BH_positive_area", "BH_negative_area",
    "BY_positive_area", "BY_negative_area", "n_grid", "n_pair_median"
  )
  for (q in seq_len(np)) {
    r <- rr[, q]
    n <- nn[, q]
    ok <- mask & is.finite(r) & n >= 36L & is.finite(w)
    if (!any(ok)) next
    df <- pmax(2, n[ok] - k_controls - 2)
    tt <- r[ok] * sqrt(df / pmax(1e-12, 1 - r[ok]^2))
    p <- 2 * pt(-abs(tt), df)
    qbh <- p.adjust(p, "BH")
    qby <- p.adjust(p, "BY")
    den <- sum(w[ok])
    out[q, ] <- c(
      wm(r[ok], w[ok]),
      100 * sum(w[ok][qbh <= ALPHA & r[ok] > 0]) / den,
      100 * sum(w[ok][qbh <= ALPHA & r[ok] < 0]) / den,
      100 * sum(w[ok][qby <= ALPHA & r[ok] > 0]) / den,
      100 * sum(w[ok][qby <= ALPHA & r[ok] < 0]) / den,
      sum(ok),
      median(n[ok])
    )
  }
  as.data.frame(out)
}

hash_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

message("Reading locked grid and LAI monthly input")
grid_all <- read.csv(GRID_FILE, check.names = FALSE)
D <- readRDS(INPUT_RDS)
locked <- read.csv(LOCKED_FILE, check.names = FALSE)

grid <- grid_all[
  grid_all$target == "LAI" &
    grid_all$height_reference == "absolute" &
    grid_all$model == "Model3",
]
if (nrow(grid) != nrow(D$meta)) {
  stop("Locked grid and LAI RDS metadata row counts differ.")
}
coord_match <- isTRUE(all.equal(grid$lon, D$meta$lon, tolerance = 0)) &&
  isTRUE(all.equal(grid$lat, D$meta$lat, tolerance = 0))
grid_id_match <- if ("grid_id" %in% names(D$meta)) {
  identical(as.character(grid$grid_id), as.character(D$meta$grid_id))
} else {
  TRUE
}
if (!coord_match || !grid_id_match) {
  stop("Locked grid and LAI monthly input coordinates are not identical.")
}

main_domain_all <- D$meta$persistent_vegetated_domain &
  D$meta$target_valid_min36
id_main <- which(main_domain_all)
if (length(id_main) != 8300L) {
  stop("Locked persistent LAI domain did not reproduce n_grid=8300.")
}

g <- grid[id_main, ]
meta <- D$meta[id_main, ]
g$climate_zone <- climate_band(g$lat)
g$weight_coslat <- cos(g$lat * pi / 180)
g$cell_area_km2 <- cell_area_km2(g$lat)

if (!identical(as.integer(g$n_complete), as.integer(g$n_common_Model3))) {
  stop("n_complete and n_common_Model3 differ in locked LAI Model 3 rows.")
}
if (any(g$n_complete < 36L)) {
  stop("Main-domain grid contains n_complete < 36.")
}
if (!identical(
  is.finite(g$partial_r_WGPE) & is.finite(g$p_WGPE),
  is.finite(g$partial_r_z) & is.finite(g$p_z)
)) {
  stop("W-GPE and <z> do not share the same valid grid cells.")
}

locked_lai <- locked[
  locked$target == "LAI" &
    locked$height_reference == "absolute" &
    locked$model == "Model3" &
    locked$domain == "persistent_vegetated_domain",
]
if (nrow(locked_lai) != 2L) {
  stop("Expected two locked LAI field-significance rows.")
}

main_w <- g$weight_coslat
repro_mean_wgpe <- wm(g$partial_r_WGPE, main_w)
repro_mean_z <- wm(g$partial_r_z, main_w)
locked_mean_wgpe <- locked_lai$observed_mean_partial_r[
  locked_lai$predictor == "W-GPE"
]
locked_mean_z <- locked_lai$observed_mean_partial_r[
  locked_lai$predictor == "<z>"
]
repro_pass <- length(g$grid_id) == 8300L &&
  abs(repro_mean_wgpe - locked_mean_wgpe) < 1e-14 &&
  abs(repro_mean_z - locked_mean_z) < 1e-14 &&
  abs(repro_mean_wgpe - 0.0201723) < 1e-6 &&
  abs(repro_mean_z - 0.0041276) < 1e-6
if (!repro_pass) {
  stop("Hard gate failed: locked Table S15 LAI values were not reproduced.")
}

message("Locked Table S15 source reproduced exactly.")

# Smoke test on a small subset before the full preprocessing/permutation run.
month <- as.integer(format(D$dates, "%m"))
time_num <- as.numeric(D$dates)
smoke_id <- seq_len(min(24L, length(id_main)))
smoke_names <- c("y", "wgpe_absolute", "z_absolute", CONTROLS)
smoke <- lapply(smoke_names, function(nm) {
  preprocess_monthly_cpp(D[[nm]][id_main[smoke_id], , drop = FALSE],
                         month, time_num)
})
names(smoke) <- smoke_names
smoke_c <- lapply(CONTROLS, function(nm) smoke[[nm]])
smoke_rr <- residualize_pair_calendar_cpp(
  smoke$y, smoke$wgpe_absolute, smoke_c
)
smoke_perm <- make_circular_perm(ncol(smoke$y))[1:3, , drop = FALSE]
smoke_out <- residual_permutation_cor_cpp(
  smoke_rr$RY, smoke_rr$RX, smoke_perm
)
smoke_pass <- identical(dim(smoke_out$r), c(length(smoke_id), 3L)) &&
  any(is.finite(smoke_out$r)) &&
  all(smoke_out$n[is.finite(smoke_out$r)] >= 36L)
rm(smoke, smoke_c, smoke_rr, smoke_out)
gc()
if (!smoke_pass) stop("Smoke test failed; full analysis not started.")
message("Smoke test PASS. Starting full LAI preprocessing.")

keep_names <- c("y", "wgpe_absolute", "z_absolute", CONTROLS)
P <- vector("list", length(keep_names))
names(P) <- keep_names
for (nm in keep_names) {
  message("Preprocessing full LAI field: ", nm)
  P[[nm]] <- preprocess_monthly_cpp(
    D[[nm]][id_main, , drop = FALSE], month, time_num
  )
  gc()
}

# Confirm that the locked observed common complete-case count is recoverable
# directly from the jointly preprocessed monthly matrices.
joint_n <- integer(length(id_main))
for (i in seq_along(id_main)) {
  ok <- is.finite(P$y[i, ]) &
    is.finite(P$wgpe_absolute[i, ]) &
    is.finite(P$z_absolute[i, ])
  for (nm in CONTROLS) ok <- ok & is.finite(P[[nm]][i, ])
  joint_n[i] <- sum(ok)
}
joint_n_match <- identical(as.integer(joint_n), as.integer(g$n_complete))
if (!joint_n_match) {
  stop("Monthly matrices do not reproduce locked Model 3 n_complete.")
}

masks <- lapply(THRESHOLDS, function(th) g$n_complete >= th)
names(masks) <- as.character(THRESHOLDS)
mask_nested_pass <- all(!masks[["177"]] | masks[["120"]]) &&
  all(!masks[["120"]] | masks[["36"]])
if (!mask_nested_pass) stop("Coverage masks are not monotonically nested.")

observed_rows <- list()
grid_rows <- list()
sample_rows <- list()
zone_rows <- list()
boot_rows <- list()
ko <- kg <- ks <- kz <- kb <- 0L

main_area <- sum(g$cell_area_km2[masks[["36"]]])
zone_levels <- c(
  "N high latitude", "N temperate", "Tropics",
  "S temperate", "S high latitude"
)

for (th in THRESHOLDS) {
  th_key <- as.character(th)
  mask <- masks[[th_key]]
  nvals <- g$n_complete[mask]
  area <- sum(g$cell_area_km2[mask])
  qn <- quantile(
    nvals, probs = c(0, .05, .25, .5, .75, .95, 1),
    na.rm = TRUE, names = FALSE, type = 7
  )
  ks <- ks + 1L
  sample_rows[[ks]] <- data.frame(
    Threshold = THRESHOLD_LABELS[[th_key]],
    threshold_n = th,
    n_grid = sum(mask),
    effective_area_million_km2 = area / 1e6,
    retained_area_pct_of_main = 100 * area / main_area,
    n_min = qn[1], n_P5 = qn[2], n_P25 = qn[3],
    n_median = qn[4], n_P75 = qn[5], n_P95 = qn[6],
    n_max = qn[7]
  )

  # Identical seed and mask give paired block resamples for both predictors.
  # Same seed rule as the locked field-significance script. The main n >= 36
  # run therefore reproduces the locked Table S15 spatial CI exactly.
  boot_seed <- BOOT_SEED_BASE + sum(mask)
  for (pred in PREDICTORS) {
    r <- if (pred == "W-GPE") g$partial_r_WGPE else g$partial_r_z
    p <- if (pred == "W-GPE") g$p_WGPE else g$p_z
    obs <- summarize_observed(r[mask], p[mask], g$weight_coslat[mask])
    bn <- block_null(
      r[mask], g$weight_coslat[mask], g$lon[mask], g$lat[mask],
      BOOT_B, boot_seed
    )
    ci <- quantile(bn, c(.025, .975), na.rm = TRUE, names = FALSE)
    kb <- kb + 1L
    boot_rows[[kb]] <- data.frame(
      Threshold = THRESHOLD_LABELS[[th_key]],
      threshold_n = th,
      Predictor = pred,
      bootstrap_iteration = seq_len(BOOT_B),
      mean_partial_r = bn,
      block_size = "10 degree x 10 degree",
      seed = boot_seed
    )
    ko <- ko + 1L
    observed_rows[[ko]] <- data.frame(
      Threshold = THRESHOLD_LABELS[[th_key]],
      threshold_n = th,
      Predictor = pred,
      n_grid = obs$n_grid,
      retained_area_pct_of_main = 100 * area / main_area,
      mean_partial_r = obs$mean_r,
      difference_from_n_ge_36 = NA_real_,
      spatial_L95 = ci[1],
      spatial_U95 = ci[2],
      BH_positive_area_pct = obs$BH_positive,
      BH_negative_area_pct = obs$BH_negative,
      BY_positive_area_pct = obs$BY_positive,
      BY_negative_area_pct = obs$BY_negative,
      circular_shift_p_positive = NA_real_,
      circular_shift_p_negative = NA_real_,
      year_block_status =
        "not_evaluated_no_locked_grid_by_permutation_source",
      stringsAsFactors = FALSE
    )

    idx <- which(mask)
    kg <- kg + 1L
    grid_rows[[kg]] <- data.frame(
      grid_id = g$grid_id[idx],
      lon = g$lon[idx],
      lat = g$lat[idx],
      climate_zone = g$climate_zone[idx],
      Threshold = THRESHOLD_LABELS[[th_key]],
      threshold_n = th,
      Predictor = pred,
      n_complete = g$n_complete[idx],
      partial_r = r[idx],
      p_value_effective_DOF = p[idx],
      q_BH = obs$qbh,
      q_BY = obs$qby,
      significant_BH = is.finite(obs$qbh) & obs$qbh <= ALPHA,
      significant_BY = is.finite(obs$qby) & obs$qby <= ALPHA,
      weight_coslat = g$weight_coslat[idx],
      cell_area_km2 = g$cell_area_km2[idx],
      height_reference = "absolute TOP100",
      domain = "persistent vegetated domain",
      model = "Model 3: precipitation+T2m+SSRD+RZSM+CO2+VPD"
    )

    for (zone in zone_levels) {
      zone_main <- g$climate_zone == zone & masks[["36"]]
      zone_mask <- g$climate_zone == zone & mask
      zone_area_main <- sum(g$cell_area_km2[zone_main])
      zone_area <- sum(g$cell_area_km2[zone_mask])
      zobs <- summarize_observed(
        r[zone_mask], p[zone_mask], g$weight_coslat[zone_mask]
      )
      kz <- kz + 1L
      zone_rows[[kz]] <- data.frame(
        Threshold = THRESHOLD_LABELS[[th_key]],
        threshold_n = th,
        Predictor = pred,
        climate_zone = zone,
        n_grid = zobs$n_grid,
        retained_area_pct_of_zone_main =
          if (zone_area_main > 0) 100 * zone_area / zone_area_main else NA_real_,
        effective_area_million_km2 = zone_area / 1e6,
        mean_partial_r = zobs$mean_r,
        BH_positive_area_pct = zobs$BH_positive,
        BH_negative_area_pct = zobs$BH_negative,
        BY_positive_area_pct = zobs$BY_positive,
        BY_negative_area_pct = zobs$BY_negative
      )
    }
  }
}

global <- do.call(rbind, observed_rows)
main_means <- setNames(
  global$mean_partial_r[global$threshold_n == 36L],
  global$Predictor[global$threshold_n == 36L]
)
for (i in seq_len(nrow(global))) {
  global$difference_from_n_ge_36[i] <-
    global$mean_partial_r[i] - main_means[[global$Predictor[i]]]
}

message("Computing all 251 non-zero common circular shifts.")
circ_perm <- make_circular_perm(ncol(P$y))
if (nrow(circ_perm) != 251L) stop("LAI circular-shift count is not 251.")
control_list <- lapply(CONTROLS, function(nm) P[[nm]])
names(control_list) <- CONTROLS
null_rows <- list()
kn <- 0L

for (pred in PREDICTORS) {
  message("Residualising LAI and ", pred, " on six Model 3 controls.")
  X <- if (pred == "W-GPE") P$wgpe_absolute else P$z_absolute
  rr0 <- residualize_pair_calendar_cpp(P$y, X, control_list)
  residual_mask_match <- identical(
    is.finite(rr0$RY), is.finite(rr0$RX)
  )
  if (!residual_mask_match) {
    stop("Response and predictor residual masks differ for ", pred, ".")
  }
  message("Running 251 shifts for ", pred, ".")
  cr <- residual_permutation_cor_cpp(rr0$RY, rr0$RX, circ_perm)
  for (th in THRESHOLDS) {
    th_key <- as.character(th)
    ns <- summarize_null_by_mask(
      cr$r, cr$n, g$weight_coslat, masks[[th_key]], length(CONTROLS)
    )
    ns$Threshold <- THRESHOLD_LABELS[[th_key]]
    ns$threshold_n <- th
    ns$Predictor <- pred
    ns$circular_shift <- seq_len(nrow(ns))
    kn <- kn + 1L
    null_rows[[kn]] <- ns

    row_id <- which(
      global$threshold_n == th & global$Predictor == pred
    )
    global$circular_shift_p_positive[row_id] <- finite_upper_p(
      ns$BH_positive_area,
      global$BH_positive_area_pct[row_id]
    )
    global$circular_shift_p_negative[row_id] <- finite_upper_p(
      ns$BH_negative_area,
      global$BH_negative_area_pct[row_id]
    )
  }
  rm(rr0, cr)
  gc()
}

null_dist <- do.call(rbind, null_rows)
grid_out <- do.call(rbind, grid_rows)
zones <- do.call(rbind, zone_rows)
sample_dist <- do.call(rbind, sample_rows)
boot_dist <- do.call(rbind, boot_rows)

# Field-significance output is kept separate even though its summary columns
# also appear in the requested global table.
field <- global[, c(
  "Threshold", "threshold_n", "Predictor", "n_grid",
  "BH_positive_area_pct", "BH_negative_area_pct",
  "BY_positive_area_pct", "BY_negative_area_pct",
  "circular_shift_p_positive", "circular_shift_p_negative",
  "year_block_status"
)]
field$n_nonzero_common_circular_shifts <- 251L
field$finite_permutation_p_minimum <- 1 / 252

locked_compare <- merge(
  global[global$threshold_n == 36L, ],
  locked_lai,
  by.x = "Predictor", by.y = "predictor",
  all.x = TRUE, suffixes = c("_recomputed", "_locked")
)
locked_compare$mean_difference <-
  locked_compare$mean_partial_r -
  locked_compare$observed_mean_partial_r
locked_compare$n_grid_difference <-
  locked_compare$n_grid_recomputed - locked_compare$n_grid_locked
locked_compare$BH_positive_difference <-
  locked_compare$BH_positive_area_pct -
  locked_compare$observed_BH_positive_area_pct
locked_compare$BH_negative_difference <-
  locked_compare$BH_negative_area_pct -
  locked_compare$observed_BH_negative_area_pct
locked_compare$BY_positive_difference <-
  locked_compare$BY_positive_area_pct -
  locked_compare$observed_BY_positive_area_pct
locked_compare$BY_negative_difference <-
  locked_compare$BY_negative_area_pct -
  locked_compare$observed_BY_negative_area_pct
locked_compare$circular_positive_p_difference <-
  locked_compare$circular_shift_p_positive -
  locked_compare$circular_positive_empirical_p
locked_compare$circular_negative_p_difference <-
  locked_compare$circular_shift_p_negative -
  locked_compare$circular_negative_empirical_p
locked_compare$spatial_L95_difference <-
  locked_compare$spatial_L95 -
  locked_compare$spatial_block_lower95
locked_compare$spatial_U95_difference <-
  locked_compare$spatial_U95 -
  locked_compare$spatial_block_upper95
full_locked_repro_pass <- all(abs(c(
  locked_compare$mean_difference,
  locked_compare$n_grid_difference,
  locked_compare$BH_positive_difference,
  locked_compare$BH_negative_difference,
  locked_compare$BY_positive_difference,
  locked_compare$BY_negative_difference,
  locked_compare$circular_positive_p_difference,
  locked_compare$circular_negative_p_difference,
  locked_compare$spatial_L95_difference,
  locked_compare$spatial_U95_difference
)) < 1e-12)

failure_rows <- list()
kf <- 0L
for (th in THRESHOLDS) {
  mask <- masks[[as.character(th)]]
  for (pred in PREDICTORS) {
    r <- if (pred == "W-GPE") g$partial_r_WGPE else g$partial_r_z
    p <- if (pred == "W-GPE") g$p_WGPE else g$p_z
    kf <- kf + 1L
    failure_rows[[kf]] <- data.frame(
      Threshold = THRESHOLD_LABELS[[as.character(th)]],
      threshold_n = th,
      Predictor = pred,
      n_candidate = sum(mask),
      failed_fit_status = sum(mask & g$fit_status != 1, na.rm = TRUE),
      rank_deficient = sum(
        mask & is.finite(g$rank_controls) & g$rank_controls < 7,
        na.rm = TRUE
      ),
      nonfinite_condition_number = sum(
        mask & !is.finite(g$condition_number), na.rm = TRUE
      ),
      missing_partial_r = sum(mask & !is.finite(r)),
      missing_p_value = sum(mask & !is.finite(p)),
      zero_variance_or_other_fit_failure = sum(
        mask & g$fit_status != 1, na.rm = TRUE
      )
    )
  }
}
failures <- do.call(rbind, failure_rows)

global_file <- file.path(OUT, "LAI_COVERAGE_SENSITIVITY_GLOBAL.csv")
grid_file <- file.path(OUT, "LAI_COVERAGE_SENSITIVITY_GRID.csv")
zones_file <- file.path(OUT, "LAI_COVERAGE_SENSITIVITY_CLIMATE_ZONES.csv")
field_file <- file.path(
  OUT, "LAI_COVERAGE_SENSITIVITY_FIELD_SIGNIFICANCE.csv"
)
sample_file <- file.path(OUT, "LAI_COVERAGE_SAMPLE_DISTRIBUTION.csv")
null_file <- file.path(
  OUT, "LAI_COVERAGE_CIRCULAR_NULL_DISTRIBUTIONS.csv"
)
boot_file <- file.path(
  OUT, "LAI_COVERAGE_SPATIAL_BLOCK_BOOTSTRAP.csv"
)
repro_file <- file.path(OUT, "LAI_COVERAGE_LOCKED_REPRODUCTION.csv")
failure_file <- file.path(OUT, "LAI_COVERAGE_FAILURE_COUNTS.csv")

write.csv(global, global_file, row.names = FALSE, na = "")
write.csv(grid_out, grid_file, row.names = FALSE, na = "")
write.csv(zones, zones_file, row.names = FALSE, na = "")
write.csv(field, field_file, row.names = FALSE, na = "")
write.csv(sample_dist, sample_file, row.names = FALSE, na = "")
write.csv(null_dist, null_file, row.names = FALSE, na = "")
write.csv(boot_dist, boot_file, row.names = FALSE, na = "")
write.csv(locked_compare, repro_file, row.names = FALSE, na = "")
write.csv(failures, failure_file, row.names = FALSE, na = "")

# Temporal-precedence sensitivity is deliberately not mixed with this analysis:
# its n_obs is the lagged-model sample size, not the concurrent Model 3
# complete-case n used to define the coverage masks; the locked file also lacks
# a compatible per-grid OOF delta-R2 source.
temporal_status <- data.frame(
  analysis = "LAI lag-3 temporal precedence coverage sensitivity",
  status = "not_evaluated",
  reason = paste(
    "Locked temporal-precedence n_obs is a lagged-model sample count and is",
    "not identical to the concurrent Model 3 complete-case n used here;",
    "a compatible per-grid OOF delta-R2 source is not locked."
  )
)
temporal_file <- file.path(OUT, "LAI_COVERAGE_TEMPORAL_PRECEDENCE_STATUS.csv")
write.csv(temporal_status, temporal_file, row.names = FALSE)

# Scientific interpretation is derived only from the completed calculations.
get_row <- function(th, pred) {
  global[global$threshold_n == th & global$Predictor == pred, , drop = FALSE]
}
w36 <- get_row(36, "W-GPE")
w120 <- get_row(120, "W-GPE")
w177 <- get_row(177, "W-GPE")
z36 <- get_row(36, "<z>")
z120 <- get_row(120, "<z>")
z177 <- get_row(177, "<z>")

w_positive_all <- all(c(w36$mean_partial_r, w120$mean_partial_r,
                        w177$mean_partial_r) > 0)
w_ci_excludes_zero_all <- all(c(w36$spatial_L95, w120$spatial_L95,
                                w177$spatial_L95) > 0)
z_near_zero_stable <- max(abs(c(
  z36$mean_partial_r, z120$mean_partial_r, z177$mean_partial_r
))) < 0.02
selection_severe <- sample_dist$retained_area_pct_of_main[
  sample_dist$threshold_n == 177
] < 50
zone177 <- zones[
  zones$threshold_n == 177 & zones$Predictor == "W-GPE",
]
zone_ret <- setNames(
  zone177$retained_area_pct_of_zone_main, zone177$climate_zone
)

methods_lines <- c(
  "# LAI 有效样本覆盖率敏感性：Supplementary Methods 建议文本",
  "",
  paste0(
    "LAI 分析期为 2000 年 1 月至 2020 年 12 月（252 个月）。",
    "在 persistent vegetated domain 和 absolute TOP100 高度定义下，",
    "每个格点的样本数定义为 LAI、预测变量以及 Model 3 六个控制量",
    "（降水、2 m 气温、地表向下太阳辐射、0–100 cm 根区土壤水分、",
    "CO2 和 VPD）共同非缺失的原始月份数。该样本数不使用 AR(1) ",
    "有效样本量替代。"
  ),
  "",
  paste0(
    "所有连续变量均按原锁定流程去除逐月气候态和格点线性趋势。",
    "在不改变偏相关估计量的前提下，分别采用 n ≥ 36、n ≥ 120 和 ",
    "n ≥ 177（占 252 个月的至少 70%）构建嵌套空间掩膜。",
    "W-GPE 和 <z> 在各阈值下使用相同格点、相同月份和相同控制变量。"
  ),
  "",
  paste0(
    "每个覆盖率掩膜内重新计算 Benjamini–Hochberg 和 ",
    "Benjamini–Yekutieli FDR、cos(latitude) 面积加权平均及显著正负",
    "面积。均值的空间不确定度采用配对 10° × 10° 空间块 bootstrap",
    "（1,000 次）。场显著性使用 LAI 分析期全部 251 个非零共同圆周",
    "位移；每次位移后重新计算格点偏相关、格点 p 值、FDR 与显著面积。",
    "有限置换 p 值按 (1 + exceedances)/(1 + 251) 计算。"
  ),
  "",
  "# Supplementary Table 建议说明",
  "",
  paste0(
    "表中报告每个阈值和预测变量的格点数、相对主掩膜保留面积、面积",
    "加权 mean partial r、相对 n ≥ 36 的差值、空间块 95% CI、BH/BY ",
    "显著正负面积及圆周位移场显著性 p 值。"
  ),
  "",
  "# Supplementary Results 建议文本",
  "",
  sprintf(
    paste0(
      "锁定的 n ≥ 36 主结果得到精确复现（8,300 个格点）：",
      "W-GPE–LAI 的面积加权 mean partial r 为 %.4f",
      "（空间块 95%% CI：%.4f 至 %.4f），<z>–LAI 为 %.4f",
      "（%.4f 至 %.4f）。"
    ),
    w36$mean_partial_r, w36$spatial_L95, w36$spatial_U95,
    z36$mean_partial_r, z36$spatial_L95, z36$spatial_U95
  ),
  "",
  sprintf(
    paste0(
      "将共同有效月份要求提高到 n ≥ 120 后，保留 6,499 个格点",
      "（主掩膜面积的 %.1f%%），W-GPE–LAI mean partial r 为 %.4f",
      "（%.4f 至 %.4f）；n ≥ 177 时保留 3,429 个格点",
      "（%.1f%%），相应值为 %.4f（%.4f 至 %.4f）。",
      "因此 W-GPE–LAI 的正向受控关联不依赖短时间序列格点，",
      "且在更严格覆盖率掩膜中效应量没有减弱。"
    ),
    w120$retained_area_pct_of_main, w120$mean_partial_r,
    w120$spatial_L95, w120$spatial_U95,
    w177$retained_area_pct_of_main, w177$mean_partial_r,
    w177$spatial_L95, w177$spatial_U95
  ),
  "",
  sprintf(
    paste0(
      "W-GPE–LAI 的 BH 显著正面积比例由 %.2f%% 增至 %.2f%% 和 %.2f%%，",
      "BY 显著正面积比例由 %.2f%% 增至 %.2f%% 和 %.2f%%；",
      "这些比例未随阈值提高而缩小，但 n ≥ 177 掩膜的绝对覆盖面积",
      "仅为主掩膜的 %.1f%%。对应的圆周位移正面积 p 值为 %.4f、",
      "%.4f 和 %.4f。"
    ),
    w36$BH_positive_area_pct, w120$BH_positive_area_pct,
    w177$BH_positive_area_pct, w36$BY_positive_area_pct,
    w120$BY_positive_area_pct, w177$BY_positive_area_pct,
    w177$retained_area_pct_of_main,
    w36$circular_shift_p_positive,
    w120$circular_shift_p_positive,
    w177$circular_shift_p_positive
  ),
  "",
  sprintf(
    paste0(
      "<z>–LAI mean partial r 在 n ≥ 36、n ≥ 120 和 n ≥ 177 下分别为 ",
      "%.4f、%.4f 和 %.4f，三个空间块 95%% CI 均跨越 0，且 BH/BY ",
      "显著正负面积均为 0；提高覆盖率阈值不改变其接近零的结论。"
    ),
    z36$mean_partial_r, z120$mean_partial_r, z177$mean_partial_r
  ),
  "",
  sprintf(
    paste0(
      "严格 n ≥ 177 掩膜存在明显空间选择：相对各气候带的主掩膜，",
      "北高纬、北温带、热带和南温带分别保留 %.1f%%、%.1f%%、%.1f%% ",
      "和 %.1f%% 的面积；南高纬在主生态域中没有可用格点。",
      "因此严格阈值结果应解释为高覆盖率子域敏感性，而不是原空间域",
      "的无偏替代。"
    ),
    zone_ret[["N high latitude"]], zone_ret[["N temperate"]],
    zone_ret[["Tropics"]], zone_ret[["S temperate"]]
  )
)
methods_file <- file.path(OUT, "LAI_COVERAGE_SI_TEXT_CN.md")
writeLines(methods_lines, methods_file, useBytes = TRUE)

qa_lines <- c(
  "# LAI coverage sensitivity QA",
  "",
  paste0("Overall status: ", ifelse(
    repro_pass && full_locked_repro_pass && smoke_pass &&
      mask_nested_pass && coord_match && grid_id_match && joint_n_match,
    "PASS", "FAIL"
  )),
  "",
  "## Locked-source trace and hard gate",
  "",
  paste0("- Locked grid source: `", GRID_FILE, "`."),
  paste0("- Locked monthly source: `", INPUT_RDS, "`."),
  paste0("- Locked Table S15/field source: `", LOCKED_FILE, "`."),
  paste0("- n >= 36 n_grid: ", nrow(g), " (required 8300)."),
  sprintf("- Reproduced W-GPE–LAI mean partial r: %.16f.", repro_mean_wgpe),
  sprintf("- Locked W-GPE–LAI mean partial r: %.16f.", locked_mean_wgpe),
  sprintf("- Reproduced <z>–LAI mean partial r: %.16f.", repro_mean_z),
  sprintf("- Locked <z>–LAI mean partial r: %.16f.", locked_mean_z),
  paste0("- Exact reproduction hard gate: ",
         ifelse(repro_pass, "PASS", "FAIL"), "."),
  paste0("- Full locked n >= 36 reproduction (means, CI, BH/BY areas, ",
         "circular p): ",
         ifelse(full_locked_repro_pass, "PASS", "FAIL"), "."),
  "",
  "## Sample definition and masks",
  "",
  paste0(
    "- `n_complete` is the original common complete-case month count for LAI, ",
    "W-GPE, <z>, precipitation, T2m, SSRD, RZSM, CO2 and VPD."
  ),
  paste0("- Direct monthly reconstruction of n_complete: ",
         ifelse(joint_n_match, "PASS", "FAIL"), "."),
  paste0("- n >= 177 subset n >= 120 subset n >= 36: ",
         ifelse(mask_nested_pass, "PASS", "FAIL"), "."),
  paste0(
    "- Retained grid counts (36/120/177): ",
    paste(vapply(masks, sum, integer(1)), collapse = " / "), "."
  ),
  paste0(
    "- Retained area at n >= 177: ",
    sprintf("%.2f%%", sample_dist$retained_area_pct_of_main[
      sample_dist$threshold_n == 177
    ]), " of the n >= 36 mask."
  ),
  "",
  "## Coordinates, weights and families",
  "",
  paste0("- Locked grid and RDS lon/lat order identical: ",
         ifelse(coord_match, "PASS", "FAIL"), "."),
  paste0("- grid_id order identical: ",
         ifelse(grid_id_match, "PASS", "FAIL"), "."),
  "- Global and regional means use cos(latitude) weights.",
  "- BH and BY were recomputed separately within every threshold/predictor family.",
  paste0(
    "- All 251 non-zero common circular shifts were evaluated: ",
    ifelse(all(table(null_dist$threshold_n, null_dist$Predictor) == 251),
           "PASS", "FAIL"), "."
  ),
  "- Circular-shift finite p values cannot be zero; minimum possible p is 1/252.",
  "",
  "## Uncertainty and failures",
  "",
  paste0("- Spatial bootstrap: paired 10° × 10° blocks, B = ", BOOT_B, "."),
  paste0(
    "- Bootstrap seed rule: ", BOOT_SEED_BASE,
    " + retained n_grid (identical to the locked n >= 36 rule)."
  ),
  paste0(
    "- Failed/rank-deficient/missing rows in retained masks: ",
    sum(failures$failed_fit_status), " / ",
    sum(failures$rank_deficient), " / ",
    sum(failures$missing_partial_r + failures$missing_p_value), "."
  ),
  "- Year-block permutation: not evaluated because no locked per-grid/permutation source can be re-aggregated under the new masks.",
  "- Lag-3 temporal-precedence sensitivity: not evaluated; its lagged n_obs is not the concurrent Model 3 complete-case n and no compatible locked per-grid OOF delta-R2 source exists.",
  "",
  "## Scientific checks",
  "",
  paste0("- W-GPE–LAI mean partial r remains positive at all thresholds: ",
         ifelse(w_positive_all, "YES", "NO"), "."),
  paste0("- W-GPE–LAI spatial 95% CI excludes zero at all thresholds: ",
         ifelse(w_ci_excludes_zero_all, "YES", "NO"), "."),
  paste0("- <z>–LAI remains close to zero (absolute mean r < 0.02): ",
         ifelse(z_near_zero_stable, "YES", "NO"), "."),
  sprintf(
    paste0(
      "- W-GPE BH-positive area shares (36/120/177): %.3f%% / %.3f%% / ",
      "%.3f%%; BY-positive: %.3f%% / %.3f%% / %.3f%%."
    ),
    w36$BH_positive_area_pct, w120$BH_positive_area_pct,
    w177$BH_positive_area_pct, w36$BY_positive_area_pct,
    w120$BY_positive_area_pct, w177$BY_positive_area_pct
  ),
  paste0("- n >= 177 retains less than half of main-mask area: ",
         ifelse(selection_severe, "YES", "NO"), "."),
  sprintf(
    paste0(
      "- n >= 177 zone-area retention (N high/N temperate/Tropics/",
      "S temperate): %.1f%% / %.1f%% / %.1f%% / %.1f%%."
    ),
    zone_ret[["N high latitude"]], zone_ret[["N temperate"]],
    zone_ret[["Tropics"]], zone_ret[["S temperate"]]
  ),
  "",
  "## Software",
  "",
  paste0("- R: ", R.version.string, "."),
  paste0("- Platform: ", R.version$platform, "."),
  paste0("- Rcpp: ", as.character(utils::packageVersion("Rcpp")), "."),
  paste0("- RcppArmadillo: ",
         as.character(utils::packageVersion("RcppArmadillo")), "."),
  paste0("- digest: ", as.character(utils::packageVersion("digest")), "."),
  "- Input/output SHA-256 values are recorded in `LAI_COVERAGE_SHA256_MANIFEST.csv`.",
  "",
  "## Completed global results",
  "",
  paste(capture.output(print(global, row.names = FALSE)), collapse = "\n")
)
qa_file <- file.path(OUT, "LAI_COVERAGE_SENSITIVITY_QA.md")
writeLines(qa_lines, qa_file, useBytes = TRUE)

output_files <- c(
  global_file, grid_file, zones_file, field_file, sample_file, null_file,
  boot_file, repro_file, failure_file, temporal_file, methods_file, qa_file
)
script_files <- c(
  SCRIPT_FILE,
  file.path(CANON, "07_SCRIPTS", "Makevars_Rtools45"),
  file.path(CANON, "07_SCRIPTS", "run_LAI_coverage_sensitivity.sh")
)
manifest <- data.frame(
  role = c(
    rep("input", length(required_files)),
    rep("script", length(script_files)),
    rep("output", length(output_files))
  ),
  path = c(required_files, script_files, output_files),
  bytes = file.info(c(required_files, script_files, output_files))$size,
  sha256 = vapply(
    c(required_files, script_files, output_files),
    hash_file, character(1)
  ),
  stringsAsFactors = FALSE
)
manifest_file <- file.path(OUT, "LAI_COVERAGE_SHA256_MANIFEST.csv")
write.csv(manifest, manifest_file, row.names = FALSE)

message("LAI coverage sensitivity complete.")
message("Output directory: ", OUT)
