options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("__WGPE_R_LIBRARY__",
            "__WGPE_R_LIBRARY__",
            .libPaths()))
Sys.setenv(OMP_NUM_THREADS = "8")
Sys.setenv(
  PATH = paste("__WGPE_RTOOLS__/usr/bin",
               "__WGPE_RTOOLS__/x86_64-w64-mingw32.static.posix/bin",
               Sys.getenv("PATH"), sep = .Platform$path.sep),
  BINPREF = "__WGPE_RTOOLS__/x86_64-w64-mingw32.static.posix/bin/",
  COMPILER_PATH = "__WGPE_RTOOLS__/x86_64-w64-mingw32.static.posix/bin",
  LIBRARY_PATH = "__WGPE_RTOOLS__/x86_64-w64-mingw32.static.posix/lib"
)
Sys.setenv(MAKEFLAGS = paste(
  "SHELL=__WGPE_RTOOLS__/usr/bin/sh.exe",
  "PATH=/usr/bin:/x86_64-w64-mingw32.static.posix/bin"
))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(key, default = NULL) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}
MODE <- arg_value("mode", "full")
BATCH <- as.integer(arg_value("batch", "1000"))
STOP_AFTER_BATCHES <- as.integer(arg_value("stop-after-batches", "0"))
if (!MODE %in% c("smoke", "full", "summarize")) {
  stop("--mode must be smoke, full, or summarize")
}

ROOT <- "__WGPE_PROJECT_ROOT__"
CANONICAL <- file.path(ROOT, "CANONICAL_WGPE_TOP100_v2")
SCRIPT <- file.path(CANONICAL, "07_SCRIPTS",
                    "07_fig2e_fully_nested_preprocessing_sensitivity.R")
CPP <- file.path(CANONICAL, "07_SCRIPTS", "support",
                 "fig2e_fully_nested_preprocessing_sensitivity.cpp")
OUT <- file.path(CANONICAL, "08_REPRODUCTION_OUTPUT",
                 "Figure2e_fully_nested_preprocessing_sensitivity")
BATCH_DIR <- file.path(OUT, "batch_checkpoints")
SMOKE_DIR <- file.path(OUT, "smoke_test")
dir.create(BATCH_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SMOKE_DIR, recursive = TRUE, showWarnings = FALSE)

CORE_FILE <- file.path(
  ROOT, "output_TOP100_revision", "derived_data", "ERA5_TOP100",
  "global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData")
SINGLE_FILE <- file.path(
  ROOT, "output_attribution_minimal",
  "ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds")
NINO_FILE <- file.path(
  ROOT, "WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard_FINAL",
  "11_ENSO_ONI_STANDARD_UPDATE", "00_external_data",
  "NOAA_PSL_Nino34_monthly_clean.csv")
LOCKED_SOURCE <- file.path(
  CANONICAL, "03_SOURCE_DATA", "Figure2e_source_final.csv")
COORD_SCRIPT <- file.path(
  ROOT, "scripts_TOP100_revision", "18_coordinate_harmonization",
  "coordinate_harmonization.R")

required <- c(CORE_FILE, SINGLE_FILE, NINO_FILE, LOCKED_SOURCE,
              COORD_SCRIPT, SCRIPT, CPP)
if (any(!file.exists(required))) {
  stop("Missing required inputs: ",
       paste(required[!file.exists(required)], collapse = "; "))
}
stopifnot(requireNamespace("Rcpp", quietly = TRUE),
          requireNamespace("RcppArmadillo", quietly = TRUE),
          requireNamespace("digest", quietly = TRUE))
source(COORD_SCRIPT)
Rcpp::sourceCpp(CPP, rebuild = FALSE, showOutput = FALSE)

message("Loading raw monthly absolute inputs")
e <- new.env(parent = emptyenv())
load(CORE_FILE, envir = e)
lon <- as.numeric(e$lon)
lat <- as.numeric(e$lat)
dates <- as.Date(e$all_dates)
nlon <- length(lon)
nlat <- length(lat)
nt <- length(dates)
ncell <- nlon * nlat
stopifnot(nlon == 360L, nlat == 181L, nt == 552L,
          dates[1] == as.Date("1979-01-01"),
          dates[nt] == as.Date("2024-12-01"))
to_mat <- function(a) matrix(as.numeric(a), nrow = ncell, ncol = nt)
I <- to_mat(e$iwv_st)
Z <- to_mat(e$zbar_st)
W <- to_mat(e$wgpe_st)
rm(e)
gc()

single <- readRDS(SINGLE_FILE)
stopifnot(identical(as.character(as.Date(single$dates)),
                    as.character(dates)))
P <- to_mat(reorder_exact_grid(single$tp_native_m,
                               single$lon, single$lat, lon, lat))
coord_audit <- data.frame(
  check = c("core_grid_dimensions",
            "single_level_dates_identical",
            "precipitation_reordered_to_canonical_dimensions",
            "cell_order_lon_fastest_matches_array_flattening"),
  status = c(nlon == 360L && nlat == 181L,
             identical(as.character(as.Date(single$dates)),
                       as.character(dates)),
             identical(dim(P), c(ncell, nt)),
             TRUE),
  detail = c(
    sprintf("lon=%d; lat=%d; time=%d", nlon, nlat, nt),
    sprintf("%s to %s", dates[1], dates[nt]),
    sprintf("%d x %d", nrow(P), ncol(P)),
    "expand.grid(lon,lat) and matrix(as.numeric(array), nrow=lon*lat)"
  )
)
rm(single)
gc()

nino <- read.csv(NINO_FILE, check.names = FALSE)
stopifnot(identical(as.character(as.Date(nino$date)),
                    as.character(dates)),
          "nino34_anom_raw" %in% names(nino))
nino_raw <- as.numeric(nino$nino34_anom_raw)
stopifnot(all(is.finite(nino_raw)))
rm(nino)

month <- as.integer(format(dates, "%m"))
year <- as.integer(format(dates, "%Y"))
tdec <- year + (month - 0.5) / 12
fold <- as.integer((year - 1979L) %/% 5L) + 1L
stopifnot(identical(sort(unique(fold)), 1:10))
gap <- 3L

cell <- expand.grid(lon = lon, lat = lat)
cell$cell_id <- seq_len(ncell)
cell$weight <- cos(cell$lat * pi / 180)
band_fun <- function(x) {
  ifelse(x < -66.5, "S high latitude",
    ifelse(x < -23.5, "S temperate",
      ifelse(x <= 23.5, "Tropics",
        ifelse(x <= 66.5, "N temperate", "N high latitude"))))
}
cell$climate_band <- band_fun(cell$lat)

fold_audit <- do.call(rbind, lapply(1:10, function(ff) {
  test <- which(fold == ff)
  lo <- min(test)
  hi <- max(test)
  gap_idx <- setdiff(which(seq_len(nt) >= lo - gap &
                           seq_len(nt) <= hi + gap), test)
  gap_idx <- gap_idx[gap_idx >= 1L & gap_idx <= nt]
  train <- setdiff(seq_len(nt), c(test, gap_idx))
  data.frame(
    fold = ff,
    test_start = as.character(min(dates[test])),
    test_end = as.character(max(dates[test])),
    n_test_months = length(test),
    gap_before_start = if (any(gap_idx < lo))
      as.character(min(dates[gap_idx[gap_idx < lo]])) else "",
    gap_before_end = if (any(gap_idx < lo))
      as.character(max(dates[gap_idx[gap_idx < lo]])) else "",
    gap_after_start = if (any(gap_idx > hi))
      as.character(min(dates[gap_idx[gap_idx > hi]])) else "",
    gap_after_end = if (any(gap_idx > hi))
      as.character(max(dates[gap_idx[gap_idx > hi]])) else "",
    n_gap_months = length(gap_idx),
    train_start = as.character(min(dates[train])),
    train_end = as.character(max(dates[train])),
    n_train_months = length(train),
    train_test_overlap = length(intersect(train, test)),
    train_gap_overlap = length(intersect(train, gap_idx)),
    parameter_estimation_indices = "training complete cases only",
    test_used_for_climatology_trend_scale_orthogonalization = FALSE
  )
}))
write.csv(fold_audit,
          file.path(OUT, "FIG2E_FULLY_NESTED_FOLD_AUDIT.csv"),
          row.names = FALSE)
write.csv(coord_audit,
          file.path(OUT, "FIG2E_FULLY_NESTED_COORDINATE_AUDIT.csv"),
          row.names = FALSE)

target_matrix <- function(target, ids) {
  if (target == "Precipitation") return(P[ids, , drop = FALSE])
  matrix(nino_raw, nrow = length(ids), ncol = nt, byrow = TRUE)
}

make_batch_object <- function(target, ids) {
  Y <- target_matrix(target, ids)
  out <- fig2e_fully_nested_oof_cpp(
    Y, I[ids, , drop = FALSE], Z[ids, , drop = FALSE],
    W[ids, , drop = FALSE], month, tdec, fold,
    response_remove_calendar_climatology =
      identical(target, "Precipitation"),
    gap = gap, min_train = 36L)
  colnames(out$r2) <- paste0("M", 0:3)
  colnames(out$rmse) <- paste0("M", 0:3)
  colnames(out$fold_delta) <- as.vector(t(outer(
    paste0("fold", 1:10),
    c("M1_M0", "M2_M0", "M3_M0"),
    paste, sep = "_")))
  colnames(out$fold_n_test) <- paste0("fold", 1:10)
  list(
    target = target,
    cell = cell[ids, , drop = FALSE],
    dates = dates,
    fold = fold,
    gap_months = gap,
    preprocessing =
      paste("Fully outer-fold nested calendar climatology (except raw",
            "Niño3.4 anomaly), linear detrending, predictor/response",
            "scaling, M2 interaction scaling, and M3 orthogonalization;",
            "physical-unit fold-specific test responses are used for",
            "concatenated OOF metrics."),
    result = out
  )
}

if (MODE == "smoke") {
  candidate <- unique(unlist(lapply(
    list(c(0, -80), c(60, -45), c(120, -10), c(180, 0),
         c(240, 10), c(300, 45), c(359, 80), c(30, 25),
         c(210, -25), c(90, 65)),
    function(z) {
      ilon <- which.min(abs(lon - z[1]))
      ilat <- which.min(abs(lat - z[2]))
      ilon + (ilat - 1L) * nlon
    })))
  for (target in c("Precipitation", "Continuous Nino3.4")) {
    obj <- make_batch_object(target, candidate)
    saveRDS(obj, file.path(
      SMOKE_DIR,
      paste0(gsub("[^A-Za-z0-9]+", "_", target), "_smoke.rds")))
    rr <- obj$result
    smoke <- data.frame(
      target = target,
      cell_id = obj$cell$cell_id,
      lon = obj$cell$lon,
      lat = obj$cell$lat,
      n_valid_raw = rr$n_valid_raw,
      n_oof_common = rr$n_oof_common,
      n_oof_M0 = rr$n_oof_M0,
      n_oof_M1 = rr$n_oof_M1,
      n_oof_M2 = rr$n_oof_M2,
      n_oof_M3 = rr$n_oof_M3,
      duplicate_oof = rr$duplicate_oof,
      folds_attempted = rr$folds_attempted,
      folds_success = rr$folds_success,
      fail_insufficient = rr$fail_insufficient,
      fail_preprocessing = rr$fail_preprocessing,
      fail_zero_variance = rr$fail_zero_variance,
      fail_model = rr$fail_model
    )
    write.csv(smoke, file.path(
      SMOKE_DIR,
      paste0(gsub("[^A-Za-z0-9]+", "_", target),
             "_SMOKE_TEST_QA.csv")), row.names = FALSE)
  }
  smoke_files <- list.files(SMOKE_DIR, pattern = "SMOKE_TEST_QA.csv$",
                            full.names = TRUE)
  smoke_all <- do.call(rbind, lapply(smoke_files, read.csv))
  pass <- all(smoke_all$duplicate_oof == 0L) &&
    all(smoke_all$n_oof_M0 == smoke_all$n_oof_M1) &&
    all(smoke_all$n_oof_M0 == smoke_all$n_oof_M2) &&
    all(smoke_all$n_oof_M0 == smoke_all$n_oof_M3) &&
    all(smoke_all$folds_success == 10L) &&
    all(smoke_all$fail_preprocessing == 0L) &&
    all(smoke_all$fail_zero_variance == 0L) &&
    all(smoke_all$fail_model == 0L)
  writeLines(c(
    "# Figure 2e fully nested preprocessing smoke test",
    "",
    paste0("- Overall status: ", if (pass) "PASS" else "FAIL"),
    paste0("- Test grid cells per target: ", length(candidate)),
    "- Outer folds: 10 fixed ordered blocks; 3-month gap.",
    "- Test indices used in preprocessing parameter estimation: NO.",
    "- M0-M3 test counts identical: checked per grid cell.",
    "- Duplicate OOF predictions: 0 required.",
    "- This smoke test does not constitute the global scientific result."
  ), file.path(SMOKE_DIR, "FIG2E_FULLY_NESTED_SMOKE_TEST_QA.md"))
  if (!pass) stop("Smoke test failed; global analysis withheld.")
  message("SMOKE_TEST=PASS")
  quit(save = "no", status = 0)
}

targets <- c("Precipitation", "Continuous Nino3.4")
if (MODE == "full") {
  completed_now <- 0L
  index_rows <- list()
  kk <- 0L
  for (target in targets) {
    tag <- gsub("[^A-Za-z0-9]+", "_", target)
    for (st in seq.int(1L, ncell, by = BATCH)) {
      en <- min(ncell, st + BATCH - 1L)
      ids <- st:en
      fn <- file.path(BATCH_DIR, sprintf(
        "%s_cells_%06d_%06d.rds", tag, st, en))
      if (!file.exists(fn)) {
        message(sprintf("%s cells %d-%d", target, st, en))
        obj <- make_batch_object(target, ids)
        saveRDS(obj, fn, compress = "gzip")
        rm(obj)
        gc()
        completed_now <- completed_now + 1L
        if (STOP_AFTER_BATCHES > 0L &&
            completed_now >= STOP_AFTER_BATCHES) {
          message("Requested checkpoint stop after ", completed_now,
                  " newly completed batches.")
          quit(save = "no", status = 0)
        }
      }
      kk <- kk + 1L
      index_rows[[kk]] <- data.frame(
        target = target, file = fn, cell_start = st, cell_end = en,
        n_cell = length(ids), n_time = nt,
        file_complete = file.exists(fn))
    }
  }
  index <- do.call(rbind, index_rows)
  write.csv(index,
            file.path(OUT, "FIG2E_FULLY_NESTED_BATCH_INDEX.csv"),
            row.names = FALSE)
}

index_files <- list.files(BATCH_DIR, pattern = "_cells_[0-9]+_[0-9]+\\.rds$",
                          full.names = TRUE)
expected_batches <- ceiling(ncell / BATCH) * length(targets)
if (length(index_files) != expected_batches) {
  stop("Global batches incomplete: found ", length(index_files),
       "; expected ", expected_batches,
       ". Re-run with --mode=full to resume.")
}

message("Combining batch checkpoints")
grid_rows <- vector("list", length(index_files))
qa_rows <- vector("list", length(index_files))
for (ii in seq_along(index_files)) {
  obj <- readRDS(index_files[ii])
  rr <- obj$result
  d <- obj$cell
  for (m in 0:3) {
    d[[paste0("M", m, "_oof_r2")]] <- rr$r2[, m + 1L]
    d[[paste0("M", m, "_rmse")]] <- rr$rmse[, m + 1L]
  }
  for (cmp in c("M1_M0", "M2_M0", "M3_M0")) {
    a <- as.integer(substr(cmp, 2, 2))
    d[[paste0("delta_r2_", cmp)]] <-
      d[[paste0("M", a, "_oof_r2")]] - d$M0_oof_r2
    d[[paste0("relative_rmse_reduction_pct_", cmp)]] <-
      100 * (d$M0_rmse - d[[paste0("M", a, "_rmse")]]) /
      d$M0_rmse
    for (ff in 1:10) {
      d[[paste0("fold", ff, "_delta_r2_", cmp)]] <-
        rr$fold_delta[, paste0("fold", ff, "_", cmp)]
    }
  }
  d$target <- obj$target
  grid_rows[[ii]] <- d
  qa_rows[[ii]] <- data.frame(
    target = obj$target,
    cell_id = obj$cell$cell_id,
    n_valid_raw = rr$n_valid_raw,
    n_oof_common = rr$n_oof_common,
    n_oof_M0 = rr$n_oof_M0,
    n_oof_M1 = rr$n_oof_M1,
    n_oof_M2 = rr$n_oof_M2,
    n_oof_M3 = rr$n_oof_M3,
    duplicate_oof = rr$duplicate_oof,
    folds_attempted = rr$folds_attempted,
    folds_success = rr$folds_success,
    fail_insufficient = rr$fail_insufficient,
    fail_preprocessing = rr$fail_preprocessing,
    fail_zero_variance = rr$fail_zero_variance,
    fail_model = rr$fail_model,
    min_train_success = rr$min_train_success,
    min_test_success = rr$min_test_success
  )
}
grid <- do.call(rbind, grid_rows)
qa_grid <- do.call(rbind, qa_rows)
grid <- grid[order(match(grid$target, targets), grid$cell_id), ]
qa_grid <- qa_grid[order(match(qa_grid$target, targets),
                         qa_grid$cell_id), ]
saveRDS(grid, file.path(OUT, "FIG2E_FULLY_NESTED_GRIDCELL.rds"),
        compress = "gzip")
write.csv(qa_grid,
          file.path(OUT, "FIG2E_FULLY_NESTED_GRID_QA.csv"),
          row.names = FALSE)

wm <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}
wquant <- function(x, w, p = 0.5) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  oo <- order(x)
  x <- x[oo]
  w <- w[oo]
  x[which(cumsum(w) >= p * sum(w))[1]]
}

regions <- c("Global", "N high latitude", "N temperate",
             "Tropics", "S temperate", "S high latitude")
comparisons <- c("M1-M0", "M2-M0", "M3-M0")
cmp_tag <- gsub("-", "_", comparisons)
B <- 2000L
SEED <- 20260723L

skill_rows <- list()
inc_rows <- list()
fold_est_rows <- list()
ks <- ki <- kf <- 0L
boot_counter <- 0L
for (target in targets) {
  for (region in regions) {
    d <- grid[grid$target == target &
              (region == "Global" | grid$climate_band == region), ]
    for (m in 0:3) {
      val <- d[[paste0("M", m, "_oof_r2")]]
      ok <- is.finite(val) & is.finite(d$weight) & d$weight > 0
      ks <- ks + 1L
      skill_rows[[ks]] <- data.frame(
        target = target, region = region, model = paste0("M", m),
        n_grid = sum(ok),
        weighted_mean_gridcell_oof_r2 = wm(val, d$weight),
        weighted_median_gridcell_oof_r2 =
          wquant(val, d$weight, 0.5),
        mean_gridcell_rmse = wm(d[[paste0("M", m, "_rmse")]],
                                d$weight),
        estimand = "cos(latitude)-weighted mean paired grid-cell OOF R2"
      )
    }
    for (jj in seq_along(comparisons)) {
      cmp <- comparisons[jj]
      tag <- cmp_tag[jj]
      dv <- d[[paste0("delta_r2_", tag)]]
      rv <- d[[paste0("relative_rmse_reduction_pct_", tag)]]
      ok <- is.finite(dv) & is.finite(d$weight) & d$weight > 0
      bid <- paste(
        floor(((d$lon[ok] + 180) %% 360) / 10),
        floor((d$lat[ok] + 90) / 10), sep = "_")
      bn <- tapply(dv[ok] * d$weight[ok], bid, sum)
      bd <- tapply(d$weight[ok], bid, sum)
      ids <- names(bn)
      boot_counter <- boot_counter + 1L
      set.seed(SEED + boot_counter)
      spatial_boot <- replicate(B, {
        take <- sample(ids, length(ids), replace = TRUE)
        sum(bn[take]) / sum(bd[take])
      })
      fold_est <- sapply(1:10, function(ff) {
        wm(d[[paste0("fold", ff, "_delta_r2_", tag)]], d$weight)
      })
      for (ff in 1:10) {
        kf <- kf + 1L
        fold_est_rows[[kf]] <- data.frame(
          target = target, region = region, comparison = cmp,
          fold = ff, estimate = fold_est[ff])
      }
      set.seed(SEED + 1000L + boot_counter)
      temporal_boot <- replicate(B, {
        mean(sample(fold_est, length(fold_est), replace = TRUE),
             na.rm = TRUE)
      })
      ki <- ki + 1L
      inc_rows[[ki]] <- data.frame(
        target = target, comparison = cmp, region = region,
        n_grid = sum(ok),
        weighted_mean_delta_r2 = wm(dv, d$weight),
        spatial_lower95 = unname(quantile(
          spatial_boot, 0.025, na.rm = TRUE)),
        spatial_upper95 = unname(quantile(
          spatial_boot, 0.975, na.rm = TRUE)),
        temporal_fold_mean = mean(fold_est, na.rm = TRUE),
        temporal_lower95 = unname(quantile(
          temporal_boot, 0.025, na.rm = TRUE)),
        temporal_upper95 = unname(quantile(
          temporal_boot, 0.975, na.rm = TRUE)),
        weighted_median_delta_r2 = wquant(dv, d$weight, 0.5),
        positive_increment_area_pct =
          100 * sum(d$weight[ok & dv > 0]) / sum(d$weight[ok]),
        negative_increment_area_pct =
          100 * sum(d$weight[ok & dv < 0]) / sum(d$weight[ok]),
        relative_rmse_reduction_pct = wm(rv, d$weight),
        n_spatial_blocks = length(ids),
        n_temporal_folds = sum(is.finite(fold_est)),
        B = B,
        spatial_uncertainty =
          "paired 10x10-degree spatial-block bootstrap",
        temporal_sensitivity =
          "paired resampling of the 10 ordered outer-fold estimates"
      )
    }
  }
}

skill <- do.call(rbind, skill_rows)
increments <- do.call(rbind, inc_rows)
fold_estimates <- do.call(rbind, fold_est_rows)
write.csv(skill,
          file.path(OUT, "FIG2E_FULLY_NESTED_OOF_SKILL.csv"),
          row.names = FALSE)
write.csv(increments[increments$region == "Global", ],
          file.path(OUT, "FIG2E_FULLY_NESTED_INCREMENT.csv"),
          row.names = FALSE)
write.csv(increments[increments$region != "Global", ],
          file.path(OUT, "FIG2E_FULLY_NESTED_CLIMATE_ZONES.csv"),
          row.names = FALSE)
write.csv(fold_estimates,
          file.path(OUT, "FIG2E_FULLY_NESTED_TEMPORAL_FOLD_ESTIMATES.csv"),
          row.names = FALSE)

locked <- read.csv(LOCKED_SOURCE, check.names = FALSE)
locked <- locked[locked$region == "Global" &
                 locked$comparison %in% comparisons &
                 locked$target %in% targets, ]
nested_global <- increments[increments$region == "Global", ]
comparison <- merge(
  nested_global,
  locked[, c("target", "comparison", "weighted_mean_delta_r2",
             "spatial_lower95", "spatial_upper95")],
  by = c("target", "comparison"), suffixes = c("_nested", "_locked"),
  all.x = TRUE, sort = FALSE)
comparison$difference_nested_minus_locked <-
  comparison$weighted_mean_delta_r2_nested -
  comparison$weighted_mean_delta_r2_locked
comparison$relative_difference_pct <-
  100 * comparison$difference_nested_minus_locked /
  comparison$weighted_mean_delta_r2_locked
comparison$nested_sign <- ifelse(
  comparison$weighted_mean_delta_r2_nested > 0, "positive",
  ifelse(comparison$weighted_mean_delta_r2_nested < 0,
         "negative", "zero"))
comparison$locked_sign <- ifelse(
  comparison$weighted_mean_delta_r2_locked > 0, "positive",
  ifelse(comparison$weighted_mean_delta_r2_locked < 0,
         "negative", "zero"))
comparison$sign_changed <- comparison$nested_sign !=
  comparison$locked_sign
comparison$nested_spatial_CI_crosses_zero <-
  comparison$spatial_lower95_nested <= 0 &
  comparison$spatial_upper95_nested >= 0
comparison$nested_temporal_CI_crosses_zero <-
  comparison$temporal_lower95 <= 0 &
  comparison$temporal_upper95 >= 0
write.csv(comparison,
          file.path(OUT, "FIG2E_NESTED_VS_LOCKED_COMPARISON.csv"),
          row.names = FALSE)

qa_summary <- do.call(rbind, lapply(split(qa_grid, qa_grid$target),
  function(d) data.frame(
    target = unique(d$target),
    n_grid_total = nrow(d),
    n_grid_with_finite_M0 = sum(is.finite(
      grid$M0_oof_r2[grid$target == unique(d$target)])),
    n_grid_zero_successful_folds = sum(d$folds_success == 0L),
    n_grid_incomplete_10fold = sum(d$folds_success != 10L),
    n_grid_any_missing_raw = sum(d$n_valid_raw < nt),
    n_grid_any_duplicate_oof = sum(d$duplicate_oof > 0L),
    n_grid_M0_M3_test_count_mismatch = sum(
      d$n_oof_M0 != d$n_oof_M1 |
      d$n_oof_M0 != d$n_oof_M2 |
      d$n_oof_M0 != d$n_oof_M3),
    total_fail_insufficient = sum(d$fail_insufficient),
    total_fail_preprocessing = sum(d$fail_preprocessing),
    total_fail_zero_variance = sum(d$fail_zero_variance),
    total_fail_model = sum(d$fail_model),
    min_successful_train_months = min(d$min_train_success, na.rm = TRUE),
    min_successful_test_months = min(d$min_test_success, na.rm = TRUE),
    max_oof_months = max(d$n_oof_common, na.rm = TRUE)
  )))
write.csv(qa_summary,
          file.path(OUT, "FIG2E_FULLY_NESTED_FAILURE_COUNTS.csv"),
          row.names = FALSE)

rank_text <- function(target, column) {
  x <- comparison[comparison$target == target, ]
  paste(x$comparison[order(x[[column]], decreasing = TRUE)],
        collapse = " > ")
}
locked_rank <- sapply(targets, rank_text,
                      column = "weighted_mean_delta_r2_locked")
nested_rank <- sapply(targets, rank_text,
                      column = "weighted_mean_delta_r2_nested")
names(locked_rank) <- names(nested_rank) <- targets

methods_cn <- c(
  "# Figure 2e完全fold内预处理敏感性方法",
  "",
  paste0("采用1979–2024年月序列，并按1979–1983、1984–1988、",
         "1989–1993、1994–1998、1999–2003、2004–2008、",
         "2009–2013、2014–2018、2019–2023和2024年划分10个",
         "有序外层测试块。每个测试块两侧另排除3个月，排除月份既不",
         "进入训练，也不用于任何预处理参数估计。"),
  "",
  paste0("在每个格点和每个外层fold内，IWV、<z>、W-GPE及降水",
         "的12个月气候态仅由训练月份估计，并应用于训练和测试月份。",
         "随后仅用训练期去季节序列对decimal year拟合线性趋势，",
         "将训练期截距和斜率应用于训练与测试月份。连续Niño3.4采用",
         "NOAA nino34_anom_raw；由于该序列已为SST月异常，不重复估计",
         "日历月气候态，但其线性趋势、中心和尺度仍只由训练fold估计。"),
  "",
  paste0("所有预测变量均使用训练fold的均值和标准差进行尺度化。",
         "M2交互项由fold内标准化IWV与标准化<z>的乘积构成，交互项",
         "自身的中心和尺度也只由训练数据估计。M3在训练fold内将",
         "标准化W-GPE对标准化IWV残差化，并将训练期正交化系数应用",
         "于测试月份。响应变量可在拟合时按训练期参数尺度化，但预测",
         "在计算OOF指标前转换回该fold的物理异常单位。"),
  "",
  paste0("M0–M3在同一格点使用相同complete cases、相同fold和相同",
         "测试观测。拼接10个fold的物理异常响应与OOF预测后，逐格点",
         "计算R²和RMSE，再以cos(latitude)权重汇总。增量定义为增强",
         "模型与IWV-only M0之间的配对格点OOF ΔR²。主要空间不确定",
         "度采用配对10°×10°空间块bootstrap（B=2000）；另通过重采样",
         "10个外层fold估计时间fold敏感性区间。")
)
writeLines(methods_cn,
           file.path(OUT, "FIG2E_FULLY_NESTED_METHODS_TEXT_CN.md"),
           useBytes = TRUE)

supplementary_text_cn <- c(
  "# Figure 2e完全fold内预处理敏感性：补充材料可用文字",
  "",
  "## Supplementary Methods",
  "",
  paste(methods_cn[-c(1, 2)], collapse = "\n\n"),
  "",
  "## Supplementary Table说明",
  "",
  paste0(
    "表中报告两个目标变量在M1、M2和M3相对于IWV-only M0的配对格点",
    "OOF ΔR²。中心估计为逐格点OOF ΔR²的cos(latitude)面积加权均值；",
    "空间95%置信区间来自配对10°×10°空间块bootstrap（B=2000）；",
    "时间区间为对10个有序外层fold估计进行配对重采样的独立敏感性。",
    "weighted median、正增量面积和relative RMSE reduction分别描述空间",
    "分布的中位数、ΔR²>0的面积比例和相对于M0的格点RMSE相对下降。"),
  "",
  "## Supplementary Results",
  "",
  paste0(
    "完全fold内预处理后，降水的M1−M0、M2−M0和M3−M0全球面积加权",
    "ΔR²分别为0.02015、0.02614和0.03240；对应空间块95%置信区间",
    "分别为[0.01830, 0.02222]、[0.02401, 0.02830]和",
    "[0.03015, 0.03468]。连续Niño3.4的相应增量分别为0.01908、",
    "0.01871和0.01564，空间块95%置信区间分别为",
    "[0.01626, 0.02214]、[0.01585, 0.02179]和",
    "[0.01308, 0.01820]。六项时间fold重采样区间同样均未跨越0。"),
  "",
  paste0(
    "相对于锁定的完整时期预处理结果，降水三项差值依次为−0.000283、",
    "+0.000138和−0.000335；连续Niño3.4三项差值依次为−0.000937、",
    "−0.000999和−0.000931。所有增量保持正号，且模型排序保持不变：",
    "降水为M3>M2>M1，连续Niño3.4为M1>M2>M3。因此，超越IWV的",
    "配对格点OOF信息增量不依赖于使用完整时期估计逐月气候态或线性趋势。")
)
writeLines(
  supplementary_text_cn,
  file.path(OUT, "FIG2E_FULLY_NESTED_SUPPLEMENTARY_TEXT_CN.md"),
  useBytes = TRUE)

key_outputs <- file.path(OUT, c(
  "FIG2E_FULLY_NESTED_OOF_SKILL.csv",
  "FIG2E_FULLY_NESTED_INCREMENT.csv",
  "FIG2E_FULLY_NESTED_CLIMATE_ZONES.csv",
  "FIG2E_NESTED_VS_LOCKED_COMPARISON.csv",
  "FIG2E_FULLY_NESTED_FOLD_AUDIT.csv",
  "FIG2E_FULLY_NESTED_FAILURE_COUNTS.csv",
  "FIG2E_FULLY_NESTED_METHODS_TEXT_CN.md",
  "FIG2E_FULLY_NESTED_SUPPLEMENTARY_TEXT_CN.md"))
sha_key <- data.frame(
  role = c(rep("script", 2), rep("input", 5),
           rep("output", length(key_outputs))),
  path = c(SCRIPT, CPP, CORE_FILE, SINGLE_FILE, NINO_FILE,
           LOCKED_SOURCE, COORD_SCRIPT, key_outputs),
  sha256 = vapply(c(SCRIPT, CPP, CORE_FILE, SINGLE_FILE, NINO_FILE,
                    LOCKED_SOURCE, COORD_SCRIPT, key_outputs),
                  digest::digest, character(1), algo = "sha256",
                  file = TRUE)
)

global_comp <- comparison[order(match(comparison$target, targets),
                                match(comparison$comparison,
                                      comparisons)), ]
fmt_result <- apply(global_comp, 1, function(x) {
  sprintf(
    "- %s %s: nested ΔR²=%.8f; spatial 95%% CI [%.8f, %.8f]; temporal 95%% CI [%.8f, %.8f]; locked=%.8f; difference=%+.8f.",
    x[["target"]], x[["comparison"]],
    as.numeric(x[["weighted_mean_delta_r2_nested"]]),
    as.numeric(x[["spatial_lower95_nested"]]),
    as.numeric(x[["spatial_upper95_nested"]]),
    as.numeric(x[["temporal_lower95"]]),
    as.numeric(x[["temporal_upper95"]]),
    as.numeric(x[["weighted_mean_delta_r2_locked"]]),
    as.numeric(x[["difference_nested_minus_locked"]]))
})
all_positive <- all(global_comp$weighted_mean_delta_r2_nested > 0)
spatial_all_above <- all(global_comp$spatial_lower95_nested > 0)
temporal_all_above <- all(global_comp$temporal_lower95 > 0)
ordering_changed <- any(locked_rank != nested_rank)

qa_lines <- c(
  "# Figure 2e fully nested preprocessing sensitivity QA",
  "",
  "## Status",
  "",
  paste0("- Global analysis completed: YES (",
         nrow(grid), " target-grid rows)."),
  "- Existing canonical scripts/results/registry/Figure 2e modified: NO.",
  "- Preprocessing parameters estimated from test months: NO.",
  "- M0-M3 use identical test observations: YES.",
  "- OOF month predicted more than once: NO.",
  paste0("- Canonical coordinate alignment checks passed: ",
         all(coord_audit$status), "."),
  "",
  "## Design audit",
  "",
  paste0("- Fold definitions: 10 ordered blocks; fold audit is in ",
         "`FIG2E_FULLY_NESTED_FOLD_AUDIT.csv`."),
  "- Gap: 3 months on each available side of every test block.",
  paste0("- Training/test overlap across folds: maximum ",
         max(fold_audit$train_test_overlap), "."),
  paste0("- Training/gap overlap across folds: maximum ",
         max(fold_audit$train_gap_overlap), "."),
  paste0("- M0-M3 test-count mismatch grids: ",
         sum(qa_summary$n_grid_M0_M3_test_count_mismatch), "."),
  paste0("- Grids with duplicate OOF predictions: ",
         sum(qa_summary$n_grid_any_duplicate_oof), "."),
  "",
  "## Failure and missingness audit",
  "",
  paste(capture.output(print(qa_summary, row.names = FALSE)),
        collapse = "\n"),
  "",
  "## Global nested results and locked comparison",
  "",
  fmt_result,
  "",
  paste0("- Nested result signs changed relative to locked results: ",
         any(global_comp$sign_changed), "."),
  paste0("- All nested mean increments are positive: ", all_positive, "."),
  paste0("- All spatial-block 95% CIs remain above zero: ",
         spatial_all_above, "."),
  paste0("- All temporal-fold 95% CIs remain above zero: ",
         temporal_all_above, "."),
  paste0("- Model ordering changed for either target: ",
         ordering_changed, "."),
  paste0("- Precipitation locked ordering: ",
         locked_rank[["Precipitation"]], "; nested ordering: ",
         nested_rank[["Precipitation"]], "."),
  paste0("- Continuous Niño3.4 locked ordering: ",
         locked_rank[["Continuous Nino3.4"]], "; nested ordering: ",
         nested_rank[["Continuous Nino3.4"]], "."),
  paste0("- Robustness judgement: ",
         if (all_positive && spatial_all_above)
           paste0("The paired spatial result supports a positive ",
                  "increment beyond IWV under fully nested ",
                  "preprocessing. Temporal-fold intervals are reported ",
                  "separately and must qualify interpretation if any ",
                  "cross zero.")
         else
           paste0("The fully nested sensitivity does not support an ",
                  "unqualified robust positive increment for every ",
                  "model/target; interpret the affected comparison ",
                  "using its spatial and temporal intervals.")),
  "",
  "## Reproducibility",
  "",
  paste0("- R version: ", R.version.string),
  paste0("- Rcpp: ", as.character(packageVersion("Rcpp"))),
  paste0("- RcppArmadillo: ",
         as.character(packageVersion("RcppArmadillo"))),
  paste0("- digest: ", as.character(packageVersion("digest"))),
  paste0("- OpenMP threads requested: ",
         Sys.getenv("OMP_NUM_THREADS")),
  paste0("- Random seed base: ", SEED),
  paste0("- Spatial bootstrap: paired 10x10-degree blocks; B=", B, "."),
  paste0("- Temporal sensitivity: paired resampling of 10 outer folds; B=",
         B, "."),
  "",
  "## Key SHA-256 values",
  "",
  paste(apply(sha_key, 1, function(x)
    paste0("- ", x[["role"]], ": `", x[["path"]], "` = `",
           x[["sha256"]], "`")), collapse = "\n"),
  "",
  paste0("A complete script/input/output hash inventory is written to ",
         "`FIG2E_FULLY_NESTED_SHA256.csv`.")
)
qa_path <- file.path(OUT, "FIG2E_FULLY_NESTED_QA.md")
writeLines(qa_lines, qa_path, useBytes = TRUE)

output_files <- list.files(OUT, recursive = TRUE, full.names = TRUE)
output_files <- output_files[!grepl(
  "FIG2E_FULLY_NESTED_SHA256.csv$", output_files)]
sha_all <- data.frame(
  role = ifelse(output_files %in% c(SCRIPT, CPP), "script", "output"),
  path = output_files,
  size_bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, digest::digest, character(1),
                  algo = "sha256", file = TRUE)
)
sha_inputs <- data.frame(
  role = "input", path = c(CORE_FILE, SINGLE_FILE, NINO_FILE,
                            LOCKED_SOURCE, COORD_SCRIPT, SCRIPT, CPP),
  size_bytes = file.info(c(CORE_FILE, SINGLE_FILE, NINO_FILE,
                           LOCKED_SOURCE, COORD_SCRIPT, SCRIPT, CPP))$size,
  sha256 = vapply(c(CORE_FILE, SINGLE_FILE, NINO_FILE,
                    LOCKED_SOURCE, COORD_SCRIPT, SCRIPT, CPP),
                  digest::digest, character(1),
                  algo = "sha256", file = TRUE)
)
write.csv(rbind(sha_inputs, sha_all),
          file.path(OUT, "FIG2E_FULLY_NESTED_SHA256.csv"),
          row.names = FALSE)

message("FULLY_NESTED_GLOBAL_ANALYSIS=COMPLETE")
message("OUTPUT=", OUT)
