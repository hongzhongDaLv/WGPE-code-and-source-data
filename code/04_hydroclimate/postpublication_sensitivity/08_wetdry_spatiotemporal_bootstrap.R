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
SCRIPT <- file.path(
  CANON, "07_SCRIPTS", "08_wetdry_spatiotemporal_bootstrap.R"
)
CPP <- file.path(
  CANON, "07_SCRIPTS", "support",
  "wetdry_spatiotemporal_bootstrap.cpp"
)
MAKEVARS_FILE <- file.path(CANON, "07_SCRIPTS", "Makevars_Rtools45")
Sys.setenv(R_MAKEVARS_USER = MAKEVARS_FILE)

OUT <- file.path(
  CANON, "08_REPRODUCTION_OUTPUT",
  "wetdry_spatiotemporal_uncertainty"
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
LOCKED_SUMMARY <- file.path(
  ROOT, "output_TOP100_revision", "results", "precip_wetdry",
  "exact_TOP100", "WET_DRY_EXACT_SUMMARY_TOP100.csv"
)

required <- c(
  SCRIPT, CPP, MAKEVARS_FILE, CORE_FILE, SINGLE_FILE, COORD_FILE,
  LOCKED_SUMMARY
)
if (!all(file.exists(required))) {
  stop(
    "Missing required file(s): ",
    paste(required[!file.exists(required)], collapse = "; ")
  )
}
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")

source(COORD_FILE)
Rcpp::sourceCpp(CPP, rebuild = TRUE, showOutput = FALSE)

G <- 9.80665
B <- 1000L
SEED <- 20260723L
TIME_BLOCK_YEARS <- 5L
N_THREADS <- 8L
BATCH <- 1000L
METRICS <- c(
  "actual", "IWV_channel", "z_channel", "interaction_channel"
)

e <- new.env(parent = emptyenv())
load(CORE_FILE, envir = e)
lon <- as.numeric(e$lon)
lat <- as.numeric(e$lat)
dates <- as.Date(e$all_dates)
nt <- length(dates)
ncell <- length(lon) * length(lat)
Imat <- matrix(as.numeric(e$iwv_st), nrow = ncell, ncol = nt)
Zmat <- matrix(as.numeric(e$zbar_st), nrow = ncell, ncol = nt)
rm(e)
gc()

s <- readRDS(SINGLE_FILE)
p_aligned <- reorder_exact_grid(
  s$tp_native_m, s$lon, s$lat, lon, lat
)
Pmat <- matrix(as.numeric(p_aligned), nrow = ncell, ncol = nt)
rm(s, p_aligned)
gc()

stopifnot(
  length(dates) == 552L,
  ncell == 65160L,
  all(dim(Imat) == c(ncell, nt)),
  all(dim(Zmat) == c(ncell, nt)),
  all(dim(Pmat) == c(ncell, nt))
)

month <- as.integer(format(dates, "%m"))
year <- as.integer(format(dates, "%Y"))
years <- sort(unique(year))
ny <- length(years)
year_id <- match(year, years)
td <- as.numeric(dates)
tc <- td - mean(td)

cell <- expand.grid(lon = lon, lat = lat)
cell$area_weight <- pmax(0, cos(cell$lat * pi / 180))
raw_block <- paste(
  floor(((cell$lon + 180) %% 360) / 10),
  floor((cell$lat + 90) / 10),
  sep = "_"
)
block_levels <- unique(raw_block)
cell$block_index <- match(raw_block, block_levels)
nblock <- length(block_levels)

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
  if (integer) matrix(0L, ncell, ny) else matrix(0, ncell, ny)
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

for (st in seq.int(1L, ncell, by = BATCH)) {
  en <- min(ncell, st + BATCH - 1L)
  id <- st:en
  message("Preparing wet/dry annual sufficient statistics: ", st, "-", en)

  I <- Imat[id, , drop = FALSE]
  Z <- Zmat[id, , drop = FALSE]
  P <- Pmat[id, , drop = FALSE]
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

  wet_count[id, ] <- annual_count(wet_mask)
  dry_count[id, ] <- annual_count(dry_mask)
  wet_actual[id, ] <- annual_sum(actual, wet_mask)
  wet_iwv[id, ] <- annual_sum(iwv_channel, wet_mask)
  wet_z[id, ] <- annual_sum(z_channel, wet_mask)
  wet_interaction[id, ] <- annual_sum(interaction_channel, wet_mask)
  dry_actual[id, ] <- annual_sum(actual, dry_mask)
  dry_iwv[id, ] <- annual_sum(iwv_channel, dry_mask)
  dry_z[id, ] <- annual_sum(z_channel, dry_mask)
  dry_interaction[id, ] <- annual_sum(interaction_channel, dry_mask)

  rm(
    I, Z, P, ai, az, ap, covariance_climatology, product_anomaly,
    actual, iwv_channel, z_channel, interaction_channel,
    wet_mask, dry_mask
  )
  gc()
}
rm(Imat, Zmat, Pmat)
gc()

make_spatial_counts <- function(B, nblock, seed) {
  set.seed(seed)
  out <- matrix(0L, B, nblock)
  for (b in seq_len(B)) {
    out[b, ] <- tabulate(
      sample.int(nblock, nblock, replace = TRUE),
      nbins = nblock
    )
  }
  out
}

make_temporal_counts <- function(B, ny, block_years, seed) {
  set.seed(seed)
  starts <- seq_len(ny - block_years + 1L)
  n_draw <- ceiling(ny / block_years)
  out <- matrix(0L, B, ny)
  for (b in seq_len(B)) {
    ss <- sample(starts, n_draw, replace = TRUE)
    selected <- unlist(
      lapply(ss, function(z) z:(z + block_years - 1L)),
      use.names = FALSE
    )[seq_len(ny)]
    out[b, ] <- tabulate(selected, nbins = ny)
  }
  out
}

spatial_random <- make_spatial_counts(B, nblock, SEED + 1L)
temporal_random <- make_temporal_counts(
  B, ny, TIME_BLOCK_YEARS, SEED + 2L
)
spatial_ones <- matrix(1L, B, nblock)
temporal_ones <- matrix(1L, B, ny)

all_temporal_counts <- rbind(
  temporal_ones,
  temporal_random,
  temporal_random
)
all_spatial_counts <- rbind(
  spatial_random,
  spatial_ones,
  spatial_random
)

draw_matrix <- wetdry_spatiotemporal_bootstrap_cpp(
  wet_count,
  wet_actual,
  wet_iwv,
  wet_z,
  wet_interaction,
  dry_count,
  dry_actual,
  dry_iwv,
  dry_z,
  dry_interaction,
  as.integer(cell$block_index),
  as.numeric(cell$area_weight),
  all_temporal_counts,
  all_spatial_counts,
  N_THREADS
)

observed_matrix <- wetdry_spatiotemporal_bootstrap_cpp(
  wet_count,
  wet_actual,
  wet_iwv,
  wet_z,
  wet_interaction,
  dry_count,
  dry_actual,
  dry_iwv,
  dry_z,
  dry_interaction,
  as.integer(cell$block_index),
  as.numeric(cell$area_weight),
  matrix(1L, 1L, ny),
  matrix(1L, 1L, nblock),
  1L
)

method <- rep(
  c("spatial_only", "temporal_5yr_only", "two_stage_spatial_x_time"),
  each = B
)
iteration <- rep(seq_len(B), 3L)
draws <- data.frame(
  method = method,
  iteration = iteration,
  as.data.frame(draw_matrix, check.names = FALSE),
  stringsAsFactors = FALSE
)
write.csv(
  draws,
  file.path(OUT, "WET_DRY_SPATIOTEMPORAL_BOOTSTRAP_DRAWS.csv"),
  row.names = FALSE
)

locked <- read.csv(LOCKED_SUMMARY, check.names = FALSE)
locked <- locked[
  locked$region == "Global" & locked$metric %in% METRICS,
]

metric_map <- expand.grid(
  event = c("Wet", "Dry"),
  metric = METRICS,
  stringsAsFactors = FALSE
)
metric_map$column <- paste(metric_map$event, metric_map$metric, sep = "_")

q2 <- function(x) {
  unname(quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
}
summary_rows <- lapply(seq_len(nrow(metric_map)), function(i) {
  ev <- metric_map$event[i]
  metric <- metric_map$metric[i]
  column <- metric_map$column[i]
  obs <- as.numeric(observed_matrix[1, column])
  qs <- q2(draws[[column]][draws$method == "spatial_only"])
  qt <- q2(draws[[column]][draws$method == "temporal_5yr_only"])
  qst <- q2(
    draws[[column]][draws$method == "two_stage_spatial_x_time"]
  )
  lock_row <- locked[locked$event == ev & locked$metric == metric,]
  data.frame(
    event = ev,
    metric = metric,
    estimate = obs,
    locked_estimate = lock_row$estimate,
    difference_from_locked = obs - lock_row$estimate,
    locked_spatial_L95 = lock_row$lower95,
    locked_spatial_U95 = lock_row$upper95,
    recomputed_spatial_L95 = qs[1],
    recomputed_spatial_U95 = qs[2],
    temporal_5yr_L95 = qt[1],
    temporal_5yr_U95 = qt[2],
    two_stage_L95 = qst[1],
    two_stage_U95 = qst[2],
    n_grid = ncell,
    n_year = ny,
    spatial_blocks = nblock,
    bootstrap_replicates = B,
    time_block_years = TIME_BLOCK_YEARS,
    seed_base = SEED,
    stringsAsFactors = FALSE
  )
})
summary_out <- do.call(rbind, summary_rows)
write.csv(
  summary_out,
  file.path(OUT, "WET_DRY_SPATIOTEMPORAL_UNCERTAINTY_GLOBAL.csv"),
  row.names = FALSE
)

wet_closure <- with(
  draws,
  Wet_actual - Wet_IWV_channel - Wet_z_channel -
    Wet_interaction_channel
)
dry_closure <- with(
  draws,
  Dry_actual - Dry_IWV_channel - Dry_z_channel -
    Dry_interaction_channel
)
max_locked_difference <- max(
  abs(summary_out$difference_from_locked),
  na.rm = TRUE
)
max_closure <- max(
  abs(c(wet_closure, dry_closure)),
  na.rm = TRUE
)
if (max_locked_difference > 1e-6) {
  stop("Locked wet/dry global estimates were not reproduced.")
}
if (max_closure > 1e-7) {
  stop("Wet/dry channel closure failed in bootstrap draws.")
}
if (any(!is.finite(as.matrix(summary_out[, 8:13])))) {
  stop("Non-finite uncertainty interval detected.")
}

methods_cn <- c(
  "# 湿干月份时间—空间两阶段bootstrap：补充材料建议文字",
  "",
  paste0(
    "为量化有限湿干月份造成的时间抽样不确定度，在锁定的格点",
    "P90/P10事件定义、TOP100诊断量和精确通道分解不变的前提下，",
    "将1979–2024年划分为46个完整日历年。时间阶段采用5年移动块",
    "bootstrap：每次有放回抽取连续5年块并截取至46年；同一次",
    "时间抽样共同用于湿月、干月、W-GPE总异常及三个分解通道。"
  ),
  "",
  paste0(
    "空间阶段以10°×10°块为单位有放回抽样，并与时间块抽样独立",
    "组合。每个重复内先在每个格点按被抽取年份中的原事件月份",
    "重新计算事件均值，再进行余弦纬度面积加权空间汇总，从而",
    "保持原格点事件复合估计量。分别报告仅空间、仅时间以及",
    "时间—空间两阶段区间；每种方案均使用1,000次重复，随机种子",
    "基值为20260723。"
  ),
  "",
  paste0(
    "所有通道在每次重复中使用完全相同的时间块和空间块抽样，",
    "W-GPE总异常与IWV、<z>和交互通道在全部重复中保持数值闭合。"
  )
)
writeLines(
  methods_cn,
  file.path(OUT, "WET_DRY_SPATIOTEMPORAL_SI_TEXT_CN.md"),
  useBytes = TRUE
)

qa_lines <- c(
  "# Wet/dry spatiotemporal uncertainty QA",
  "",
  "## Status",
  "",
  "- Global computation completed: YES.",
  "- Existing canonical wet/dry results modified: NO.",
  paste0("- Grid cells: ", ncell, "."),
  paste0("- Calendar years: ", ny, " (", min(years), "-", max(years), ")."),
  paste0("- Spatial blocks: ", nblock, "."),
  paste0("- Replicates per method: ", B, "."),
  paste0("- Time block: ", TIME_BLOCK_YEARS, " calendar years."),
  paste0("- Random seed base: ", SEED, "."),
  "",
  "## Reproduction and closure",
  "",
  paste0(
    "- Maximum absolute difference from locked global estimates: ",
    format(max_locked_difference, scientific = TRUE, digits = 8), "."
  ),
  paste0(
    "- Maximum absolute channel-closure residual across all bootstrap draws: ",
    format(max_closure, scientific = TRUE, digits = 8), "."
  ),
  "- Locked event masks and percentile definitions were preserved.",
  "- Wet, dry and all exact channels used paired temporal and spatial draws.",
  "",
  "## Interpretation",
  "",
  paste0(
    "- Spatial-only intervals reproduce the original uncertainty target; ",
    "temporal-only intervals quantify finite event-year sampling; ",
    "two-stage intervals combine both sources."
  )
)
qa_path <- file.path(OUT, "WET_DRY_SPATIOTEMPORAL_QA.md")
writeLines(qa_lines, qa_path, useBytes = TRUE)

output_files <- list.files(OUT, full.names = TRUE)
manifest_path <- file.path(OUT, "WET_DRY_SPATIOTEMPORAL_SHA256.csv")
output_files <- setdiff(output_files, manifest_path)
manifest_inputs <- c(
  CORE_FILE, SINGLE_FILE, COORD_FILE, LOCKED_SUMMARY,
  SCRIPT, CPP, MAKEVARS_FILE
)
manifest <- data.frame(
  role = c(
    rep("input_or_script", length(manifest_inputs)),
    rep("output", length(output_files))
  ),
  path = c(manifest_inputs, output_files),
  size_bytes = file.info(c(manifest_inputs, output_files))$size,
  sha256 = vapply(
    c(manifest_inputs, output_files),
    digest::digest,
    character(1),
    algo = "sha256",
    file = TRUE
  ),
  stringsAsFactors = FALSE
)
write.csv(manifest, manifest_path, row.names = FALSE)

message("WET_DRY_SPATIOTEMPORAL_ANALYSIS=COMPLETE")
message("OUTPUT=", OUT)
