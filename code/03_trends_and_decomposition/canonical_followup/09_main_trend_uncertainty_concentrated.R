options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("__WGPE_R_LIBRARY__", .libPaths()))

ROOT <- "__WGPE_PROJECT_ROOT__"
CANON <- file.path(ROOT, "CANONICAL_WGPE_TOP100_v2")
SCRIPT <- file.path(
  CANON, "07_SCRIPTS", "09_main_trend_uncertainty_concentrated.R"
)
OUT <- file.path(
  CANON, "08_REPRODUCTION_OUTPUT", "main_trend_uncertainty"
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

CORE_FILE <- file.path(
  ROOT, "output_TOP100_revision", "derived_data", "ERA5_TOP100",
  "global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData"
)
LOCKED_TREND_FILE <- file.path(
  ROOT, "output_TOP100_revision", "results", "core",
  "Trend_results_FINAL_v3.csv"
)

if (!all(file.exists(c(SCRIPT, CORE_FILE, LOCKED_TREND_FILE)))) {
  stop("Required main-trend input or script is missing.")
}
if (!requireNamespace("sandwich", quietly = TRUE)) {
  stop("sandwich is required for Newey-West covariance.")
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required for SHA-256.")
}

B <- 2000L
SEED <- 20260724L
TIME_BLOCK_MONTHS <- 60L
START_YEARS <- 1979:1989

decimal_year <- function(d) {
  y <- as.integer(format(d, "%Y"))
  y0 <- as.Date(sprintf("%04d-01-01", y))
  y1 <- as.Date(sprintf("%04d-01-01", y + 1L))
  y + as.numeric(d - y0) / as.numeric(y1 - y0)
}

deseason_vector <- function(y, month) {
  y - ave(y, month, FUN = function(z) mean(z, na.rm = TRUE))
}

deseason_rows <- function(y, month) {
  out <- y
  for (mm in 1:12) {
    jj <- which(month == mm)
    out[, jj] <- sweep(
      y[, jj, drop = FALSE],
      1,
      rowMeans(y[, jj, drop = FALSE], na.rm = TRUE),
      "-"
    )
  }
  out
}

slope_fast <- function(y, x) {
  xc <- x - mean(x)
  sum(xc * (y - mean(y))) / sum(xc^2)
}

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

make_time_indices <- function(B, n, block_months, seed) {
  set.seed(seed)
  starts <- seq_len(n - block_months + 1L)
  n_draw <- ceiling(n / block_months)
  out <- matrix(0L, B, n)
  for (b in seq_len(B)) {
    ss <- sample(starts, n_draw, replace = TRUE)
    out[b, ] <- unlist(
      lapply(ss, function(z) z:(z + block_months - 1L)),
      use.names = FALSE
    )[seq_len(n)]
  }
  out
}

residual_block_slopes <- function(y, x, time_indices) {
  b0 <- slope_fast(y, x)
  a0 <- mean(y) - b0 * mean(x)
  fitted <- a0 + b0 * x
  residual <- y - fitted
  residual <- residual - mean(residual)
  xc <- x - mean(x)
  denom <- sum(xc^2)
  apply(
    time_indices,
    1,
    function(ii) {
      yy <- fitted + residual[ii]
      sum(xc * (yy - mean(yy))) / denom
    }
  )
}

e <- new.env(parent = emptyenv())
load(CORE_FILE, envir = e)
lon <- as.numeric(e$lon)
lat <- as.numeric(e$lat)
dates <- as.Date(e$all_dates)
arrays <- list(
  iwv = e$iwv_st,
  zbar = e$zbar_st,
  wgpe = e$wgpe_st
)
rm(e)
gc()

nt <- length(dates)
ncell <- length(lon) * length(lat)
stopifnot(nt == 552L, ncell == 65160L)
month <- as.integer(format(dates, "%m"))
year <- as.integer(format(dates, "%Y"))
x <- decimal_year(dates)
xc <- x - mean(x)
denom_x <- sum(xc^2)

cell <- expand.grid(lon = lon, lat = lat)
cell$area_weight <- pmax(0, cos(cell$lat * pi / 180))
raw_block <- paste(
  floor(((cell$lon + 180) %% 360) / 10),
  floor((cell$lat + 90) / 10),
  sep = "_"
)
block_factor <- factor(raw_block, levels = unique(raw_block))
nblock <- nlevels(block_factor)

spatial_counts <- make_spatial_counts(B, nblock, SEED + 1L)
time_indices <- make_time_indices(
  B, nt, TIME_BLOCK_MONTHS, SEED + 2L
)

locked <- read.csv(LOCKED_TREND_FILE, check.names = FALSE)
locked <- locked[
  locked$height_reference == "absolute" & locked$region == "Global",
]
locked_value <- c(
  iwv = locked$iwv_trend_kg_m2_yr,
  zbar = locked$z_trend_m_yr,
  wgpe = locked$wgpe_trend_J_m2_yr
)
unit <- c(
  iwv = "kg m-2 yr-1",
  zbar = "m yr-1",
  wgpe = "J m-2 yr-1"
)

summary_rows <- list()
start_rows <- list()
draw_rows <- list()

for (variable in names(arrays)) {
  message("Main trend uncertainty: ", variable)
  m <- matrix(
    as.numeric(arrays[[variable]]),
    nrow = ncell,
    ncol = nt
  )
  ok <- is.finite(m)
  m0 <- m
  m0[!ok] <- 0

  numerator <- rowsum(
    m0 * cell$area_weight,
    block_factor,
    reorder = FALSE
  )
  denominator <- rowsum(
    ok * cell$area_weight,
    block_factor,
    reorder = FALSE
  )

  global_raw <- colSums(numerator) / colSums(denominator)
  global_anomaly <- deseason_vector(global_raw, month)
  estimate <- slope_fast(global_anomaly, x)

  fit <- lm(global_anomaly ~ x)
  ols_ci <- unname(confint(fit, "x", level = 0.95))
  nw_lag <- floor(4 * (nt / 100)^(2 / 9))
  nw_v <- sandwich::NeweyWest(
    fit,
    lag = nw_lag,
    prewhite = FALSE,
    adjust = TRUE
  )
  nw_se <- sqrt(nw_v["x", "x"])
  nw_ci <- estimate + c(-1, 1) * qt(0.975, df.residual(fit)) * nw_se

  time_slopes <- residual_block_slopes(
    global_anomaly,
    x,
    time_indices
  )

  spatial_num <- spatial_counts %*% numerator
  spatial_den <- spatial_counts %*% denominator
  spatial_series <- spatial_num / spatial_den
  spatial_anomaly <- deseason_rows(spatial_series, month)
  spatial_slopes <- as.numeric(spatial_anomaly %*% xc / denom_x)

  combined_slopes <- numeric(B)
  for (b in seq_len(B)) {
    yb <- spatial_anomaly[b, ]
    bb <- slope_fast(yb, x)
    aa <- mean(yb) - bb * mean(x)
    fitted <- aa + bb * x
    residual <- yb - fitted
    residual <- residual - mean(residual)
    yy <- fitted + residual[time_indices[b, ]]
    combined_slopes[b] <- slope_fast(yy, x)
  }

  q_time <- unname(quantile(time_slopes, c(0.025, 0.975)))
  q_spatial <- unname(quantile(spatial_slopes, c(0.025, 0.975)))
  q_combined <- unname(quantile(combined_slopes, c(0.025, 0.975)))

  summary_rows[[variable]] <- data.frame(
    variable = variable,
    unit = unit[[variable]],
    estimate = estimate,
    locked_estimate = locked_value[[variable]],
    difference_from_locked = estimate - locked_value[[variable]],
    relative_difference_from_locked = (
      estimate - locked_value[[variable]]
    ) / locked_value[[variable]],
    OLS_L95 = ols_ci[1],
    OLS_U95 = ols_ci[2],
    Newey_West_L95 = nw_ci[1],
    Newey_West_U95 = nw_ci[2],
    time_block_5yr_L95 = q_time[1],
    time_block_5yr_U95 = q_time[2],
    spatial_10deg_L95 = q_spatial[1],
    spatial_10deg_U95 = q_spatial[2],
    spatiotemporal_L95 = q_combined[1],
    spatiotemporal_U95 = q_combined[2],
    all_robust_CI_exclude_zero = all(
      c(nw_ci[1], q_time[1], q_spatial[1], q_combined[1]) > 0
    ),
    n_month = nt,
    n_grid = ncell,
    n_spatial_blocks = nblock,
    bootstrap_replicates = B,
    time_block_months = TIME_BLOCK_MONTHS,
    Newey_West_lag = nw_lag,
    seed_base = SEED,
    stringsAsFactors = FALSE
  )

  for (start_year in START_YEARS) {
    jj <- which(year >= start_year)
    yy <- global_raw[jj]
    mm <- month[jj]
    xx <- x[jj]
    an <- deseason_vector(yy, mm)
    start_rows[[length(start_rows) + 1L]] <- data.frame(
      variable = variable,
      unit = unit[[variable]],
      start_year = start_year,
      end_year = max(year),
      n_month = length(jj),
      slope = slope_fast(an, xx),
      stringsAsFactors = FALSE
    )
  }

  draw_rows[[variable]] <- rbind(
    data.frame(
      variable = variable,
      method = "time_block_5yr",
      iteration = seq_len(B),
      slope = time_slopes
    ),
    data.frame(
      variable = variable,
      method = "spatial_10deg",
      iteration = seq_len(B),
      slope = spatial_slopes
    ),
    data.frame(
      variable = variable,
      method = "spatiotemporal",
      iteration = seq_len(B),
      slope = combined_slopes
    )
  )

  rm(
    m, ok, m0, numerator, denominator, spatial_num, spatial_den,
    spatial_series, spatial_anomaly
  )
  gc()
}

summary_out <- do.call(rbind, summary_rows)
start_out <- do.call(rbind, start_rows)
draws_out <- do.call(rbind, draw_rows)

start_summary <- do.call(
  rbind,
  lapply(split(start_out, start_out$variable), function(d) {
    full <- d$slope[d$start_year == 1979]
    data.frame(
      variable = d$variable[1],
      unit = d$unit[1],
      start_year_min = min(d$start_year),
      start_year_max = max(d$start_year),
      minimum_slope = min(d$slope),
      maximum_slope = max(d$slope),
      full_period_slope = full,
      maximum_absolute_change_from_full = max(abs(d$slope - full)),
      all_positive = all(d$slope > 0),
      stringsAsFactors = FALSE
    )
  })
)

max_locked_difference <- max(
  abs(summary_out$difference_from_locked),
  na.rm = TRUE
)
max_relative_locked_difference <- max(
  abs(summary_out$relative_difference_from_locked),
  na.rm = TRUE
)
if (max_relative_locked_difference > 1e-6) {
  stop("Main global trend estimates do not reproduce the locked values.")
}
if (!all(summary_out$all_robust_CI_exclude_zero)) {
  stop("At least one principal robust uncertainty interval crosses zero.")
}
if (!all(start_summary$all_positive)) {
  stop("At least one start-year sensitivity slope changes sign.")
}

write.csv(
  summary_out,
  file.path(OUT, "MAIN_TREND_UNCERTAINTY_CONCENTRATED.csv"),
  row.names = FALSE
)
write.csv(
  start_out,
  file.path(OUT, "MAIN_TREND_START_YEAR_SENSITIVITY.csv"),
  row.names = FALSE
)
write.csv(
  start_summary,
  file.path(OUT, "MAIN_TREND_START_YEAR_SUMMARY.csv"),
  row.names = FALSE
)
write.csv(
  draws_out,
  file.path(OUT, "MAIN_TREND_BOOTSTRAP_DRAWS.csv"),
  row.names = FALSE
)

si_text <- c(
  "# 主趋势不确定度与起始年份敏感性：补充材料建议文字",
  "",
  paste0(
    "对1979–2024年全球TOP100 IWV、<z>和W-GPE趋势，集中报告",
    "普通最小二乘区间、Newey–West区间、5年残差移动块bootstrap、",
    "10°×10°空间块bootstrap和空间—时间联合bootstrap。时间块、",
    "空间块和联合方案均使用2,000次重复；联合方案在每次重复中",
    "先重采样空间块形成全球月序列，再对该序列的趋势残差实施",
    "60个月移动块重采样。"
  ),
  "",
  paste0(
    "起始年份敏感性保持2024年终点不变，将起始年份从1979年逐年",
    "移动至1989年，并在每个时间窗口内重新估计逐月气候态和线性",
    "趋势。该分析检验主趋势是否由最初若干年份决定。"
  ),
  "",
  paste0(
    "三项主趋势在Newey–West、时间块、空间块和空间—时间联合",
    "95%区间下均保持正值；1979–1989各起始年份窗口的趋势符号",
    "同样保持为正。"
  )
)
writeLines(
  si_text,
  file.path(OUT, "MAIN_TREND_UNCERTAINTY_SI_TEXT_CN.md"),
  useBytes = TRUE
)

qa_lines <- c(
  "# Main trend uncertainty QA",
  "",
  "## Status",
  "",
  "- Concentrated global trend uncertainty table completed: YES.",
  "- Existing canonical trends, registry and figures modified: NO.",
  paste0("- Variables: ", paste(names(arrays), collapse = ", "), "."),
  paste0("- Months: ", nt, "; grid cells: ", ncell, "."),
  paste0("- Spatial blocks: ", nblock, "."),
  paste0("- Bootstrap replicates: ", B, "."),
  paste0("- Time-block length: ", TIME_BLOCK_MONTHS, " months."),
  paste0("- Start-year range: ", min(START_YEARS), "-", max(START_YEARS), "."),
  paste0("- Random seed base: ", SEED, "."),
  "",
  "## Reproduction and robustness",
  "",
  paste0(
    "- Maximum absolute difference from locked global trend estimates: ",
    format(max_locked_difference, scientific = TRUE, digits = 8), "."
  ),
  paste0(
    "- Maximum relative difference from locked global trend estimates: ",
    format(max_relative_locked_difference, scientific = TRUE, digits = 8),
    "."
  ),
  paste0(
    "- All Newey-West, time-block, spatial-block and combined intervals ",
    "exclude zero: ", all(summary_out$all_robust_CI_exclude_zero), "."
  ),
  paste0(
    "- All 1979-1989 start-year sensitivity slopes remain positive: ",
    all(start_summary$all_positive), "."
  )
)
qa_path <- file.path(OUT, "MAIN_TREND_UNCERTAINTY_QA.md")
writeLines(qa_lines, qa_path, useBytes = TRUE)

manifest_path <- file.path(OUT, "MAIN_TREND_UNCERTAINTY_SHA256.csv")
output_files <- setdiff(list.files(OUT, full.names = TRUE), manifest_path)
all_files <- c(CORE_FILE, LOCKED_TREND_FILE, SCRIPT, output_files)
manifest <- data.frame(
  role = c(
    "input", "input", "script",
    rep("output", length(output_files))
  ),
  path = all_files,
  size_bytes = file.info(all_files)$size,
  sha256 = vapply(
    all_files,
    digest::digest,
    character(1),
    algo = "sha256",
    file = TRUE
  ),
  stringsAsFactors = FALSE
)
write.csv(manifest, manifest_path, row.names = FALSE)

message("MAIN_TREND_UNCERTAINTY=COMPLETE")
message("OUTPUT=", OUT)
