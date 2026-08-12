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
  CANON, "08_REPRODUCTION_OUTPUT",
  "wetdry_spatiotemporal_uncertainty"
)
CHUNK_DIR <- file.path(OUT, "annual_chunks")
SCRIPT <- file.path(
  CANON, "07_SCRIPTS",
  "08b_finalize_wetdry_spatiotemporal_bootstrap.R"
)
PREP_SCRIPT <- file.path(
  CANON, "07_SCRIPTS", "08a_prepare_wetdry_annual_chunk.R"
)
CPP <- file.path(
  CANON, "07_SCRIPTS", "support",
  "wetdry_spatiotemporal_bootstrap.cpp"
)
MAKEVARS_FILE <- file.path(CANON, "07_SCRIPTS", "Makevars_Rtools45")
Sys.setenv(R_MAKEVARS_USER = MAKEVARS_FILE)

LOCKED_SUMMARY <- file.path(
  ROOT, "output_TOP100_revision", "results", "precip_wetdry",
  "exact_TOP100", "WET_DRY_EXACT_SUMMARY_TOP100.csv"
)
chunk_files <- file.path(
  CHUNK_DIR, sprintf("WETDRY_ANNUAL_CHUNK_%02d.rds", 1:7)
)
required <- c(
  SCRIPT, PREP_SCRIPT, CPP, MAKEVARS_FILE, LOCKED_SUMMARY,
  chunk_files
)
if (!all(file.exists(required))) {
  stop(
    "Missing required file(s): ",
    paste(required[!file.exists(required)], collapse = "; ")
  )
}
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required.")
Rcpp::sourceCpp(CPP, rebuild = TRUE, showOutput = FALSE)

B <- 1000L
SEED <- 20260723L
TIME_BLOCK_YEARS <- 5L
N_THREADS <- 8L
METRICS <- c(
  "actual", "IWV_channel", "z_channel", "interaction_channel"
)

chunks <- lapply(chunk_files, readRDS)
cell <- do.call(rbind, lapply(chunks, `[[`, "cell"))
years <- chunks[[1]]$years
dates <- chunks[[1]]$dates
matrix_names <- c(
  "wet_count", "wet_actual", "wet_iwv", "wet_z", "wet_interaction",
  "dry_count", "dry_actual", "dry_iwv", "dry_z", "dry_interaction"
)
annual <- lapply(
  matrix_names,
  function(nm) do.call(rbind, lapply(chunks, `[[`, nm))
)
names(annual) <- matrix_names
rm(chunks)
gc()

stopifnot(
  nrow(cell) == 65160L,
  identical(cell$cell_id, seq_len(nrow(cell))),
  length(years) == 46L,
  all(vapply(annual, nrow, integer(1)) == nrow(cell)),
  all(vapply(annual, ncol, integer(1)) == length(years))
)

block_levels <- unique(cell$raw_block)
cell$block_index <- match(cell$raw_block, block_levels)
nblock <- length(block_levels)
ny <- length(years)

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

run_bootstrap <- function(temporal_counts, spatial_counts) {
  wetdry_spatiotemporal_bootstrap_cpp(
    annual$wet_count,
    annual$wet_actual,
    annual$wet_iwv,
    annual$wet_z,
    annual$wet_interaction,
    annual$dry_count,
    annual$dry_actual,
    annual$dry_iwv,
    annual$dry_z,
    annual$dry_interaction,
    as.integer(cell$block_index),
    as.numeric(cell$area_weight),
    temporal_counts,
    spatial_counts,
    N_THREADS
  )
}

message("Running spatial-only bootstrap.")
draw_spatial <- run_bootstrap(
  matrix(1L, B, ny),
  spatial_random
)
message("Running temporal-only bootstrap.")
draw_temporal <- run_bootstrap(
  temporal_random,
  matrix(1L, B, nblock)
)
message("Running paired two-stage bootstrap.")
draw_two_stage <- run_bootstrap(
  temporal_random,
  spatial_random
)
observed <- run_bootstrap(
  matrix(1L, 1L, ny),
  matrix(1L, 1L, nblock)
)

draws <- rbind(
  data.frame(
    method = "spatial_only",
    iteration = seq_len(B),
    as.data.frame(draw_spatial, check.names = FALSE)
  ),
  data.frame(
    method = "temporal_5yr_only",
    iteration = seq_len(B),
    as.data.frame(draw_temporal, check.names = FALSE)
  ),
  data.frame(
    method = "two_stage_spatial_x_time",
    iteration = seq_len(B),
    as.data.frame(draw_two_stage, check.names = FALSE)
  )
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
  obs <- as.numeric(observed[1, column])
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
    n_grid = nrow(cell),
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
  paste0("- Grid cells: ", nrow(cell), "."),
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
    "- Maximum absolute channel-closure residual across all draws: ",
    format(max_closure, scientific = TRUE, digits = 8), "."
  ),
  "- Locked event masks and percentile definitions were preserved.",
  "- Wet, dry and all exact channels used paired temporal and spatial draws."
)
qa_path <- file.path(OUT, "WET_DRY_SPATIOTEMPORAL_QA.md")
writeLines(qa_lines, qa_path, useBytes = TRUE)

manifest_path <- file.path(OUT, "WET_DRY_SPATIOTEMPORAL_SHA256.csv")
output_files <- list.files(
  OUT,
  full.names = TRUE,
  recursive = TRUE
)
output_files <- setdiff(
  output_files,
  c(
    manifest_path,
    file.path(OUT, "run_stdout.log"),
    file.path(OUT, "run_stderr.log")
  )
)
inputs_and_scripts <- c(
  LOCKED_SUMMARY, SCRIPT, PREP_SCRIPT, CPP, MAKEVARS_FILE
)
all_files <- c(inputs_and_scripts, output_files)
manifest <- data.frame(
  role = c(
    rep("input_or_script", length(inputs_and_scripts)),
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

message("WET_DRY_SPATIOTEMPORAL_ANALYSIS=COMPLETE")
message("OUTPUT=", OUT)
