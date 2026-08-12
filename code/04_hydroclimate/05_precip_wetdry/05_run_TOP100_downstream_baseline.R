# Purpose: Recompute the audited hydroclimate/ecology downstream baseline using
#          the formal TOP100 absolute-height core while reusing only immutable
#          preprocessed non-WGPE controls from the earlier audited workflow.
# Inputs:  formal TOP100 core RData; ERA5 single-level monthly arrays; GPP/LAI;
#          processed controls; land-cover summaries; audited legacy runner source.
# Outputs: output_TOP100_revision/temp/TOP100_DOWNSTREAM_BASELINE/*
# Parameters: 1979-2024 hydroclimate; product-specific ecology periods; seed 20260710.
# Dependencies: R >= 4.5, ncdf4, ggplot2, cowplot, scales.
# Overwrite policy: writes only inside the TOP100 revision tree; legacy outputs are untouched.
# Important: the legacy event-year ENSO branch emitted by this baseline is NOT formal.
#            Formal ENSO values are replaced by the ONI-standard TOP100 module.

options(stringsAsFactors = FALSE, warn = 1)
ROOT <- "__WGPE_PROJECT_ROOT__"
SRC <- file.path(ROOT, "output_attribution_minimal", "FULL_RERUN_FINAL_absolute_relative_v3",
                 "13_scripts", "run_FULL_RERUN_FINAL_v3.R")
CORE <- file.path(ROOT, "output_TOP100_revision", "derived_data", "ERA5_TOP100",
                  "global_wgpe_1deg_surface_truncated_100hPa_TOP100.RData")
OUT <- Sys.getenv("TOP100_CLOSURE_BASELINE_OUT", unset = file.path(ROOT, "output_TOP100_revision", "temp", "TOP100_DOWNSTREAM_BASELINE"))
CONTROLS <- file.path(ROOT, "output_attribution_minimal", "FULL_RERUN_FINAL_absolute_relative_v3",
                      "01_inputs_and_processing", "processed_controls_FINAL_v3.rds")
PATCHED <- file.path(OUT, "run_FULL_RERUN_FINAL_v3_PATCHED_TOP100.R")

stopifnot(file.exists(SRC), file.exists(CORE), file.exists(CONTROLS))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
x <- readLines(SRC, warn = FALSE, encoding = "UTF-8")

replace_one <- function(pattern, replacement, fixed = TRUE) {
  hit <- grep(pattern, x, fixed = fixed)
  if (length(hit) != 1L) stop("Expected one match for: ", pattern, "; found ", length(hit))
  x[hit] <<- replacement
}

replace_one('OUT <- file.path(ROOT, "output_attribution_minimal", "FULL_RERUN_FINAL_absolute_relative_v3")',
            sprintf('OUT <- "%s"', gsub("\\\\", "/", OUT)))
replace_one('core_file <- file.path(ROOT, "output_attribution_minimal", "global_wgpe_1deg_surface_truncated_300hPa_FIXED_FULL.RData")',
            sprintf('core_file <- "%s"', gsub("\\\\", "/", CORE)))
replace_one('process_file <- file.path(OUT, "01_inputs_and_processing", "processed_controls_FINAL_v3.rds")',
            sprintf('process_file <- "%s"', gsub("\\\\", "/", CONTROLS)))

# Reconcile the legacy -179--180 external grid to the TOP100 0--359 grid in
# the generated execution copy.  The original v3 runner paired raw array
# positions, which shifts cross-dataset fields by about 180 degrees.
insert_after <- function(pattern, lines) {
  hit <- grep(pattern, x, fixed = TRUE)
  if (length(hit) != 1L) stop("Expected one insertion point for: ", pattern)
  x <<- append(x, lines, after = hit)
}
insert_after('single <- readRDS(single_file)', c(
  'single_lon <- as.numeric(single$lon); single_lat <- as.numeric(single$lat)',
  'single_mapping <- exact_grid_mapping(single_lon, single_lat, lon, lat)',
  'single_lon_idx <- single_mapping$lon_index; single_lat_idx <- single_mapping$lat_index'
))
replace_one('tp_arr <- single$tp_native_m', 'tp_arr <- single$tp_native_m[single_lon_idx, single_lat_idx, , drop = FALSE]')
replace_one('t2m_arr <- single$t2m_K', 't2m_arr <- single$t2m_K[single_lon_idx, single_lat_idx, , drop = FALSE]')
replace_one('ssrd_arr <- single$ssrd_native_J_m2', 'ssrd_arr <- single$ssrd_native_J_m2[single_lon_idx, single_lat_idx, , drop = FALSE]')
insert_after('controls <- readRDS(process_file)', c(
  'controls$surface_height <- reorder_exact_grid(controls$surface_height, single_lon, single_lat, lon, lat)',
  'controls$vpd_arr <- reorder_exact_grid(controls$vpd_arr, single_lon, single_lat, lon, lat)',
  'controls$rzsm_arr <- reorder_exact_grid(controls$rzsm_arr, single_lon, single_lat, lon, lat)'
))
replace_one('  ii <- match(nearest_to(round(panel$lon_val), lon), lon)',
            '  ii <- vapply(panel$lon_val, function(v) which.min(abs(((lon - v + 180) %% 360) - 180)), integer(1))')
replace_one('  jj <- match(nearest_to(round(panel$lat_val), lat), lat)',
            '  jj <- vapply(panel$lat_val, function(v) which.min(abs(lat - v)), integer(1))')
insert_after('write.csv(tp_rel$summary, file.path(OUT, "03_surface_relative_results", "Precip_temporal_precedence_FINAL_v3.csv"), row.names = FALSE)', c(
  'if (identical(Sys.getenv("TOP100_HYDRO_ONLY", unset="0"), "1")) {',
  '  writeLines(c("coordinate_harmonization=unified_module", "grid_id=TOP100_CANONICAL_GRID", "scope=hydroclimate_and_surface_relative"), file.path(OUT,"HYDRO_COORDINATE_SCOPE.txt"))',
  '  message("Hydro-only coordinate-fixed run complete")',
  '  quit(save="no",status=0)',
  '}'
))

# The old script's diagnostic figures are intentionally suppressed. Formal figures
# are rebuilt only after all TOP100 result branches have been replaced and audited.
i0 <- grep('message("Generating figures")', x, fixed = TRUE)
i1 <- grep('message("Writing source tables")', x, fixed = TRUE)
if (length(i0) != 1L || length(i1) != 1L || i1 <= i0) stop("Could not isolate old figure block")
x <- c(x[seq_len(i0 - 1L)],
       'message("Skipping legacy figure block; formal TOP100 plots are built later")',
       x[i1:length(x)])

header <- c(
  '# AUTO-PATCHED FORMAL TOP100 EXECUTION COPY',
  '# Generated by 05_run_TOP100_downstream_baseline.R.',
  '# The fixed-year ENSO rows are exploratory and must be replaced by ONI-standard outputs.',
  sprintf('source("%s")', gsub("\\\\", "/", file.path(ROOT, "scripts_TOP100_revision", "18_coordinate_harmonization", "coordinate_harmonization.R"))),
  ''
)
writeLines(c(header, x), PATCHED, useBytes = TRUE)
status <- system2("__WGPE_RSCRIPT__", PATCHED)
if (!identical(status, 0L)) stop("TOP100 downstream baseline failed with exit status ", status)

writeLines(c(
  "TOP100 downstream baseline completed.",
  paste("Formal core:", CORE),
  paste("Output:", OUT),
  "Fixed-year ENSO rows are not formal; use the ONI-standard TOP100 replacement module."
), file.path(OUT, "TOP100_BASELINE_SCOPE.txt"))
