# Purpose: one-command orchestration of the formal TOP100 revision.
# Inputs/parameters: 00_config/config_TOP100_revision.json and registered project inputs.
# Outputs: output_TOP100_revision and its main/SI figures, registry, workbooks and QA.
# Dependencies: R 4.5+, bundled Python, Node artifact-tool runtime, ncdf4 and plotting packages.
# Overwrite policy: refreshes only files under output_TOP100_revision; legacy project outputs are untouched.

options(stringsAsFactors = FALSE)
ROOT <- normalizePath("__WGPE_PROJECT_ROOT__", winslash = "/", mustWork = TRUE)
S <- file.path(ROOT, "scripts_TOP100_revision")
O <- file.path(ROOT, "output_TOP100_revision")
dir.create(file.path(O, "logs"), recursive = TRUE, showWarnings = FALSE)
logfile <- file.path(O, "logs", "run_all_TOP100_revision.log")
PY <- "__WGPE_PYTHON__"
NODE <- "__WGPE_NODE__"
RS <- file.path(R.home("bin"), "Rscript.exe")

run <- function(label, exe, script) {
  msg <- sprintf("[%s] START %s :: %s", format(Sys.time()), label, script)
  cat(msg, "\n"); cat(msg, "\n", file = logfile, append = TRUE)
  status <- system2(exe, normalizePath(file.path(S, script), winslash = "/", mustWork = TRUE),
                    stdout = logfile, stderr = logfile)
  end <- sprintf("[%s] END %s status=%s", format(Sys.time()), label, status)
  cat(end, "\n"); cat(end, "\n", file = logfile, append = TRUE)
  if (!identical(status, 0L)) stop("Stage failed: ", label, call. = FALSE)
}

run_env <- function(label, exe, script, env) {
  msg <- sprintf("[%s] START %s :: %s", format(Sys.time()), label, script)
  cat(msg, "\n"); cat(msg, "\n", file = logfile, append = TRUE)
  status <- system2(exe, normalizePath(file.path(S, script), winslash = "/", mustWork = TRUE),
                    stdout = logfile, stderr = logfile, env = env)
  end <- sprintf("[%s] END %s status=%s", format(Sys.time()), label, status)
  cat(end, "\n"); cat(end, "\n", file = logfile, append = TRUE)
  if (!identical(status, 0L)) stop("Stage failed: ", label, call. = FALSE)
}

stages <- list(
  c("input audit", PY, "01_input_audit/01_freeze_and_audit_TOP100.py"),
  c("TOP100/200/300 monthly core", PY, "02_TOP100_core_recalculation/04_compute_ERA5_top_boundary_monthly.py"),
  c("trend and top sensitivity", PY, "03_trend_and_layers/05_analyze_top_boundary_results.py"),
  c("six-layer contribution", PY, "03_trend_and_layers/09_compute_top100_layer_contributions.py"),
  c("top-boundary domain reconciliation", PY, "03_trend_and_layers/15_reconcile_top_boundary_domains.py"),
  c("export TOP100 R core", RS, "02_TOP100_core_recalculation/02_export_TOP100_core_RData.R"),
  c("downstream absolute/relative baseline", RS, "05_precip_wetdry/05_run_TOP100_downstream_baseline.R"),
  c("Figure 2a deseasonalized-only association", RS, "05_precip_wetdry/05c_compute_Fig2a_deseasoned_only_TOP100.R"),
  c("exact wet/dry", RS, "05_precip_wetdry/05b_exact_wetdry_decomposition_TOP100.R"),
  c("ONI-standard ENSO", RS, "06_ENSO/06_run_ONIstandard_TOP100.R"),
  c("conditional temporal precedence", RS, "08_temporal_precedence/08_conditional_temporal_precedence_lag1_4_TOP100.R"),
  c("precip temporal effect sizes", RS, "08_temporal_precedence/08a_precip_temporal_effects_TOP100.R"),
  c("ecology temporal effect sizes", RS, "08_temporal_precedence/08b_ecology_temporal_effects_TOP100.R"),
  c("ecology temporal effect resummary", RS, "08_temporal_precedence/08c_resummarize_ecology_temporal_effects_TOP100.R"),
  c("JRA55 native TOP100", PY, "09_JRA55/09_compute_JRA55_TOP100_native_1982_2023.py"),
  c("ERA5 common-period export", RS, "09_JRA55/09b_export_ERA5_TOP100_1982_2023.R"),
  c("JRA55 ERA5 comparison", PY, "09_JRA55/09c_compare_JRA55_ERA5_TOP100.py"),
  c("ERA5 common-period layers", PY, "09_JRA55/09d_compute_ERA5_TOP100_layers_1982_2023.py"),
  c("JRA55 six-layer comparison", PY, "09_JRA55/09e_compare_six_layer_contributions_TOP100.py"),
  c("surface endpoint audit", RS, "10_surface_truncation/16_audit_surface_endpoint_inputs.R"),
  c("surface endpoint methods", RS, "10_surface_truncation/17_compute_surface_endpoint_methods.R"),
  c("surface endpoint summaries", RS, "10_surface_truncation/18_analyze_surface_endpoint_methods.R"),
  c("top-boundary spatial blocks", PY, "11_uncertainty/06_build_top_boundary_10deg_blocks.py"),
  c("top-boundary zone monthly", PY, "11_uncertainty/07_build_top_boundary_zone_monthly.py"),
  c("top-boundary uncertainty", RS, "11_uncertainty/08_top_boundary_uncertainty.R"),
  c("latitude-profile uncertainty", RS, "11_uncertainty/11_build_lat5_monthly_and_uncertainty.R"),
  c("trend effective n and FDR", RS, "12_field_significance/10_top_boundary_effective_n_FDR.R"),
  c("precip monthly NetCDF export", RS, "12_field_significance/12a_export_TOP100_precip_monthly_nc.R"),
  c("precip field significance", PY, "12_field_significance/12b_precip_BH_BY_field_significance_TOP100.py"),
  c("ecology field significance", PY, "12_field_significance/12c_ecology_BH_BY_TOP100.py"),
  c("hydroclimate incremental value", RS, "13_incremental_value/13a_precip_wetdry_increment_TOP100.R"),
  c("ecology incremental value", RS, "13_incremental_value/13b_ecology_increment_TOP100.R"),
  c("ENSO incremental value", RS, "13_incremental_value/13c_ENSO_increment_TOP100.R"),
  c("increment robustness classes", PY, "13_incremental_value/13d_classify_incremental_value_TOP100.py"),
  c("Figure 1d coordinate-aligned spatial correlations", RS, "14_plot_main_figures/14e_Figure1d_spatial_correlation_uncertainty_TOP100.R"),
  c("Figure 1 adapter", RS, "14_plot_main_figures/14b_build_Figure1_TOP100_adapter.R"),
  c("Figure 2 adapter", RS, "14_plot_main_figures/14f_build_Figure2_TOP100_adapter.R"),
  c("Figure 3 adapter", RS, "14_plot_main_figures/14g_build_Figure3_TOP100_adapter.R"),
  c("Figure 1", RS, "14_plot_main_figures/build_Figure1_TOP100_FINAL.R"),
  c("Figure 2", RS, "14_plot_main_figures/build_Figure2_TOP100_FINAL.R"),
  c("Figure 2d north-to-south standalone", RS, "14_plot_main_figures/14h_build_Figure2d_NS_order_TOP100.R"),
  c("Figure 3", RS, "14_plot_main_figures/build_Figure3_TOP100_FINAL.R"),
  c("Figure 1 no-panel-labels", RS, "14_plot_main_figures/build_Figure1_TOP100_FINAL.R", "TOP100_NO_PANEL_LABELS=1"),
  c("Figure 2 no-panel-labels", RS, "14_plot_main_figures/build_Figure2_TOP100_FINAL.R", "TOP100_NO_PANEL_LABELS=1"),
  c("Figure 3 no-panel-labels", RS, "14_plot_main_figures/build_Figure3_TOP100_FINAL.R", "TOP100_NO_PANEL_LABELS=1"),
  c("SI figures", RS, "15_plot_SI_figures/build_SI_TOP100_FINAL.R"),
  c("result registry", PY, "16_build_workbook/16a_build_result_registry_TOP100.py"),
  c("data product workbook", NODE, "16_build_workbook/16d_build_DATA_PRODUCTS_AND_PROCESSING.mjs"),
  c("SI indexes and tables", NODE, "16_build_workbook/16e_build_SI_INDEX_AND_TABLES.mjs"),
  c("preliminary manifests", PY, "16_build_workbook/16c_finalize_manifests_TOP100.py"),
  c("preliminary master workbook", NODE, "16_build_workbook/16b_build_ALL_RESULTS_AND_METHODS_DETAIL.mjs"),
  c("final QA", PY, "16_build_workbook/16f_build_final_QA_TOP100.py"),
  c("post-QA manifests", PY, "16_build_workbook/16c_finalize_manifests_TOP100.py"),
  c("final master workbook", NODE, "16_build_workbook/16b_build_ALL_RESULTS_AND_METHODS_DETAIL.mjs"),
  c("final QA verification", PY, "16_build_workbook/16f_build_final_QA_TOP100.py"),
  c("longitude-alignment audit", RS, "QA/audit_longitude_alignment_TOP100.R"),
  c("final checksums", PY, "16_build_workbook/16c_finalize_manifests_TOP100.py")
)

for (x in stages) {
  if (length(x) >= 4L) run_env(x[1], x[2], x[3], x[4]) else run(x[1], x[2], x[3])
}
cat("TOP100 revision pipeline completed.\n")
