# September 2026 five-control ecology and figure update

This development update is based on v1.0.0. No new release tag or DOI is assigned. Manuscript and SI Word/PDF files are not distributed here.

## Current ecology results

The primary concurrent analysis controls precipitation, T2m, SSRD, root-zone soil moisture and CO2. All time series are deseasonalized and linearly detrended. Both models use the same complete-case months established with all six variables (including VPD) and the same persistent-vegetation domain. Thus removing VPD changes the adjustment set, not the missing-data selection.

| Target | Predictor | Valid cells | Mean partial r |
|---|---|---:|---:|
| GPP | W-GPE | 12,009 | 0.1464207755 |
| GPP | <z> | 12,009 | -0.0339226698 |
| LAI | W-GPE | 8,300 | 0.0669563850 |
| LAI | <z> | 8,300 | -0.0154267647 |

Global and climate-zone means use cos(latitude) weights over **all valid cells**, including nonsignificant map cells. Missing/invalid cells are excluded. BH/BY significant percentages use the same valid-area denominator. The supplied spatial confidence intervals use 1,000 resamples of 10-degree blocks and seed 20260904. Exact results and paired six-control comparisons are in `source_data/ecology_20260904/`.

`Fig3a_source.csv` through `Fig3d_source.csv` retain their historical latitude-profile schema; full map grids are in `GRID_RESULTS.csv`. `Fig3e_source.csv` now contains the five-control regional summaries. The original files are retained under `source_data/historical/v1.0.0/`. Current registry entries have new `ECO5_20260904_` identifiers. Historical result IDs keep their original estimates and are not relabelled as five-control results.

## Temporal classification

The lag-3 analysis jointly tests lags 1–3 with calendar lags preserved. It controls the five covariates' lagged values and target autoregression, uses BH correction separately in each direction, and holds six-control valid lag windows fixed for the paired comparison. Both tests must be valid for classification. GPP has 12,009 tested cells; LAI has 5,788, different from the concurrent LAI domain.

The f panel plots **grid_pct**, the fraction of tested grid cells. `area_pct` is also provided as a distinct statistic; it is not interchangeable with grid_pct. New temporal-class confidence intervals or ecological OOF scores were not computed. The retained v1.0.0 field-significance, coverage, HAC, OOF and other sensitivity files describe their original six-control analyses and original domains. They are not updated five-control tests.

## Approved figure outputs and uncertainty

Reference PNGs under `figures/approved_20260916/` are copies of the approved originals, without resizing or rerendering. Historical code labels Figure 3 refer to the ecology figure delivered as Supplementary Figure S1.

- Figure 1: `assemble_channel_dots_ci.R`.
- Figure 2: `align_fig2_cde_titles.R` (the aligned c/d/e version).
- Supplementary Figure S1: `fix_latitude_alignment_v7.R`, built on `layout_bottom_v5.R`.

The archived dependency chain and entrypoint map are in `code/07_figures_and_tables/approved_20260916/`. Panel-e stars indicate at least 5% BH-significant area under the existing plotting rule, not a circular-shift p value or a confidence interval excluding zero.

`source_data/channel_intervals_20260916/Fig1g_values.csv` contains channel-specific 10-degree spatial-block intervals (1,000 replicates) from the Figure 1 adapter's `U_CORE`. `Fig2d_values.csv` contains channel-specific **episode-cluster** ENSO intervals from the Figure 2 adapter's `enso_uncertainty`; estimates and endpoints are divided by 1,000 for the plotted units. These are distinct uncertainty designs. Total-anomaly intervals were not substituted for channel intervals.

## Reproduction levels

1. **Included-data verification:** run `run_source_data_workflow.ps1` or `.sh`. No ERA5 download or R installation is required. The additional test reconstructs means, area fractions and temporal class percentages from supplied grids.
2. **Ecological recalculation from existing inputs:** scripts in `code/05_ecology/september_five_control/` require the external canonical RDS inputs below. Set `WGPE_PROJECT_ROOT` to their parent project and optionally `WGPE_ECOLOGY_OUTPUT` to a new writable output directory. `run_no_vpd.R` uses base R and ggplot2 (openxlsx optional); `compute_f_no_vpd.R` uses base R. Analysis algorithms are retained; the changes are path configuration and removal of local R-library/toolchain overrides. The inherited concurrent AR(1) calculation uses complete-case residual sequences, whereas the temporal models explicitly preserve calendar lag windows.
3. **Exact full-figure rendering:** the approved script chain is a project-layout snapshot, not a standalone figure API. Its project-relative paths, evaluation order and expression-based wrappers are retained. Set `WGPE_PROJECT_ROOT`, install the required R packages, and place the snapshot scripts at their original project-relative paths in a separate reproduction workspace. Adapters, legacy grid inputs, coastline data and other files named by those scripts are still required. The PNGs provide reference outputs. This update has not rerun the entire raw-data pipeline or every full-figure dependency in a clean environment.

External inputs for level 2:

```
output_TOP100_revision/FINAL_COORDINATE_FIXED_CLOSURE/
  derived_data/ecology/GPP_CANONICAL_INPUTS.rds
  derived_data/ecology/LAI_CANONICAL_INPUTS.rds
  results/ecology/ECOLOGY_PARTIAL_SUMMARY_CANONICAL.csv
  results/ecology/ECOLOGY_TEMPORAL_GRID_CANONICAL.csv
```

The included gridded CSVs support audit and plotting; they do not replace the monthly RDS inputs for refitting models. Existing `run_full_analysis.ps1` retains its original workflow; the two September calculation entrypoints above are explicit additional steps.
