"""
Purpose: Build the formal long-form TOP100 result registry and concise writing-source tables.
Inputs: Completed TOP100 summary CSV files under output_TOP100_revision/results.
Outputs: manifests/RESULT_REGISTRY_EXPORT.csv, QA/RESULT_ID_AUDIT.csv, and
         derived_data/figure_source_data/TOP100_WRITING_SOURCE_VALUES.csv.
Parameters: Formal top=100 hPa; absolute height is main; TOP200/TOP300 and
            surface-relative values are sensitivity records.
Dependencies: Python standard library only.
Overwrite policy: Writes only inside output_TOP100_revision and replaces only
                  same-named TOP100 revision products.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

ROOT = Path(r"__WGPE_PROJECT_ROOT__")
OUT = ROOT / "output_TOP100_revision"
RESULTS = OUT / "results"
REGISTRY = OUT / "manifests" / "RESULT_REGISTRY_EXPORT.csv"
AUDIT = OUT / "QA" / "RESULT_ID_AUDIT.csv"
WRITING = OUT / "derived_data" / "figure_source_data" / "TOP100_WRITING_SOURCE_VALUES.csv"

FIELDS = [
    "result_id", "module", "analysis", "dataset", "variable", "target", "predictor",
    "metric", "estimate", "lower95", "upper95", "unit", "unit_caption", "period",
    "region", "height_reference", "top_hPa", "domain", "method", "uncertainty_method",
    "n_grid", "result_status", "approved_for_main_TRUE_FALSE", "main_usage_role",
    "source_file", "source_row", "notes",
]

rows: list[dict[str, str]] = []
counter = 0


def clean(v):
    if v is None:
        return ""
    s = str(v).strip()
    return "" if s.lower() in {"na", "nan", "none"} else s


def fnum(v):
    s = clean(v)
    if not s:
        return ""
    try:
        x = float(s)
        if not math.isfinite(x):
            return ""
        return format(x, ".15g")
    except ValueError:
        return s


def caption_unit(unit: str) -> str:
    u = clean(unit)
    mapping = {
        "m yr-1": "m yr⁻¹",
        "m/yr": "m yr⁻¹",
        "kg m-2 yr-1": "kg m⁻² yr⁻¹",
        "J m-2 yr-1": "J m⁻² yr⁻¹",
        "J m-2": "J m⁻²",
        "kg m-2": "kg m⁻²",
        "m": "m",
        "%": "%",
        "Pearson r": "Pearson r",
        "Partial r": "Partial r",
        "R2": "R²",
        "delta R2": "ΔR²",
        "RMSE": "RMSE",
        "AUC": "AUC",
        "fraction": "fraction",
    }
    return mapping.get(u, u)


def add(*, module, analysis, variable="", metric, estimate, source, source_row,
        dataset="ERA5", target="", predictor="", lower95="", upper95="", unit="",
        period="", region="Global", height_reference="absolute", top_hPa="100",
        domain="all valid grid cells", method="", uncertainty_method="", n_grid="",
        status="main_approved", approved=True, role="main_core", notes=""):
    global counter
    counter += 1
    estimate_clean = fnum(estimate)
    if estimate_clean == "":
        status = "insufficient_sample"
        approved = False
        role = "extended_data"
    rid = f"TOP100R{counter:06d}"
    row = {k: "" for k in FIELDS}
    row.update({
        "result_id": rid,
        "module": module,
        "analysis": analysis,
        "dataset": dataset,
        "variable": variable,
        "target": target,
        "predictor": predictor,
        "metric": metric,
        "estimate": estimate_clean,
        "lower95": fnum(lower95),
        "upper95": fnum(upper95),
        "unit": clean(unit),
        "unit_caption": caption_unit(unit),
        "period": clean(period),
        "region": clean(region),
        "height_reference": clean(height_reference),
        "top_hPa": clean(top_hPa),
        "domain": clean(domain),
        "method": clean(method),
        "uncertainty_method": clean(uncertainty_method),
        "n_grid": clean(n_grid),
        "result_status": status,
        "approved_for_main_TRUE_FALSE": "TRUE" if approved else "FALSE",
        "main_usage_role": role,
        "source_file": str(Path(source).relative_to(ROOT)).replace("\\", "/"),
        "source_row": str(source_row),
        "notes": clean(notes),
    })
    rows.append(row)


def read_csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def iterate(path: Path):
    for i, r in enumerate(read_csv(path), start=2):
        yield i, r


# Core trend and climatological values: TOP100 absolute main and relative sensitivity.
p = RESULTS / "core" / "Trend_results_FINAL_v3.csv"
for i, r in iterate(p):
    main = r["height_reference"] == "absolute"
    common = dict(module="core", analysis="TOP100 global and climate-band trends", source=p,
                  source_row=i, region=r["region"], height_reference=r["height_reference"],
                  period="1979-2024", n_grid=r["n_grid"],
                  domain="all-valid common TOP100 grid",
                  method="calendar-month anomalies; grid-cell OLS; cos(latitude) regional mean",
                  status="main_approved" if main else "sensitivity_only",
                  approved=main, role="main_core" if main else "sensitivity")
    for variable, col, unit in [
        ("zbar", "z_trend_m_yr", "m yr-1"),
        ("W-GPE", "wgpe_trend_J_m2_yr", "J m-2 yr-1"),
        ("IWV", "iwv_trend_kg_m2_yr", "kg m-2 yr-1"),
        ("zbar", "mean_z_m", "m"),
        ("IWV", "mean_iwv_kg_m2", "kg m-2"),
    ]:
        add(variable=variable, metric=col, estimate=r[col], unit=unit, **common)

# Exact centered-product W-GPE trend decomposition.
p = RESULTS / "core" / "WGPE_decomposition_FINAL_v3.csv"
for i, r in iterate(p):
    main = r["height_reference"] == "absolute"
    for metric in ["actual", "IWV", "z", "interaction", "sum_channels", "residual"]:
        add(module="core", analysis="exact W-GPE trend product decomposition", variable="W-GPE",
            metric=metric, estimate=r[metric], unit="J m-2 yr-1", source=p, source_row=i,
            region=r["region"], height_reference=r["height_reference"], period="1979-2024",
            n_grid=r["n_grid"], domain="all-valid common TOP100 grid",
            method="grid-cell exact centered-product decomposition, then cos(latitude) aggregation",
            status="main_approved" if main else "sensitivity_only", approved=main,
            role="main_core" if main else "sensitivity")

# Top-boundary comparison and six-layer vertical attribution.
p = RESULTS / "top_boundary" / "TOP_BOUNDARY_GLOBAL_RESULTS.csv"
for i, r in iterate(p):
    main = r["top_hPa"] == "100"
    for metric, col, unit in [
        ("climatological_mean", "global_climatological_mean", r["mean_unit"]),
        ("trend_deseasonalized", "global_trend_deseasonalized", r["trend_unit"]),
        ("trend_raw", "global_trend_raw", r["trend_unit"]),
    ]:
        add(module="top_boundary", analysis="TOP100/TOP200/TOP300 boundary comparison",
            variable=r["variable"], metric=metric, estimate=r[col], unit=unit, source=p, source_row=i,
            period=r["period"], region="Global", height_reference=r["height_reference"],
            top_hPa=r["top_hPa"], n_grid=r["n_grid"], domain=r["mask"], method=r["trend_method"],
            status="main_approved" if main else "sensitivity_only", approved=main,
            role="main_core" if main else "sensitivity",
            notes="TOP200 is convergence sensitivity; TOP300 is truncation sensitivity")

p = RESULTS / "top_boundary" / "TOP_BOUNDARY_LAYER_CONTRIBUTIONS.csv"
for i, r in iterate(p):
    for metric, col, unit in [
        ("layer_mass_trend", "layer_mass_trend_kg_m2_yr", "kg m-2 yr-1"),
        ("layer_height_trend", "layer_height_trend_m_yr", "m yr-1"),
        ("layer_WGPE_trend", "layer_WGPE_trend_J_m2_yr", "J m-2 yr-1"),
        ("mass_contribution", "mass_contribution_m_yr", "m yr-1"),
        ("height_coordinate_contribution", "height_contribution_m_yr", "m yr-1"),
        ("total_first_order_contribution", "total_first_order_contribution_m_yr", "m yr-1"),
    ]:
        add(module="vertical_layers", analysis="six-layer first-order zbar trend attribution",
            variable="zbar", metric=f"{r['layer']}:{metric}", estimate=r[col], unit=unit,
            source=p, source_row=i, region=r["region"], period="1979-2024",
            n_grid=r["n_grid"], domain="layer-valid TOP100 grid", method="first-order quotient decomposition",
            role="main_core" if r["region"] == "Global" else "main_companion_descriptive")

# Surface endpoint methods are sensitivity-only except P0, which reproduces the formal core.
p = RESULTS / "surface_truncation" / "SURFACE_METHOD_COMPARISON.csv"
for i, r in iterate(p):
    main = r["method"].startswith("P0_")
    add(module="surface_truncation", analysis="surface endpoint method comparison",
        variable=r["variable"], metric="trend", estimate=r["trend"], unit=r["unit"], source=p,
        source_row=i, region=r["region"], period=r["period"], n_grid=r["n_grid"],
        method=r["method"], status="main_approved" if main else "sensitivity_only",
        approved=main, role="main_core" if main else "sensitivity")

# JRA-55 validation.
p = RESULTS / "JRA55" / "JRA55_TOP100_GLOBAL_TRENDS.csv"
for i, r in iterate(p):
    add(module="JRA55", analysis="cross-reanalysis TOP100 trend validation", dataset=r["dataset"],
        variable=r["variable"], metric="trend", estimate=r["trend"], lower95=r["newey_west_ci_low"],
        upper95=r["newey_west_ci_high"], unit=r["unit"], source=p, source_row=i,
        region=r["region"], period=r["period"], n_grid=r["n_grid"],
        method="native-grid vertical integration followed by harmonization to common 1-degree grid",
        uncertainty_method=f"Newey-West lag {r['newey_west_lag']}; 5-year block bootstrap",
        role="extended_data")

# Precipitation association and BH/BY field significance.
p = RESULTS / "precip_wetdry" / "Precipitation_correlation_summary_FINAL_v3.csv"
for i, r in iterate(p):
    main = r["height_reference"] == "absolute"
    for metric, col, unit in [
        ("mean_Pearson_r", "mean_r_area", "Pearson r"),
        ("BH_positive_area", "positive_significant_area_pct", "%"),
        ("BH_negative_area", "negative_significant_area_pct", "%"),
        ("BH_total_area", "total_significant_area_pct", "%"),
    ]:
        add(module="precipitation", analysis="W-GPE precipitation association", variable="W-GPE",
            target="precipitation", predictor="W-GPE", metric=metric, estimate=r[col], unit=unit,
            source=p, source_row=i, region=r["region"], height_reference=r["height_reference"],
            period="1979-2024", n_grid=r["n_grid"], method="calendar-month climatology removed and each grid-cell series linearly detrended; direct Pearson r (IWV not controlled); AR(1) effective n; BH-FDR",
            status="main_approved" if main else "sensitivity_only", approved=main,
            role="main_core" if main else "sensitivity")

p = RESULTS / "correlation_significance" / "WGPE_precip_significance_summary_TOP100.csv"
for i, r in iterate(p):
    for metric, col, unit in [
        ("mean_Pearson_r", "mean_r_coslat", "Pearson r"),
        ("raw:p05_area", "raw_p05_area_pct", "%"),
        ("BH:positive_significant_area", "BH_positive_area_pct", "%"),
        ("BH:negative_significant_area", "BH_negative_area_pct", "%"),
        ("BH:total_significant_area", "BH_q05_area_pct", "%"),
        ("BY:positive_significant_area", "BY_positive_area_pct", "%"),
        ("BY:negative_significant_area", "BY_negative_area_pct", "%"),
        ("BY:total_significant_area", "BY_q05_area_pct", "%"),
        ("effective_spatial_dof_variance_ratio", "effective_spatial_dof_variance_ratio", ""),
    ]:
        add(module="field_significance", analysis="precipitation BH/BY sensitivity", variable="W-GPE",
            target="precipitation", predictor="W-GPE", metric=metric,
            estimate=r[col], unit=unit, source=p, source_row=i, period="1979-2024", n_grid=r["n_valid_grid"],
            method="AR(1) effective n; BH/BY-FDR; common circular shift and latitude-stratified spatial blocks",
            role="extended_data")

# Exact wet/dry composite decomposition with spatial-block uncertainty.
p = RESULTS / "precip_wetdry" / "exact_TOP100" / "WET_DRY_EXACT_SUMMARY_TOP100.csv"
for i, r in iterate(p):
    add(module="wet_dry", analysis="exact grid-month wet/dry W-GPE decomposition",
        variable="W-GPE", target=r["event"], predictor="precipitation event", metric=r["metric"],
        estimate=r["estimate"], lower95=r["lower95"], upper95=r["upper95"],
        unit="m" if r["metric"] == "precip_anomaly" else "J m-2",
        source=p, source_row=i, region=r["region"], period="1979-2024", n_grid=r["n_grid"],
        method="grid-cell >90th/<10th percentile of detrended deseasoned precipitation; exact monthly-climatology product decomposition",
        uncertainty_method=r["uncertainty_method"], role="main_core")

# ONI-standard ENSO exact decomposition.
p = RESULTS / "ENSO" / "ONI_standard_TOP100" / "03_absolute_results" / "ENSO_ONIstandard_absolute_channel_decomposition.csv"
for i, r in iterate(p):
    add(module="ENSO", analysis="ONI-standard El Nino minus La Nina exact decomposition",
        variable="W-GPE", target="ENSO phase", predictor="ONI", metric=r["metric"],
        estimate=r["estimate"], lower95=r["lower95"], upper95=r["upper95"], unit="J m-2",
        source=p, source_row=i, region=r["region"], period="1979-2024", n_grid="",
        method="Oct-Mar ONI-standard event-window composite; exact grid-month decomposition",
        uncertainty_method="episode-cluster x 10-degree spatial-block bootstrap, 1000 resamples",
        notes=f"percent contribution={clean(r.get('percent_contribution'))}; unstable={clean(r.get('unstable_percentage'))}")

# Ecology partial correlations and BH/BY sensitivity.
p = RESULTS / "ecology" / "Ecology_partial_correlations_summary_FINAL_v3.csv"
for i, r in iterate(p):
    main = r["height_reference"] == "absolute" and r["model"] == "Model3"
    for metric, col, unit in [
        ("mean_partial_r", "mean_partial_r_area", "Partial r"),
        ("BH_positive_area", "positive_significant_area_pct", "%"),
        ("BH_negative_area", "negative_significant_area_pct", "%"),
    ]:
        add(module="ecology", analysis="ecological partial correlation", variable=r["predictor"],
            target=r["target"], predictor=r["predictor"], metric=f"{r['model']}:{metric}",
            estimate=r[col], unit=unit, source=p, source_row=i, region=r["region"],
            height_reference=r["height_reference"], period="product-specific ecology period", n_grid=r["n_grid"],
            domain="persistent vegetated land", method=f"partial correlation controlling {r['controls_used']}; AR(1) effective n; BH-FDR",
            status="main_approved" if main else "sensitivity_only", approved=main,
            role="main_core" if main else "sensitivity")

p = RESULTS / "correlation_significance" / "ECOLOGY_PARTIAL_CORRELATION_BH_BY_SUMMARY_TOP100.csv"
for i, r in iterate(p):
    for metric, col, unit in [
        ("mean_partial_r", "mean_partial_r_coslat", "Partial r"),
        ("positive_significant_area", "positive_significant_area_pct", "%"),
        ("negative_significant_area", "negative_significant_area_pct", "%"),
        ("total_significant_area", "total_significant_area_pct", "%"),
    ]:
        add(module="field_significance", analysis="ecology BH/BY sensitivity", variable=r["predictor"],
            target=r["target"], predictor=r["predictor"], metric=f"{r['model']}:{r['FDR_method']}:{metric}",
            estimate=r[col], unit=unit, source=p, source_row=i, n_grid=r["n_grid"],
            domain="persistent vegetated land", method=f"AR(1) effective n and {r['FDR_method']}-FDR", role="extended_data")

# Four-class conditional temporal precedence summaries.
p = RESULTS / "temporal_precedence" / "TOP100_conditional" / "conditional_temporal_precedence_summary_TOP100.csv"
for i, r in iterate(p):
    est_col = "area_fraction_pct" if "area_fraction_pct" in r else ("area_pct" if "area_pct" in r else "fraction_pct")
    add(module="temporal_precedence", analysis="conditional four-class temporal precedence",
        variable="W-GPE", target=r.get("target", ""), predictor=r.get("predictor", "W-GPE"),
        metric=f"lag{r.get('lag','')}:{r.get('direction_class',r.get('class',''))}", estimate=r.get(est_col, ""),
        unit="%", source=p, source_row=i, region=r.get("region", "Global"),
        height_reference=r.get("height_reference", "absolute"), period=r.get("period", "product-specific"),
        n_grid=r.get("n_grid", ""), domain=r.get("domain", "persistent vegetated land"),
        method="nested F test, conditional controls, BH-FDR separately by target/predictor/lag/direction",
        role="main_companion_descriptive" if r.get("lag") == "3" and r.get("height_reference") == "absolute" else "extended_data")

# Temporal effect sizes.
for name in ["PRECIP_TEMPORAL_EFFECT_SIZE_SUMMARY_TOP100.csv", "ECOLOGY_TEMPORAL_EFFECT_SIZE_SUMMARY_TOP100.csv"]:
    p = RESULTS / "temporal_precedence" / "TOP100_effect_sizes" / name
    for i, r in iterate(p):
        unit = "delta R2" if "R2" in r["metric"] else "RMSE"
        add(module="temporal_precedence", analysis="temporal precedence effect size",
            variable="W-GPE", target=r["target"], predictor=r["predictor"], metric=f"lag{r['lag']}:{r['metric']}",
            estimate=r["estimate"], lower95=r["lower95"], upper95=r["upper95"], unit=unit,
            source=p, source_row=i, region=r["region"], n_grid=r["n_grid"],
            domain=r.get("domain", "all valid grid cells"), method="nested lag model and blocked cross-validation",
            uncertainty_method=r["uncertainty_method"], role="main_companion_descriptive" if r["lag"] == "3" else "extended_data")

# Incremental value robustness registry.
p = RESULTS / "incremental_value" / "INCREMENTAL_VALUE_ROBUSTNESS_TOP100.csv"
for i, r in iterate(p):
    add(module="incremental_value", analysis="increment beyond IWV baseline", variable=r["increment"],
        target=r["target"], predictor="IWV + added diagnostic", metric=r["primary_metric"],
        estimate=r["global_estimate"], lower95=r["global_lower95"], upper95=r["global_upper95"],
        unit="delta R2" if "R2" in r["primary_metric"] else "AUC", source=p, source_row=i,
        period="target-specific", n_grid="", domain="shared-fold valid domain",
        method=f"blocked cross-validation, gap={r['gap_months']} months; {r['classification_rule']}",
        role="main_companion_descriptive", notes=f"robustness={r['robustness_class']}; positive area={r['positive_increment_area_pct']}%")

# Write outputs and audit uniqueness/status consistency.
REGISTRY.parent.mkdir(parents=True, exist_ok=True)
WRITING.parent.mkdir(parents=True, exist_ok=True)
AUDIT.parent.mkdir(parents=True, exist_ok=True)
with REGISTRY.open("w", encoding="utf-8-sig", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    w.writeheader(); w.writerows(rows)

writing_rows = [r for r in rows if r["approved_for_main_TRUE_FALSE"] == "TRUE" and r["main_usage_role"] in {"main_core", "main_companion_descriptive"}]
with WRITING.open("w", encoding="utf-8-sig", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    w.writeheader(); w.writerows(writing_rows)

ids = [r["result_id"] for r in rows]
bad_approval = [r for r in rows if r["approved_for_main_TRUE_FALSE"] == "TRUE" and r["result_status"] in {"pending", "legacy_excluded", "rejected", "sensitivity_only"}]
audit_rows = [
    {"check": "total_registry_rows", "value": len(rows), "status": "PASS"},
    {"check": "unique_result_ids", "value": len(set(ids)), "status": "PASS" if len(ids) == len(set(ids)) else "FAIL"},
    {"check": "duplicate_result_ids", "value": len(ids) - len(set(ids)), "status": "PASS" if len(ids) == len(set(ids)) else "FAIL"},
    {"check": "main_approved_writing_rows", "value": len(writing_rows), "status": "PASS"},
    {"check": "invalid_approved_status_rows", "value": len(bad_approval), "status": "PASS" if not bad_approval else "FAIL"},
]
with AUDIT.open("w", encoding="utf-8-sig", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["check", "value", "status"])
    w.writeheader(); w.writerows(audit_rows)

print(f"registry rows={len(rows)} unique={len(set(ids))} writing={len(writing_rows)}")
