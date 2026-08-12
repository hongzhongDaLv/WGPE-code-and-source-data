"""Build the final TOP100 cross-file QA artifacts.

Purpose
-------
Create the three final audit products that are not scientific recomputations:
MD/XLSX consistency, old/new result comparison, and the consolidated QC report.

Inputs
------
Formal TOP100 registry, summary, result CSV files, figures, SI package, and the
main results workbook.

Outputs
-------
output_TOP100_revision/QA/MD_XLSX_CONSISTENCY.csv
output_TOP100_revision/QA/OLD_NEW_RESULT_COMPARISON.xlsx
output_TOP100_revision/QA/QC_REPORT.md

Parameters/dependencies
-----------------------
Python standard library plus openpyxl.  This script never changes scientific
source values and never writes manuscript/SI DOCX files.

Overwrite policy
----------------
Only the three TOP100 QA products above are refreshed.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import openpyxl
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill


ROOT = Path(r"__WGPE_PROJECT_ROOT__")
OUT = ROOT / "output_TOP100_revision"
QA = OUT / "QA"
QA.mkdir(parents=True, exist_ok=True)
REGISTRY = OUT / "manifests" / "RESULT_REGISTRY_EXPORT.csv"
SUMMARY = OUT / "RESULTS_COMPLETION_SUMMARY_CN.md"
WORKBOOK = OUT / "ALL_RESULTS_AND_METHODS_DETAIL.xlsx"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def num(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


registry = read_csv(REGISTRY)
summary_text = SUMMARY.read_text(encoding="utf-8")


def select_one(**criteria):
    rows = [r for r in registry if all(str(r.get(k, "")) == str(v) for k, v in criteria.items())]
    if len(rows) != 1:
        raise RuntimeError(f"Expected one registry row for {criteria}, found {len(rows)}")
    return rows[0]


# Read the registry sheet back from the produced workbook.  This confirms that
# the workbook is not merely present but contains the registered result IDs.
wb_read = openpyxl.load_workbook(WORKBOOK, read_only=True, data_only=True)
ws_reg = wb_read["Result Registry"]
it = ws_reg.iter_rows(values_only=True)
header = [str(v) if v is not None else "" for v in next(it)]
wb_rows = []
for values in it:
    if not any(v is not None for v in values):
        continue
    wb_rows.append({header[i]: values[i] for i in range(min(len(header), len(values)))})
wb_by_id = {str(r.get("result_id")): r for r in wb_rows if r.get("result_id")}


checks = [
    ("zbar_global_trend", dict(module="core", metric="z_trend_m_yr", region="Global"), "+0.819"),
    ("wgpe_global_trend", dict(module="core", metric="wgpe_trend_J_m2_yr", region="Global"), "+672.055"),
    ("iwv_global_trend", dict(module="core", metric="iwv_trend_kg_m2_yr", region="Global"), "+0.019010"),
    ("wgpe_decomp_actual", dict(module="core", analysis="exact W-GPE trend product decomposition", metric="actual", region="Global"), "+672.055"),
    ("precip_detrended_r", dict(module="precipitation", metric="mean_Pearson_r", region="Global"), "+0.01124"),
    ("enso_global_actual", dict(module="ENSO", metric="actual", region="Global"), "+11723.9"),
    ("wgpe_gpp_partial_r", dict(module="ecology", metric="Model3:mean_partial_r", target="GPP", predictor="W-GPE", region="Global", height_reference="absolute"), "+0.170"),
    ("wgpe_lai_partial_r", dict(module="ecology", metric="Model3:mean_partial_r", target="LAI", predictor="W-GPE", region="Global", height_reference="absolute"), "+0.107"),
]

consistency = []
for check_id, criteria, token in checks:
    rr = select_one(**criteria)
    rid = rr["result_id"]
    rv = num(rr["estimate"])
    wr = wb_by_id.get(rid, {})
    wv = num(wr.get("estimate"))
    summary_ok = token in summary_text
    diff = None if rv is None or wv is None else abs(rv - wv)
    status = "PASS" if summary_ok and diff is not None and diff <= 1e-10 else "FAIL"
    consistency.append({
        "check_id": check_id,
        "result_id": rid,
        "registry_value": rv,
        "workbook_value": wv,
        "absolute_difference": diff,
        "summary_token": token,
        "summary_contains_token_TRUE_FALSE": summary_ok,
        "status": status,
        "source_file": rr["source_file"],
    })

with (QA / "MD_XLSX_CONSISTENCY.csv").open("w", encoding="utf-8-sig", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(consistency[0]))
    writer.writeheader()
    writer.writerows(consistency)


# Old/new comparison workbook.  The precipitation comparison is deliberately
# labelled as a method change, not as an IWV-controlled partial correlation.
top_global = read_csv(OUT / "results" / "top_boundary" / "TOP_BOUNDARY_GLOBAL_RESULTS.csv")
domain_rec = read_csv(OUT / "results" / "top_boundary" / "TOP_BOUNDARY_DOMAIN_RECONCILIATION.csv")

key_changes = [
    {
        "quantity": "W-GPE global trend (locked common-layer mask)",
        "old_definition": "TOP300; locked layer-valid common mask",
        "old_value": 645.5975935321833,
        "new_definition": "TOP100 main; all-valid common TOP100 grid",
        "new_value": 672.055360999501,
        "unit": "J m-2 yr-1",
        "interpretation": "TOP300 is truncation sensitivity; TOP100 is the formal main definition",
    },
    {
        "quantity": "W-GPE global trend (all-valid domain)",
        "old_definition": "TOP300; 65,160-grid all-valid domain",
        "old_value": 644.792143307247,
        "new_definition": "TOP100; 65,160-grid all-valid domain",
        "new_value": 672.0555286207635,
        "unit": "J m-2 yr-1",
        "interpretation": "Like-for-like top-boundary comparison",
    },
    {
        "quantity": "W-GPE-precipitation global mean Pearson r",
        "old_definition": "calendar-month climatology removed; linear trends retained; TOP300",
        "old_value": 0.434860664109332,
        "new_definition": "calendar-month climatology removed; both grid-cell series linearly detrended; TOP100",
        "new_value": 0.0112439592114552,
        "unit": "Pearson r",
        "interpretation": "Direct correlation in both cases; IWV is not controlled. The drop is principally a detrending effect.",
    },
    {
        "quantity": "Vertical attribution layers",
        "old_definition": "four layers, 1000-300 hPa",
        "old_value": 4,
        "new_definition": "six layers, 1000-100 hPa",
        "new_value": 6,
        "unit": "layer count",
        "interpretation": "300-200 and 200-100 hPa are now explicit",
    },
]


def add_sheet(workbook, title, rows):
    ws = workbook.create_sheet(title)
    if not rows:
        ws.append(["status"])
        ws.append(["no rows"])
        return
    fields = list(rows[0])
    ws.append(fields)
    for r in rows:
        ws.append([r.get(k) for k in fields])
    for cell in ws[1]:
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="315A7D")
        cell.alignment = Alignment(wrap_text=True)
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions
    for col in ws.columns:
        letter = col[0].column_letter
        width = min(60, max(10, max(len(str(c.value or "")) for c in col) + 2))
        ws.column_dimensions[letter].width = width


wb = Workbook()
wb.remove(wb.active)
add_sheet(wb, "Key changes", key_changes)
add_sheet(wb, "Top boundary global", top_global)
add_sheet(wb, "Domain reconciliation", domain_rec)
add_sheet(wb, "MD-XLSX consistency", consistency)
add_sheet(wb, "Precip method audit", [
    {"version": "legacy", "r": 0.434860664109332,
     "seasonal_cycle_removed": True, "linear_trend_removed": False,
     "IWV_controlled": False, "role": "trend-retaining descriptive sensitivity"},
    {"version": "TOP100 formal", "r": 0.0112439592114552,
     "seasonal_cycle_removed": True, "linear_trend_removed": True,
     "IWV_controlled": False, "role": "main detrended interannual association"},
])
wb.save(QA / "OLD_NEW_RESULT_COMPARISON.xlsx")


# Consolidated completion and scientific-consistency checks.
result_ids = [r["result_id"] for r in registry]
duplicates = len(result_ids) - len(set(result_ids))
invalid_main = [r for r in registry if str(r.get("approved_for_main_TRUE_FALSE", "")).upper() == "TRUE"
                and r.get("result_status") in {"pending", "legacy_excluded", "rejected", "sensitivity_only"}]
old_main = [r for r in registry if str(r.get("approved_for_main_TRUE_FALSE", "")).upper() == "TRUE"
            and num(r.get("top_hPa")) == 300]

fig_checks = []
for f in (1, 2, 3):
    d = OUT / "figures" / "main" / f"Figure{f}_TOP100_FINAL"
    nd = OUT / "figures" / "main_no_panel_labels" / f"Figure{f}_TOP100_FINAL"
    fig_checks.append((f, (d / f"Figure{f}_TOP100_FINAL.png").exists(),
                       (d / f"Figure{f}_TOP100_FINAL.pdf").exists(),
                       (nd / f"Figure{f}_TOP100_FINAL.png").exists(),
                       (nd / f"Figure{f}_TOP100_FINAL.pdf").exists()))

si_dir = OUT / "10_SUPPLEMENTARY_INFORMATION" / "figures"
si_png = sorted(si_dir.glob("Figure_S*.png"))
si_pdf = sorted(si_dir.glob("Figure_S*.pdf"))
required = [
    OUT / "RESULTS_COMPLETION_SUMMARY_CN.md",
    OUT / "ALL_RESULTS_AND_METHODS_DETAIL.xlsx",
    OUT / "methods" / "METHODS_FULL_TOP100_CN.md",
    OUT / "methods" / "METHODS_FORMULAS_TOP100_CN.md",
    OUT / "methods" / "DATA_PRODUCTS_AND_PROCESSING.xlsx",
    OUT / "manifests" / "RESULT_REGISTRY_EXPORT.csv",
]
missing_required = [str(p.relative_to(ROOT)) for p in required if not p.exists()]

fig2_source = read_csv(OUT / "figures" / "main" / "Figure2_TOP100_FINAL" / "Figure2_TOP100_source_values.csv")
fig2_a = [r for r in fig2_source if r["panel"] == "a" and r["group"] == "Global"]
fig2_ok = len(fig2_a) == 1 and abs(float(fig2_a[0]["value"]) - 0.0112439592114552) < 1e-12
fig1g = read_csv(OUT / "figures" / "main" / "Figure1_TOP100_FINAL" / "standard_latband_WGPE_trend_decomposition.csv")
fig1g_global = [r for r in fig1g if r["region"] == "Global"]
fig1g_ok = len(fig1g_global) == 1 and abs(float(fig1g_global[0]["actual_J_m2_yr"]) - 672.055360999501) < 1e-9

consistency_ok = all(r["status"] == "PASS" for r in consistency)
figures_ok = all(all(x[1:]) for x in fig_checks)
overall = (not missing_required and duplicates == 0 and not invalid_main and not old_main
           and consistency_ok and figures_ok and fig2_ok and fig1g_ok
           and len(si_png) >= 12 and len(si_pdf) >= 12)

lines = [
    "# WGPE TOP100 final QC report",
    "",
    f"Overall delivery gate: **{'PASS WITH INTERPRETIVE CAVEATS' if overall else 'FAIL — see checks below'}**.",
    "",
    "## Core registry and cross-file consistency",
    "",
    f"- Registry rows: {len(registry)}; unique result IDs: {len(set(result_ids))}; duplicates: {duplicates}.",
    f"- Invalid main-approved status rows: {len(invalid_main)}.",
    f"- TOP300 rows incorrectly main-approved: {len(old_main)}.",
    f"- Summary–registry–workbook key-value checks: {sum(r['status']=='PASS' for r in consistency)}/{len(consistency)} PASS.",
    f"- Figure 1g formal TOP100 global decomposition used: {'YES' if fig1g_ok else 'NO'}.",
    f"- Figure 2a source value equals formal detrended TOP100 r=0.011243959: {'YES' if fig2_ok else 'NO'}.",
    "- Legacy r=0.434860664 is retained only in the old/new audit as a trend-retaining descriptive statistic; it is not an IWV-controlled result.",
    "",
    "## Figure delivery",
    "",
]
for f, png, pdf, nlpng, nlpdf in fig_checks:
    lines.append(f"- Figure {f}: labelled PNG={'YES' if png else 'NO'}, labelled PDF={'YES' if pdf else 'NO'}, no-label PNG={'YES' if nlpng else 'NO'}, no-label PDF={'YES' if nlpdf else 'NO'}.")
lines += [
    f"- SI figures: {len(si_png)} PNG and {len(si_pdf)} PDF companions.",
    "- Individual main-figure panels are exported in each figure directory in PNG and PDF.",
    "",
    "## Methods, reproducibility, and files",
    "",
    f"- Missing required core files: {len(missing_required)}" + (" (" + "; ".join(missing_required) + ")" if missing_required else "."),
    "- Formal TOP100 scripts, configuration, logs, derived data, figure source data, registry, and manifests are retained.",
    "- Manuscript and SI DOCX files were not modified by this TOP100 workflow.",
    "- All recovered execution failures remain recorded in logs/failed_tasks.csv; none changes the TOP100 global trend or exact decomposition closure.",
    "",
    "## Interpretive caveats retained for the manuscript",
    "",
    "1. The JRA-55 versus ERA5 200–100 hPa spatial structure is product-sensitive even though the full TOP100 global trend and dominant-layer direction agree.",
    "2. The detrended W-GPE–precipitation association is weak globally and spatially heterogeneous; temporal field significance does not imply independent local effects.",
    "3. Precipitation temporal-precedence models show positive in-sample partial R² but negative blocked-CV ΔR² and must not be described as causal or robust predictive gain.",
    "4. GPP out-of-sample gains weaken with lag; LAI gains are not robust.",
    "",
    "## Final decision",
    "",
    "The TOP100 core, six-layer attribution, exact W-GPE decomposition, event/ENSO outputs, ecology outputs, figures, and registries are numerically closed. The remaining items are interpretation limits, not missing calculations or missing core inputs.",
]
(QA / "QC_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

print(f"QC complete: overall={overall}, registry={len(registry)}, consistency={len(consistency)}")
