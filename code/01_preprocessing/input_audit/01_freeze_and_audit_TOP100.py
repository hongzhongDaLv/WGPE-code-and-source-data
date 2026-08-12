"""
Purpose: freeze the legacy TOP300/manuscript packages and register all core TOP100 inputs.
Inputs: project inputs and locked legacy outputs listed below.
Outputs: OLD_RESULTS_SHA256.csv, INPUT_INVENTORY.csv, PROJECT_AUDIT.md.
Parameters: no scientific parameters.
Dependencies: Python standard library.
Overwrite policy: writes only under output_TOP100_revision/temp; never edits source files.
"""
from __future__ import annotations
import csv, hashlib, os
from datetime import datetime
from pathlib import Path

ROOT=Path(r"__WGPE_PROJECT_ROOT__")
OUT=ROOT/"output_TOP100_revision"/"temp"
OUT.mkdir(parents=True,exist_ok=True)

def sha256(p:Path)->str:
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()

critical=[
 ROOT/"WGPE_manuscript_20250529.docx",
 ROOT/"WGPE_SI 20250529.docx",
 ROOT/"output_attribution_minimal/global_wgpe_1deg_surface_truncated_300hPa_FIXED_FULL.RData",
 ROOT/"output_attribution_minimal/ERA5_surface_pressure_1deg_1979_2024_FIXED_FULL.rds",
 ROOT/"output_attribution_minimal/FULL_RERUN_FINAL_absolute_relative_v3/11_registry_combined/Writing_source_values_FINAL_v3.csv",
 ROOT/"output_core_validation/ALL_CORE_VALIDATION_RESULTS.xlsx",
 ROOT/"WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard_FINAL/07_QA_AND_AUDIT/FINAL_VALUE_LOCK_v6_ONIstandard.md",
 ROOT/"WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard_FINAL/11_ENSO_ONI_STANDARD_UPDATE/11_QA/ENSO_event_count_lock_FINAL.csv",
 ROOT/"WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard_FINAL/01_MAIN_FIGURES/Figure1_TOP100_FINAL.png", # expected absent; records package state
]
# Add the actual final legacy main figures without assuming exact filenames.
figdir=ROOT/"WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard_FINAL"/"01_MAIN_FIGURES"
if figdir.exists(): critical.extend(sorted(p for p in figdir.rglob("*.png") if "text_plus2pt" not in str(p))[:12])

rows=[]
for p in dict.fromkeys(critical):
    exists=p.exists() and p.is_file()
    rows.append({"absolute_path":str(p),"exists":exists,"size_bytes":p.stat().st_size if exists else "",
      "last_modified":datetime.fromtimestamp(p.stat().st_mtime).isoformat(timespec="seconds") if exists else "",
      "sha256":sha256(p) if exists else "","role":"legacy_locked_input_or_output","modified_by_TOP100_revision":False})
with (OUT/"OLD_RESULTS_SHA256.csv").open("w",newline="",encoding="utf-8-sig") as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

inputs=[]
def add(category,p,variable,period,resolution,role):
    p=Path(p);inputs.append({"category":category,"absolute_path":str(p),"exists":p.exists(),
      "size_bytes":p.stat().st_size if p.exists() and p.is_file() else "","variable":variable,
      "period":period,"resolution":resolution,"analysis_role":role,"read_only":True})
for p in sorted((ROOT/"era5_plev_ascii").glob("*.nc")):
    add("ERA5 pressure levels",p,"q;geopotential","1979-2024","1 degree; 37 pressure levels","TOP100/TOP200/TOP300 core")
for cat,p,var,per,res,role in [
 ("ERA5 surface pressure",ROOT/"output_attribution_minimal/ERA5_surface_pressure_1deg_1979_2024_FIXED_FULL.rds","sp","1979-2024","1 degree","surface truncation"),
 ("ERA5 surface geopotential",ROOT/"geopotential.nc","z","1979","0.25 degree","surface endpoint/relative height"),
 ("ERA5 dewpoint",ROOT/"2m dewpoint temperature.nc","d2m","1979-2026","0.25 degree","surface q and VPD"),
 ("ERA5 single levels",ROOT/"output_attribution_minimal/ERA5_single_levels_T2m_tp_ssrd_1deg_1979_2024_FIXED_FULL.rds","tp;t2m;ssrd","1979-2024","1 degree","hydroclimate controls"),
 ("ERA5-Land soil water",ROOT/"soil water.nc","swvl1-4","1979-2024","0.1/0.25 degree local file","ecology controls"),
 ("GPP",ROOT/"output 1deg v3/step5_greening/step5_fixest.RData","panel_gpp","2000-03 to 2024-12","1 degree panel","ecology"),
 ("LAI",ROOT/"output 1deg v3/step6_lai/step6_fixest.RData","panel_lai","2000-01 to 2020-12","1 degree panel","ecology"),
 ("Land cover",ROOT/"output_attribution_minimal/FULL_RERUN_input_inventory_rescue_v2/landcover_aligned_cell_summary_v2.csv","IGBP modal/stability masks","2001-2024","1 degree","categorical domain"),
 ("CO2",ROOT/"co2_mm_gl.csv","monthly mean CO2","analysis overlap","global monthly","ecology control"),
 ("ONI",ROOT/"WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard_FINAL/11_ENSO_ONI_STANDARD_UPDATE/00_external_data/NOAA_CPC_ONI_historical_clean.csv","ONI","1979-2024 overlap","monthly/overlapping seasons","ENSO event definition"),
 ("Nino3.4",ROOT/"WGPE_FINAL_MANUSCRIPT_PACKAGE_v6_ONIstandard_FINAL/11_ENSO_ONI_STANDARD_UPDATE/00_external_data/NOAA_PSL_Nino34_monthly_clean.csv","Nino3.4","1979-2024 overlap","monthly","continuous ENSO"),
]: add(cat,p,var,per,res,role)
with (OUT/"INPUT_INVENTORY.csv").open("w",newline="",encoding="utf-8-sig") as f:
    w=csv.DictWriter(f,fieldnames=inputs[0].keys());w.writeheader();w.writerows(inputs)

missing=[r["absolute_path"] for r in inputs if not r["exists"]]
text=f"""# TOP100 project audit

- Formal definition: local effective lower boundary to 100 hPa.
- Convergence reference: 200 hPa.
- Truncation sensitivity: 300 hPa.
- Legacy files modified: **NO**.
- Frozen legacy paths: {len(rows)}; existing: {sum(bool(x['exists']) for x in rows)}.
- Registered input paths: {len(inputs)}; missing: {len(missing)}.
- SHA-256 ledger: `OLD_RESULTS_SHA256.csv`.
- Input inventory: `INPUT_INVENTORY.csv`.

## Missing registered inputs
{os.linesep.join('- '+x for x in missing) if missing else '- None'}

The deliberately absent `Figure1_TOP100_FINAL.png` entry records that the legacy package predates the TOP100 switch; it is not an input failure.
"""
(OUT/"PROJECT_AUDIT.md").write_text(text,encoding="utf-8")
print(f"Wrote TOP100 Stage 0 audit to {OUT}")

