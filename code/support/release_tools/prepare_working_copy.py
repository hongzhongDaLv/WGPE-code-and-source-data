#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, shutil
from pathlib import Path

REQUIRED = [
 "project_root","raw_data_root","derived_data_root","output_root","temporary_root",
 "era5_root","era5_land_root","jra55_root","rharm_root","gosif_gpp_root",
 "gimms_lai_root","mcd12c1_root","noaa_indices_root","co2_root",
 "era5_pressure_level_dir","era5_surface_pressure_file","era5_surface_geopotential_file",
 "era5_dewpoint_file","era5_single_levels_file","era5_land_soil_water_file",
 "gpp_array_file","gpp_panel_file","lai_panel_file","landcover_summary_file",
 "co2_file","oni_file","nino34_file","oni_winter_inventory_file"
]
PATH_FIELDS = [x for x in REQUIRED if x not in {"derived_data_root","output_root","temporary_root"}]

def load_simple_yaml(path):
    data={}
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line=raw.strip()
        if not line or line.startswith("#") or ":" not in line: continue
        k,v=line.split(":",1); data[k.strip()]=v.strip().strip('"').strip("'")
    return data

def resolve(value, root):
    p=Path(value).expanduser()
    return p if p.is_absolute() else (root/p).resolve()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); ap.add_argument("--prepare",action="store_true"); a=ap.parse_args()
    cfg=load_simple_yaml(a.config); missing_fields=[k for k in REQUIRED if not cfg.get(k)]
    if missing_fields: raise SystemExit("Missing config fields: "+", ".join(missing_fields))
    repo=Path(__file__).resolve().parents[3]
    missing=[]
    for k in PATH_FIELDS:
        v=cfg[k]
        if v.startswith("<") or not resolve(v, repo).exists(): missing.append(f"{k}: {v}")
    print("Configuration fields: PASS")
    print("External inputs present:", "YES" if not missing else "NO")
    for m in missing: print("MISSING",m)
    if not a.prepare: return 0 if not missing else 3
    if missing: raise SystemExit("Refusing to prepare an incomplete full-analysis run.")
    work=resolve(cfg["project_root"],repo)
    work.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader((repo/"config"/"SCRIPT_COMPATIBILITY_MAP.csv").open(encoding="utf-8-sig")))
    tokens={"__WGPE_PROJECT_ROOT__":str(work).replace("\\","/"),"__WGPE_PYTHON__":cfg.get("python_executable","python"),"__WGPE_RSCRIPT__":cfg.get("rscript_executable","Rscript"),"__WGPE_NODE__":cfg.get("node_executable","node"),"__WGPE_PYTHON_SITE_PACKAGES__":cfg.get("python_site_packages",""),"__WGPE_R_LIBRARY__":cfg.get("r_library",""),"__WGPE_RTOOLS__":cfg.get("rtools_root","")}
    for row in rows:
        src=repo/row["public_script"]; dst=work/row["compatibility_relpath"]
        dst.parent.mkdir(parents=True,exist_ok=True); text=src.read_text(encoding="utf-8")
        for old,new in tokens.items(): text=text.replace(old,new)
        dst.write_text(text,encoding="utf-8",newline="\n")
    inputs={
      "ERA5_pressure_level_dir":str(resolve(cfg["era5_pressure_level_dir"],repo)),
      "ERA5_surface_pressure":str(resolve(cfg["era5_surface_pressure_file"],repo)),
      "ERA5_surface_geopotential":str(resolve(cfg["era5_surface_geopotential_file"],repo)),
      "ERA5_dewpoint_2m":str(resolve(cfg["era5_dewpoint_file"],repo)),
      "ERA5_single_levels":str(resolve(cfg["era5_single_levels_file"],repo)),
      "ERA5_land_soil_water":str(resolve(cfg["era5_land_soil_water_file"],repo)),
      "GPP_array":str(resolve(cfg["gpp_array_file"],repo)),
      "GPP_panel":str(resolve(cfg["gpp_panel_file"],repo)),
      "LAI_panel":str(resolve(cfg["lai_panel_file"],repo)),
      "landcover_summary":str(resolve(cfg["landcover_summary_file"],repo)),
      "CO2":str(resolve(cfg["co2_file"],repo)),
      "ONI":str(resolve(cfg["oni_file"],repo)),
      "Nino34":str(resolve(cfg["nino34_file"],repo)),
      "ONI_winter_inventory":str(resolve(cfg["oni_winter_inventory_file"],repo))
    }
    analysis={
      "project_root":str(work),"scripts_root":str(work/"scripts_TOP100_revision"),
      "output_root":str(resolve(cfg["output_root"],repo)),
      "main_pressure_top_hPa":100,"convergence_pressure_top_hPa":200,
      "truncation_sensitivity_top_hPa":300,
      "effective_lower_boundary":"min(local monthly surface pressure, 1000 hPa)",
      "analysis_period":["1979-01","2024-12"],
      "native_pressure_levels_hPa":[1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100],
      "standard_layers_hPa":[[1000,850],[850,700],[700,500],[500,300],[300,200],[200,100]],
      "latitude_profile_bin_degrees":5,
      "climate_bands":{"S high latitude":[-90,-66.5],"S temperate":[-66.5,-23.5],"Tropics":[-23.5,23.5],"N temperate":[23.5,66.5],"N high latitude":[66.5,90]},
      "g_m_s2":9.80665,"bootstrap_replicates":1000,"time_block_months":60,
      "spatial_block_degrees":10,"fdr_alpha":0.05,"field_significance_permutations":1000,
      "circular_shift_min_months":24,"cv_folds":5,"cv_gap_months":[0,3],
      "random_seed":20260713,"figure_dpi":600,"figure_main_width_mm":180,
      "figure_main_height_mm":220,"overwrite_policy":"Write only under the configured working and output roots.",
      "inputs":inputs
    }
    cp=work/"scripts_TOP100_revision"/"00_config"/"config_TOP100_revision.json"
    cp.parent.mkdir(parents=True,exist_ok=True); cp.write_text(json.dumps(analysis,indent=2),encoding="utf-8",newline="\n")
    print("Prepared compatibility tree:",work)
    return 0
if __name__=="__main__": raise SystemExit(main())
