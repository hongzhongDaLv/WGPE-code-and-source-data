import csv
from datetime import datetime
from pathlib import Path

import netCDF4
import numpy as np

ROOT = Path(r"__WGPE_PROJECT_ROOT__")
OUT = ROOT / "output_reviewer_P0_resolution"
TOPS = [300, 200, 100]
VARS = ["iwv", "zbar", "wgpe"]
ZONES = [
    ("Global", -90, 90),
    ("S high latitude", -90, -66.5),
    ("S temperate", -66.5, -23.5),
    ("Tropics", -23.5, 23.5),
    ("N temperate", 23.5, 66.5),
    ("N high latitude", 66.5, 90),
]
UNITS = {"iwv":"kg m-2", "zbar":"m", "wgpe":"J m-2"}
dates = [datetime(y,m,1) for y in range(1979,2025) for m in range(1,13)]
rows = []

for top in TOPS:
    ds = netCDF4.Dataset(OUT / f"ERA5_TOP{top}_monthly.nc")
    lat = np.asarray(ds.variables["latitude"][:], float)
    wlat = np.cos(np.deg2rad(lat))
    for var in VARS:
        a = np.asarray(ds.variables[var][:], np.float32)
        for region, lo, hi in ZONES:
            lm = (lat >= lo) & (lat <= hi)
            if region != "Global" and hi != 90:
                lm &= lat < hi
            aa = a[:, lm, :]
            ww = wlat[lm][None, :, None]
            s = np.nansum(aa * ww, axis=(1,2)) / np.sum(np.isfinite(aa) * ww, axis=(1,2))
            clim = np.asarray([s[m::12].mean() for m in range(12)])
            an = s - clim[np.arange(552)%12]
            for i,d in enumerate(dates):
                rows.append({"date":d.strftime("%Y-%m-%d"),"top_hPa":top,"variable":var,"region":region,
                             "area_mean_raw":float(s[i]),"area_mean_deseasonalized":float(an[i]),"unit":UNITS[var]})
    ds.close()

with (OUT / "TOP_BOUNDARY_ZONE_MONTHLY_SERIES.csv").open("w", newline="", encoding="utf-8-sig") as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
print(len(rows))
