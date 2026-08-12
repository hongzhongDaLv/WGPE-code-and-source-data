from __future__ import annotations

import csv
from datetime import datetime
from pathlib import Path

import netCDF4
import numpy as np


ROOT = Path(r"__WGPE_PROJECT_ROOT__")
OUT = ROOT / "output_reviewer_P0_resolution"
TOPS = [300, 200, 100]
VARS = ["iwv", "zbar", "wgpe"]
UNITS = {"iwv": "kg m-2", "zbar": "m", "wgpe": "J m-2"}
TREND_UNITS = {"iwv": "kg m-2 yr-1", "zbar": "m yr-1", "wgpe": "J m-2 yr-1"}
ZONE_SPECS = [
    ("Global", -90, 90),
    ("S high latitude", -90, -66.5),
    ("S temperate", -66.5, -23.5),
    ("Tropics", -23.5, 23.5),
    ("N temperate", 23.5, 66.5),
    ("N high latitude", 66.5, 90),
]


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


def decimal_years():
    out = []
    for y in range(1979, 2025):
        y0 = datetime(y, 1, 1); y1 = datetime(y + 1, 1, 1)
        for m in range(1, 13):
            d = datetime(y, m, 1)
            out.append(y + (d - y0).days / (y1 - y0).days)
    return np.asarray(out, dtype=np.float64)


def slope_field(a, x):
    # a: time, lat, lon. Full 552-month coverage is required by Stage 1 gate.
    clim = np.stack([np.nanmean(a[m::12], axis=0) for m in range(12)])
    anom = a - clim[np.arange(a.shape[0]) % 12]
    xc = x - x.mean()
    slope = np.nansum(anom * xc[:, None, None], axis=0) / np.sum(xc * xc)
    raw = np.nansum((a - np.nanmean(a, axis=0)) * xc[:, None, None], axis=0) / np.sum(xc * xc)
    return clim.mean(axis=0), slope, raw, anom


def weighted_mean(field, weight, mask=None):
    ok = np.isfinite(field)
    if mask is not None:
        ok &= mask
    ww = np.where(ok, weight, 0.0)
    return float(np.nansum(field * ww) / np.sum(ww))


def weighted_corr(x, y, w):
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    xx, yy, ww = x[ok], y[ok], w[ok]
    ww = ww / ww.sum()
    mx = np.sum(ww * xx); my = np.sum(ww * yy)
    return float(np.sum(ww * (xx - mx) * (yy - my)) /
                 np.sqrt(np.sum(ww * (xx - mx) ** 2) * np.sum(ww * (yy - my) ** 2)))


def unweighted_corr(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    return float(np.corrcoef(x[ok], y[ok])[0, 1])


x = decimal_years()
sample = netCDF4.Dataset(OUT / "ERA5_TOP300_monthly.nc")
lat = np.asarray(sample.variables["latitude"][:], dtype=float)
lon = np.asarray(sample.variables["longitude"][:], dtype=float)
sample.close()
w2 = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], (len(lat), len(lon))).copy()

means = {}; trends = {}; raws = {}; global_series = {}
global_rows = []; zone_rows = []; lat_rows = []
rid = 1

trend_nc = netCDF4.Dataset(OUT / "TOP_BOUNDARY_GRIDCELL_FIELDS.nc", "w", format="NETCDF4")
trend_nc.createDimension("top", 3); trend_nc.createDimension("latitude", len(lat)); trend_nc.createDimension("longitude", len(lon))
trend_nc.createVariable("top_hPa", "i4", ("top",))[:] = TOPS
trend_nc.createVariable("latitude", "f4", ("latitude",))[:] = lat
trend_nc.createVariable("longitude", "f4", ("longitude",))[:] = lon
for var in VARS:
    trend_nc.createVariable(f"{var}_climatology", "f4", ("top", "latitude", "longitude"), zlib=True, complevel=2)
    trend_nc.createVariable(f"{var}_trend", "f4", ("top", "latitude", "longitude"), zlib=True, complevel=2)
    trend_nc.variables[f"{var}_trend"].units = TREND_UNITS[var]

for top_i, top in enumerate(TOPS):
    ds = netCDF4.Dataset(OUT / f"ERA5_TOP{top}_monthly.nc")
    for var in VARS:
        a = np.asarray(ds.variables[var][:], dtype=np.float32)
        clim, tr, raw, anom = slope_field(a, x)
        means[(top, var)] = clim
        trends[(top, var)] = tr
        raws[(top, var)] = raw
        trend_nc.variables[f"{var}_climatology"][top_i] = clim.astype(np.float32)
        trend_nc.variables[f"{var}_trend"][top_i] = tr.astype(np.float32)

        denom = np.sum(np.isfinite(a) * w2[None, :, :], axis=(1, 2))
        series = np.nansum(a * w2[None, :, :], axis=(1, 2)) / denom
        series_anom = series - np.asarray([series[m::12].mean() for m in range(12)])[np.arange(552) % 12]
        global_series[(top, var)] = (series, series_anom)
        global_rows.append({
            "result_id": f"P0TB{rid:05d}", "top_hPa": top, "variable": var,
            "period": "1979-01 to 2024-12", "n_grid": int(np.isfinite(tr).sum()),
            "global_climatological_mean": weighted_mean(clim, w2), "mean_unit": UNITS[var],
            "global_trend_deseasonalized": weighted_mean(tr, w2), "global_trend_raw": weighted_mean(raw, w2),
            "trend_unit": TREND_UNITS[var], "height_reference": "absolute geopotential height",
            "mask": "common all-valid TOP100/TOP200/TOP300 grid", "area_weighting": "cos(latitude)",
            "trend_method": "calendar-month climatology removed; grid-cell OLS first; area mean second",
            "input_file": str(OUT / f"ERA5_TOP{top}_monthly.nc"),
            "script": str(Path(__file__).resolve()),
        }); rid += 1

        for name, lo, hi in ZONE_SPECS:
            lmask = (lat >= lo) & (lat <= hi)
            if name != "Global" and hi != 90:
                lmask &= lat < hi
            m2 = np.broadcast_to(lmask[:, None], tr.shape)
            zone_rows.append({
                "result_id": f"P0TB{rid:05d}", "top_hPa": top, "variable": var,
                "region": name, "lat_min": lo, "lat_max": hi,
                "n_grid": int(np.sum(m2 & np.isfinite(tr))),
                "climatological_mean": weighted_mean(clim, w2, m2), "mean_unit": UNITS[var],
                "trend": weighted_mean(tr, w2, m2), "trend_unit": TREND_UNITS[var],
                "area_weighting": "cos(latitude)", "mask": "common all-valid top-boundary mask",
            }); rid += 1

        for lo in range(-90, 90, 5):
            hi = lo + 5
            lmask = (lat >= lo) & ((lat < hi) if hi < 90 else (lat <= hi))
            m2 = np.broadcast_to(lmask[:, None], tr.shape)
            lat_rows.append({
                "result_id": f"P0TB{rid:05d}", "top_hPa": top, "variable": var,
                "latitude_bin_start": lo, "latitude_bin_end": hi, "latitude_bin_mid": lo + 2.5,
                "n_grid": int(np.sum(m2 & np.isfinite(tr))),
                "trend": weighted_mean(tr, w2, m2), "trend_unit": TREND_UNITS[var],
                "climatological_mean": weighted_mean(clim, w2, m2), "mean_unit": UNITS[var],
                "area_weighting": "cos(latitude)",
            }); rid += 1
        del a, anom
    ds.close()

trend_nc.source_script = str(Path(__file__).resolve())
trend_nc.common_mask = "finite across TOP300, TOP200, TOP100; all 65160 grid cells in monthly fields"
trend_nc.close()

# Month-level area means for HAC and paired time-block bootstrap.
series_rows = []
dates = [f"{y:04d}-{m:02d}-01" for y in range(1979, 2025) for m in range(1, 13)]
for top in TOPS:
    for var in VARS:
        raw_s, an_s = global_series[(top, var)]
        for i, d in enumerate(dates):
            series_rows.append({"date": d, "decimal_year": x[i], "top_hPa": top, "variable": var,
                                "area_mean_raw": raw_s[i], "area_mean_deseasonalized": an_s[i], "unit": UNITS[var]})

# Pairwise fair comparisons.
pair_rows = []; pattern_rows = []
pairs = [(300, 100), (200, 100), (300, 200)]
for a_top, b_top in pairs:
    for var in VARS:
        ta, tb = trends[(a_top, var)], trends[(b_top, var)]
        ma, mb = means[(a_top, var)], means[(b_top, var)]
        common = np.isfinite(ta) & np.isfinite(tb)
        ga, gb = weighted_mean(ta, w2, common), weighted_mean(tb, w2, common)
        ca, cb = weighted_mean(ma, w2, common), weighted_mean(mb, w2, common)
        pair_rows.append({
            "result_id": f"P0TB{rid:05d}", "comparison": f"TOP{a_top}-TOP{b_top}", "variable": var,
            "n_grid_common": int(common.sum()), "mean_difference": ca - cb, "mean_unit": UNITS[var],
            "trend_difference": ga - gb, "trend_unit": TREND_UNITS[var],
            "relative_trend_difference_pct_of_reference": 100 * (ga - gb) / gb if gb != 0 else np.nan,
            "reference": f"TOP{b_top}", "mask": "pairwise common finite; identical to TOP100 common all-valid mask",
        }); rid += 1

        sign_agree = np.mean(np.sign(ta[common]) == np.sign(tb[common])) * 100
        profile_a = np.asarray([r["trend"] for r in lat_rows if r["top_hPa"] == a_top and r["variable"] == var])
        profile_b = np.asarray([r["trend"] for r in lat_rows if r["top_hPa"] == b_top and r["variable"] == var])
        pattern_rows.append({
            "result_id": f"P0TB{rid:05d}", "comparison": f"TOP{a_top}-TOP{b_top}", "variable": var,
            "n_grid_common": int(common.sum()), "trend_map_r_equal_grid": unweighted_corr(ta, tb),
            "trend_map_r_coslat": weighted_corr(ta, tb, w2),
            "trend_sign_agreement_pct": sign_agree,
            "latitude_profile_r": float(np.corrcoef(profile_a, profile_b)[0, 1]),
            "significant_trend_sign_agreement_pct": "PENDING_EFFECTIVE_N_FDR",
            "script": str(Path(__file__).resolve()),
        }); rid += 1

# Omitted fractions relative to TOP100.
omitted_rows = []
for top in [300, 200]:
    for var in ["iwv", "wgpe"]:
        mean_full = weighted_mean(means[(100, var)], w2)
        mean_top = weighted_mean(means[(top, var)], w2)
        trend_full = weighted_mean(trends[(100, var)], w2)
        trend_top = weighted_mean(trends[(top, var)], w2)
        omitted_rows.append({
            "result_id": f"P0TB{rid:05d}", "top_hPa": top, "reference_top_hPa": 100,
            "variable": var, "omitted_mean_fraction_pct": 100 * (mean_full - mean_top) / mean_full,
            "omitted_trend_fraction_pct": 100 * (trend_full - trend_top) / trend_full,
            "note": "For W-GPE this is identical to omitted geopotential-moment fraction because W-GPE=g times moment.",
        }); rid += 1

write_csv(OUT / "TOP_BOUNDARY_GLOBAL_RESULTS.csv", global_rows)
write_csv(OUT / "TOP_BOUNDARY_ZONE_RESULTS.csv", zone_rows)
write_csv(OUT / "TOP_BOUNDARY_LAT5_PROFILES.csv", lat_rows)
write_csv(OUT / "TOP_BOUNDARY_DIFFERENCES.csv", pair_rows)
write_csv(OUT / "TOP_BOUNDARY_PATTERN_METRICS.csv", pattern_rows)
write_csv(OUT / "TOP_BOUNDARY_OMITTED_FRACTIONS.csv", omitted_rows)
write_csv(OUT / "TOP_BOUNDARY_GLOBAL_MONTHLY_SERIES.csv", series_rows)

print("global rows")
for r in global_rows:
    print(r["top_hPa"], r["variable"], r["global_trend_deseasonalized"])
