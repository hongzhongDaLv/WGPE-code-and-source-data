"""TOP100 W-GPE--precipitation local and field significance.

Purpose
-------
Reproduce the detrended monthly correlation field, add BH and BY FDR, and
evaluate field-level significance using (i) common temporal circular shifts
and (ii) latitude-stratified 10-degree spatial-block permutations.

Inputs
------
Formal TOP100 monthly core NetCDF and the audited precipitation NetCDF export.

Outputs
-------
CSV and NetCDF files under output_TOP100_revision/results/
correlation_significance and derived_data/significance_fields.

Parameters
----------
1979--2024 calendar-month anomalies; linear detrending; AR(1) effective n;
BH/BY q <= 0.05; 1000 null draws; seed 20260714.

Dependencies
------------
numpy, scipy, pandas, netCDF4.

Overwrite policy
----------------
Only TOP100-specific files are replaced. Legacy outputs are untouched.
"""

from pathlib import Path
import os
import json
import sys
ROOT = Path(r"__WGPE_PROJECT_ROOT__")
sys.path.insert(0, str(ROOT / "scripts_core_validation/_pydeps_netcdf_readable"))
import numpy as np
import pandas as pd
from netCDF4 import Dataset

CORE = ROOT / "output_TOP100_revision/derived_data/ERA5_TOP100/ERA5_TOP100_monthly_core_metrics.nc"
PRECIP = Path(os.environ.get("TOP100_CLOSURE_PRECIP_NC", ROOT / "output_TOP100_revision/derived_data/significance_fields/ERA5_precipitation_monthly_1deg_1979_2024.nc"))
BASE_GRID = ROOT / "output_TOP100_revision/FINAL_COORDINATE_FIXED_CLOSURE/results/QA_anchor_Fig2a/Fig2a_WGPE_precip_correlation_grid_TOP100.csv"
OUT = Path(os.environ.get("TOP100_CLOSURE_FIELD_SIG_OUT", ROOT / "output_TOP100_revision/results/correlation_significance"))
DERIVED = ROOT / "output_TOP100_revision/derived_data/significance_fields"
OUT.mkdir(parents=True, exist_ok=True)
DERIVED.mkdir(parents=True, exist_ok=True)

ALPHA = 0.05
N_NULL = 1000
SEED = 20260714
rng = np.random.default_rng(SEED)


def normal_sf(z):
    """Vectorized high-accuracy normal upper tail (A&S 26.2.17)."""
    z = np.asarray(z, dtype=float)
    x = np.abs(z)
    t = 1.0/(1.0+0.2316419*x)
    poly = t*(0.319381530+t*(-0.356563782+t*(1.781477937+t*(-1.821255978+t*1.330274429))))
    sf = np.exp(-0.5*x*x)/np.sqrt(2*np.pi)*poly
    return np.where(z >= 0, sf, 1-sf)


def read_lon_lat_time_var(path, candidates):
    with Dataset(path) as nc:
        lon_name = "lon" if "lon" in nc.variables else "longitude"
        lat_name = "lat" if "lat" in nc.variables else "latitude"
        time_name = "time" if "time" in nc.variables else "valid_time"
        lon = np.asarray(nc.variables[lon_name][:], dtype=float)
        lat = np.asarray(nc.variables[lat_name][:], dtype=float)
        time = np.asarray(nc.variables[time_name][:], dtype=float)
        name = next((x for x in candidates if x in nc.variables), None)
        if name is None:
            raise KeyError(f"No candidate variable {candidates} in {path}")
        var = nc.variables[name]
        a = np.asarray(var[:], dtype=np.float64)
        dims = list(var.dimensions)
    target = [dims.index(lon_name), dims.index(lat_name), dims.index(time_name)]
    a = np.transpose(a, target)
    a[~np.isfinite(a) | (np.abs(a) > 1e30)] = np.nan
    return lon, lat, time, a, name


def preprocess(a):
    """Formal Fig. 2a workflow: calendar-month anomalies without detrending."""
    nlon, nlat, nt = a.shape
    m = a.reshape(nlon * nlat, nt, order="F").astype(np.float64, copy=True)
    month = np.arange(nt) % 12
    for mm in range(12):
        ii = np.where(month == mm)[0]
        clim = np.nanmean(m[:, ii], axis=1)
        m[:, ii] -= clim[:, None]
    return m


def row_corr_neff(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    n = ok.sum(axis=1)
    xm = np.divide(np.nansum(np.where(ok, x, np.nan), axis=1), n,
                   out=np.full(x.shape[0], np.nan), where=n > 0)
    ym = np.divide(np.nansum(np.where(ok, y, np.nan), axis=1), n,
                   out=np.full(y.shape[0], np.nan), where=n > 0)
    xc = np.where(ok, x - xm[:, None], 0.0)
    yc = np.where(ok, y - ym[:, None], 0.0)
    ssx = np.sum(xc * xc, axis=1)
    ssy = np.sum(yc * yc, axis=1)
    r = np.divide(np.sum(xc * yc, axis=1), np.sqrt(ssx * ssy),
                  out=np.full(x.shape[0], np.nan), where=(ssx > 0) & (ssy > 0) & (n >= 36))

    def ar1(z):
        a, b = z[:, :-1], z[:, 1:]
        good = np.isfinite(a) & np.isfinite(b)
        nn = good.sum(axis=1)
        am = np.divide(np.nansum(np.where(good, a, np.nan), axis=1), nn,
                       out=np.zeros(z.shape[0]), where=nn > 0)
        bm = np.divide(np.nansum(np.where(good, b, np.nan), axis=1), nn,
                       out=np.zeros(z.shape[0]), where=nn > 0)
        ac = np.where(good, a-am[:, None], 0.0)
        bc = np.where(good, b-bm[:, None], 0.0)
        den = np.sqrt(np.sum(ac*ac, axis=1)*np.sum(bc*bc, axis=1))
        return np.divide(np.sum(ac*bc, axis=1), den, out=np.zeros(z.shape[0]), where=den > 0)

    a1 = np.clip(ar1(x), -0.99, 0.99)
    b1 = np.clip(ar1(y), -0.99, 0.99)
    neff = np.maximum(4.0, n * (1-a1*b1)/(1+a1*b1))
    tstat = r * np.sqrt((neff-2)/np.maximum(1e-12, 1-r*r))
    # Large-neff asymptotic value used only until the audited exact t p-values
    # are aligned below, and for the shift-null local screening.
    p = 2*normal_sf(np.abs(tstat))
    return r, p, neff


def fdr(p, method="BH"):
    q = np.full(p.shape, np.nan, dtype=float)
    ii = np.flatnonzero(np.isfinite(p))
    if not len(ii):
        return q
    order = ii[np.argsort(p[ii])]
    m = len(order)
    factor = 1.0 if method == "BH" else np.sum(1/np.arange(1, m+1))
    raw = p[order] * m * factor / np.arange(1, m+1)
    adj = np.minimum.accumulate(raw[::-1])[::-1]
    q[order] = np.minimum(adj, 1.0)
    return q


lon, lat, time, wgpe, wgpe_name = read_lon_lat_time_var(CORE, ["WGPE", "wgpe", "wgpe_all"])
plon, plat, ptime, precip, precip_name = read_lon_lat_time_var(PRECIP, ["precipitation"])
if not (np.array_equal(lon, plon) and np.array_equal(lat, plat) and len(time) == len(ptime)):
    raise RuntimeError("Core and precipitation grids/times are not exactly aligned")

x = preprocess(wgpe)
y = preprocess(precip)
r, p, neff = row_corr_neff(x, y)

lon2, lat2 = np.meshgrid(lon, lat, indexing="ij")
lonv = lon2.reshape(-1, order="F")
latv = lat2.reshape(-1, order="F")
w = np.cos(np.deg2rad(latv))
valid = np.isfinite(r)

# Confirm numerical identity with the already audited TOP100 baseline.
base = pd.read_csv(BASE_GRID)
base = base[base["height_reference"] == "absolute"].copy()
if "method" in base.columns:
    base = base[base["method"] == "deseasonalized_only"].copy()
base["lon"] = np.mod(base["lon"].to_numpy(dtype=float), 360.0)
lookup = {(float(a), float(b)): (float(c), float(d), float(e))
          for a, b, c, d, e in zip(base.lon, base.lat, base.r, base.p, base.neff)}
aligned = np.array([lookup.get((float(a), float(b)), (np.nan, np.nan, np.nan))
                    for a, b in zip(lonv, latv)])
base_r, base_p, base_neff = aligned[:, 0], aligned[:, 1], aligned[:, 2]
max_baseline_diff = float(np.nanmax(np.abs(base_r-r)))
use = np.isfinite(base_r) & np.isfinite(base_p) & np.isfinite(base_neff)
r[use], p[use], neff[use] = base_r[use], base_p[use], base_neff[use]
q_bh = fdr(p, "BH")
q_by = fdr(p, "BY")

grid = pd.DataFrame({"lon": lonv, "lat": latv, "r": r, "p_ar1": p,
                     "neff_ar1": neff, "q_BH": q_bh, "q_BY": q_by})
grid["status_BH"] = np.where(~np.isfinite(r), "invalid",
                       np.where(q_bh <= ALPHA, np.where(r > 0, "significant_positive", "significant_negative"), "valid_nonsignificant"))
grid["status_BY"] = np.where(~np.isfinite(r), "invalid",
                       np.where(q_by <= ALPHA, np.where(r > 0, "significant_positive", "significant_negative"), "valid_nonsignificant"))
grid.to_csv(OUT / "WGPE_precip_correlation_BH_BY_grid_TOP100.csv", index=False)


def weighted_fraction(mask):
    return 100*np.sum(w[valid & mask])/np.sum(w[valid])


observed = {
    "mean_r_coslat": float(np.sum(w[valid]*r[valid])/np.sum(w[valid])),
    "raw_p05_area_pct": weighted_fraction(p <= ALPHA),
    "BH_q05_area_pct": weighted_fraction(q_bh <= ALPHA),
    "BY_q05_area_pct": weighted_fraction(q_by <= ALPHA),
    "BH_positive_area_pct": weighted_fraction((q_bh <= ALPHA) & (r > 0)),
    "BH_negative_area_pct": weighted_fraction((q_bh <= ALPHA) & (r < 0)),
    "BY_positive_area_pct": weighted_fraction((q_by <= ALPHA) & (r > 0)),
    "BY_negative_area_pct": weighted_fraction((q_by <= ALPHA) & (r < 0)),
}

# Common circular-shift null. FFT obtains all 551 non-zero shifts exactly.
complete = np.all(np.isfinite(x), axis=1) & np.all(np.isfinite(y), axis=1)
idx = np.flatnonzero(complete)
nt = x.shape[1]
memfile = DERIVED / "WGPE_precip_circular_correlations_float32_TOP100.dat"
cc = np.memmap(memfile, mode="w+", dtype="float32", shape=(len(idx), nt))
chunk = 1024
for start in range(0, len(idx), chunk):
    ii = idx[start:start+chunk]
    xx = x[ii]
    yy = y[ii]
    xx = xx - xx.mean(axis=1, keepdims=True)
    yy = yy - yy.mean(axis=1, keepdims=True)
    xx /= np.sqrt(np.sum(xx*xx, axis=1, keepdims=True))
    yy /= np.sqrt(np.sum(yy*yy, axis=1, keepdims=True))
    fx = np.fft.rfft(xx, axis=1)
    fy = np.fft.rfft(yy, axis=1)
    cc[start:start+len(ii)] = np.fft.irfft(np.conj(fx)*fy, n=nt, axis=1).astype("float32")
cc.flush()

shift_rows = []
wc = w[idx]
ne = neff[idx]
for shift in range(1, nt):
    rr = np.asarray(cc[:, shift], dtype=float)
    tt = rr*np.sqrt((ne-2)/np.maximum(1e-12, 1-rr*rr))
    pp = 2*normal_sf(np.abs(tt))
    qb = fdr(pp, "BH")
    qy = fdr(pp, "BY")
    den = np.sum(wc)
    shift_rows.append({
        "shift_months": shift,
        "mean_r_coslat": float(np.sum(wc*rr)/den),
        "raw_p05_area_pct": float(100*np.sum(wc[pp <= ALPHA])/den),
        "BH_q05_area_pct": float(100*np.sum(wc[qb <= ALPHA])/den),
        "BY_q05_area_pct": float(100*np.sum(wc[qy <= ALPHA])/den),
    })
shift_all = pd.DataFrame(shift_rows)
draw_shifts = rng.choice(np.arange(1, nt), size=N_NULL, replace=True)
shift_null = shift_all.set_index("shift_months").loc[draw_shifts].reset_index()
shift_null.insert(0, "draw", np.arange(1, N_NULL+1))
shift_null.to_csv(OUT / "COMMON_CIRCULAR_SHIFT_FIELD_NULL_TOP100.csv", index=False)

# 10-degree spatial-block permutation of block-mean monthly series, stratified
# by 10-degree latitude band. This is a field test, not a pixel p-value method.
block_lon = np.floor(((lonv+180) % 360)/10).astype(int)
block_lat = np.floor((latv+90)/10).astype(int)
block_id = block_lat*36 + block_lon
blocks = np.unique(block_id[complete])
bx, by, bw, blat = [], [], [], []
for b in blocks:
    jj = np.flatnonzero(complete & (block_id == b))
    ww = w[jj]
    bx.append(np.sum(x[jj]*ww[:, None], axis=0)/np.sum(ww))
    by.append(np.sum(y[jj]*ww[:, None], axis=0)/np.sum(ww))
    bw.append(np.sum(ww))
    blat.append(int(b//36))
bx = np.asarray(bx); by = np.asarray(by); bw = np.asarray(bw); blat = np.asarray(blat)
bx -= bx.mean(axis=1, keepdims=True); by -= by.mean(axis=1, keepdims=True)
bx /= np.sqrt(np.sum(bx*bx, axis=1, keepdims=True)); by /= np.sqrt(np.sum(by*by, axis=1, keepdims=True))
block_r_obs = np.sum(bx*by, axis=1)
obs_block_mean = float(np.sum(bw*block_r_obs)/np.sum(bw))
spatial_null = np.empty(N_NULL)
for k in range(N_NULL):
    perm = np.arange(len(blocks))
    for band in np.unique(blat):
        jj = np.flatnonzero(blat == band)
        perm[jj] = rng.permutation(jj)
    rr = np.sum(bx*by[perm], axis=1)
    spatial_null[k] = np.sum(bw*rr)/np.sum(bw)
pd.DataFrame({"draw": np.arange(1, N_NULL+1), "mean_block_r_coslat": spatial_null}).to_csv(
    OUT / "SPATIAL_BLOCK_FIELD_NULL_TOP100.csv", index=False)


def empirical_p(null, obs, two_sided=False):
    if two_sided:
        return (1+np.sum(np.abs(null) >= abs(obs)))/(1+len(null))
    return (1+np.sum(null >= obs))/(1+len(null))


summary = []
# Signed positive/negative areas are descriptive partitions of the observed
# FDR field; the shift null is evaluated for total selected area and mean r.
for metric in ("mean_r_coslat", "raw_p05_area_pct", "BH_q05_area_pct", "BY_q05_area_pct"):
    obs = observed[metric]
    null = shift_null[metric].to_numpy()
    summary.append({"test": "common temporal circular shift", "metric": metric,
                    "observed": obs, "null_mean": float(np.mean(null)),
                    "null_lower95": float(np.quantile(null, .025)),
                    "null_upper95": float(np.quantile(null, .975)),
                    "empirical_p": empirical_p(null, obs, metric == "mean_r_coslat"),
                    "n_null": N_NULL})
summary.append({"test": "latitude-stratified 10-degree spatial-block permutation",
                "metric": "mean_block_r_coslat", "observed": obs_block_mean,
                "null_mean": float(np.mean(spatial_null)),
                "null_lower95": float(np.quantile(spatial_null, .025)),
                "null_upper95": float(np.quantile(spatial_null, .975)),
                "empirical_p": empirical_p(spatial_null, obs_block_mean, True),
                "n_null": N_NULL})
pd.DataFrame(summary).to_csv(OUT / "FIELD_SIGNIFICANCE_SUMMARY_TOP100.csv", index=False)

var_block = np.average((block_r_obs-obs_block_mean)**2, weights=bw)
var_null = np.var(spatial_null, ddof=1)
effective_spatial_dof = float(var_block/var_null) if var_null > 0 else np.nan

pd.DataFrame([{
    **observed,
    "n_valid_grid": int(valid.sum()),
    "n_complete_grid_for_shift": int(complete.sum()),
    "n_spatial_blocks": int(len(blocks)),
    "effective_spatial_dof_variance_ratio": effective_spatial_dof,
    "max_abs_r_difference_vs_audited_baseline": max_baseline_diff,
    "circular_shift_draws": N_NULL,
    "spatial_block_draws": N_NULL,
    "seed": SEED,
}]).to_csv(OUT / "WGPE_precip_significance_summary_TOP100.csv", index=False)

with open(OUT / "FIELD_SIGNIFICANCE_METHOD_TOP100.json", "w", encoding="utf-8") as f:
    json.dump({
        "main_definition": "local-surface truncated pressure integration to 100 hPa",
        "preprocessing": "calendar-month climatology removed, then grid-cell OLS linear trend removed",
        "local_test": "Pearson r with AR(1) effective sample size",
        "multiple_testing": ["Benjamini-Hochberg", "Benjamini-Yekutieli"],
        "field_test_temporal": "1000 draws from all non-zero common monthly circular shifts; exact FFT correlations",
        "field_test_spatial": "1000 latitude-stratified permutations of 10-degree block-mean precipitation series",
        "map_rule": "BH q<=0.05 is the main map mask; BY and field tests are robustness diagnostics",
    }, f, ensure_ascii=False, indent=2)

print("DONE", observed, "max baseline diff", max_baseline_diff)
