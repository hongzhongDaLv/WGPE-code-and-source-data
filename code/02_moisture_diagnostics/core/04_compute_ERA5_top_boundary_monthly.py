from __future__ import annotations

import glob
import json
import os
import time
from datetime import datetime
from pathlib import Path

import netCDF4
import numpy as np


ROOT = Path(r"__WGPE_PROJECT_ROOT__")
PLEV = ROOT / "era5_plev_ascii"
OUT = ROOT / "output_reviewer_P0_resolution"
SPFILE = OUT / "ERA5_surface_pressure_aligned_0_359.nc"
G = 9.80665
TOPS = (300, 200, 100)


def decode_dates(ds, var="valid_time"):
    t = ds.variables[var]
    return [datetime(x.year, x.month, 1) for x in netCDF4.num2date(t[:], t.units, only_use_cftime_datetimes=False)]


def create_output(path: Path, lon, lat, time_days, top):
    ds = netCDF4.Dataset(path, "w", format="NETCDF4")
    ds.createDimension("time", len(time_days))
    ds.createDimension("latitude", len(lat))
    ds.createDimension("longitude", len(lon))
    tv = ds.createVariable("time", "f8", ("time",))
    tv.units = "days since 1970-01-01"; tv.calendar = "standard"; tv[:] = time_days
    lav = ds.createVariable("latitude", "f4", ("latitude",)); lav.units = "degrees_north"; lav[:] = lat
    lov = ds.createVariable("longitude", "f4", ("longitude",)); lov.units = "degrees_east"; lov[:] = lon
    # Keep the formal month fields uncompressed. On this Windows/OneDrive
    # workspace HDF5 compression was the dominant cost and did not affect the
    # science. Fixed-size uncompressed arrays are fast and fully auditable.
    kwargs = dict(zlib=False, contiguous=True, fill_value=np.float32(-9999.0))
    iwv = ds.createVariable("iwv", "f4", ("time", "latitude", "longitude"), **kwargs)
    zbar = ds.createVariable("zbar", "f4", ("time", "latitude", "longitude"), **kwargs)
    wgpe = ds.createVariable("wgpe", "f4", ("time", "latitude", "longitude"), **kwargs)
    iwv.units = "kg m-2"; zbar.units = "m"; wgpe.units = "J m-2"
    iwv.long_name = f"Integrated water vapour from local effective lower boundary to {top} hPa"
    zbar.long_name = f"Water-vapour mass-weighted geopotential height, top boundary {top} hPa"
    wgpe.long_name = f"Water-vapour gravitational potential energy, top boundary {top} hPa"
    ds.analysis_period = "1979-01 to 2024-12"
    ds.top_boundary_hPa = top
    ds.lower_boundary = "min(monthly local surface pressure, 1000 hPa)"
    ds.crossing_interval = "pressure-linear interpolation of q and geopotential at local surface pressure"
    ds.integration = "trapezoidal in pressure coordinates"
    ds.geopotential_conversion = f"geopotential divided by g={G} m s-2"
    ds.below_ground_rule = "below-ground pressure-level values excluded except endpoint used by archived crossing-interval interpolation"
    ds.result_id_prefix = f"P0_TOP{top}"
    ds.source_pressure_levels = str(PLEV)
    ds.source_surface_pressure = str(SPFILE)
    return ds


# Build a complete month -> source-file/local-index lookup.
lookup = {}
level_reference = None
lat = lon = None
for f in sorted(PLEV.glob("*.nc")):
    ds = netCDF4.Dataset(f)
    dates = decode_dates(ds)
    levels = np.asarray(ds.variables["pressure_level"][:], dtype=float)
    if level_reference is None:
        level_reference = levels
        lat = np.asarray(ds.variables["latitude"][:], dtype=float)
        lon = np.asarray(ds.variables["longitude"][:], dtype=float)
    else:
        if not np.array_equal(level_reference, levels):
            raise RuntimeError(f"Pressure levels differ in {f}")
    for i, d in enumerate(dates):
        if d in lookup:
            raise RuntimeError(f"Duplicate month {d}")
        lookup[d] = (f, i)
    ds.close()

dates = [datetime(y, m, 1) for y in range(1979, 2025) for m in range(1, 13)]
missing = [d for d in dates if d not in lookup]
if missing:
    raise RuntimeError(f"Missing source months: {missing}")

levels = level_reference
use_idx = np.where((levels <= 1000) & (levels >= 100))[0]
plev = levels[use_idx]
if not all(t in plev for t in TOPS):
    raise RuntimeError(f"Required top boundaries absent: {TOPS}; actual={plev}")

spds = netCDF4.Dataset(SPFILE)
spdates = decode_dates(spds, "time")
sp_lookup = {d: i for i, d in enumerate(spdates)}
if any(d not in sp_lookup for d in dates):
    raise RuntimeError("Surface-pressure months incomplete")
if spds.variables["sp"].shape != (552, 181, 360):
    raise RuntimeError(f"Unexpected sp shape {spds.variables['sp'].shape}")

epoch = datetime(1970, 1, 1)
time_days = np.asarray([(d - epoch).days for d in dates], dtype=float)
outs = {top: create_output(OUT / f"ERA5_TOP{top}_monthly.nc", lon, lat, time_days, top) for top in TOPS}

t0 = time.time()
valid_counts = {t: 0 for t in TOPS}
try:
    # Source files are internally chunked in blocks of 40 months. Reading the
    # same block once (instead of once per month) avoids repeated decompression.
    source_info = []
    for path in sorted(PLEV.glob("*.nc")):
        ds0 = netCDF4.Dataset(path)
        dd = decode_dates(ds0)
        ds0.close()
        source_info.append((dd[0], path, dd))
    source_info.sort()

    for _, path, file_dates in source_info:
        src = netCDF4.Dataset(path)
        for local_start in range(0, len(file_dates), 40):
            local_stop = min(local_start + 40, len(file_dates))
            block_dates = file_dates[local_start:local_stop]
            global_start = dates.index(block_dates[0])
            global_stop = global_start + len(block_dates)

            q = np.asarray(src.variables["q"][local_start:local_stop, use_idx, :, :], dtype=np.float32)
            zh = np.asarray(src.variables["z"][local_start:local_stop, use_idx, :, :], dtype=np.float32) / np.float32(G)
            # Match the archived corrected-core validity convention exactly.
            q[(~np.isfinite(q)) | (q < 0)] = np.nan
            zh[(~np.isfinite(zh)) | (zh < 0)] = np.nan
            sp = np.asarray(spds.variables["sp"][global_start:global_stop, :, :], dtype=np.float32)
            mass = np.zeros(sp.shape, dtype=np.float64)
            moment = np.zeros(sp.shape, dtype=np.float64)
            valid = np.isfinite(sp)
            block_out = {
                top: {
                    "iwv": np.full(sp.shape, np.nan, dtype=np.float32),
                    "zbar": np.full(sp.shape, np.nan, dtype=np.float32),
                    "wgpe": np.full(sp.shape, np.nan, dtype=np.float32),
                } for top in TOPS
            }

            for k in range(len(plev) - 1):
                p0, p1 = float(plev[k]), float(plev[k + 1])
                q0, q1 = q[:, k], q[:, k + 1]
                z0, z1 = zh[:, k], zh[:, k + 1]
                endpoints_ok = np.isfinite(q0) & np.isfinite(q1) & np.isfinite(z0) & np.isfinite(z1)
                full = sp >= p0 * 100.0
                partial = (sp < p0 * 100.0) & (sp > p1 * 100.0)
                active = (full | partial) & endpoints_ok & valid
                start_hpa = np.where(full, p0, sp / 100.0)
                alpha = np.clip((p0 - start_hpa) / (p0 - p1), 0.0, 1.0)
                qs = q0 + (q1 - q0) * alpha
                zs = z0 + (z1 - z0) * alpha
                dp_pa = np.maximum(start_hpa - p1, 0.0) * 100.0
                mass += np.where(active, 0.5 * (qs + q1) * dp_pa / G, 0.0)
                moment += np.where(active, 0.5 * (qs * zs + q1 * z1) * dp_pa / G, 0.0)

                if int(p1) in TOPS:
                    top = int(p1)
                    good = valid & (mass > 0) & np.isfinite(mass) & np.isfinite(moment)
                    block_out[top]["iwv"][:] = np.where(good, mass, np.nan)
                    block_out[top]["zbar"][:] = np.where(good, moment / mass, np.nan)
                    block_out[top]["wgpe"][:] = np.where(good, G * moment, np.nan)
                    valid_counts[top] += int(good.sum())

            for top in TOPS:
                for name in ("iwv", "zbar", "wgpe"):
                    outs[top].variables[name][global_start:global_stop] = block_out[top][name]
            elapsed = time.time() - t0
            print(f"{global_stop:03d}/552 {block_dates[-1]:%Y-%m} elapsed={elapsed:.1f}s", flush=True)
            del q, zh, sp, mass, moment, block_out
        src.close()
finally:
    spds.close()
    for ds in outs.values():
        ds.close()

manifest = {
    "script": str(Path(__file__).resolve()),
    "created": datetime.now().isoformat(timespec="seconds"),
    "tops_hPa": TOPS,
    "months": len(dates),
    "pressure_levels_used_hPa": plev.tolist(),
    "valid_cell_month_counts": valid_counts,
    "g_m_s2": G,
    "elapsed_seconds": time.time() - t0,
}
(OUT / "TOP_BOUNDARY_MONTHLY_BUILD_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print(json.dumps(manifest, indent=2), flush=True)
