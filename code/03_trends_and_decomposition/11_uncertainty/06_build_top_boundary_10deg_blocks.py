from pathlib import Path

import netCDF4
import numpy as np

ROOT = Path(r"__WGPE_PROJECT_ROOT__")
OUT = ROOT / "output_reviewer_P0_resolution"
TOPS = [300, 200, 100]
VARS = ["iwv", "zbar", "wgpe"]

src0 = netCDF4.Dataset(OUT / "ERA5_TOP300_monthly.nc")
lat = np.asarray(src0.variables["latitude"][:], float)
lon = np.asarray(src0.variables["longitude"][:], float)
time = np.asarray(src0.variables["time"][:], float)
time_units = src0.variables["time"].units
src0.close()

lat_edges = np.arange(-90, 91, 10)
lon_edges = np.arange(0, 361, 10)
nlatb, nlonb = len(lat_edges)-1, len(lon_edges)-1

nc = netCDF4.Dataset(OUT / "TOP_BOUNDARY_10DEG_BLOCK_MONTHLY.nc", "w", format="NETCDF4")
nc.createDimension("top", 3); nc.createDimension("variable", 3); nc.createDimension("time", 552)
nc.createDimension("lat_block", nlatb); nc.createDimension("lon_block", nlonb)
nc.createVariable("top_hPa", "i4", ("top",))[:] = TOPS
nc.createVariable("variable_code", str, ("variable",))[:] = np.asarray(VARS, dtype=object)
tv = nc.createVariable("time", "f8", ("time",)); tv.units = time_units; tv[:] = time
nc.createVariable("lat_block_mid", "f4", ("lat_block",))[:] = (lat_edges[:-1]+lat_edges[1:])/2
nc.createVariable("lon_block_mid", "f4", ("lon_block",))[:] = (lon_edges[:-1]+lon_edges[1:])/2
v = nc.createVariable("block_area_mean", "f4", ("top", "variable", "time", "lat_block", "lon_block"), zlib=True, complevel=2)
aw = nc.createVariable("block_area_weight", "f8", ("lat_block", "lon_block"))

block_weight = np.zeros((nlatb, nlonb), float)
for ib in range(nlatb):
    lm = (lat >= lat_edges[ib]) & ((lat < lat_edges[ib+1]) if ib < nlatb-1 else (lat <= lat_edges[ib+1]))
    block_weight[ib, :] = np.cos(np.deg2rad(lat[lm])).sum() * 10
aw[:] = block_weight

for it, top in enumerate(TOPS):
    ds = netCDF4.Dataset(OUT / f"ERA5_TOP{top}_monthly.nc")
    for iv, var in enumerate(VARS):
        a = np.asarray(ds.variables[var][:], np.float32)
        out = np.full((552, nlatb, nlonb), np.nan, np.float32)
        for ib in range(nlatb):
            lm = (lat >= lat_edges[ib]) & ((lat < lat_edges[ib+1]) if ib < nlatb-1 else (lat <= lat_edges[ib+1]))
            lw = np.cos(np.deg2rad(lat[lm]))[:, None]
            for jb in range(nlonb):
                xm = (lon >= lon_edges[jb]) & ((lon < lon_edges[jb+1]) if jb < nlonb-1 else (lon <= lon_edges[jb+1]))
                aa = a[:, lm, :][:, :, xm]
                den = np.sum(np.isfinite(aa) * lw[None, :, :], axis=(1,2))
                out[:, ib, jb] = np.nansum(aa * lw[None, :, :], axis=(1,2)) / den
        v[it, iv] = out
        print(top, var, flush=True)
    ds.close()

nc.block_definition = "fixed 10 degree latitude by 10 degree longitude blocks; cos(latitude)-weighted cell means"
nc.analysis_role = "spatial-block and combined spatiotemporal bootstrap"
nc.close()
