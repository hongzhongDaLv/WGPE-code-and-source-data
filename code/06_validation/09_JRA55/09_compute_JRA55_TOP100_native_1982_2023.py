# -*- coding: utf-8 -*-
"""Compute the formal JRA-55 100--1000 hPa water-vapour core.

The calculation is deliberately native-grid first.  Pressure-layer q and gh
are linearly interpolated to the local surface pressure for the one interval
crossing the surface.  Derived fields are then bilinearly mapped to the locked
ERA5 1-degree grid.  Categorical QC fields use nearest-neighbour mapping.
"""
from __future__ import annotations

import csv
import os
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(r"__WGPE_PROJECT_ROOT__")
OUT = ROOT / "output_TOP100_revision" / "results" / "JRA55"
OUT.mkdir(parents=True, exist_ok=True)
INV = ROOT / "output_core_validation" / "JRA55_GRIB_FILE_INVENTORY_V2.csv"
sys.path.insert(0, str(ROOT / "scripts_core_validation" / "_pydeps_grib_readable"))
sys.path.insert(0, str(ROOT / "scripts_core_validation" / "_pydeps_netcdf_readable"))
import eccodes as ec  # noqa: E402
import numpy as np  # noqa: E402
from netCDF4 import Dataset, date2num  # noqa: E402

G = 9.80665
LEVELS = np.array([1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100], dtype=float)
BANDS = [(1000,850),(850,700),(700,500),(500,300),(300,200),(200,100)]
DATES = [datetime(y,m,1) for y in range(1982,2024) for m in range(1,13)]
NATIVE = OUT / "JRA55_TOP100_WGPE_native_1p25deg_1982_2023.nc"
ONEDEG = OUT / "JRA55_TOP100_WGPE_1deg_1982_2023.nc"
QCNC = OUT / "JRA55_TOP100_surface_truncation_QC.nc"
LAYNC = OUT / "JRA55_TOP100_z_layer_contributions_1982_2023.nc"
QCCSV = OUT / "JRA55_TOP100_surface_truncation_summary.csv"


def inventory_files():
    out = {"q":{}, "gh":{}, "sp":{}}
    with open(INV, encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            v = r["variable_short_name"]
            if v not in out: continue
            a, b = r["first_time"][:7], r["last_time"][:7]
            if a[:4] == b[:4] and 1982 <= int(a[:4]) <= 2023:
                out[v][int(a[:4])] = Path(r["absolute_path"])
    for v in out:
        miss = sorted(set(range(1982,2024))-set(out[v]))
        if miss: raise RuntimeError(f"Missing {v} annual files: {miss}")
    return out


def read_year(path, wanted, levels=None):
    records = {}
    with open(path, "rb") as f:
        while True:
            gid = ec.codes_grib_new_from_file(f)
            if gid is None: break
            short = str(ec.codes_get(gid, "shortName"))
            if short == wanted:
                date = int(ec.codes_get(gid, "validityDate"))
                month = (date // 100) % 100
                lev = float(ec.codes_get(gid, "level")) if wanted != "sp" else 0.0
                if levels is None or lev in levels:
                    ni, nj = int(ec.codes_get(gid,"Ni")), int(ec.codes_get(gid,"Nj"))
                    a = np.asarray(ec.codes_get_array(gid,"values"), dtype=np.float32).reshape(nj,ni)
                    records[(month,lev)] = a
            ec.codes_release(gid)
    if wanted == "sp":
        return np.stack([records[(m,0.0)] for m in range(1,13)])
    return np.stack([[records[(m,float(p))] for p in levels] for m in range(1,13)])


def make_core_nc(path, lats, lons, native=True):
    ds=Dataset(path,"w",format="NETCDF4")
    ds.createDimension("time",len(DATES)); ds.createDimension("lat",len(lats)); ds.createDimension("lon",len(lons))
    t=ds.createVariable("time","f8",("time",)); t.units="days since 1982-01-01 00:00:00"; t.calendar="standard"; t[:]=date2num(DATES,t.units,t.calendar)
    ds.createVariable("lat","f4",("lat",))[:]=lats; ds.createVariable("lon","f4",("lon",))[:]=lons
    kw=dict(zlib=True,complevel=4,shuffle=True,fill_value=np.float32(9.96921e36),chunksizes=(1,min(60,len(lats)),min(120,len(lons))))
    for n,u,long in [("iwv","kg m-2","integrated water vapour"),("zbar","m","water-vapour mass-weighted height"),("wgpe","J m-2","water-vapour gravitational potential energy")]:
        v=ds.createVariable(n,"f4",("time","lat","lon"),**kw); v.units=u; v.long_name=long
    if native:
        for n,u in [("effective_dp_pa","Pa"),("surface_pressure_pa","Pa")]:
            v=ds.createVariable(n,"f4",("time","lat","lon"),**kw); v.units=u
        for n in ["n_valid_pressure_levels","excluded_below_ground_levels","surface_interpolation_flag"]:
            v=ds.createVariable(n,"i1",("time","lat","lon"),zlib=True,complevel=4,fill_value=np.int8(-127),chunksizes=(1,min(60,len(lats)),min(120,len(lons))))
    ds.integration="native 1.25-degree pressure integration before any horizontal remapping"
    ds.vertical_domain="100 hPa to min(local surface pressure, 1000 hPa); crossing interval truncated by linear q/gh interpolation"
    return ds


def native_to_1deg(a, target_lat, target_lon, categorical=False):
    # Native latitude is 90..-90; longitude is 0..358.75 and periodic.
    y=(90.0-target_lat)/1.25
    x=np.mod(target_lon,360.0)/1.25
    if categorical:
        yi=np.clip(np.rint(y).astype(int),0,a.shape[-2]-1); xi=np.mod(np.rint(x).astype(int),a.shape[-1])
        return a[np.ix_(yi,xi)]
    y0=np.floor(y).astype(int); y1=np.clip(y0+1,0,a.shape[-2]-1); wy=y-y0
    x0=np.floor(x).astype(int)%a.shape[-1]; x1=(x0+1)%a.shape[-1]; wx=x-np.floor(x)
    z00=a[np.ix_(y0,x0)]; z01=a[np.ix_(y0,x1)]; z10=a[np.ix_(y1,x0)]; z11=a[np.ix_(y1,x1)]
    return ((1-wy)[:,None]*((1-wx)[None,:]*z00+wx[None,:]*z01)+wy[:,None]*((1-wx)[None,:]*z10+wx[None,:]*z11)).astype(np.float32)


def main():
    files=inventory_files()
    nlat,nlon=145,288; lats=np.linspace(90,-90,nlat,dtype=np.float32); lons=np.arange(nlon,dtype=np.float32)*1.25
    tlats=np.arange(90,-91,-1,dtype=np.float32); tlons=np.arange(-179,181,1,dtype=np.float32)
    dn=make_core_nc(NATIVE,lats,lons,True); d1=make_core_nc(ONEDEG,tlats,tlons,False)
    # Sufficient statistics for native-grid first-order layer decomposition.
    shp=(len(BANDS),nlat,nlon); hcount=np.zeros(shp,np.int32); hst=np.zeros(shp); hst2=np.zeros(shp); sm=np.zeros(shp); sh=np.zeros(shp); stm=np.zeros(shp); sth=np.zeros(shp)
    total_count=np.zeros((nlat,nlon),np.int32); total_z=np.zeros((nlat,nlon)); total_tz=np.zeros((nlat,nlon))
    qc=[]; ti=0
    for year in range(1982,2024):
        q=read_year(files["q"][year],"q",LEVELS); gh=read_year(files["gh"][year],"gh",LEVELS); sp=read_year(files["sp"][year],"sp")
        if q.shape != (12,len(LEVELS),nlat,nlon) or gh.shape != q.shape or sp.shape != (12,nlat,nlon): raise RuntimeError(f"Shape mismatch {year}: {q.shape} {gh.shape} {sp.shape}")
        for mm in range(12):
            s=sp[mm].astype(np.float64); mass=np.zeros((nlat,nlon)); mz=np.zeros_like(mass)
            bm=np.zeros(shp); bz=np.zeros(shp); interp=np.zeros((nlat,nlon),bool)
            for ii in range(len(LEVELS)-1):
                pl,pu=LEVELS[ii],LEVELS[ii+1]; ql=q[mm,ii].astype(float); qu=q[mm,ii+1].astype(float); zl=gh[mm,ii].astype(float); zu=gh[mm,ii+1].astype(float)
                sph=s/100.; full=sph>=pl; part=(sph>pu)&(sph<pl)
                dp=np.zeros_like(s); dm=np.zeros_like(s); dz=np.zeros_like(s)
                ok=full & np.isfinite(ql)&np.isfinite(qu)&np.isfinite(zl)&np.isfinite(zu)
                dp[ok]=(pl-pu)*100.; dm[ok]=0.5*(ql[ok]+qu[ok])*dp[ok]/G; dz[ok]=0.5*(ql[ok]*zl[ok]+qu[ok]*zu[ok])*dp[ok]/G
                okp=part & np.isfinite(ql)&np.isfinite(qu)&np.isfinite(zl)&np.isfinite(zu)
                f=(pl-sph)/(pl-pu); qs=ql+f*(qu-ql); zs=zl+f*(zu-zl)
                dp[okp]=(sph[okp]-pu)*100.; dm[okp]=0.5*(qs[okp]+qu[okp])*dp[okp]/G; dz[okp]=0.5*(qs[okp]*zs[okp]+qu[okp]*zu[okp])*dp[okp]/G; interp|=okp
                k=next(k for k,(bl,bu) in enumerate(BANDS) if pl<=bl and pu>=bu)
                bm[k]+=dm; bz[k]+=dz; mass+=dm; mz+=dz
            valid=(mass>0)&np.isfinite(mass)&np.isfinite(mz); z=np.full_like(mass,np.nan); z[valid]=mz[valid]/mass[valid]; w=G*mass*z
            eff=np.maximum(0,np.minimum(s,100000.)-10000.); nvalid=np.sum(LEVELS[:,None,None]*100<=s[None,:,:],axis=0).astype(np.int8); nexcl=(len(LEVELS)-nvalid).astype(np.int8)
            for n,a in [("iwv",mass),("zbar",z),("wgpe",w),("effective_dp_pa",eff),("surface_pressure_pa",s),("n_valid_pressure_levels",nvalid),("excluded_below_ground_levels",nexcl),("surface_interpolation_flag",interp.astype(np.int8))]: dn[n][ti]=a
            for n,a in [("iwv",mass),("zbar",z),("wgpe",w)]: d1[n][ti]=native_to_1deg(a,tlats,tlons)
            t=(year-1982)+(mm+0.5)/12.; bh=np.divide(bz,bm,out=np.full_like(bz,np.nan),where=bm>0)
            ok=np.isfinite(bh); hcount+=ok; hst+=ok*t; hst2+=ok*t*t
            # A completely below-ground band has physically zero atmospheric
            # mass and must remain in the mass time series.  Its representative
            # height is undefined and is therefore excluded only from height.
            sm+=bm; stm+=t*bm; sh+=np.where(ok,bh,0); sth+=np.where(ok,t*bh,0)
            total_count+=valid; total_z+=np.where(valid,z,0); total_tz+=np.where(valid,t*z,0)
            qc.append(dict(date=f"{year}-{mm+1:02d}",sp_min_pa=float(np.nanmin(s)),sp_max_pa=float(np.nanmax(s)),fraction_surface_interpolated=float(interp.mean()),fraction_sp_below_1000hPa=float((s<100000).mean()),fraction_invalid_column=float((~valid).mean()),mean_excluded_below_ground_levels=float(nexcl.mean()),mean_effective_dp_pa=float(eff.mean())))
            ti+=1
        print(f"completed {year}",flush=True)
    dn.close(); d1.close()
    # OLS slopes from sufficient statistics; t sums are grid-independent for valid values but retained explicitly.
    times=np.array([(d.year-1982)+(d.month-0.5)/12 for d in DATES]); st=times.sum(); st2=(times**2).sum()
    def slope(sy,sty,n,st_local=None,st2_local=None):
        sl=st if st_local is None else st_local; s2=st2 if st2_local is None else st2_local
        den=n*s2-sl*sl; return np.divide(n*sty-sl*sy,den,out=np.full_like(sy,np.nan),where=(n>2)&(den!=0))
    mcount=np.full(shp,len(DATES),dtype=np.int32)
    mt=sm/len(DATES); ht=np.divide(sh,hcount,out=np.full_like(sh,np.nan),where=hcount>0)
    mtrend=slope(sm,stm,mcount); htrend=slope(sh,sth,hcount,hst,hst2); zmean=np.divide(total_z,total_count,out=np.full_like(total_z,np.nan),where=total_count>0); ztrend=slope(total_z,total_tz,total_count)
    # Total mean column mass is the sum of the four layer mean masses.  It is
    # not the mean over all layer observations (which would divide by four).
    M=np.nansum(mt,axis=0); massc=(ht-zmean[None,:,:])/M[None,:,:]*mtrend; heightc=mt/M[None,:,:]*htrend
    absent=(mt==0)&np.isfinite(M[None,:,:])&(M[None,:,:]>0);massc[absent]=0;heightc[absent]=0
    totalc=massc+heightc; residual=ztrend-np.sum(totalc,axis=0)
    dl=Dataset(LAYNC,"w",format="NETCDF4"); dl.createDimension("layer",len(BANDS)); dl.createDimension("lat",nlat); dl.createDimension("lon",nlon)
    dl.createVariable("lat","f4",("lat",))[:]=lats; dl.createVariable("lon","f4",("lon",))[:]=lons
    lv=dl.createVariable("layer","i1",("layer",)); lv[:]=np.arange(1,len(BANDS)+1); lv.layer_bounds="1000-850,850-700,700-500,500-300,300-200,200-100 hPa"
    for n,a in [("mass_redistribution_contribution",massc),("height_coordinate_contribution",heightc),("total_first_order_contribution",totalc)]: v=dl.createVariable(n,"f4",("layer","lat","lon"),zlib=True,complevel=4); v.units="m yr-1"; v[:]=a
    for n,a in [("actual_zbar_trend",ztrend),("closure_residual",residual)]: v=dl.createVariable(n,"f4",("lat","lon"),zlib=True,complevel=4); v.units="m yr-1"; v[:]=a
    dl.method="native-grid first-order centered ratio decomposition; complete balanced monthly record, for which raw and calendar-month-anomaly OLS slopes are algebraically identical"; dl.close()
    # Separate compact QC NetCDF requested by the protocol.
    src=Dataset(NATIVE); dq=Dataset(QCNC,"w",format="NETCDF4")
    for d,n in src.dimensions.items(): dq.createDimension(d,len(n))
    for n in ["time","lat","lon","effective_dp_pa","surface_pressure_pa","n_valid_pressure_levels","excluded_below_ground_levels","surface_interpolation_flag"]:
        sv=src[n]; dv=dq.createVariable(n,sv.dtype,sv.dimensions,zlib=(len(sv.dimensions)>1),complevel=4 if len(sv.dimensions)>1 else 0); dv[:]=sv[:]
        for a in sv.ncattrs():
            if a != "_FillValue": setattr(dv,a,getattr(sv,a))
    dq.close(); src.close()
    with open(QCCSV,"w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=list(qc[0])); w.writeheader(); w.writerows(qc)
    print(f"DONE months={ti} native={NATIVE} one_degree={ONEDEG}")


if __name__ == "__main__": main()
