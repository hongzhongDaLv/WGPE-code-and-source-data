# -*- coding: utf-8 -*-
"""Common-period ERA5/JRA-55 trends, uncertainty, profiles and sensitivities."""
from __future__ import annotations
import csv, sys
from pathlib import Path
import numpy as np
ROOT=Path(r"__WGPE_PROJECT_ROOT__"); OUT=ROOT/"output_TOP100_revision"/"results"/"JRA55"
sys.path.insert(0,str(ROOT/"scripts_core_validation"/"_pydeps_netcdf_readable"))
from netCDF4 import Dataset  # noqa
JFILE=OUT/"JRA55_TOP100_WGPE_1deg_1982_2023.nc"; EFILE=OUT/"ERA5_TOP100_WGPE_1deg_1982_2023.nc"
VARS=[("iwv","kg m-2 yr-1"),("zbar","m yr-1"),("wgpe","J m-2 yr-1")]
REGIONS=[("Global",-90,90),("S high latitude",-90,-66.5),("S temperate",-66.5,-23.5),("Tropics",-23.5,23.5),("N temperate",23.5,66.5),("N high latitude",66.5,90)]
PERIODS=[(1982,2023),(1985,2023),(1994,2023),(2000,2023)]
RNG=np.random.default_rng(20260713)

def write(path,rows):
    with open(path,"w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else []);w.writeheader();w.writerows(rows)
def get_data(ds,n):
    a=np.asarray(ds[n][:],float)
    if ds[n].dimensions==( "lon","lat","time"): a=np.transpose(a,(2,1,0))
    return a
def deseason(a,months):
    b=a.copy()
    for m in range(12):
        ix=months==m;b[ix]-=np.nanmean(b[ix],axis=0)
    return b
def slope_map(a,t):
    ok=np.isfinite(a); n=ok.sum(0); sy=np.nansum(a,0); st=np.sum(ok*t[:,None,None],0); sty=np.nansum(a*t[:,None,None],0); st2=np.sum(ok*t[:,None,None]**2,0)
    den=n*st2-st*st
    return np.divide(n*sty-st*sy,den,out=np.full(a.shape[1:],np.nan),where=(n>24)&(den>0))
def wmean_map(a,lat,mask=None):
    w=np.cos(np.deg2rad(lat))[:,None]*np.ones((1,a.shape[-1]));ok=np.isfinite(a)
    if mask is not None: ok&=mask
    return np.nansum(np.where(ok,a*w,0))/np.sum(np.where(ok,w,0))
def spatial_series(a,lat,mask):
    w=np.cos(np.deg2rad(lat))[:,None]*np.ones(a.shape[1:]); ok=np.isfinite(a)&mask[None,:,:]
    return np.sum(np.where(ok,a*w,0),(1,2))/np.sum(np.where(ok,w,0),(1,2))
def ols(y,t):
    ok=np.isfinite(y);x=t[ok];z=y[ok];return np.sum((x-x.mean())*(z-z.mean()))/np.sum((x-x.mean())**2)
def nw(y,t,lag=None):
    ok=np.isfinite(y);y=y[ok];t=t[ok];X=np.column_stack((np.ones(len(t)),t)); b=np.linalg.solve(X.T@X,X.T@y);u=y-X@b
    if lag is None: lag=int(np.floor(4*(len(y)/100)**(2/9)))
    S=np.zeros((2,2))
    for i in range(len(y)): S+=u[i]**2*np.outer(X[i],X[i])
    for L in range(1,lag+1):
        q=1-L/(lag+1);g=np.zeros((2,2))
        for i in range(L,len(y)): g+=u[i]*u[i-L]*np.outer(X[i],X[i-L])
        S+=q*(g+g.T)
    V=np.linalg.inv(X.T@X)@S@np.linalg.inv(X.T@X);return b[1],np.sqrt(max(V[1,1],0)),lag
def block_ci(y,t,reps=1000,block=60):
    n=len(y); vals=[]; b0=ols(y,t); fitted=np.nanmean(y)+b0*(t-np.mean(t)); resid=y-fitted
    starts=np.arange(0,n-block+1)
    for _ in range(reps):
        ids=[]
        while len(ids)<n:
            s=int(RNG.choice(starts));ids.extend(range(s,s+block))
        rr=resid[np.array(ids[:n])];yy=fitted+rr;vals.append(ols(yy,t))
    return np.percentile(vals,[2.5,97.5])
def corr(a,b,w=None):
    ok=np.isfinite(a)&np.isfinite(b)
    x=a[ok];y=b[ok]
    if w is None:return np.corrcoef(x,y)[0,1]
    ww=w[ok];mx=np.sum(ww*x)/sum(ww);my=np.sum(ww*y)/sum(ww);return np.sum(ww*(x-mx)*(y-my))/np.sqrt(np.sum(ww*(x-mx)**2)*np.sum(ww*(y-my)**2))

def main():
    J=Dataset(JFILE);E=Dataset(EFILE);lat=np.asarray(J["lat"][:]);lon=np.asarray(J["lon"][:])
    assert np.allclose(lat,E["lat"][:])
    # JRA target longitude is -179..180 whereas ERA5 is 0..359.  Reorder ERA5
    # by exact modulo-360 coordinate matching; no horizontal interpolation.
    elon=np.asarray(E["lon"][:]); eorder=np.array([np.argmin(np.abs(((elon-x+180)%360)-180)) for x in lon])
    assert np.allclose(np.mod(elon[eorder],360),np.mod(lon,360))
    months=np.tile(np.arange(12),42); years=1982+np.arange(504)/12.; w2=np.cos(np.deg2rad(lat))[:,None]*np.ones((181,360))
    # Mean JRA surface pressure supplies reproducible common terrain masks.
    Q=Dataset(OUT/"JRA55_TOP100_surface_truncation_QC.nc"); spn=np.nanmean(np.asarray(Q["surface_pressure_pa"][:],float),axis=0);Q.close()
    # Nearest-neighbour terrain classification from the native 1.25-degree
    # surface-pressure grid; no interpolation of the mask itself.
    yi=np.clip(np.rint((90-lat)/1.25).astype(int),0,spn.shape[0]-1)
    xi=np.mod(np.rint(np.mod(lon,360)/1.25).astype(int),spn.shape[1])
    sp=spn[np.ix_(yi,xi)]
    terrain={"all-valid":np.ones((181,360),bool),"exclude_mean_sp_lt_850hPa":sp>=85000,"exclude_mean_sp_lt_700hPa":sp>=70000}
    global_rows=[];zone_rows=[];prof_rows=[];pattern_rows=[];sub_rows=[];subpattern_rows=[];terr_rows=[];maps={}
    for var,unit in VARS:
        data={"JRA55":get_data(J,var),"ERA5":get_data(E,var)[...,eorder]}
        trends={}; submaps={}
        for ds,a in data.items():
            da=deseason(a,months); tr=slope_map(da,years);trends[ds]=tr;maps[(ds,var)]=tr
            common=np.isfinite(trends["JRA55"] if "JRA55" in trends else tr)&np.isfinite(tr)
            for rn,lo,hi in REGIONS:
                rm=(lat[:,None]>=lo)&(lat[:,None]<=hi)&np.ones((1,360),bool); mask=rm&common
                est=wmean_map(tr,lat,mask); ser=spatial_series(da,lat,mask); b,se,L=nw(ser,years);ci=block_ci(ser,years)
                row=dict(dataset=ds,variable=var,region=rn,period="1982-2023",trend=est,unit=unit,newey_west_se=se,newey_west_ci_low=b-1.96*se,newey_west_ci_high=b+1.96*se,newey_west_lag=L,block5yr_ci_low=ci[0],block5yr_ci_high=ci[1],bootstrap_repetitions=1000,n_grid=int(mask.sum()))
                (global_rows if rn=="Global" else zone_rows).append(row)
            for s in range(-90,90,5):
                m=(lat[:,None]>=s)&(lat[:,None]<(s+5))&np.isfinite(tr);prof_rows.append(dict(dataset=ds,variable=var,latitude_bin_start=s,latitude_bin_end=s+5,latitude_bin_mid=s+2.5,trend=wmean_map(tr,lat,m),unit=unit,n_grid=int(m.sum())))
            for p0,p1 in PERIODS:
                ix=(years>=p0)&(years<p1+1);dm=deseason(a[ix],months[ix]);tm=years[ix];tmap=slope_map(dm,tm);submaps[(ds,p0)]=tmap
                for rn in ["Global","Tropics","N high latitude"]:
                    lo,hi=[(x[1],x[2]) for x in REGIONS if x[0]==rn][0];m=(lat[:,None]>=lo)&(lat[:,None]<=hi)&np.isfinite(tmap)
                    sub_rows.append(dict(dataset=ds,variable=var,region=rn,period=f"{p0}-{p1}",trend=wmean_map(tmap,lat,m),unit=unit,n_grid=int(m.sum())))
        for p0,p1 in PERIODS:
            aa=submaps[("ERA5",p0)];bb=submaps[("JRA55",p0)];cm=np.isfinite(aa)&np.isfinite(bb)
            for rn,lo,hi in [("Global",-90,90),("Tropics",-23.5,23.5),("N high latitude",66.5,90)]:
                mm=cm&(lat[:,None]>=lo)&(lat[:,None]<=hi);subpattern_rows.append(dict(variable=var,region=rn,period=f"{p0}-{p1}",coslat_weighted_pattern_r=corr(aa[mm],bb[mm],w2[mm]),sign_agreement_fraction=float(np.mean(np.sign(aa[mm])==np.sign(bb[mm]))),n_grid=int(mm.sum())))
        common=np.isfinite(trends["JRA55"])&np.isfinite(trends["ERA5"])
        eq=corr(trends["ERA5"],trends["JRA55"]);cw=corr(trends["ERA5"],trends["JRA55"],w2)
        for rn,lo,hi in REGIONS:
            m=common&(lat[:,None]>=lo)&(lat[:,None]<=hi);sg=np.mean(np.sign(trends["ERA5"][m])==np.sign(trends["JRA55"][m]));
            pattern_rows.append(dict(variable=var,region=rn,period="1982-2023",equal_grid_pattern_r=corr(trends["ERA5"][m],trends["JRA55"][m]),coslat_weighted_pattern_r=corr(trends["ERA5"][m],trends["JRA55"][m],w2[m]),sign_agreement_fraction=sg,n_grid=int(m.sum())))
        for tn,tmask in terrain.items():
            m=common&tmask
            for ds in ["ERA5","JRA55"]:terr_rows.append(dict(dataset=ds,variable=var,terrain_mask=tn,trend=wmean_map(trends[ds],lat,m),unit=unit,n_grid=int(m.sum()),pattern_r_common=corr(trends["ERA5"][m],trends["JRA55"][m],w2[m]),sign_agreement_fraction=float(np.mean(np.sign(trends["ERA5"][m])==np.sign(trends["JRA55"][m])))))
    J.close();E.close()
    # Profile correlations.
    for var,_ in VARS:
        a=np.array([r["trend"] for r in prof_rows if r["variable"]==var and r["dataset"]=="ERA5"]);b=np.array([r["trend"] for r in prof_rows if r["variable"]==var and r["dataset"]=="JRA55"])
        pattern_rows.append(dict(variable=var,region="5deg_latitude_profile",period="1982-2023",equal_grid_pattern_r=np.corrcoef(a,b)[0,1],coslat_weighted_pattern_r=np.nan,sign_agreement_fraction=float(np.mean(np.sign(a)==np.sign(b))),n_grid=36))
    write(OUT/"JRA55_TOP100_GLOBAL_TRENDS.csv",global_rows);write(OUT/"JRA55_TOP100_ZONE_TRENDS.csv",zone_rows);write(OUT/"JRA55_TOP100_LAT5_PROFILES.csv",prof_rows);write(OUT/"ERA5_JRA55_TOP100_PATTERN_METRICS.csv",pattern_rows);write(OUT/"TOP100_SUBPERIOD_SENSITIVITY.csv",sub_rows);write(OUT/"TOP100_SUBPERIOD_PATTERN_METRICS.csv",subpattern_rows);write(OUT/"TOP100_HIGH_TERRAIN_SENSITIVITY.csv",terr_rows)
    # Trend map NetCDF.
    d=Dataset(OUT/"ERA5_JRA55_TOP100_TREND_MAPS_1982_2023.nc","w",format="NETCDF4");d.createDimension("lat",181);d.createDimension("lon",360);d.createVariable("lat","f4",("lat",))[:]=lat;d.createVariable("lon","f4",("lon",))[:]=lon
    for (ds,v),a in maps.items():x=d.createVariable(f"{ds}_{v}_trend","f4",("lat","lon"),zlib=True,complevel=4);x[:]=a
    for v,_ in VARS:x=d.createVariable(f"ERA5_minus_JRA55_{v}_trend","f4",("lat","lon"),zlib=True,complevel=4);x[:]=maps[("ERA5",v)]-maps[("JRA55",v)]
    d.close();print("comparison complete")
if __name__=="__main__":main()
