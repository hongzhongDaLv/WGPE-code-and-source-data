# Purpose: Compare ERA5 and JRA-55 TOP100 six-layer <z>-trend contributions
#          on the common 1982-2023 period.
# Inputs: JRA55_TOP100_z_layer_contributions_1982_2023.nc and
#         ERA5_TOP100_LAYER_GRIDCELL_FIELDS_1982_2023.nc.
# Outputs: JRA55_ERA5_TOP100_SIX_LAYER_SUMMARY.csv and pattern metrics CSV.
# Parameters: six locked layers; JRA native-first then bilinear field remapping;
#             cos(latitude) regional weighting.
# Dependencies: Python 3.12, numpy, netCDF4.
# Overwrite policy: TOP100-specific files only.

from pathlib import Path
import csv, sys
ROOT=Path(r"__WGPE_PROJECT_ROOT__")
OUT=ROOT/"output_TOP100_revision"/"results"/"JRA55"
sys.path.insert(0,str(ROOT/"scripts_core_validation"/"_pydeps_netcdf_readable"))
import numpy as np
from netCDF4 import Dataset

LAYERS=["1000-850","850-700","700-500","500-300","300-200","200-100"]
REGIONS=[("Global",-90,90),("S high latitude",-90,-66.5),("S temperate",-66.5,-23.5),
         ("Tropics",-23.5,23.5),("N temperate",23.5,66.5),("N high latitude",66.5,90)]
VARS=["mass_redistribution_contribution","height_coordinate_contribution","total_first_order_contribution"]

def write(path,rows):
    with open(path,"w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def native_to_1deg(a,target_lat,target_lon):
    y=(90.0-target_lat)/1.25; x=np.mod(target_lon,360.0)/1.25
    y0=np.floor(y).astype(int);y1=np.clip(y0+1,0,a.shape[-2]-1);wy=y-y0
    x0=np.floor(x).astype(int)%a.shape[-1];x1=(x0+1)%a.shape[-1];wx=x-np.floor(x)
    z00=a[np.ix_(y0,x0)];z01=a[np.ix_(y0,x1)];z10=a[np.ix_(y1,x0)];z11=a[np.ix_(y1,x1)]
    return (1-wy)[:,None]*((1-wx)[None,:]*z00+wx[None,:]*z01)+wy[:,None]*((1-wx)[None,:]*z10+wx[None,:]*z11)

def wm(a,lat,mask):
    w=np.cos(np.deg2rad(lat))[:,None];ok=np.isfinite(a)&mask
    return float(np.nansum(np.where(ok,a*w,0))/np.sum(np.where(ok,w,0)))

def wcorr(a,b,lat,mask):
    w=np.broadcast_to(np.cos(np.deg2rad(lat))[:,None],a.shape);ok=np.isfinite(a)&np.isfinite(b)&mask
    x=a[ok];y=b[ok];ww=w[ok];mx=np.sum(ww*x)/sum(ww);my=np.sum(ww*y)/sum(ww)
    return float(np.sum(ww*(x-mx)*(y-my))/np.sqrt(np.sum(ww*(x-mx)**2)*np.sum(ww*(y-my)**2)))

def main():
    je=Dataset(OUT/"JRA55_TOP100_z_layer_contributions_1982_2023.nc")
    ee=Dataset(OUT/"ERA5_TOP100_LAYER_GRIDCELL_FIELDS_1982_2023.nc")
    lat=np.asarray(ee["latitude"][:],float);lon=np.asarray(ee["longitude"][:],float)
    ef={v:np.asarray(ee[v][:],float) for v in VARS}
    jf={v:np.stack([native_to_1deg(np.asarray(je[v][k],float),lat,lon) for k in range(6)]) for v in VARS}
    rows=[];pat=[]
    for ri,(region,lo,hi) in enumerate(REGIONS):
        lm=(lat>=lo)&(lat<=hi)
        if region!="Global" and hi!=90: lm&=lat<hi
        mask=np.broadcast_to(lm[:,None],(181,360))
        eglob=[];jglob=[]
        for k,layer in enumerate(LAYERS):
            et=wm(ef["total_first_order_contribution"][k],lat,mask);jt=wm(jf["total_first_order_contribution"][k],lat,mask)
            eglob.append(et);jglob.append(jt)
            for v in VARS:
                ev=wm(ef[v][k],lat,mask);jv=wm(jf[v][k],lat,mask)
                rows.append(dict(region=region,layer=layer,component=v,dataset="ERA5",estimate_m_yr=ev,period="1982-2023",height_reference="absolute",area_weighting="cos(latitude)"))
                rows.append(dict(region=region,layer=layer,component=v,dataset="JRA55",estimate_m_yr=jv,period="1982-2023",height_reference="absolute",area_weighting="cos(latitude)"))
                rows.append(dict(region=region,layer=layer,component=v,dataset="JRA55_minus_ERA5",estimate_m_yr=jv-ev,period="1982-2023",height_reference="absolute",area_weighting="cos(latitude)"))
            a=ef["total_first_order_contribution"][k];b=jf["total_first_order_contribution"][k];ok=np.isfinite(a)&np.isfinite(b)&mask
            pat.append(dict(region=region,layer=layer,coslat_pattern_r=wcorr(a,b,lat,mask),sign_agreement_fraction=float(np.mean(np.sign(a[ok])==np.sign(b[ok]))),n_grid=int(ok.sum())))
        erank=np.argsort(-np.abs(eglob));jrank=np.argsort(-np.abs(jglob))
        for rank,k in enumerate(erank,1): rows.append(dict(region=region,layer=LAYERS[k],component="absolute_rank",dataset="ERA5",estimate_m_yr=rank,period="1982-2023",height_reference="absolute",area_weighting="cos(latitude)"))
        for rank,k in enumerate(jrank,1): rows.append(dict(region=region,layer=LAYERS[k],component="absolute_rank",dataset="JRA55",estimate_m_yr=rank,period="1982-2023",height_reference="absolute",area_weighting="cos(latitude)"))
    je.close();ee.close()
    write(OUT/"JRA55_ERA5_TOP100_SIX_LAYER_SUMMARY.csv",rows)
    write(OUT/"JRA55_ERA5_TOP100_SIX_LAYER_PATTERN_METRICS.csv",pat)
    print("six-layer comparison complete")

if __name__=="__main__": main()
