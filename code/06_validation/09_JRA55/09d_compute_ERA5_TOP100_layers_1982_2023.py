from __future__ import annotations

import csv
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(r"__WGPE_PROJECT_ROOT__") / "scripts_core_validation" / "_pydeps_netcdf_readable"))
import netCDF4
import numpy as np

ROOT = Path(r"__WGPE_PROJECT_ROOT__")
PLEV = ROOT / "era5_plev_ascii"
OUT = ROOT / "output_TOP100_revision" / "results" / "JRA55"
SPFILE = ROOT / "output_reviewer_P0_resolution" / "ERA5_surface_pressure_aligned_0_359.nc"
G = 9.80665
LAYERS = [(1000,850),(850,700),(700,500),(500,300),(300,200),(200,100)]
LNAMES = [f"{a}-{b}" for a,b in LAYERS]
ZONES = [("Global",-90,90),("S high latitude",-90,-66.5),("S temperate",-66.5,-23.5),
         ("Tropics",-23.5,23.5),("N temperate",23.5,66.5),("N high latitude",66.5,90)]

def dates_for(ds,name="valid_time"):
    v=ds.variables[name]
    return [datetime(x.year,x.month,1) for x in netCDF4.num2date(v[:],v.units,only_use_cftime_datetimes=False)]

def decyear(d):
    y0=datetime(d.year,1,1); y1=datetime(d.year+1,1,1)
    return d.year+(d-y0).days/(y1-y0).days

files=[]; levels0=None
for p in sorted(PLEV.glob("*.nc")):
    d=netCDF4.Dataset(p); dd=dates_for(d); lev=np.asarray(d.variables["pressure_level"][:],float)
    if levels0 is None:
        levels0=lev; lat=np.asarray(d.variables["latitude"][:],float); lon=np.asarray(d.variables["longitude"][:],float)
    sel=np.asarray([i for i,x0 in enumerate(dd) if datetime(1982,1,1)<=x0<=datetime(2023,12,1)],int)
    if len(sel): files.append((dd[sel[0]],p,[dd[i] for i in sel],sel))
    d.close()
files.sort()
use=np.where((levels0<=1000)&(levels0>=100))[0]; plev=levels0[use]
dates=[datetime(y,m,1) for y in range(1982,2024) for m in range(1,13)]; x=np.asarray([d.year+(d.month-1)/12 for d in dates])
date_index={d:i for i,d in enumerate(dates)}
spds=netCDF4.Dataset(SPFILE)

shape=(12,6,len(lat),len(lon))
msum=np.zeros(shape,np.float64); mcount=np.zeros(shape,np.uint16)
hsum=np.zeros(shape,np.float64); hcount=np.zeros(shape,np.uint16)
qzsum=np.zeros(shape,np.float64); qzcount=np.zeros(shape,np.uint16)
hn=np.zeros((6,len(lat),len(lon)),np.uint16); hst=np.zeros((6,len(lat),len(lon)),np.float64); hst2=np.zeros_like(hst)

def integrate_block(src,ls,le,gs,ge):
    q=np.asarray(src.variables['q'][ls:le,use,:,:],np.float32)
    zh=np.asarray(src.variables['z'][ls:le,use,:,:],np.float32)/np.float32(G)
    q[(~np.isfinite(q))|(q<0)]=np.nan; zh[(~np.isfinite(zh))|(zh<0)]=np.nan
    sp=np.asarray(spds.variables['sp'][gs:ge],np.float32)
    nb=ge-gs; mass=np.zeros((nb,6,len(lat),len(lon)),np.float64); mom=np.zeros_like(mass)
    has=np.zeros((nb,6,len(lat),len(lon)),bool)
    for k in range(len(plev)-1):
        p0,p1=float(plev[k]),float(plev[k+1]); li=next(i for i,(lo,hi) in enumerate(LAYERS) if p0<=lo and p1>=hi)
        q0,q1=q[:,k],q[:,k+1]; z0,z1=zh[:,k],zh[:,k+1]
        ok=np.isfinite(q0)&np.isfinite(q1)&np.isfinite(z0)&np.isfinite(z1)&np.isfinite(sp)
        full=ok&(sp>=p0*100); partial=ok&(sp<p0*100)&(sp>p1*100); active=full|partial
        start=np.where(full,p0,sp/100); alpha=np.clip((p0-start)/(p0-p1),0,1)
        qs=q0+(q1-q0)*alpha; zs=z0+(z1-z0)*alpha; dp=np.maximum(start-p1,0)*100
        mass[:,li]+=np.where(active,.5*(qs+q1)*dp/G,0)
        mom[:,li]+=np.where(active,.5*(qs*zs+q1*z1)*dp/G,0)
        has[:,li]|=active
    height=np.where(mass>0,mom/mass,np.nan)
    return mass,height,mom

# Pass 1: monthly climatologies and valid-time geometry.
for _,p,dd,srcidx in files:
    src=netCDF4.Dataset(p)
    for ls in range(0,len(dd),40):
        le=min(ls+40,len(dd)); ss=int(srcidx[ls]); se=int(srcidx[le-1])+1
        spgs=(dd[ls].year-1979)*12+(dd[ls].month-1); spge=spgs+(le-ls)
        mass,height,mom=integrate_block(src,ss,se,spgs,spge)
        for bi in range(le-ls):
            mm=dd[ls+bi].month-1; tt=x[date_index[dd[ls+bi]]]
            for arr,su,co in [(mass,msum,mcount),(height,hsum,hcount),(mom,qzsum,qzcount)]:
                yy=arr[bi]; ok=np.isfinite(yy); su[mm]+=np.where(ok,yy,0); co[mm]+=ok
            hok=np.isfinite(height[bi]); hn+=hok; hst+=hok*tt; hst2+=hok*tt*tt
        print('pass1',dd[ls].strftime('%Y-%m'),dd[le-1].strftime('%Y-%m'),flush=True)
    src.close()

mclim=np.divide(msum,mcount,where=mcount>0,out=np.full_like(msum,np.nan))
hclim=np.divide(hsum,hcount,where=hcount>0,out=np.full_like(hsum,np.nan))
qzclim=np.divide(qzsum,qzcount,where=qzcount>0,out=np.full_like(qzsum,np.nan))
mmean=np.nansum(msum,axis=0)/np.maximum(np.sum(mcount,axis=0),1)
hmean=np.nansum(hsum,axis=0)/np.maximum(np.sum(hcount,axis=0),1)
qzmean=np.nansum(qzsum,axis=0)/np.maximum(np.sum(qzcount,axis=0),1)

# Pass 2: trend numerators after removal of each calendar-month climatology.
m_sy=np.zeros_like(mmean); m_sty=np.zeros_like(mmean)
h_sy=np.zeros_like(hmean); h_sty=np.zeros_like(hmean)
qz_sy=np.zeros_like(qzmean); qz_sty=np.zeros_like(qzmean)
for _,p,dd,srcidx in files:
    src=netCDF4.Dataset(p)
    for ls in range(0,len(dd),40):
        le=min(ls+40,len(dd)); ss=int(srcidx[ls]); se=int(srcidx[le-1])+1
        spgs=(dd[ls].year-1979)*12+(dd[ls].month-1); spge=spgs+(le-ls)
        mass,height,mom=integrate_block(src,ss,se,spgs,spge)
        for bi in range(le-ls):
            mm=dd[ls+bi].month-1; tt=x[date_index[dd[ls+bi]]]
            for arr,cl,sy,sty in [(mass,mclim,m_sy,m_sty),(height,hclim,h_sy,h_sty),(mom,qzclim,qz_sy,qz_sty)]:
                yy=arr[bi]-cl[mm]; ok=np.isfinite(yy); sy+=np.where(ok,yy,0); sty+=np.where(ok,yy*tt,0)
        print('pass2',dd[ls].strftime('%Y-%m'),dd[le-1].strftime('%Y-%m'),flush=True)
    src.close()
spds.close()

n=float(len(x)); st=np.sum(x); st2=np.sum(x*x); den=st2-st*st/n
mtrend=(m_sty-st*m_sy/n)/den; qztrend=(qz_sty-st*qz_sy/n)/den
hden=hst2-hst*hst/np.maximum(hn,1); htrend=(h_sty-hst*h_sy/np.maximum(hn,1))/hden
htrend[(hn<24)|(~np.isfinite(hden))|(hden<=0)]=np.nan

tf=netCDF4.Dataset(OUT/'ERA5_TOP100_WGPE_1deg_1982_2023.nc')
za=np.asarray(tf.variables['zbar'][:],float); tf.close()
zmean=np.nanmean(za,axis=0)
zclim=np.stack([np.nanmean(za[m::12],axis=0) for m in range(12)])
zc=za-zclim[np.arange(len(za))%12]
ztrend=np.nansum(zc*(x-x.mean())[:,None,None],axis=0)/np.sum((x-x.mean())**2)
Mmean=np.sum(mmean,axis=0)
mass_contrib=(hmean-zmean[None,:,:])/Mmean[None,:,:]*mtrend
height_contrib=mmean/Mmean[None,:,:]*htrend
zero=mmean==0; mass_contrib[zero]=0; height_contrib[zero]=0
total_contrib=mass_contrib+height_contrib

nc=netCDF4.Dataset(OUT/'ERA5_TOP100_LAYER_GRIDCELL_FIELDS_1982_2023.nc','w',format='NETCDF4')
nc.createDimension('layer',6);nc.createDimension('latitude',len(lat));nc.createDimension('longitude',len(lon))
nc.createVariable('layer_name',str,('layer',))[:]=np.asarray(LNAMES,dtype=object)
nc.createVariable('latitude','f4',('latitude',))[:]=lat;nc.createVariable('longitude','f4',('longitude',))[:]=lon
for name,data,unit in [('layer_mass_mean',mmean,'kg m-2'),('layer_height_mean',hmean,'m'),('layer_moment_mean',qzmean,'kg m-1'),
                       ('layer_mass_trend',mtrend,'kg m-2 yr-1'),('layer_height_trend',htrend,'m yr-1'),('layer_moment_trend',qztrend,'kg m-1 yr-1'),
                       ('mass_redistribution_contribution',mass_contrib,'m yr-1'),('height_coordinate_contribution',height_contrib,'m yr-1'),
                       ('total_first_order_contribution',total_contrib,'m yr-1')]:
    vv=nc.createVariable(name,'f4',('layer','latitude','longitude'),zlib=True,complevel=2);vv.units=unit;vv[:]=data.astype(np.float32)
nc.top_boundary_hPa=100;nc.method='same local-surface truncation and first-order diagnostic as archived FIXED_FULL, extended with 300-200 and 200-100 hPa'
nc.close()

w=np.broadcast_to(np.cos(np.deg2rad(lat))[:,None],zmean.shape)
def wm(a,mask):
    ok=np.isfinite(a)&mask; ww=np.where(ok,w,0); return float(np.nansum(a*ww)/ww.sum())
rows=[];latrows=[];rid=1
for region,lo,hi in ZONES:
    lm=(lat>=lo)&(lat<=hi)
    if region!='Global' and hi!=90: lm&=lat<hi
    mask=np.broadcast_to(lm[:,None],zmean.shape)
    actual=wm(ztrend,mask)
    for il,name in enumerate(LNAMES):
        rows.append({'result_id':f'P0TL{rid:05d}','region':region,'layer':name,'p_bottom_hPa':LAYERS[il][0],'p_top_hPa':LAYERS[il][1],
                     'n_grid':int(np.sum(mask&np.isfinite(total_contrib[il]))),'mean_layer_mass_kg_m2':wm(mmean[il],mask),
                     'mean_layer_height_m':wm(hmean[il],mask),'layer_mass_trend_kg_m2_yr':wm(mtrend[il],mask),
                     'layer_height_trend_m_yr':wm(htrend[il],mask),'layer_moment_trend_kg_m1_yr':wm(qztrend[il],mask),
                     'layer_WGPE_trend_J_m2_yr':G*wm(qztrend[il],mask),'mass_contribution_m_yr':wm(mass_contrib[il],mask),
                     'height_contribution_m_yr':wm(height_contrib[il],mask),'total_first_order_contribution_m_yr':wm(total_contrib[il],mask),
                     'actual_TOP100_zbar_trend_m_yr':actual,'height_reference':'absolute','area_weighting':'cos(latitude)'});rid+=1
for lo in range(-90,90,5):
    hi=lo+5;lm=(lat>=lo)&((lat<hi) if hi<90 else (lat<=hi));mask=np.broadcast_to(lm[:,None],zmean.shape)
    for il,name in enumerate(LNAMES):
        latrows.append({'result_id':f'P0TL{rid:05d}','latitude_bin_start':lo,'latitude_bin_end':hi,'latitude_bin_mid':lo+2.5,'layer':name,
                        'mass_contribution_m_yr':wm(mass_contrib[il],mask),'height_contribution_m_yr':wm(height_contrib[il],mask),
                        'total_first_order_contribution_m_yr':wm(total_contrib[il],mask),'layer_mass_mean_kg_m2':wm(mmean[il],mask),
                        'layer_moment_mean_kg_m1':wm(qzmean[il],mask)});rid+=1

for fn,rr in [('ERA5_TOP100_LAYER_CONTRIBUTIONS_1982_2023.csv',rows),('ERA5_TOP100_LAYER_LAT5_1982_2023.csv',latrows)]:
    with (OUT/fn).open('w',newline='',encoding='utf-8-sig') as f:
        dw=csv.DictWriter(f,fieldnames=rr[0].keys());dw.writeheader();dw.writerows(rr)

upper_mass=wm(np.sum(mmean[4:6],axis=0),np.ones(zmean.shape,bool));total_mass=wm(Mmean,np.ones(zmean.shape,bool))
upper_mom=wm(np.sum(qzmean[4:6],axis=0),np.ones(zmean.shape,bool));total_mom=wm(np.sum(qzmean,axis=0),np.ones(zmean.shape,bool))
summary=[{'metric':'above_300_IWV_fraction_pct','value':100*upper_mass/total_mass,'unit':'%'},
         {'metric':'above_300_geopotential_moment_fraction_pct','value':100*upper_mom/total_mom,'unit':'%'},
         {'metric':'above_300_first_order_ztrend_contribution_m_yr','value':sum(wm(total_contrib[i],np.ones(zmean.shape,bool)) for i in [4,5]),'unit':'m yr-1'},
         {'metric':'TOP100_actual_ztrend_m_yr','value':wm(ztrend,np.ones(zmean.shape,bool)),'unit':'m yr-1'}]
with (OUT/'ERA5_TOP100_UPPER_LAYER_SUMMARY_1982_2023.csv').open('w',newline='',encoding='utf-8-sig') as f:
    dw=csv.DictWriter(f,fieldnames=summary[0].keys());dw.writeheader();dw.writerows(summary)
print(summary)
