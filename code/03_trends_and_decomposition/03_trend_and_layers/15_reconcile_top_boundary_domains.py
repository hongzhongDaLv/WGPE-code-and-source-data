import csv
from pathlib import Path

import netCDF4
import numpy as np

ROOT=Path(r"__WGPE_PROJECT_ROOT__");OUT=ROOT/"output_reviewer_P0_resolution"
rows=list(csv.DictReader(open(OUT/"LOCKED_64529_COMMON_MASK.csv",encoding="utf-8-sig")))
mask_lookup={(int(round(float(r['lon'])))%360,int(round(float(r['lat'])))):r['locked_common_mask'].upper()=='TRUE' for r in rows}
nc=netCDF4.Dataset(OUT/"TOP_BOUNDARY_GRIDCELL_FIELDS.nc");lon=np.asarray(nc.variables['longitude'][:]);lat=np.asarray(nc.variables['latitude'][:]);w=np.broadcast_to(np.cos(np.deg2rad(lat))[:,None],(181,360));mask=np.asarray([[mask_lookup[(int(round(x))%360,int(round(y)))] for x in lon] for y in lat])
out=[];rid=1
for it,top in enumerate([300,200,100]):
 for var in ['iwv','zbar','wgpe']:
  a=np.asarray(nc.variables[f'{var}_trend'][it],float)
  for domain,m in [('all_valid_65160',np.isfinite(a)),('locked_layer_common_64529',mask&np.isfinite(a))]:
   val=float(np.nansum(a*np.where(m,w,0))/np.sum(np.where(m,w,0)))
   out.append({'result_id':f'P0TD{rid:05d}','top_hPa':top,'variable':var,'domain':domain,'n_grid':int(m.sum()),'trend':val,'unit':{'iwv':'kg m-2 yr-1','zbar':'m yr-1','wgpe':'J m-2 yr-1'}[var]});rid+=1
nc.close()
with open(OUT/'TOP_BOUNDARY_DOMAIN_RECONCILIATION.csv','w',newline='',encoding='utf-8-sig') as f:
 wri=csv.DictWriter(f,fieldnames=out[0].keys());wri.writeheader();wri.writerows(out)
print([r for r in out if r['variable']=='wgpe'])
