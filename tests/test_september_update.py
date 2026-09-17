"""Numerical/semantic checks of the published September source data (no raw rerun)."""
from pathlib import Path
import csv, math
from collections import defaultdict

ROOT=Path(__file__).resolve().parents[1]
E=ROOT/'source_data/ecology_20260904'
checks=[]
def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def check(name,condition):
    checks.append((name,bool(condition)))
def close(a,b,tol=1e-9):return math.isclose(float(a),float(b),abs_tol=tol,rel_tol=tol)
def mean(g):
    w=[math.cos(float(r['lat'])*math.pi/180) for r in g]
    return sum(float(r['r'])*v for r,v in zip(g,w))/sum(w)

grid=rows(E/'GRID_RESULTS.csv');groups=defaultdict(list)
for r in grid:groups[(r['target'],r['predictor'],r['model'])].append(r)
summary=rows(E/'SUMMARY_RESULTS.csv')
for r in summary:
    key=(r['target'],r['predictor'],r['model'])
    g=[x for x in groups[key] if r['region']=='Global' or x['region']==r['region']]
    label='_'.join(key+(r['region'],))
    check(label+'_n',len(g)==int(r['n_grid']))
    check(label+'_mean_all_valid',close(mean(g),r['mean_r']))
    w=[math.cos(float(x['lat'])*math.pi/180) for x in g]
    for field,q,sign in [('BH_positive_pct','q_BH',1),('BH_negative_pct','q_BH',-1),('BY_positive_pct','q_BY',1)]:
        p=100*sum(v for x,v in zip(g,w) if float(x[q])<=.05 and sign*float(x['r'])>0)/sum(w)
        check(label+'_'+field,close(p,r[field]))
    check(label+'_CI_finite_ordered',math.isfinite(float(r['lower95'])) and float(r['lower95'])<=float(r['upper95']))
expected={('GPP','W-GPE'):(12009,.146420775501812),('GPP','<z>'):(12009,-.0339226698380928),('LAI','W-GPE'):(8300,.066956385007),('LAI','<z>'):(8300,-.0154267647)}
for (target,pred),(n,value) in expected.items():
    a=groups[(target,pred,'Without_VPD')];b=groups[(target,pred,'With_VPD')]
    check(f'{target}_{pred}_frozen_samples',[(r['lon'],r['lat'],r['n']) for r in a]==[(r['lon'],r['lat'],r['n']) for r in b])
    check(f'{target}_{pred}_headline',len(a)==n and close(mean(a),value,1e-8))
    check(f'{target}_{pred}_grey_cells_in_mean',any(float(r['q_BH'])>.05 for r in a))

classes=rows(E/'temporal_lag3/F_CLASS_SUMMARY.csv')
for target,n in [('GPP',12009),('LAI',5788)]:
    gr=rows(E/f'temporal_lag3/{target}_GRID.csv')
    for model in ['Without_VPD','With_VPD']:
        g=[r for r in gr if r['model']==model and r['class'] in ['forward-only','reverse-only','bidirectional','none']]
        ss=[r for r in classes if r['target']==target and r['model']==model]
        check(target+model+'_valid_tested_n',len(g)==n)
        check(target+model+'_grid_pct_sum',close(sum(float(r['grid_pct']) for r in ss),100))
        check(target+model+'_area_pct_sum',close(sum(float(r['area_pct']) for r in ss),100))
        weights=[math.cos(float(r['lat'])*math.pi/180) for r in g]
        for r in ss:
            count=sum(x['class']==r['class'] for x in g)
            area=100*sum(w for x,w in zip(g,weights) if x['class']==r['class'])/sum(weights)
            check(target+model+r['class']+'_count',count==int(r['n_class']) and int(r['n_tested'])==len(g))
            check(target+model+r['class']+'_grid_pct',close(100*count/len(g),r['grid_pct']))
            check(target+model+r['class']+'_area_pct',close(area,r['area_pct']))
    a=[r for r in gr if r['model']=='Without_VPD'];b=[r for r in gr if r['model']=='With_VPD']
    check(target+'_paired_lag_windows',[(r['lon'],r['lat'],r['n_forward'],r['n_reverse']) for r in a]==[(r['lon'],r['lat'],r['n_forward'],r['n_reverse']) for r in b])
for r in rows(ROOT/'source_data/final/Figure3f_source_final.csv'):
    check(r['target']+r['class']+'_plotted_denominator',close(r['estimate'],r['grid_pct']) and r['lower95']==r['upper95']=='')

registry=rows(ROOT/'source_data/registry/WGPE_PUBLIC_RESULT_REGISTRY.csv')
check('registry_unique',len(registry)==len({r['result_id'] for r in registry}))
lookup={r['result_id']:r for r in registry}
for r in rows(ROOT/'source_data/final/Fig3e_source.csv'):
    check(r['result_id']+'_registered',close(r['estimate'],lookup[r['result_id']]['estimate']) and r['n_grid']==lookup[r['result_id']]['n_grid'])
for name in ['Fig1g_values.csv','Fig2d_values.csv']:
    for i,r in enumerate(rows(ROOT/'source_data/channel_intervals_20260916'/name)):
        check(name+str(i)+'_channel_CI',float(r['lower95'])<=float(r['value'])<=float(r['upper95']))

out=ROOT/'tests/output';out.mkdir(exist_ok=True)
with (out/'SEPTEMBER_UPDATE_CHECKS.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f);w.writerow(['check','status']);w.writerows((n,'PASS' if ok else 'FAIL') for n,ok in checks)
fail=[n for n,ok in checks if not ok]
(out/'SEPTEMBER_UPDATE_REPORT.md').write_text(f'# September update source-data QA\n\n{len(checks)-len(fail)}/{len(checks)} checks passed.\n\nChecks recompute all-valid-cell cosine-weighted means and significant area percentages from supplied grids, paired sample membership, temporal class grid/area denominators, current registry identity and channel-CI ordering. They do not rerun raw-data preprocessing, bootstrap resampling or full figure rendering.\n\nFailures: '+(', '.join(fail) or 'none')+'\n',encoding='utf-8')
print(f'{len(checks)-len(fail)}/{len(checks)} September checks passed')
if fail: print('\n'.join(fail))
raise SystemExit(bool(fail))
