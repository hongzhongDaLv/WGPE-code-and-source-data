from __future__ import annotations
import csv, math
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"source_data"/"final"
OUT=ROOT/"tests"/"output"; OUT.mkdir(parents=True,exist_ok=True)

def rows(name):
    with (SRC/name).open(encoding="utf-8-sig",newline="") as f: return list(csv.DictReader(f))
def num(x): return float(x)
def close(a,b,tol=1e-10): return math.isfinite(a) and math.isfinite(b) and abs(a-b)<=tol*max(1,abs(a),abs(b))

registry={r["result_id"]:r for r in csv.DictReader((ROOT/"source_data"/"registry"/"WGPE_PUBLIC_RESULT_REGISTRY.csv").open(encoding="utf-8-sig",newline=""))}
out=[]
z=rows("Fig1a_source.csv")[0]; zv=num(z["estimate"]); zr=num(registry[z["result_id"]]["estimate"])
out.append({"check":"global_z_trend_registry_match","value":zv,"reference":zr,"tolerance":"1e-10 relative","status":"PASS" if close(zv,zr) else "FAIL"})
w=rows("Fig1b_source.csv"); actual=next(num(r["estimate"]) for r in w if r["metric"]=="actual"); channels=sum(next(num(r["estimate"]) for r in w if r["metric"]==k) for k in ["IWV","z","interaction"])
out.append({"check":"global_WGPE_channel_closure","value":channels,"reference":actual,"tolerance":"1e-6 J m-2 yr-1","status":"PASS" if abs(channels-actual)<1e-6 else "FAIL"})
wg=next(r for r in w if r["metric"]=="wgpe_trend_J_m2_yr"); wr=num(registry[wg["result_id"]]["estimate"])
out.append({"check":"global_WGPE_trend_registry_match","value":num(wg["estimate"]),"reference":wr,"tolerance":"1e-10 relative","status":"PASS" if close(num(wg["estimate"]),wr) else "FAIL"})
f2=rows("Figure2e_source_final.csv"); p=[r for r in f2 if r["target"]=="Precipitation"]
for r in p:
    model={"M1-M0":"M1: IWV + <z>","M2-M0":"M2: IWV + <z> + IWV×<z>","M3-M0":"M3: IWV + W-GPE"}[r["comparison"]]
    rr=next(x for x in registry.values() if x.get("figure_panel")=="Figure 2e" and x.get("model")==model and x.get("metric")=="weighted_mean_delta_r2")
    v=num(r["weighted_mean_delta_r2"]); ref=num(rr["estimate"])
    out.append({"check":f"Figure2e_Precipitation_{r['comparison']}","value":v,"reference":ref,"tolerance":"1e-10 relative","status":"PASS" if close(v,ref) else "FAIL"})
eco=rows("Ecology_field_significance_final.csv"); ev=num(next(r for r in eco if r["target"]=="GPP" and r["predictor"]=="W-GPE")["observed_mean_partial_r"]); er=num(registry["FINAL_ECO_01"]["estimate"])
out.append({"check":"historical_six_control_GPP_WGPE_mean_partial_r","value":ev,"reference":er,"tolerance":"1e-10 relative","status":"PASS" if close(ev,er) else "FAIL"})

with (OUT/"SMOKE_TEST_RESULTS.csv").open("w",encoding="utf-8-sig",newline="") as f:
    wri=csv.DictWriter(f,fieldnames=list(out[0])); wri.writeheader(); wri.writerows(out)

layers=rows("Figure1e_source_final.csv")
img=Image.new("RGB",(1000,650),"white"); d=ImageDraw.Draw(img); font=ImageFont.load_default(); d.text((35,20),"WGPE Source Data smoke-test panel: vertical contributions",fill="black",font=font)
x0=450; scale=700
for i,r in enumerate(layers):
    y=80+i*100; v=num(r["exact_layer_total_m_yr"]); x1=x0+int(v*scale); color=(32,99,170) if v>=0 else (151,87,12)
    d.text((35,y+12),r["layer"]+" hPa",fill="black",font=font); d.rectangle((min(x0,x1),y,max(x0,x1),y+45),fill=color,outline="black"); d.line((x0,65,x0,575),fill=(80,80,80),width=2)
    d.text((max(x0,x1)+8 if v>=0 else min(x0,x1)-75,y+12),f"{v:+.3f}",fill="black",font=font)
d.text((35,610),"Author-generated test rendering; units: m yr-1",fill="black",font=font)
img.save(OUT/"SMOKE_TEST_VERTICAL_CONTRIBUTION_PANEL.png",dpi=(150,150))

status="PASS" if all(r["status"]=="PASS" for r in out) else "FAIL"
report=f"""# Smoke test report\n\nStatus: **{status}**\n\n- Checks run: {len(out)}\n- Global <z> trend: {zv:.12f} m yr-1.\n- Global W-GPE trend: {num(wg['estimate']):.12f} J m-2 yr-1.\n- Three-channel closure difference: {channels-actual:.3e} J m-2 yr-1.\n- Figure 2e precipitation increments checked: {len(p)}.\n- GPP-W-GPE partial r: {ev:.12f}.\n- A low-cost vertical-contribution test panel was regenerated from included Source Data.\n\nThis test validates compact Source Data and registry consistency; it is not a raw-data rerun.\n"""
(OUT/"SMOKE_TEST_REPORT.md").write_text(report,encoding="utf-8",newline="\n")
raise SystemExit(0 if status=="PASS" else 1)
