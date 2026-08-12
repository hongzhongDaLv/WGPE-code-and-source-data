"""Add BY-FDR sensitivity and signed-area summaries to formal TOP100 ecology maps.

Inputs are the TOP100-recomputed grid-level partial correlations. The script
does not recompute the monthly regression models and does not alter legacy data.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(r"__WGPE_PROJECT_ROOT__")
INFILE = ROOT / "output_TOP100_revision/results/ecology/Ecology_partial_correlations_grid_FINAL_v3.csv"
OUT = ROOT / "output_TOP100_revision/results/correlation_significance"
OUT.mkdir(parents=True, exist_ok=True)
ALPHA = 0.05


def fdr(p, method):
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ii = np.flatnonzero(np.isfinite(p))
    if not len(ii):
        return q
    oo = ii[np.argsort(p[ii])]
    m = len(oo)
    cm = 1.0 if method == "BH" else np.sum(1/np.arange(1, m+1))
    raw = p[oo]*m*cm/np.arange(1, m+1)
    q[oo] = np.minimum(1, np.minimum.accumulate(raw[::-1])[::-1])
    return q


d = pd.read_csv(INFILE)
rows = []
out = []
keys = ["height_reference", "model", "target", "predictor"]
for key, x in d.groupby(keys, dropna=False, sort=False):
    x = x.copy()
    x["q_BH_recomputed"] = fdr(x.p_value, "BH")
    x["q_BY"] = fdr(x.p_value, "BY")
    x["q_BH_input"] = x.q_value
    w = np.cos(np.deg2rad(x.lat.to_numpy()))
    r = x.partial_r.to_numpy()
    valid = np.isfinite(r)
    den = np.sum(w[valid])
    for method, qcol in [("BH", "q_BH_recomputed"), ("BY", "q_BY")]:
        q = x[qcol].to_numpy()
        sig = valid & np.isfinite(q) & (q <= ALPHA)
        rows.append({
            **dict(zip(keys, key)), "FDR_method": method,
            "mean_partial_r_coslat": np.sum(w[valid]*r[valid])/den,
            "positive_significant_area_pct": 100*np.sum(w[sig & (r > 0)])/den,
            "negative_significant_area_pct": 100*np.sum(w[sig & (r < 0)])/den,
            "total_significant_area_pct": 100*np.sum(w[sig])/den,
            "n_grid": int(valid.sum()),
        })
    out.append(x)

allgrid = pd.concat(out, ignore_index=True)
allgrid.to_csv(OUT / "ECOLOGY_PARTIAL_CORRELATION_BH_BY_GRID_TOP100.csv", index=False)
pd.DataFrame(rows).to_csv(OUT / "ECOLOGY_PARTIAL_CORRELATION_BH_BY_SUMMARY_TOP100.csv", index=False)
maxdiff = np.nanmax(np.abs(allgrid.q_BH_input-allgrid.q_BH_recomputed))
pd.DataFrame([{"max_abs_input_vs_recomputed_BH_q": maxdiff,
               "n_rows": len(allgrid), "alpha": ALPHA}]).to_csv(
    OUT / "ECOLOGY_BH_BY_QA_TOP100.csv", index=False)
print("DONE rows", len(allgrid), "max BH difference", maxdiff)
