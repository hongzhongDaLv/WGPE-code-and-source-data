"""Create the formal TOP100 incremental-value robustness register.

The classification is deliberately conservative and restricted to the five
labels requested for the manuscript: robust, spatially selective, weak,
unstable, and no increment.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(r"__WGPE_PROJECT_ROOT__")
BASE = ROOT / "output_TOP100_revision/results/incremental_value"
OUT = BASE / "INCREMENTAL_VALUE_ROBUSTNESS_TOP100.csv"


def load_standard(path, target, metric, gap=3, extra=None):
    d = pd.read_csv(path)
    if "gap_months" in d:
        d = d[d.gap_months == gap]
    if extra:
        for k, v in extra.items():
            d = d[d[k] == v]
    est_col = metric
    area_col = next((c for c in ["positive_delta_CV_R2_area_pct", "positive_delta_AUC_area_pct", "positive_increment_area_pct"] if c in d), None)
    out = d[["increment", "region", est_col, "lower95", "upper95"] + ([area_col] if area_col else [])].copy()
    out = out.rename(columns={est_col: "estimate", area_col: "positive_area_pct"} if area_col else {est_col: "estimate"})
    out["target"] = target
    out["primary_metric"] = metric
    out["primary_gap_months"] = gap
    return out


parts = []
h = BASE / "hydroclimate_TOP100"
parts += [load_standard(h/"PRECIP_INCREMENT_SUMMARY.csv", "precipitation", "delta_CV_R2_mean")]
parts += [load_standard(h/"WET_INCREMENT_SUMMARY.csv", "wet_event", "delta_AUC")]
parts += [load_standard(h/"DRY_INCREMENT_SUMMARY.csv", "dry_event", "delta_AUC")]
e = BASE / "ENSO_TOP100"
parts += [load_standard(e/"NINO34_INCREMENT_SUMMARY_TOP100.csv", "Nino3.4", "delta_CV_R2")]
parts += [load_standard(e/"ENSO_WINTER_PRECIP_INCREMENT_SUMMARY_TOP100.csv", "ENSO_winter_precipitation", "delta_CV_R2")]
parts += [load_standard(e/"ENSO_PHASE_INCREMENT_SUMMARY_TOP100.csv", "ENSO_phase", "delta_AUC")]
eco = BASE / "ecology_TOP100"
parts += [load_standard(eco/"GPP_INCREMENT_BY_CONTROL.csv", "GPP", "delta_CV_R2", extra={"control_set":"Control3", "domain":"persistent_vegetated_domain"})]
parts += [load_standard(eco/"LAI_INCREMENT_BY_CONTROL.csv", "LAI", "delta_CV_R2", extra={"control_set":"Control3", "domain":"persistent_vegetated_domain"})]
d = pd.concat(parts, ignore_index=True)

rows = []
for (target, inc), x in d.groupby(["target", "increment"], sort=False):
    g = x[x.region == "Global"].iloc[0]
    r = x[x.region != "Global"]
    n_region = int(np.isfinite(r.estimate).sum())
    n_pos = int((r.estimate > 0).sum())
    n_ci_pos = int((r.lower95 > 0).sum())
    n_ci_neg = int((r.upper95 < 0).sum())
    area = float(g.positive_area_pct) if "positive_area_pct" in x and np.isfinite(g.positive_area_pct) else np.nan
    est, lo, hi = float(g.estimate), float(g.lower95), float(g.upper95)

    if hi <= 0 and (not np.isfinite(area) or area < 35):
        cls = "no increment"
    elif lo > 0 and n_pos >= max(3, n_region-1) and (not np.isfinite(area) or area >= 50):
        cls = "robust"
    elif lo > 0 and (n_pos >= 2 or (np.isfinite(area) and area >= 30)):
        cls = "spatially selective"
    elif lo <= 0 <= hi and est > 0 and n_ci_neg == 0:
        cls = "weak"
    elif n_ci_pos > 0 and n_ci_neg > 0:
        cls = "unstable"
    elif hi <= 0:
        cls = "no increment"
    else:
        cls = "unstable"

    rows.append({
        "target": target, "increment": inc, "primary_metric": g.primary_metric,
        "gap_months": int(g.primary_gap_months), "global_estimate": est,
        "global_lower95": lo, "global_upper95": hi,
        "positive_increment_area_pct": area,
        "n_finite_climate_zones": n_region, "n_positive_climate_zones": n_pos,
        "n_CI_positive_climate_zones": n_ci_pos, "n_CI_negative_climate_zones": n_ci_neg,
        "robustness_class": cls,
        "classification_rule": "global 95% CI + signed area + five standard-climate-zone consistency",
    })

out = pd.DataFrame(rows)
out.to_csv(OUT, index=False)
print(out.to_string(index=False))
