"""Check that the critical hours used in the bisection (60 hours) are enough: at the found hosting capacity no hour among the
600 highest-scoring hours may violate a limit, and 5% above it at least one must (three profiles)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pvregional import data, feeder, pv  # noqa: E402

hc = pd.read_csv(Path(__file__).resolve().parents[1] / "results" / "hosting_capacity.csv", index_col=0)
f = feeder.Feeder()
rows = []
for a, pick in (("de", "mean"), ("de", "point"), ("no", "mean")):
    pts, w = data.load_points(a)
    cf = pv.matrix(pv.all_points(pts, w))
    name = f"{a}:mean" if pick == "mean" else hc.index[hc.index.str.startswith(a + ":") & (hc.index != a + ":mean")][5]
    prof = cf.mean(axis=1).to_numpy() if pick == "mean" else cf[name.split(":")[1]].to_numpy()
    cap = float(hc.loc[name, "kwp_per_house"])
    hrs = f.critical_hours(prof, k=600)
    at = sum(f.violated(int(h), cap / 1000, float(prof[h])) for h in hrs)
    above = sum(f.violated(int(h), 1.05 * cap / 1000, float(prof[h])) for h in hrs)
    rows.append({"profile": name, "capacity_kwp": cap, "violating_hours_at_capacity_of_600": at, "violating_hours_at_105pct_of_600": above})
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(Path(__file__).resolve().parents[1] / "results" / "hosting_check.csv", index=False)
