"""Spatial smoothing of PV production: ramps of the aggregate against the number of points, correlation against distance."""
from __future__ import annotations

import numpy as np
import pandas as pd


def haversine_km(lat1, lon1, lat2, lon2):
    p = np.pi / 180
    a = np.sin((lat2 - lat1) * p / 2) ** 2 + np.cos(lat1 * p) * np.cos(lat2 * p) * np.sin((lon2 - lon1) * p / 2) ** 2
    return 12742.0 * np.arcsin(np.sqrt(a))


def distance_matrix(pts: pd.DataFrame) -> np.ndarray:
    la, lo = pts["lat"].to_numpy()[:, None], pts["lon"].to_numpy()[:, None]
    return haversine_km(la, lo, la.T, lo.T)


def daytime(clear: pd.DataFrame, min_mean: float = 0.05) -> np.ndarray:
    """Hours in which the mean clear-sky production over the points is above `min_mean` (and the hour before):
    ratios and ramps are taken only when the sun is clearly up."""
    ok = clear.mean(axis=1).to_numpy() > min_mean
    return ok & np.roll(ok, 1)


def index_series(cf: np.ndarray, clear: np.ndarray) -> np.ndarray:
    """Clear-sky index of an aggregate: production over clear-sky production (0 where the denominator is 0)."""
    den = clear.sum(axis=1)
    return np.divide(cf.sum(axis=1), den, out=np.zeros(len(den)), where=den > 1e-9)


def ramp_stats(series: np.ndarray, mask: np.ndarray) -> dict:
    r = np.diff(series, prepend=series[0])[mask]
    return {"ramp_sd": float(r.std()), "ramp_p99": float(np.quantile(np.abs(r), 0.99)), "ramp_max": float(np.abs(r).max())}


def aggregation_curve(cf: pd.DataFrame, clear: pd.DataFrame, sizes, draws: int = 200, seed: int = 0) -> pd.DataFrame:
    """For n = sizes: random subsets of n points; clear-sky index of the aggregate; hourly ramps of that index (sd, 99th
    percentile of the absolute ramp, maximum) and the sd of its level, averaged over the draws and divided by the single point."""
    rng = np.random.default_rng(seed)
    X, Cl = cf.to_numpy(), clear.to_numpy()
    mask = daytime(clear)
    rows = []
    for n in sizes:
        st = []
        for _ in range(1 if n == X.shape[1] else draws):
            idx = rng.choice(X.shape[1], size=n, replace=False)
            k = index_series(X[:, idx], Cl[:, idx])
            st.append({**ramp_stats(k, mask), "level_sd": float(k[mask].std())})
        rows.append({"n": n, **pd.DataFrame(st).mean().to_dict()})
    out = pd.DataFrame(rows)
    one = out.iloc[0]
    for c in ("ramp_sd", "ramp_p99", "ramp_max", "level_sd"):
        out[c + "_rel"] = out[c] / one[c]
    return out


def correlogram(cf: pd.DataFrame, clear: pd.DataFrame, pts: pd.DataFrame, bins) -> pd.DataFrame:
    """Correlation of the hourly clear-sky index (and of its hourly ramps) between pairs of points against their distance."""
    mask = daytime(clear)
    K = np.divide(cf.to_numpy(), clear.to_numpy(), out=np.zeros(cf.shape), where=clear.to_numpy() > 0.05)[mask]
    mk = mask & np.roll(mask, 1)
    Kfull = np.divide(cf.to_numpy(), clear.to_numpy(), out=np.zeros(cf.shape), where=clear.to_numpy() > 0.05)
    R = np.diff(Kfull, axis=0, prepend=Kfull[:1])[mk]
    C_r, C_l = np.corrcoef(R.T), np.corrcoef(K.T)
    D = distance_matrix(pts.loc[cf.columns])
    iu = np.triu_indices(len(D), 1)
    d, cr, cl = D[iu], C_r[iu], C_l[iu]
    rows = []
    for a, b in zip(bins[:-1], bins[1:]):
        s = (d >= a) & (d < b)
        if s.sum():
            rows.append({"km_from": a, "km_to": b, "pairs": int(s.sum()), "corr_ramps": float(cr[s].mean()), "corr_levels": float(cl[s].mean())})
    return pd.DataFrame(rows)


def decorrelation_distance(corr: pd.DataFrame, col: str = "corr_ramps", level: float = 0.5) -> float:
    """Distance (bin centre, km) at which the mean correlation first falls below `level` (linear interpolation); nan if never."""
    x = ((corr["km_from"] + corr["km_to"]) / 2).to_numpy()
    y = corr[col].to_numpy()
    for i in range(1, len(y)):
        if y[i] < level <= y[i - 1]:
            return float(x[i - 1] + (x[i] - x[i - 1]) * (y[i - 1] - level) / (y[i - 1] - y[i]))
    return float("nan") if (y >= level).all() else float(x[0])
