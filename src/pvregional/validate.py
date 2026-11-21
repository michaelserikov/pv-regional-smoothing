"""Check of the modelled regional production against ENTSO-E solar generation of the German control areas (2024)."""
from __future__ import annotations

import numpy as np
import pandas as pd


def zone_of(lat: float, lon: float) -> str:
    """Rough assignment of a grid point to a control area (approximate, by rectangles; not the real boundaries)."""
    if 47.0 <= lat < 49.8 and 7.5 <= lon < 10.2:
        return "transnetbw"
    if lon >= 11.0 and lat >= 50.5:
        return "50hertz"
    if lon < 8.5 and lat < 52.5:
        return "amprion"
    return "tennet"


def hourly(s: pd.Series) -> pd.Series:
    """15-minute MW (start of interval) to hourly mean MW, timestamp = end of the hour as in the weather data."""
    h = s.resample("1h").mean()
    h.index = h.index + pd.Timedelta(hours=1)
    return h


def compare(model_cf: pd.Series, meas_mw: pd.Series, clear: pd.Series) -> dict:
    """Scale the model to the measurement (one factor, installed capacity unknown) and score the shape.

    Level: correlation and normalised RMSE of the hourly series. Variability: standard deviation of the hourly ramps of the
    production divided by the clear-sky shape (clear-sky index), measured against modelled, in daytime hours."""
    idx = model_cf.index.intersection(meas_mw.index)
    m, y, c = model_cf.loc[idx].to_numpy(), meas_mw.loc[idx].to_numpy(), clear.loc[idx].to_numpy()
    a = float((m * y).sum() / (m * m).sum())
    ok = c > 0.05 * c.max() if False else (c > 0.05)
    mask = ok & np.roll(ok, 1)
    cap_clear = a * c                                   # the clear-sky shape in MW with the same scale as the model
    km = np.divide(m * a, cap_clear, out=np.zeros(len(m)), where=cap_clear > 1e-9)
    ky = np.divide(y, cap_clear, out=np.zeros(len(y)), where=cap_clear > 1e-9)
    rm, ry = np.diff(km, prepend=km[0])[mask], np.diff(ky, prepend=ky[0])[mask]
    return {"scale_mw_per_kwp_unit": a, "corr": float(np.corrcoef(m, y)[0, 1]),
            "nrmse": float(np.sqrt(((m * a - y) ** 2).mean()) / y.mean()), "bias_pct": float(100 * (m.sum() * a - y.sum()) / y.sum()),
            "ramp_sd_model": float(rm.std()), "ramp_sd_measured": float(ry.std()),
            "ramp_ratio_model_over_measured": float(rm.std() / ry.std()),
            "corr_ramps": float(np.corrcoef(rm, ry)[0, 1]), "hours": int(mask.sum())}
