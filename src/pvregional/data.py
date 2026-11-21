"""Weather points (ERA5 via Open-Meteo), the German and Norwegian point sets, and ENTSO-E solar generation by control area."""
from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.path import Path as Polygon

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"

# Coarse outline of Germany (lat, lon), hand-made from the border towns; used only to keep grid points that lie inside.
GERMANY = [(47.6, 7.6), (48.0, 7.55), (49.0, 8.2), (49.1, 6.4), (49.5, 6.1), (50.2, 6.1), (51.0, 5.9), (51.9, 6.9), (52.5, 7.0),
           (53.3, 7.1), (53.6, 6.7), (53.9, 8.7), (54.8, 8.6), (54.9, 9.5), (54.5, 10.2), (54.3, 11.0), (54.1, 12.0), (54.6, 13.4),
           (53.9, 14.2), (53.0, 14.3), (52.3, 14.6), (51.5, 14.8), (51.0, 15.0), (50.8, 14.3), (50.4, 12.2), (49.6, 12.5),
           (49.0, 13.8), (48.6, 13.9), (47.7, 13.0), (47.5, 12.2), (47.4, 10.4), (47.55, 9.7), (47.65, 9.2)]
VARS = ["shortwave_radiation", "direct_radiation", "diffuse_radiation", "direct_normal_irradiance", "temperature_2m",
        "wind_speed_10m", "snow_depth", "snowfall", "cloud_cover"]


def inside_germany(lat, lon) -> np.ndarray:
    poly = Polygon(np.array(GERMANY))
    return poly.contains_points(np.column_stack([lat, lon]))


def norway_land(lat, lon, elev) -> np.ndarray:
    """Southern Norway: land cells (elevation above 5 m) west of the Swedish border (about 11.5E at 59-60N, 12.0E further north)."""
    lat, lon, elev = map(np.asarray, (lat, lon, elev))
    return (elev > 5) & (lon <= np.where(lat < 60.3, 11.5, 12.0))


def load_points(area: str, folder: Path = RAW / "weather") -> tuple[pd.DataFrame, dict]:
    """Returns the point table (lat, lon, elevation, kept) and {variable: DataFrame hours x kept points} in UTC."""
    meta, frames = [], {v: {} for v in VARS}
    for f in sorted(folder.glob(f"{area}_*.json")):
        d = json.loads(f.read_text(encoding="utf-8"))
        meta.append({"id": f.stem, "lat": float(f.stem.split("_")[1]), "lon": float(f.stem.split("_")[2]),
                     "elevation": d.get("elevation", np.nan)})
        h = d["hourly"]
        idx = pd.to_datetime(h["time"], utc=True)
        for v in VARS:
            frames[v][f.stem] = pd.Series(h[v], index=idx, dtype=float)
    pts = pd.DataFrame(meta).set_index("id")
    keep = inside_germany(pts.lat, pts.lon) if area == "de" else norway_land(pts.lat, pts.lon, pts.elevation)
    pts["kept"] = keep
    ids = pts.index[keep]
    return pts, {v: pd.DataFrame({i: frames[v][i] for i in ids}) for v in VARS}


def parse_entsoe(xml: bytes) -> pd.Series:
    """ENTSO-E generation document -> MW indexed by UTC start of the interval. Curve type A03 omits positions whose value
    equals the previous one, so every period is expanded to all its positions by carrying the last value forward."""
    root = ET.fromstring(xml)
    out = {}
    for ts in (e for e in root.iter() if e.tag.endswith("TimeSeries")):
        if any(c.tag.endswith("outBiddingZone_Domain.mRID") for c in ts):
            continue                 # consumption series of the same unit (storage or own use): generation only
        for per in (e for e in ts if e.tag.endswith("Period")):
            ti = next(c for c in per if c.tag.endswith("timeInterval"))
            start = pd.Timestamp(next(c.text for c in ti if c.tag.endswith("start")))
            end = pd.Timestamp(next(c.text for c in ti if c.tag.endswith("end")))
            res = next(c.text for c in per if c.tag.endswith("resolution"))
            step = {"PT15M": 15, "PT30M": 30, "PT60M": 60}[res]
            n = int((end - start) / pd.Timedelta(minutes=step))
            vals = np.full(n, np.nan)
            for pt in (c for c in per if c.tag.endswith("Point")):
                pos = int(next(c.text for c in pt if c.tag.endswith("position")))
                if 1 <= pos <= n:
                    vals[pos - 1] = float(next(c.text for c in pt if c.tag.endswith("quantity")))
            idx = start + pd.to_timedelta(np.arange(n) * step, unit="min")
            for t, v in zip(idx, pd.Series(vals).ffill().to_numpy()):
                out[t] = v
    return pd.Series(out, dtype=float).sort_index()


def entsoe_zone(zone: str, folder: Path = RAW / "entsoe_zone") -> pd.Series:
    parts = []
    for f in sorted(folder.glob(f"{zone}_2024*.xml")):
        b = f.read_bytes()
        if b.lstrip().startswith(b"<?xml") and b"GL_MarketDocument" in b[:400]:
            parts.append(parse_entsoe(b))
    s = pd.concat(parts).sort_index()
    return s[~s.index.duplicated()]
