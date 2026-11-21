"""ERA5 hourly weather at a grid of points in Germany (1 degree) and southern Norway (0.5 degree, land points only) for 2024,
via the Open-Meteo archive API (model era5, CC BY 4.0). One file per point in data/raw/weather/.

    python scripts/download_weather.py
"""
import json
import time
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError

RAW = Path(__file__).resolve().parents[1] / "data" / "raw" / "weather"
VARS = "shortwave_radiation,direct_radiation,diffuse_radiation,direct_normal_irradiance,temperature_2m,wind_speed_10m,snow_depth,snowfall,cloud_cover"
DE = [(round(la, 2), round(lo, 2)) for la in [47.5 + i for i in range(8)] for lo in [6.5 + j for j in range(9)]]
NO = [(round(la, 2), round(lo, 2)) for la in [58.75 + 0.5 * i for i in range(7)] for lo in [7.0 + 0.5 * j for j in range(12)]]


def fetch(area, la, lo):
    out = RAW / f"{area}_{la:.2f}_{lo:.2f}.json"
    if out.exists():
        return
    url = (f"https://archive-api.open-meteo.com/v1/archive?latitude={la}&longitude={lo}&start_date=2024-01-01&end_date=2024-12-31"
           f"&hourly={VARS}&models=era5&timezone=UTC&wind_speed_unit=ms")
    for k in range(6):
        try:
            out.write_text(urlopen(url, timeout=60).read().decode(), encoding="utf-8")
            return
        except HTTPError as e:
            print("HTTP", e.code, area, la, lo, flush=True)
            time.sleep(30 if e.code == 429 else 5)
        except Exception as e:
            print("err", e, flush=True)
            time.sleep(5)


if __name__ == "__main__":
    RAW.mkdir(parents=True, exist_ok=True)
    for area, pts in (("de", DE), ("no", NO)):
        for la, lo in pts:
            fetch(area, la, lo)
            time.sleep(0.3)
        print(area, "done", flush=True)
