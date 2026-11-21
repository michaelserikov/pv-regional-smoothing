"""ENTSO-E actual solar generation (A75/A16, B16) of the four German control areas and of DE-LU, 2024, one file per month.
The API key is read from ENTSOE_API_KEY or ~/.entsoe_token.

    python scripts/download_entsoe_zones.py
"""
import os
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import urlopen

import pandas as pd

CACHE = Path(__file__).resolve().parents[1] / "data" / "raw" / "entsoe_zone"
ZONES = {"50hertz": "10YDE-VE-------2", "amprion": "10YDE-RWENET---I", "tennet": "10YDE-EON------1",
         "transnetbw": "10YDE-ENBW-----N", "delu": "10Y1001A1001A82H"}


def main():
    key = os.environ.get("ENTSOE_API_KEY") or (Path.home() / ".entsoe_token").read_text(encoding="utf-8").strip()
    CACHE.mkdir(parents=True, exist_ok=True)
    months = pd.date_range("2024-01-01", "2025-01-01", freq="MS")
    for name, eic in ZONES.items():
        for a, b in zip(months[:-1], months[1:]):
            out = CACHE / f"{name}_{a:%Y%m}.xml"
            if out.exists():
                continue
            q = dict(securityToken=key, documentType="A75", processType="A16", in_Domain=eic, psrType="B16",
                     periodStart=a.strftime("%Y%m%d%H%M"), periodEnd=b.strftime("%Y%m%d%H%M"))
            for k in range(4):
                try:
                    out.write_bytes(urlopen("https://web-api.tp.entsoe.eu/api?" + urlencode(q), timeout=120).read())
                    break
                except HTTPError as e:
                    body = e.read()[:300]
                    if e.code in (429, 503):
                        time.sleep(20 * (k + 1)); continue
                    out.write_bytes(body); print(name, a, e.code, body[:120]); break
            time.sleep(0.5)
        print(name, "done", flush=True)


if __name__ == "__main__":
    main()
