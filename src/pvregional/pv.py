"""PV capacity factor per weather point: transposition to a tilted plane, cell temperature, DC to AC, and a snow cover model.

Plane-of-array irradiance from ERA5 global, direct and diffuse radiation (Hay-Davies), cell temperature by the Faiman model,
PVWatts-type DC power with a temperature coefficient and system losses, inverter clipping at a DC/AC ratio of 1.15. The
snow loss follows the NREL model in pvlib (coverage from snowfall, with sliding and melting by temperature and irradiance)
and is the only part that is not checked against measurements of any kind. Hourly ERA5 values are means of the hour ending
at the timestamp, so the sun position is taken half an hour earlier.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import pvlib


@dataclass(frozen=True)
class System:
    tilt: float = 35.0
    azimuth: float = 180.0          # south
    gamma: float = -0.0035          # power temperature coefficient per K
    losses: float = 0.14            # soiling, wiring, mismatch, availability
    dc_ac: float = 1.15
    albedo: float = 0.2


def poa(lat, lon, ghi, dni, dhi, system: System, snow_albedo: pd.Series | None = None) -> pd.DataFrame:
    times = ghi.index
    mid = times - pd.Timedelta(minutes=30)
    sp = pvlib.solarposition.get_solarposition(mid, lat, lon)
    sp.index = times
    dni_extra = pvlib.irradiance.get_extra_radiation(mid).set_axis(times)
    am = pvlib.atmosphere.get_relative_airmass(sp["apparent_zenith"]).fillna(0)
    alb = system.albedo if snow_albedo is None else snow_albedo
    t = pvlib.irradiance.get_total_irradiance(system.tilt, system.azimuth, sp["apparent_zenith"], sp["azimuth"],
                                              dni.clip(lower=0), ghi.clip(lower=0), dhi.clip(lower=0), dni_extra=dni_extra,
                                              airmass=am, albedo=alb, model="haydavies")
    t["zenith"] = sp["apparent_zenith"]
    return t


def capacity_factor(lat, lon, w: dict, col: str, system: System = System(), snow: bool = True) -> pd.DataFrame:
    """AC power per kWp for one point; columns: cf, cf_nosnow, poa, snow_cover."""
    ghi, dni, dhi = w["shortwave_radiation"][col], w["direct_normal_irradiance"][col], w["direct_radiation"][col]
    dhi = w["diffuse_radiation"][col]
    # ground albedo: 0.2 bare, 0.6 with more than 5 cm of snow
    albedo = pd.Series(np.where(w["snow_depth"][col] > 0.05, 0.6, system.albedo), index=ghi.index)
    p = poa(lat, lon, ghi, dni, dhi, system, albedo)
    temp, wind = w["temperature_2m"][col], w["wind_speed_10m"][col]
    tcell = pvlib.temperature.faiman(p["poa_global"], temp, wind)
    dc = p["poa_global"] / 1000.0 * (1 + system.gamma * (tcell - 25.0)) * (1 - system.losses)
    ac = np.minimum(dc, 1.0 / system.dc_ac)
    out = pd.DataFrame({"cf_nosnow": ac.clip(lower=0), "poa": p["poa_global"], "tcell": tcell})
    cover = pd.Series(0.0, index=ghi.index)
    cf = out["cf_nosnow"]
    if snow:
        cover = pvlib.snow.coverage_nrel(w["snowfall"][col].fillna(0).clip(lower=0), p["poa_global"].fillna(0), temp,
                                         system.tilt, initial_coverage=0.0)
        loss = pvlib.snow.dc_loss_nrel(cover, num_strings=1)
        cf = cf * (1 - loss)
    out["cf"] = cf
    out["snow_cover"] = cover
    out["cf_clear"] = clear_sky_cf(lat, lon, ghi.index, temp, wind, system)
    return out


def clear_sky_cf(lat, lon, times, temp, wind, system: System = System()) -> pd.Series:
    """The same chain with clear-sky irradiance (Ineichen, Linke turbidity from the pvlib lookup), no snow: the daily shape."""
    mid = times - pd.Timedelta(minutes=30)
    loc = pvlib.location.Location(lat, lon)
    cs = loc.get_clearsky(mid, model="ineichen")
    cs.index = times
    p = poa(lat, lon, cs["ghi"], cs["dni"], cs["dhi"], system)
    tcell = pvlib.temperature.faiman(p["poa_global"], temp, wind)
    dc = p["poa_global"] / 1000.0 * (1 + system.gamma * (tcell - 25.0)) * (1 - system.losses)
    return np.minimum(dc, 1.0 / system.dc_ac).clip(lower=0)


def all_points(pts: pd.DataFrame, w: dict, system: System = System(), snow: bool = True) -> dict:
    """{point id: DataFrame} for the kept points."""
    return {i: capacity_factor(pts.loc[i, "lat"], pts.loc[i, "lon"], w, i, system, snow) for i in w["temperature_2m"].columns}


def matrix(res: dict, col: str = "cf") -> pd.DataFrame:
    return pd.DataFrame({i: r[col] for i, r in res.items()})
