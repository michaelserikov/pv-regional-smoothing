"""Hosting capacity of a low-voltage feeder for PV, with the production profile of one weather point against a regional mean.

Network: SimBench 1-LV-rural1--0-sw (13 household loads, one 160 kVA transformer, radial feeders). Every load bus gets the same
PV size (kWp). The hosting capacity is the largest size for which, in every hour of a set of critical hours, no bus is above
1.05 pu and no line or the transformer is above 100% loading. Loads: the SimBench 15-minute profile of the net, averaged to
hours and mapped by position in the year (2016 profile, 2024 weather: both leap years, weekdays differ, which is ignored).
The critical hours are the ones with the highest PV output minus load in a reference case, not all hours.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

V_MAX = 1.05
CODE = "1-LV-rural1--0-sw"


class Feeder:
    def __init__(self):
        import pandapower as pp
        import simbench as sb
        self.pp = pp
        net = sb.get_simbench_net(CODE)
        prof = sb.get_absolute_values(net, profiles_instead_of_study_cases=True)
        lp = prof[("load", "p_mw")]
        n = len(lp) // 4
        self.load_h = lp.to_numpy().reshape(n, 4, -1).mean(axis=1)           # hours x loads, MW
        self.load_q = (net.load.q_mvar / net.load.p_mw).to_numpy()
        net.sgen.drop(net.sgen.index, inplace=True)
        for b in net.load.bus:
            pp.create_sgen(net, b, p_mw=0.0, q_mvar=0.0)
        self.net = net
        self.n_loads = len(net.load)
        self.hours = n

    def violated(self, hour: int, pv_mw_per_bus: float, cf: float) -> bool:
        net = self.net
        p = self.load_h[hour]
        net.load["p_mw"] = p
        net.load["q_mvar"] = p * self.load_q
        net.sgen["p_mw"] = cf * pv_mw_per_bus
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.pp.runpp(net, numba=False)
        return bool(net.res_bus.vm_pu.max() > V_MAX or net.res_line.loading_percent.max() > 100 or net.res_trafo.loading_percent.max() > 100)

    def critical_hours(self, cf: np.ndarray, k: int = 100, ref_mw: float = 0.13) -> np.ndarray:
        n = min(len(cf), self.hours)
        score = cf[:n] * ref_mw - self.load_h[:n].sum(axis=1)
        return np.argsort(score)[::-1][:k]

    def hosting_capacity(self, cf: np.ndarray, k: int = 100, hi_kwp: float = 40.0, tol: float = 0.1) -> dict:
        """kWp per household bus; bisection on the size, violation tested in the k critical hours."""
        hrs = self.critical_hours(cf, k)
        lo, hi = 0.0, hi_kwp
        if not any(self.violated(int(h), hi / 1000, float(cf[h])) for h in hrs):
            return {"kwp_per_house": hi_kwp, "bound": "above search range", "peak_cf": float(cf[hrs].max())}
        while hi - lo > tol:
            mid = (lo + hi) / 2
            if any(self.violated(int(h), mid / 1000, float(cf[h])) for h in hrs):
                hi = mid
            else:
                lo = mid
        return {"kwp_per_house": float(lo), "peak_cf": float(cf[hrs].max())}
