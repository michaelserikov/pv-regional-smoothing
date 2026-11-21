"""All runs of the study: python -m pvregional.run   (about 10 minutes on 8 cores, hosting capacity is most of it)."""
from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor

import matplotlib
import numpy as np
import pandas as pd

from . import data, feeder, pv, smoothing as sm, validate as va

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RES, FIG = data.ROOT / "results", data.ROOT / "figures"
BINS = [0, 50, 100, 150, 200, 300, 400, 600, 800, 1200]
NAME = {"de": "Germany", "no": "South Norway"}


def build(area: str):
    pts, w = data.load_points(area)
    res = pv.all_points(pts, w)
    return pts, w, pv.matrix(res), pv.matrix(res, "cf_nosnow"), pv.matrix(res, "cf_clear"), pv.matrix(res, "snow_cover")


def _hc(args):
    f = feeder.Feeder()
    name, cf, k = args
    return name, f.hosting_capacity(np.asarray(cf), k=k)


def main():
    RES.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    out = {}
    areas = {a: build(a) for a in ("de", "no")}
    agg, cor = {}, {}
    for a, (pts, w, cf, cf0, cl, cover) in areas.items():
        n = cf.shape[1]
        sizes = sorted({1, 2, 4, 8, 16, 32, n})
        agg[a] = sm.aggregation_curve(cf0, cl, sizes)
        agg[a].to_csv(RES / f"aggregation_{a}.csv", index=False, float_format="%.4f")
        cor[a] = sm.correlogram(cf0, cl, pts, BINS)
        cor[a].to_csv(RES / f"correlogram_{a}.csv", index=False, float_format="%.4f")
        out[a] = {"points": int(n), "mean_cf": float(cf0.mean().mean()), "yield_kwh_per_kwp": float(cf0.mean().mean() * 8784),
                  "decorrelation_km_ramps": sm.decorrelation_distance(cor[a], "corr_ramps"),
                  "decorrelation_km_levels": sm.decorrelation_distance(cor[a], "corr_levels"),
                  "ramp_sd_rel_all_points": float(agg[a]["ramp_sd_rel"].iloc[-1]),
                  "ramp_p99_rel_all_points": float(agg[a]["ramp_p99_rel"].iloc[-1]),
                  "ramp_max_rel_all_points": float(agg[a]["ramp_max_rel"].iloc[-1]),
                  "ramp_sd_rel_8_points": float(agg[a].loc[agg[a].n == 8, "ramp_sd_rel"].iloc[0])}
    snow = {}
    for a, (pts, w, cf, cf0, cl, cover) in areas.items():
        loss_pt = 1 - cf.sum() / cf0.sum()
        m = 1 - cf.sum(axis=1).groupby(cf.index.month).sum() / cf0.sum(axis=1).groupby(cf0.index.month).sum()
        snow[a] = {"annual_loss_mean_pct": float(100 * (1 - cf.to_numpy().sum() / cf0.to_numpy().sum())),
                   "annual_loss_point_max_pct": float(100 * loss_pt.max()), "annual_loss_point_min_pct": float(100 * loss_pt.min()),
                   "monthly_loss_pct": {int(k): float(100 * v) for k, v in m.items()},
                   "corr_loss_elevation": float(np.corrcoef(loss_pt.to_numpy(), pts.loc[loss_pt.index, "elevation"].to_numpy())[0, 1])}
        pd.DataFrame({"lat": pts.loc[loss_pt.index, "lat"], "lon": pts.loc[loss_pt.index, "lon"],
                      "elevation": pts.loc[loss_pt.index, "elevation"], "snow_loss_pct": 100 * loss_pt}).to_csv(
            RES / f"snow_loss_{a}.csv", float_format="%.3f")
    out["snow"] = snow
    seas = pd.DataFrame({a: areas[a][3].mean(axis=1).groupby(areas[a][3].index.month).mean() for a in areas})
    seas.to_csv(RES / "monthly_cf.csv", float_format="%.4f")
    pts, w, cf, cf0, cl, cover = areas["de"]
    zones = pd.Series({i: va.zone_of(pts.loc[i, "lat"], pts.loc[i, "lon"]) for i in cf.columns})
    rows = {}
    for zone in ["50hertz", "amprion", "tennet", "transnetbw", "all"]:
        cols = cf.columns if zone == "all" else zones.index[zones == zone]
        meas = data.entsoe_zone("delu" if zone == "all" else zone)
        r = va.compare(cf0[cols].mean(axis=1), va.hourly(meas), cl[cols].mean(axis=1))
        r["points"] = int(len(cols))
        rows[zone] = r
    pd.DataFrame(rows).T.to_csv(RES / "validation_zones.csv", float_format="%.4f")
    out["validation"] = rows
    meas = data.entsoe_zone("delu")
    ref = meas.quantile(0.999)
    out["measured_ramps_national"] = {"note": "ENTSO-E DE-LU, ramps divided by the 99.9th percentile of production",
                                      "sd_15min_ramp": float((meas.diff().dropna() / ref).std()),
                                      "sd_60min_ramp_of_hourly_mean": float((va.hourly(meas).diff().dropna() / ref).std())}
    tasks = []
    for a, (pts, w, cf, cf0, cl, cover) in areas.items():
        tasks.append((f"{a}:mean", cf.mean(axis=1).to_numpy(), 60))
        for c in cf.columns:
            tasks.append((f"{a}:{c}", cf[c].to_numpy(), 60))
    hc = {}
    with ProcessPoolExecutor(8) as ex:
        for name, r in ex.map(_hc, tasks, chunksize=2):
            hc[name] = r
    hct = pd.DataFrame(hc).T
    hct.to_csv(RES / "hosting_capacity.csv", float_format="%.3f")
    for a in areas:
        pt = hct[hct.index.str.startswith(a + ":") & ~(hct.index == a + ":mean")]["kwp_per_house"].astype(float)
        mean = float(hct.loc[a + ":mean", "kwp_per_house"])
        out[a]["hosting_kwp_per_house"] = {"mean_profile": mean, "single_point_mean": float(pt.mean()),
                                           "single_point_min": float(pt.min()), "single_point_max": float(pt.max()),
                                           "mean_profile_vs_single_pct": float(100 * (mean / pt.mean() - 1)),
                                           "mean_profile_vs_lowest_pct": float(100 * (mean / pt.min() - 1))}
    (RES / "summary.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    figures(agg, cor, seas, hct)
    print(json.dumps(out, indent=1, default=float)[:7000])


def figures(agg, cor, seas, hct):
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    for a, c in (("de", "#1f77b4"), ("no", "#d62728")):
        ax[0].plot(agg[a]["n"], agg[a]["ramp_sd_rel"], "o-", color=c, label=NAME[a])
        x = (cor[a]["km_from"] + cor[a]["km_to"]) / 2
        ax[1].plot(x, cor[a]["corr_ramps"], "o-", color=c)
        ax[1].plot(x, cor[a]["corr_levels"], "s--", color=c, alpha=0.6)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("points aggregated")
    ax[0].set_ylabel("hourly ramp sd of clear-sky index,\nrelative to one point")
    ax[0].legend(fontsize=8)
    ax[0].set_title("Smoothing by aggregation")
    ax[1].set_xlabel("distance between points (km)")
    ax[1].set_ylabel("correlation")
    ax[1].set_title("Solid: hourly ramps, dashed: levels")
    ax[1].axhline(0.5, color="grey", lw=0.5)
    ax[2].plot(seas.index, seas["de"], "o-", color="#1f77b4", label="Germany")
    ax[2].plot(seas.index, seas["no"], "o-", color="#d62728", label="South Norway")
    ax[2].set_xlabel("month 2024")
    ax[2].set_ylabel("mean capacity factor")
    ax[2].legend(fontsize=8)
    ax[2].set_title("Seasonality (35 deg south roof, no snow)")
    fig.tight_layout()
    fig.savefig(FIG / "smoothing_and_seasons.png", dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    for a, c in (("de", "#1f77b4"), ("no", "#d62728")):
        v = hct[hct.index.str.startswith(a + ":") & ~(hct.index == a + ":mean")]["kwp_per_house"].astype(float)
        ax.hist(v, bins=12, alpha=0.5, color=c, label=f"{NAME[a]}: single points")
        ax.axvline(float(hct.loc[a + ":mean", "kwp_per_house"]), color=c, lw=2, ls="--")
    ax.set_xlabel("hosting capacity (kWp per household connection)")
    ax.set_ylabel("points")
    ax.legend(fontsize=8)
    ax.set_title("Dashed: regional mean profile")
    fig.tight_layout()
    fig.savefig(FIG / "hosting_capacity.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
