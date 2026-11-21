"""Small cases with known answers for the PV chain, smoothing statistics, ENTSO-E parsing and the point sets."""
import numpy as np
import pandas as pd
import pytest

from pvregional import data, pv, smoothing as sm, validate as va


def test_inside_germany_and_norway_filters():
    assert data.inside_germany([52.5], [13.4])[0]          # Berlin
    assert data.inside_germany([48.1], [11.6])[0]          # Munich
    assert not data.inside_germany([48.9], [2.3])[0]       # Paris
    assert not data.inside_germany([54.5], [6.5])[0]       # North Sea
    assert data.norway_land([59.9], [10.7], [20.0])[0]     # Oslo
    assert not data.norway_land([59.9], [10.7], [0.0])[0]  # sea cell
    assert not data.norway_land([59.5], [12.2], [100.0])[0]  # Swedish side


def test_parse_entsoe_expands_a03_and_skips_consumption():
    xml = b"""<GL_MarketDocument xmlns="urn:x"><TimeSeries><inBiddingZone_Domain.mRID>A</inBiddingZone_Domain.mRID>
    <Period><timeInterval><start>2024-06-01T00:00Z</start><end>2024-06-01T01:00Z</end></timeInterval><resolution>PT15M</resolution>
    <Point><position>1</position><quantity>5</quantity></Point><Point><position>3</position><quantity>7</quantity></Point></Period></TimeSeries>
    <TimeSeries><outBiddingZone_Domain.mRID>A</outBiddingZone_Domain.mRID>
    <Period><timeInterval><start>2024-06-01T00:00Z</start><end>2024-06-01T01:00Z</end></timeInterval><resolution>PT15M</resolution>
    <Point><position>1</position><quantity>99</quantity></Point></Period></TimeSeries></GL_MarketDocument>"""
    s = data.parse_entsoe(xml)
    assert s.tolist() == [5.0, 5.0, 7.0, 7.0]             # positions 2 and 4 carry the last value forward


def test_clear_sky_noon_in_summer_is_high_and_night_is_zero():
    times = pd.date_range("2024-06-21 00:00", periods=24, freq="h", tz="UTC")
    temp, wind = pd.Series(20.0, index=times), pd.Series(2.0, index=times)
    cf = pv.clear_sky_cf(50.0, 10.0, times, temp, wind)
    assert cf.iloc[0] == 0.0 and cf.iloc[23] == 0.0
    assert 0.5 < cf.max() <= 1 / 1.15 + 1e-9
    assert cf.idxmax().hour in (11, 12, 13)


def test_snow_reduces_output_only_when_snow_falls():
    times = pd.date_range("2024-01-10 00:00", periods=96, freq="h", tz="UTC")
    z = pd.Series(0.0, index=times)
    ghi = pd.Series(np.where((times.hour >= 9) & (times.hour <= 15), 300.0, 0.0), index=times)
    w = {"shortwave_radiation": pd.DataFrame({"p": ghi}), "direct_normal_irradiance": pd.DataFrame({"p": ghi * 1.5}),
         "direct_radiation": pd.DataFrame({"p": ghi * 0.7}), "diffuse_radiation": pd.DataFrame({"p": ghi * 0.3}),
         "temperature_2m": pd.DataFrame({"p": z - 8.0}), "wind_speed_10m": pd.DataFrame({"p": z + 2.0}),
         "snow_depth": pd.DataFrame({"p": z}), "snowfall": pd.DataFrame({"p": z})}
    dry = pv.capacity_factor(60.0, 10.0, w, "p")
    assert dry["cf"].sum() == pytest.approx(dry["cf_nosnow"].sum())
    w["snowfall"] = pd.DataFrame({"p": pd.Series(np.where(times < times[6], 2.0, 0.0), index=times)})
    snowy = pv.capacity_factor(60.0, 10.0, w, "p")
    assert snowy["cf"].sum() < snowy["cf_nosnow"].sum()


def test_aggregation_smooths_independent_points_by_sqrt_n():
    rng = np.random.default_rng(0)
    T, n = 4000, 16
    clear = pd.DataFrame(np.ones((T, n)) * 0.5)
    cf = clear * (1 + 0.3 * rng.standard_normal((T, n)))
    c = sm.aggregation_curve(cf, clear, [1, 4, 16], draws=30)
    rel = dict(zip(c["n"], c["ramp_sd_rel"]))
    assert rel[4] == pytest.approx(0.5, abs=0.1) and rel[16] == pytest.approx(0.25, abs=0.07)


def test_identical_points_do_not_smooth_and_correlogram_is_one():
    rng = np.random.default_rng(1)
    base = 0.5 * (1 + 0.3 * rng.standard_normal(2000))
    cf = pd.DataFrame({"a": base, "b": base, "c": base})
    clear = pd.DataFrame({k: np.full(2000, 0.5) for k in "abc"})
    c = sm.aggregation_curve(cf, clear, [1, 3], draws=10)
    assert c["ramp_sd_rel"].iloc[1] == pytest.approx(1.0, abs=1e-6)
    pts = pd.DataFrame({"lat": [50.0, 50.1, 50.2], "lon": [10.0, 10.0, 10.0]}, index=list("abc"))
    g = sm.correlogram(cf, clear, pts, [0, 100])
    assert g["corr_levels"].iloc[0] == pytest.approx(1.0, abs=1e-6)


def test_haversine_one_degree_of_latitude():
    assert sm.haversine_km(50.0, 10.0, 51.0, 10.0) == pytest.approx(111.2, abs=0.3)


def test_decorrelation_distance_interpolates():
    corr = pd.DataFrame({"km_from": [0, 100, 200], "km_to": [100, 200, 300], "corr_ramps": [0.9, 0.6, 0.2]})
    assert sm.decorrelation_distance(corr) == pytest.approx(150 + 100 * (0.6 - 0.5) / (0.6 - 0.2), abs=1e-6)


def test_zone_assignment_examples():
    assert va.zone_of(48.7, 9.2) == "transnetbw"           # Stuttgart
    assert va.zone_of(52.5, 13.4) == "50hertz"             # Berlin
    assert va.zone_of(51.2, 6.8) == "amprion"              # Duesseldorf
    assert va.zone_of(48.1, 11.6) == "tennet"              # Munich


def test_compare_recovers_scale_and_perfect_correlation():
    idx = pd.date_range("2024-06-01", periods=96, freq="h", tz="UTC")
    clear = pd.Series(np.clip(np.sin((idx.hour - 6) / 12 * np.pi), 0, None) * 0.7, index=idx)
    model = clear * (0.8 + 0.2 * np.sin(np.arange(96)))
    meas = model * 1000.0
    r = va.compare(model, meas, clear)
    assert r["scale_mw_per_kwp_unit"] == pytest.approx(1000.0)
    assert r["corr"] == pytest.approx(1.0) and r["ramp_ratio_model_over_measured"] == pytest.approx(1.0)
