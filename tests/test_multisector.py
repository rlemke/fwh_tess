import numpy as np

from tesstransit import _lightcurve as lc
from tesstransit import _search as sr
from tesstransit import multisector as ms

from .test_tess import _box


def _curve(sector: int, start: float, planet, seed: int) -> lc.LightCurve:
    """One ~27 d sector of white noise with its own slow trend and the planet."""
    rng = np.random.default_rng(seed)
    t = np.arange(start, start + 27.0, 2.0 / 1440.0)
    f = (
        1.0
        + 0.002 * np.sin(2 * np.pi * t / (7.0 + sector))
        + rng.normal(0, 4e-4, len(t))
    )
    f = f + _box(t, *planet)
    n = len(t)
    return lc.LightCurve(sector, sector, t, f, np.full(n, 4e-4), n, 0, 0)


def test_parse_sectors():
    assert ms.parse_sectors("1-3,7") == [1, 2, 3, 7]


def test_bin_curve_means_within_bins():
    t = np.arange(0, 1, 2.0 / 1440.0)
    tb, fb = ms.bin_curve(t, np.ones_like(t), 10.0)
    assert len(tb) == 144 and np.allclose(fb, 1.0)


def test_stitched_search_finds_a_period_no_single_sector_shows_twice():
    # 35 d period: at most one transit in any 27 d sector, two across three.
    planet = (35.0, 1340.0, 4.0, 2500.0)
    curves = [_curve(s, 1325.0 + 28.0 * (s - 1), planet, s) for s in (1, 2, 3)]
    cfg = sr.SearchConfig(prewhiten_max_period_days=1.5)
    t, f, per = ms.stitch(curves, cfg, 10.0)
    assert [p["sector"] for p in per] == [1, 2, 3]
    scfg = ms.stitched_config(cfg, float(t[-1] - t[0]), 5.0)
    assert scfg.detrend_method == "none" and scfg.max_period_days > 35.0
    res = sr.search(t, f, 1, 0, scfg)
    assert res.signals and abs(res.signals[0].period_days - 35.0) < 0.05


def test_eligible_counts_sectors_per_star():
    m = {1: {10: "a", 11: "b"}, 2: {10: "c"}, 3: {10: "d", 12: "e"}}
    assert ms.eligible(m, 3) == {10: {1: "a", 2: "c", 3: "d"}}


def test_paired_sample_regroups_against_the_new_table_and_counts_drops():
    from tesstransit import _match as mt

    tces = [mt.TCE("1-01", 30, 1, 5.0, 0.0, 2.0, 500.0, 9.0, 3)]
    stars = {10: {1: "a"}, 30: {1: "b"}, 40: {1: "c"}}
    # 30 was quiet in the earlier run but has a TCE now; 20 lost its data.
    tce, quiet, dropped = ms.paired_sample(frozenset({10, 20, 30, 40}), tces, stars)
    assert (tce, quiet, dropped) == ([30], [10, 40], [20])


def test_stitched_defaults_differ_from_single_sector_only_where_declared():
    single = sr.SearchConfig()
    cfg = ms.base_config()
    assert cfg.sde_threshold == 13.0 and cfg.min_snr_without_two_strongest == 3.0
    assert single.min_snr_without_two_strongest == 0.0  # single sector unchanged
    # --set still wins over the stitched defaults
    assert cfg.override(["sde_threshold=11"]).sde_threshold == 11.0
