"""fwh_tess — light-curve loading, the BLS search, and matching against TCEs.

⚠️ The tests that matter are about what gets REPORTED: a dropped cadence is
counted, a TCE on a star we never searched is a miss rather than absent, an
alias is labelled as one, a TCE beyond the period grid is "outside the grid"
rather than "not detected", and the TCE table's duration is read in HOURS.

No test touches the network: light curves are synthetic and S3 is stubbed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tesstransit import _lightcurve as lc
from tesstransit import _match as mt
from tesstransit import _search as sr

FAST = sr.SearchConfig(n_periods=5000)


# --- a synthetic sector ------------------------------------------------------


def _box(
    t: np.ndarray, period: float, epoch: float, dur_h: float, depth_ppm: float
) -> np.ndarray:
    phase = (t - epoch + 0.5 * period) % period - 0.5 * period
    return np.where(np.abs(phase) < dur_h / 48.0, -depth_ppm * 1e-6, 0.0)


def _sector(planets=(), noise_ppm=400.0, seed=1) -> tuple[np.ndarray, np.ndarray]:
    """~27 d at 2-min cadence, a mid-sector downlink gap, a slow trend, noise."""
    rng = np.random.default_rng(seed)
    t = np.arange(1325.3, 1352.6, 2.0 / 1440.0)
    t = t[(t < 1338.5) | (t > 1339.6)]
    f = (
        1.0
        + 0.002 * np.sin(2 * np.pi * t / 9.0)
        + rng.normal(0, noise_ppm * 1e-6, len(t))
    )
    for p in planets:
        f = f + _box(t, *p)
    return t, f


# --- locating and loading a light curve -------------------------------------


def test_target_prefix_splits_padded_tic():
    assert (
        lc.target_prefix(261136679, 1) == "tess/public/tid/s0001/0000/0002/6113/6679/"
    )


@pytest.mark.parametrize("tic,sector", [(0, 1), (5, 0)])
def test_target_prefix_rejects_nonpositive(tic, sector):
    with pytest.raises(ValueError):
        lc.target_prefix(tic, sector)


class _StubS3:
    def __init__(self, keys: dict[str, bytes]):
        self.keys = keys
        self.downloads = 0

    def list_objects_v2(self, Bucket, Prefix):
        hits = [{"Key": k} for k in self.keys if k.startswith(Prefix)]
        return {"Contents": hits} if hits else {}

    def head_object(self, Bucket, Key):
        return {"ContentLength": len(self.keys[Key])}

    def download_file(self, Bucket, Key, Filename):
        self.downloads += 1
        Path(Filename).write_bytes(self.keys[Key])


def test_find_lc_key_ignores_other_products():
    pre = lc.target_prefix(42, 3)
    s3 = _StubS3(
        {pre + "x-s_tp.fits": b"", pre + "x-s_lc.fits": b"", pre + "x_dvr.pdf": b""}
    )
    assert lc.find_lc_key(s3, 42, 3) == pre + "x-s_lc.fits"


def test_find_lc_key_missing_is_an_error_not_none():
    with pytest.raises(FileNotFoundError, match="TIC 42 sector 3"):
        lc.find_lc_key(_StubS3({}), 42, 3)


def test_fetch_reuses_complete_file_and_refetches_truncated(tmp_path):
    key = lc.target_prefix(42, 3) + "x-s_lc.fits"
    s3 = _StubS3({key: b"0123456789"})
    p = lc.fetch_lc(s3, 42, 3, tmp_path)
    lc.fetch_lc(s3, 42, 3, tmp_path)
    assert s3.downloads == 1
    p.write_bytes(b"0123")  # a truncated earlier download
    lc.fetch_lc(s3, 42, 3, tmp_path)
    assert s3.downloads == 2 and p.read_bytes() == b"0123456789"


def _write_fits(path: Path, quality: np.ndarray, flux: np.ndarray) -> None:
    from astropy.io import fits

    n = len(flux)
    cols = fits.ColDefs(
        [
            fits.Column("TIME", "D", array=np.linspace(1325.0, 1326.0, n)),
            fits.Column("PDCSAP_FLUX", "E", array=flux),
            fits.Column("PDCSAP_FLUX_ERR", "E", array=np.full(n, 2.0)),
            fits.Column("QUALITY", "J", array=quality),
        ]
    )
    primary = fits.PrimaryHDU()
    primary.header["TICID"] = 42
    primary.header["SECTOR"] = 3
    fits.HDUList([primary, fits.BinTableHDU.from_columns(cols)]).writeto(path)


def test_load_lc_reports_what_it_dropped(tmp_path):
    flux = np.full(10, 1000.0)
    flux[3] = np.nan
    quality = np.zeros(10, dtype=np.int32)
    quality[[0, 1]] = 128  # flagged (manual exclude)
    quality[3] = 128  # flagged AND nan: counted once, as flagged
    _write_fits(tmp_path / "a.fits", quality, flux)
    got = lc.load_lc(tmp_path / "a.fits")
    assert (got.tic_id, got.sector) == (42, 3)
    assert (got.n_cadences, got.n_flagged, got.n_nonfinite, got.n_kept) == (10, 3, 0, 7)
    assert np.allclose(got.flux, 1.0) and np.allclose(got.flux_err, 0.002)


def test_load_lc_quality_bitmask_is_a_parameter(tmp_path):
    quality = np.zeros(10, dtype=np.int32)
    quality[:4] = 128
    _write_fits(tmp_path / "a.fits", quality, np.full(10, 1000.0))
    assert lc.load_lc(tmp_path / "a.fits", quality_bitmask=0).n_kept == 10


def test_load_lc_all_flagged_raises(tmp_path):
    _write_fits(
        tmp_path / "a.fits", np.full(10, 1, dtype=np.int32), np.full(10, 1000.0)
    )
    with pytest.raises(ValueError, match="no usable cadences"):
        lc.load_lc(tmp_path / "a.fits")


# --- the search --------------------------------------------------------------


def test_config_refuses_window_that_would_flatten_transits():
    with pytest.raises(ValueError, match="flatten"):
        sr.SearchConfig(detrend_method="median", detrend_window_days=0.5).validate()


def test_detrend_removes_slow_trend_keeps_transit():
    t, f = _sector(planets=[(3.7, 1326.1, 3.0, 2000.0)], noise_ppm=1.0)
    flat = sr.detrend(t, f, 1.0)
    out = np.abs((t - 1326.1 + 1.85) % 3.7 - 1.85) > 0.2
    assert np.std(flat[out]) < 5e-5
    assert np.min(flat) < 1 - 1500e-6


def test_search_recovers_injected_planet():
    t, f = _sector(planets=[(3.7, 1326.1, 3.0, 1500.0)])
    res = sr.search(t, f, 42, 1, FAST)
    s = res.signals[0]
    assert abs(s.period_days - 3.7) < 0.01
    assert abs((s.epoch_btjd - 1326.1 + 1.85) % 3.7 - 1.85) < 0.05
    assert s.snr > 20 and s.n_transits >= 6
    assert res.config["n_periods"] == 5000


def test_search_finds_second_planet_after_masking_first():
    t, f = _sector(planets=[(3.7, 1326.1, 3.0, 1500.0), (8.9, 1328.4, 4.0, 1200.0)])
    periods = sorted(s.period_days for s in sr.search(t, f, 42, 1, FAST).signals[:2])
    assert abs(periods[0] - 3.7) < 0.01 and abs(periods[1] - 8.9) < 0.03


@pytest.mark.parametrize("seed", [1, 4, 10])  # the seeds whose noise peak beat SNR 7.1
def test_search_on_noise_reports_why_it_stopped(seed):
    # The DEFAULT grid: SDE thresholds are calibrated on it, not on FAST.
    cfg = sr.SearchConfig(detrend_method="median")
    t, f = _sector(seed=seed)
    res = sr.search(t, f, 42, 1, cfg)
    assert res.signals == []
    assert res.stop_reason == "below_threshold"
    # SNR alone would have accepted this noise peak; the SDE gate is what refuses it.
    assert res.final_snr > cfg.snr_threshold and res.final_sde < cfg.sde_threshold
    # The spline default reports nothing on the same noise either.
    assert sr.search(t, f, 42, 1, sr.SearchConfig()).signals == []


# --- the TCE table and matching ---------------------------------------------

TCE_CSV = """# TESS DV TCE statistics
# second comment line
tceid,ticid,tce_plnt_num,tce_period,tce_time0bt,tce_duration,tce_depth,tce_model_snr,tce_num_transits
00261136679-01,261136679,1,6.2675,1325.5,2.97585,301.8,28.29,5
00000000042-01,42,1,3.7,1326.1,3.0,1500,40,7
00000000042-02,42,2,30.5,1330.0,5.0,900,9,1
00000000077-01,77,1,2.0,1325.9,2.0,800,12,13
"""


def test_parse_tce_csv_skips_comments_and_keeps_hours():
    tces = mt.parse_tce_csv(TCE_CSV)
    assert len(tces) == 4
    pi_men = tces[0]
    assert (pi_men.tic_id, pi_men.planet_num) == (261136679, 1)
    assert pi_men.duration_hours == pytest.approx(2.97585)  # HOURS, not days


def test_parse_tce_csv_missing_period_raises():
    bad = (
        "tceid,ticid,tce_plnt_num,tce_period,tce_time0bt,tce_duration\nx,1,1,,1325,2\n"
    )
    with pytest.raises(ValueError, match="lacks"):
        mt.parse_tce_csv(bad)


def _result(tic, signals, stop="below_threshold"):
    return sr.SearchResult(
        tic_id=tic,
        sector=1,
        config=FAST.to_dict(),
        n_points=1,
        scatter_ppm=1.0,
        signals=signals,
        stop_reason=stop,
    )


def _sig(rank, period, epoch, dur_h=3.0):
    return sr.Signal(rank, period, epoch, dur_h, 1000.0, 20.0, 10.0, 5)


def test_match_buckets_every_tce_and_signal():
    tces = mt.parse_tce_csv(TCE_CSV)
    results = {
        261136679: _result(261136679, [_sig(1, 6.2678, 1325.505)]),
        42: _result(42, [_sig(1, 7.4, 1326.1), _sig(2, 1.13, 1327.0)]),
    }
    rows = mt.match_sector(results, tces)
    by = {(r.tic_id, r.tce.planet_num if r.tce else None, r.bucket): r for r in rows}
    assert (261136679, 1, "recovered") in by
    alias = by[(42, 1, "recovered_alias")]
    assert alias.alias == 2.0  # we found twice the period
    assert by[(42, 2, "missed")].cause == "outside_search_grid"  # 30.5 d > 26 d
    assert by[(77, 1, "missed")].cause == "not_searched"
    assert (42, None, "unmatched") in by
    assert mt.summarise(rows) == {
        "recovered": 1,
        "recovered_alias": 1,
        "missed": 2,
        "missed:outside_search_grid": 1,
        "missed:not_searched": 1,
        "unmatched": 1,
    }


def test_match_right_period_wrong_epoch_is_not_a_match():
    tce = mt.parse_tce_csv(TCE_CSV)[1]  # P=3.7, t0=1326.1
    rows = mt.match_star(42, _result(42, [_sig(1, 3.7, 1326.1 + 1.0)]), [tce])
    assert [r.bucket for r in rows] == ["missed", "unmatched"]
    assert rows[0].cause == "not_detected"


def test_match_half_period_alias_checks_epoch_on_shorter_period():
    tce = mt.parse_tce_csv(TCE_CSV)[1]  # P=3.7, t0=1326.1
    # We found P/2 with an epoch at an odd half-cycle: still every TCE transit is ours.
    rows = mt.match_star(42, _result(42, [_sig(1, 1.85, 1326.1 + 3 * 1.85)]), [tce])
    assert rows[0].bucket == "recovered_alias" and rows[0].alias == 0.5


def test_match_reference_epoch_far_from_sector():
    """A TOI epoch fitted ~350 orbits after the sector (measured: TOI 105.01)."""
    from dataclasses import replace

    tce = replace(
        mt.parse_tce_csv(TCE_CSV)[1], period_days=2.184667, epoch_btjd=2084.58501
    )
    rows = mt.match_star(42, _result(42, [_sig(1, 2.18447, 1326.5086, 2.0)]), [tce])
    assert rows[0].bucket == "recovered"


def test_match_prefers_unaliased_signal():
    tce = mt.parse_tce_csv(TCE_CSV)[1]
    res = _result(42, [_sig(1, 7.4, 1326.1), _sig(2, 3.7, 1326.1)])
    rows = mt.match_star(42, res, [tce])
    assert rows[0].bucket == "recovered" and rows[0].signal.rank == 2
    assert rows[1].bucket == "unmatched" and rows[1].signal.rank == 1


def test_miss_after_max_signals_is_search_exhausted():
    tce = mt.parse_tce_csv(TCE_CSV)[1]
    res = _result(42, [_sig(1, 1.1, 1325.4)], stop="max_signals")
    assert mt.match_star(42, res, [tce])[0].cause == "search_exhausted"


def test_direct_match_wins_over_an_earlier_tce_alias():
    # Measured on TIC 350479101: TCE 1 at 11.92 d, TCE 2 at 5.97 d, one signal at 5.97 d.
    t1 = mt.TCE("x-01", 9, 1, 11.92, 1326.0, 4.0, None, None, None)
    t2 = mt.TCE("x-02", 9, 2, 5.97, 1326.0, 4.0, None, None, None)
    rows = mt.match_star(9, _result(9, [_sig(1, 5.97, 1326.0, 4.0)]), [t1, t2])
    assert [(r.tce.tce_id, r.bucket) for r in rows] == [
        ("x-01", "missed"),
        ("x-02", "recovered"),
    ]


def test_duration_cap_follows_keplers_law():
    # 0.1 rho_sun: ~3.9 h at 1 d, ~9.2 h at 13 d.
    assert sr.max_duration_days(1.0, 0.1) * 24 == pytest.approx(3.9, abs=0.1)
    assert sr.max_duration_days(13.0, 0.1) * 24 == pytest.approx(9.2, abs=0.1)


def test_chunks_never_try_a_duration_longer_than_allowed():
    periods = np.exp(np.linspace(np.log(0.25), np.log(13.0), 2000))
    durations = np.array(FAST.durations_hours) / 24.0
    chunks = sr._period_chunks(periods, durations, 0.1)
    assert sum(int(m.sum()) for m, _ in chunks) == len(periods)
    for m, d in chunks:
        p = periods[m]
        assert d.max() <= sr.max_duration_days(p.min(), 0.1) + 1e-12
        assert d.max() < p.min()


def test_wide_box_at_short_period_is_not_reported():
    # Variability shaped like an 8 h dip every 0.8 d (in "transit" 42% of the
    # time): no plausible transit is that long at that period.
    t, f = _sector(planets=[(0.8, 1325.5, 8.0, 400.0)])
    for s in sr.search(t, f, 42, 1, FAST).signals:
        assert s.duration_hours <= sr.max_duration_days(s.period_days, 0.1) * 24 + 1e-9


def test_long_period_two_transits_recovered():
    # 18 d: transits at 1328.0 and 1346.0, both outside the downlink gap.
    t, f = _sector(planets=[(18.0, 1328.0, 5.0, 1500.0)])
    s = sr.search(t, f, 42, 1, FAST).signals[0]
    assert abs(s.period_days - 18.0) < 0.1 and s.n_transits == 2


def test_single_transit_is_never_reported_as_periodic():
    # One 1500 ppm dip, its second transit would fall after the sector ends.
    t, f = _sector(planets=[(40.0, 1335.0, 5.0, 1500.0)])
    for s in sr.search(t, f, 42, 1, FAST).signals:
        assert s.n_transits >= 2


# --- per-transit consistency -------------------------------------------------


def test_binned_scatter_sees_red_noise():
    rng = np.random.default_rng(3)
    t = np.arange(0, 20, 2 / 1440)
    white = rng.normal(0, 1e-3, len(t))
    red = white + 5e-4 * np.sin(2 * np.pi * t / 0.3)  # structure on transit timescales
    w = sr.binned_scatter(t, white, 3 / 24)
    assert w == pytest.approx(1e-3 / np.sqrt(90), rel=0.25)
    assert sr.binned_scatter(t, red, 3 / 24) > 3 * w


def test_event_stats_equal_events_both_carry_the_signal():
    t, f = _sector(planets=[(18.0, 1328.0, 5.0, 1500.0)], noise_ppm=300)
    flat = sr.detrend(t, f, 1.0)
    ev = sr.event_stats(t, flat, 18.0, 1328.0, 5 / 24)
    assert len(ev.depths) == 2
    assert ev.snr_without_strongest > 0.5 * ev.snr_red / np.sqrt(2)


def test_one_strong_dip_folded_with_noise_is_vetoed():
    # One 2000 ppm, 5 h dip at 1330.0 and nothing at 1347.0: BLS at 17 d folds
    # the two into a "2-transit" signal. The second event alone shows nothing.
    t, f = _sector(noise_ppm=300)
    f = f + _box(t, 1000.0, 1330.0, 5.0, 2000.0)
    res = sr.search(t, f, 42, 1, FAST)
    assert res.vetoes, "the single-event fold should have been vetoed"
    v = res.vetoes[0]
    assert abs(v.strongest_event_btjd - 1330.0) < 0.1
    assert v.snr_without_strongest < FAST.min_snr_without_strongest
    for s in res.signals:  # whatever remains, it is not that dip
        assert s.snr_without_strongest >= FAST.min_snr_without_strongest


def test_real_two_transit_planet_survives_the_veto():
    t, f = _sector(planets=[(18.0, 1328.0, 5.0, 1500.0)])
    res = sr.search(t, f, 42, 1, FAST)
    assert abs(res.signals[0].period_days - 18.0) < 0.1
    assert not any(abs(v.period_days - 18.0) < 0.1 for v in res.vetoes)


def test_final_peak_is_located():
    t, f = _sector(seed=1)
    res = sr.search(t, f, 42, 1, sr.SearchConfig())
    assert res.final_period_days is not None and res.final_epoch_btjd is not None


def test_event_depth_is_measured_from_the_out_of_transit_level():
    # Flux sitting at 1.01 out of transit with 1000 ppm dips: depth 1000 ppm,
    # not -9000 (what measuring from 1.0 gives).
    t, f = _sector(planets=[(3.7, 1326.1, 3.0, 1000.0)], noise_ppm=100)
    flat = sr.detrend(t, f, 1.0) + 0.01
    ev = sr.event_stats(t, flat, 3.7, 1326.1, 3 / 24)
    assert np.median(ev.depths) == pytest.approx(1000e-6, rel=0.2)


# --- detrending variants -----------------------------------------------------


@pytest.mark.parametrize("method", ["median", "biweight"])
@pytest.mark.parametrize("per_duration", [False, True])
def test_every_detrend_variant_recovers_an_injected_planet(method, per_duration):
    cfg = sr.SearchConfig(
        n_periods=5000, detrend_method=method, detrend_per_duration=per_duration
    )
    t, f = _sector(planets=[(3.7, 1326.1, 3.0, 1500.0)])
    s = sr.search(t, f, 42, 1, cfg).signals[0]
    assert abs(s.period_days - 3.7) < 0.01


def test_per_duration_window_removes_short_period_variability():
    # 0.6 d stellar variability survives a 1-day window; a 3x-duration window
    # for a 2 h transit (6 h) removes most of it.
    t, _ = _sector(noise_ppm=50)
    f = 1.0 + 1e-3 * np.sin(2 * np.pi * t / 0.6)
    one_day = np.std(sr.detrend(t, f, 1.0))
    short = np.std(sr.detrend(t, f, 0.25))
    assert short < 0.5 * one_day


def test_config_refuses_per_duration_factor_that_would_flatten():
    with pytest.raises(ValueError, match="flatten"):
        sr.SearchConfig(detrend_per_duration=True, detrend_window_factor=2.0).validate()


def test_holdout_sample_shares_no_stars_and_keeps_snr_coverage():
    from tesstransit import pilot

    tces = [
        mt.TCE(f"{i}-01", i, 1, 3.0, 1326.0, 2.0, None, float(i), None)
        for i in range(1, 201)
    ]
    targets = {i: f"k{i}" for i in range(1, 401)}
    dev, dev_q, _ = pilot.select_sample(tces, targets, 20, 10, seed=1)
    ho, ho_q, _ = pilot.select_sample(
        tces, targets, 20, 10, seed=2, exclude=frozenset(dev + dev_q)
    )
    assert not set(ho + ho_q) & set(dev + dev_q)
    assert len(ho) == 20 and min(ho) < 20 and max(ho) > 180  # faint to strong


# --- the planet (TOI) sample -------------------------------------------------

_TOI_HEAD = (
    "# TESS Objects of Interest\n"
    "TIC,Full TOI ID,TOI Disposition,Orbital Epoch Value,Orbital Period (days) Value,"
    "Transit Duration (hours) Value,Transit Depth Value,Sectors,Detection Pipeline(s)\n"
)


def test_parse_tois_sectors_bjd_and_single_transit():
    from tesstransit import planets as pl

    text = _TOI_HEAD + (
        '42,101.01,CP,1326.1,3.7,3.0,1500,"[1, 27]",SPOC\n'
        "42,101.02,PC,2458330.0,0,5.0,800,[1],QLP\n"
    )
    a, b = pl.parse_tois(text)
    assert a.sectors == (1, 27) and a.epoch_btjd == 1326.1
    assert b.epoch_btjd == pytest.approx(1330.0)  # BJD converted, not mis-phased
    assert b.period_days is None  # single transit: not searchable
    assert a.as_tce().planet_num == 1 and b.as_tce().planet_num == 2


def test_forced_stats_separates_absent_from_present():
    from tesstransit import pilot

    t, f = _sector(planets=[(3.7, 1326.1, 3.0, 1500.0)])
    out = pilot.forced_stats(
        t, f, [("real", 3.7, 1326.1, 3.0), ("absent", 5.3, 1327.0, 3.0),
               ("long", 60.0, 1420.0, 3.0)], sr.SearchConfig()
    )  # fmt: skip
    assert out["real"]["snr_red"] > 20 and out["real"]["n_transits"] >= 6
    assert abs(out["absent"]["snr_red"]) < 5
    assert out["long"]["n_transits"] == 0 and out["long"]["snr_red"] is None


def test_classify_does_not_blame_the_search_for_absent_signals():
    from tesstransit import planets as pl

    cfg = sr.SearchConfig()
    toi = pl.TOI("1.01", 42, "PC", 5.0, 1326.0, 3.0, 500.0, (1,), "SPOC")
    miss = mt.MatchRow(42, "missed", tce=toi.as_tce(), cause="not_detected")
    assert pl.classify(toi, None, None, cfg) == "not_searched"
    single = pl.TOI("1.02", 42, "PC", None, 1326.0, 3.0, 500.0, (1,), "SPOC")
    assert pl.classify(single, None, None, cfg) == "single_transit"
    assert pl.classify(toi, miss, {"n_transits": 1, "snr_red": 30.0}, cfg) == (
        "weak_in_sector:few_transits"
    )
    assert pl.classify(toi, miss, {"n_transits": 5, "snr_red": 4.0}, cfg) == (
        "weak_in_sector:low_forced_snr"
    )
    # SPOC saw it in this sector: the signal is in the data and WE lost it.
    assert pl.classify(toi, miss, {"n_transits": 5, "snr_red": 4.0}, cfg, 25.0) == (
        "missed:weak_after_detrend"
    )
    assert pl.classify(toi, miss, {"n_transits": 5, "snr_red": 12.0}, cfg) == (
        "missed:not_detected"
    )
    hit = mt.MatchRow(42, "recovered", tce=toi.as_tce())
    assert pl.classify(toi, hit, None, cfg) == "recovered"


def test_local_ephemeris_keeps_toi_period_takes_tce_epoch():
    from tesstransit import planets as pl

    toi = pl.TOI("121.01", 42, "PC", 14.77411, 2074.329, 5.5, 15000.0, (1,), "SPOC",
                 period_err=0.00057, epoch_err=0.013)  # fmt: skip
    tce = mt.TCE("x-01", 42, 1, 2 * 14.7748, 1336.14, 5.4, 15000.0, 150.0, 2)
    loc = pl.local_ephemeris(toi, tce)
    assert loc.period_days == 14.77411 and loc.epoch_btjd == 1336.14
    assert loc.epoch_source == "tce" and pl.local_ephemeris(toi, None) is toi
    # 50 orbits back: sqrt(0.013^2 + (50 * 0.00057)^2) d
    assert pl.epoch_sigma_hours(toi, 1338.0) == pytest.approx(0.7516, rel=1e-3)


# --- spline detrending --------------------------------------------------------


def _spotted(seed=3):
    """A DS Tuc-like star: 2% rotation at 2.85 d, plus a 3 h, 3000 ppm transit."""
    t, f = _sector(planets=[(8.14, 1332.3, 3.0, 3000.0)], noise_ppm=300.0, seed=seed)
    return t, f + 0.02 * np.sin(2 * np.pi * t / 2.85)


def test_spline_follows_rotation_the_median_cannot():
    t, f = _spotted()
    eph = [("p", 8.14, 1332.3, 3.0)]
    from tesstransit import pilot

    med = pilot.forced_stats(t, f, eph, sr.SearchConfig(detrend_method="median"))
    med = med["p"]["snr_red"]
    spl = pilot.forced_stats(t, f, eph, sr.SearchConfig(detrend_method="spline"))
    assert spl["p"]["snr_red"] > 3 * med and spl["p"]["snr_red"] > 15


def test_spline_bic_picks_dense_knots_for_active_sparse_for_quiet():
    t, f = _spotted()
    knots = (0.5, 1.0, 3.0)
    assert sr.spline_trend(t, f, knots)[1] == 0.5
    tq, fq = _sector(noise_ppm=300.0)
    fq = fq - 0.002 * np.sin(2 * np.pi * tq / 9.0)  # remove _sector's slow trend
    assert sr.spline_trend(tq, fq, knots)[1] == 3.0


def test_spline_knots_shorter_than_transits_are_refused():
    with pytest.raises(ValueError, match="knot spacing"):
        sr.SearchConfig(
            detrend_method="spline", spline_knot_days=(0.25, 1.0)
        ).validate()


def test_search_records_the_knot_spacing_it_used():
    t, f = _spotted()
    res = sr.search(
        t, f, 42, 1, sr.SearchConfig(detrend_method="spline", n_periods=5000)
    )
    assert res.knot_days == [0.5]
    assert abs(res.signals[0].period_days - 8.14) < 0.05


# --- the stand-alone results page ---------------------------------------------


def test_site_page_embeds_data_safely():
    from tesstransit import site

    page = site.render({"tois": [{"toi": "</script><b>x"}], "tce_sample": [], "variants": []},
                       "2026-11-05")  # fmt: skip
    assert "/*__DATA__*/" not in page and "2026-11-05" in page
    # A value containing "</script>" must not close the data script early.
    assert page.count("</script>") == page.count("<script")
    assert "<\\/script><b>x" in page


def test_tce_backing_pairs_each_toi_with_its_own_tce():
    """L 98-59, sector 2: d's period is 0.94% from 2x c's; d has its own TCE."""
    from tesstransit import planets as pl

    def toi(n, p):
        return pl.TOI(f"175.0{n}", 7, "CP", p, 1355.0, 1.0, 1000.0, (2,), "SPOC")

    tois = [toi(1, 3.6907), toi(2, 7.4507), toi(3, 2.2531)]
    tces = [mt.TCE(f"x-0{i}", 7, i, p, 1355.0, 1.0, 1000.0, 10.0, 3)
            for i, p in ((1, 3.6905), (2, 7.451), (3, 2.2531))]  # fmt: skip
    got = {k: v.tce_id for k, v in pl.tce_backing(tois, tces).items()}
    assert got == {"175.01": "x-01", "175.02": "x-02", "175.03": "x-03"}
    # A TCE at an alias still backs a TOI with no direct TCE, but only one TOI.
    one = mt.TCE("z", 7, 1, 2 * 3.6907, 1355.0, 1.0, 1.0, 9.0, 2)
    shared = pl.tce_backing([toi(1, 3.6907), toi(4, 3.6907)], [one])
    assert [v.tce_id for v in shared.values() if v] == ["z"]
    halved = pl.tce_backing(
        [toi(2, 7.4507)], [mt.TCE("y", 7, 1, 3.72535, 1355.0, 1.0, 1.0, 9.0, 4)]
    )
    assert halved["175.02"].tce_id == "y"


# --- pre-whitening fast stellar variability -----------------------------------

PW = sr.SearchConfig(prewhiten_max_period_days=1.5)


def test_prewhiten_removes_fast_variability_and_restores_the_transit():
    from tesstransit import pilot

    t, f = _sector(planets=[(4.6, 1327.8, 3.0, 3000.0)], noise_ppm=300.0)
    f = f + 0.004 * np.sin(2 * np.pi * t / 0.63) + 0.0015 * np.sin(4 * np.pi * t / 0.63)
    _, removed = sr.prewhiten(t, f, PW)
    assert removed and abs(removed[0].period_days - 0.63) < 0.005
    eph = [("p", 4.6, 1327.8, 3.0)]
    plain = pilot.forced_stats(t, f, eph, sr.SearchConfig())["p"]["snr_red"]
    white = pilot.forced_stats(t, f, eph, PW)["p"]["snr_red"]
    assert white > 2 * plain and white > 15


def test_prewhiten_leaves_an_ultra_short_period_planet_alone():
    t, f = _sector(planets=[(0.75, 1325.6, 1.0, 2000.0)], noise_ppm=300.0)
    _, removed = sr.prewhiten(t, f, PW)
    assert removed == []  # one-sided: a transit only dims the star


def test_prewhiten_leaves_a_quiet_star_alone_and_is_recorded():
    t, f = _sector(noise_ppm=300.0)
    assert sr.prewhiten(t, f, PW)[1] == []
    tv, fv = _sector(planets=[(4.6, 1327.8, 3.0, 3000.0)], noise_ppm=300.0)
    fv = fv + 0.004 * np.sin(2 * np.pi * tv / 0.63)
    res = sr.search(
        tv, fv, 42, 1, sr.SearchConfig(prewhiten_max_period_days=1.5, n_periods=5000)
    )
    assert res.prewhitened and abs(res.signals[0].period_days - 4.6) < 0.02


def test_prewhiten_needs_the_spline():
    with pytest.raises(ValueError, match="spline"):
        sr.SearchConfig(
            detrend_method="median", prewhiten_max_period_days=1.5
        ).validate()
