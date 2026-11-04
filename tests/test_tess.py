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
        sr.SearchConfig(detrend_window_days=0.5).validate()


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
    cfg = sr.SearchConfig()
    t, f = _sector(seed=seed)
    res = sr.search(t, f, 42, 1, cfg)
    assert res.signals == []
    assert res.stop_reason == "below_threshold"
    # SNR alone would have accepted this noise peak; the SDE gate is what refuses it.
    assert res.final_snr > cfg.snr_threshold and res.final_sde < cfg.sde_threshold


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
    assert sum(len(p) for p, _ in chunks) == len(periods)
    for p, d in chunks:
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
