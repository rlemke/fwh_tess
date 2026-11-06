"""Injection-recovery: the injected transit's geometry and the recovery verdict."""

from __future__ import annotations

import numpy as np
import pytest

from tesstransit import _search as sr
from tesstransit import inject as ij

from .test_tess import FAST, _sector


def test_solar_duration_matches_earth_and_hot_jupiter():
    # Earth across the Sun is ~13 h; a 3-day orbit ~2.6 h. Grazing shortens it.
    assert ij.solar_duration_hours(365.25, 0.0) == pytest.approx(13.0, abs=0.1)
    assert ij.solar_duration_hours(3.0, 0.0) == pytest.approx(2.62, abs=0.05)
    assert ij.solar_duration_hours(3.0, 0.8) < 0.65 * ij.solar_duration_hours(3.0, 0.0)


def test_trapezoid_flat_bottom_and_ingress():
    t = np.linspace(-0.5, 0.5, 100001)
    loss = ij.trapezoid(t, 10.0, 0.0, 2.4, 0.24, 1e-3)  # 0.1 d, ingress 0.01 d
    assert loss.max() == pytest.approx(1e-3)
    assert np.ptp(t[loss > 0]) == pytest.approx(0.1, abs=1e-4)
    assert np.ptp(t[loss >= 1e-3]) == pytest.approx(0.08, abs=1e-4)
    # Periodic: the next transit, one period later, is the same shape.
    assert ij.trapezoid(np.array([10.0]), 10.0, 0.0, 2.4, 0.24, 1e-3)[0] == 1e-3


def test_draw_sets_depth_from_the_expected_snr():
    t, _ = _sector()
    rng = np.random.default_rng(3)
    for i in range(20):
        inj = ij.draw(rng, t, 4e-4, i)
        assert ij.MIN_PERIOD_DAYS <= inj.period_days <= ij.MAX_PERIOD_DAYS
        assert inj.n_transits >= 2
        ing = ij._ingress_hours(inj.duration_hours, inj.impact, inj.radius_ratio)
        flat = ij.trapezoid(
            t, inj.period_days, inj.epoch_btjd, inj.duration_hours, ing, 1.0
        )
        n_in = int((flat >= 1.0).sum())
        # The depth is from the pre-update duration, so allow a small mismatch.
        assert inj.depth_ppm * 1e-6 * np.sqrt(n_in) / 4e-4 == pytest.approx(
            inj.snr_expected, rel=0.1
        )


def test_strong_injection_is_recovered_and_a_weak_one_is_not():
    t, f = _sector(noise_ppm=400.0)
    rng = np.random.default_rng(7)
    ij_strong = ij.draw(rng, t, 4e-4, 0)
    found = {}
    for snr in (40.0, 3.0):
        inj = ij.Injection(**{**ij_strong.__dict__, "snr_expected": snr})
        ing = ij._ingress_hours(inj.duration_hours, inj.impact, inj.radius_ratio)
        n_in = (
            ij.trapezoid(
                t, inj.period_days, inj.epoch_btjd, inj.duration_hours, ing, 1.0
            )
            >= 1
        ).sum()
        depth = snr * 4e-4 / np.sqrt(n_in)
        loss = ij.trapezoid(
            t, inj.period_days, inj.epoch_btjd, inj.duration_hours, ing, depth
        )
        res = sr.search(t, f * (1 - loss), 42, 1, FAST)
        found[snr] = ij.recovered(res, inj)[0]
    assert found[40.0] in ("recovered", "recovered_alias")
    assert found[3.0] == "missed"
