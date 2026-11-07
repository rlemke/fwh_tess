import numpy as np
import pytest

from tesstransit import vet as vt

from .test_tess import _box


def _flat(seed: int = 1) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    t = np.arange(1325.3, 1352.6, 2.0 / 1440.0)
    return t, 1.0 + rng.normal(0, 4e-4, len(t))


def test_planet_passes_all_three():
    t, f = _flat()
    f = f + _box(t, 3.0, 1326.0, 3.0, 3000.0)
    v = vt.vet(t, f, 3.0, 1326.0, 3.0 / 24)
    assert v.odd_even_sigma < 3 and abs(v.secondary_snr) < 3
    assert v.shape_ratio > 0.85
    assert vt.eb_like(v, 3000.0) == []


def test_binary_at_half_period_fails_odd_even():
    t, f = _flat()
    # primary 8000 ppm and secondary 3000 ppm, 3 d apart: found at P = 3 d
    f = f + _box(t, 6.0, 1326.0, 3.0, 8000.0) + _box(t, 6.0, 1329.0, 3.0, 3000.0)
    v = vt.vet(t, f, 3.0, 1326.0, 3.0 / 24)
    assert v.odd_even_sigma > 10
    assert {round(v.depth_odd_ppm, -3), round(v.depth_even_ppm, -3)} == {8000, 3000}
    assert "odd_even" in vt.eb_like(v, 5500.0)


def test_secondary_eclipse_is_measured_at_phase_half():
    t, f = _flat()
    f = f + _box(t, 4.0, 1326.0, 3.0, 8000.0) + _box(t, 4.0, 1328.0, 3.0, 1500.0)
    v = vt.vet(t, f, 4.0, 1326.0, 3.0 / 24)
    assert v.secondary_snr > 10
    assert v.secondary_ppm == pytest.approx(1500, rel=0.2)
    assert vt.eb_like(v, 8000.0) == ["secondary"]


def test_v_shape_has_a_lower_shape_ratio():
    t, f = _flat()
    phase = (t - 1326.0 + 1.5) % 3.0 - 1.5
    half = 1.5 / 24
    f = f - np.where(np.abs(phase) < half, 8000e-6 * (1 - np.abs(phase) / half), 0.0)
    v = vt.vet(t, f, 3.0, 1326.0, 3.0 / 24)
    assert v.shape_ratio == pytest.approx(2 / 3, abs=0.05)
