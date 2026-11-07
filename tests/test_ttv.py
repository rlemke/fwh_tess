import numpy as np

from tesstransit import ttv


def _curve(amp_min: float, seed: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """80 d at 10-min cadence, a 10 d planet whose times swing by ``amp_min``."""
    rng = np.random.default_rng(seed)
    t = np.arange(1325.0, 1405.0, 10.0 / 1440.0)
    f = 1.0 + rng.normal(0, 3e-4, len(t))
    for n in range(8):
        tc = 1330.0 + 10.0 * n + amp_min / 1440.0 * np.sin(2 * np.pi * n / 6.0)
        f += ttv._trapezoid(t - tc, 3.0 / 24, 2000e-6)
    return t, f


def test_steady_planet_fits_a_fixed_period():
    t, f = _curve(0.0)
    fit = ttv.fit_ttv(ttv.time_transits(t, f, 10.0, 1330.0, 3.0 / 24, 2000e-6))
    assert fit.n_transits == 8 and fit.chi2_dof < 3
    assert abs(fit.period_days - 10.0) < 1e-3


def test_swinging_times_are_flagged():
    t, f = _curve(20.0)
    fit = ttv.fit_ttv(ttv.time_transits(t, f, 10.0, 1330.0, 3.0 / 24, 2000e-6))
    assert fit.chi2_dof > 10 and fit.max_abs_oc_min > 12
