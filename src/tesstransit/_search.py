"""Detrend a light curve, then search it for periodic transits with BLS.

Everything that shapes the result is a field of :class:`SearchConfig`, and
the config travels with the result, so a reported signal can always be traced
to the exact search that produced it.

The search is iterative: find the strongest periodic box, record it if it
clears the threshold, mask its transits, search again. It stops at the first
signal below threshold or at ``max_signals``.

⚠️ Declared limits, each a reason a real transit can be missed:

- **Detrending window.** A running median over ``detrend_window_days``
  flattens anything slower than it, transits included. The window must stay
  well above the longest duration searched (checked in :meth:`validate`).
- **Period grid.** ``max_period_days`` defaults to 13 d because one sector
  spans ~27 d: anything longer shows fewer than two transits, which BLS cannot
  phase. A TCE outside the grid is reported as *outside the search grid*, not
  as a miss.
- **Box model.** BLS fits a flat-bottomed box. Depths come out shallower than
  a limb-darkened fit like SPOC's (measured on π Men c: 259 vs 302 ppm), so
  depth disagreement alone is a method difference, not a defect.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields, replace
from itertools import pairwise
from typing import Any

import numpy as np


@dataclass(frozen=True)
class SearchConfig:
    # How the light curve is flattened before the search.
    #   detrend_method: "median" (running median), "none" (the flux is
    #     already flattened -- the multi-sector search flattens each sector
    #     before stitching) or "biweight" (Tukey's
    #     biweight location via wotan; robust to the transit itself pulling
    #     the trend down -- needs the [wotan] extra).
    #   detrend_per_duration: False = ONE flattening with detrend_window_days
    #     for every trial duration; True = a flattening per trial duration with
    #     window = detrend_window_factor x that duration (floored at
    #     detrend_min_window_days). A short window removes sub-day stellar
    #     variability that a 1-day window leaves (measured: most no-TCE pilot
    #     "signals" were 0.36-1.1 d periodicities), but would flatten a long
    #     transit -- so each duration gets the shortest window that keeps it.
    #   "spline": an iteratively sigma-clipped cubic spline whose knot spacing
    #     is chosen PER STAR by BIC from spline_knot_days (Vanderburg &
    #     Johnson 2014's K2 recipe; Shallue & Vanderburg 2018). A running
    #     median cannot follow curvature: on a 2.85 d, 2% rotation it leaves
    #     ~4000 ppm residuals, and DS Tuc Ab (SPOC SNR 25) measured SNR 0.9
    #     on our flux. Dense knots for an active star, sparse for a quiet one.
    detrend_method: str = "spline"
    spline_knot_days: tuple[float, ...] = (0.5, 0.75, 1.0, 1.5, 2.0, 3.0)
    # Pre-whitening (spline only; 0 = off): before the spline, remove periodic
    # stellar variability faster than prewhiten_max_period_days, which a 0.5 d
    # knot spacing cannot follow -- about 10 of 31 planet misses over sectors
    # 1-13 were stars varying in under a day (TOI 1066.01: 13,793 ppm, forced
    # SNR 2.7 beside a 0.63 d variable). Up to prewhiten_max_terms Lomb-Scargle
    # peaks, each fitted with prewhiten_harmonics harmonics and removed only if
    # significant (false-alarm probability below prewhiten_fap) AND two-sided:
    # a transit or eclipse only dims the star, so a fitted shape whose
    # brightening is under prewhiten_min_symmetry of its dimming is left in --
    # ultra-short-period planets and eclipsing binaries are what we search for.
    # Ten terms down to 0.01 d (was three down to 0.05 d): a delta Scuti star
    # pulsates in many modes, some under an hour, and three terms left most of
    # it in. Over sectors 1-13 that found 112 of 229 injected transits at
    # expected SNR >= 15 on pre-whitened stars (was 88), with the same 682 of
    # 707 planets and no change on other stars (docs/experiments.md).
    prewhiten_max_period_days: float = 0.0
    prewhiten_min_period_days: float = 0.01
    prewhiten_max_terms: int = 10
    prewhiten_harmonics: int = 3
    prewhiten_fap: float = 1e-4
    prewhiten_min_symmetry: float = 0.3
    # Gate and transit protection (added after the first 13-sector run lost
    # 29 planets for 10 gained): only pre-whiten while the spline residual,
    # averaged over an hour, scatters at least prewhiten_min_red_excess times
    # what white noise would (0 = no gate) -- on TOI 396's quiet star 40-50 ppm
    # terms were removed and two shallow planets went under threshold. And
    # with prewhiten_mask_dips, stretches dipping below 3 sigma in 1-hour
    # averages are left out of the fit, so it cannot absorb a transit
    # (TOI 407.01: forced SNR 14 -> 5.8 without it).
    prewhiten_min_red_excess: float = 2.0
    prewhiten_mask_dips: bool = True
    detrend_per_duration: bool = False
    detrend_window_days: float = 1.0
    detrend_window_factor: float = 3.0
    detrend_min_window_days: float = 0.25
    min_period_days: float = 0.25
    max_period_days: float = 26.0
    n_periods: int = 24000  # same log density as 20000 over the earlier 0.25-13 d
    # A periodic search needs at least two transits WITH DATA: one transit
    # fits every period that puts the second in a gap. Long periods are only
    # admitted on these terms (measured: SPOC reports 2-transit TCEs out to
    # ~24 d; the 13 d cap excluded 33 of the 87 pilot TCEs).
    min_transits: int = 2
    max_peaks_tried: int = 500  # candidates checked for min_transits, best first
    # Per-transit consistency. A real transit has a similar depth in every
    # event; an instrumental dip is usually ONE strong event plus noise, which
    # BLS still folds into a "periodic" signal. So the events other than the
    # strongest must, on their own, show the dip at >= this SNR -- with each
    # event's uncertainty from the scatter of the flux binned to the transit
    # duration, which includes red noise (Pont et al. 2006). A vetoed peak's
    # strongest event is masked and the star searched again.
    min_snr_without_strongest: float = 3.0
    # The same test with the TWO strongest events left out (0 = off). The
    # one-event test cannot catch two coincident noise dips: each carries
    # SNR ~5-7, so either alone still "shows" the signal. With 355,000 trial
    # periods over a year of stitched data such pairs, plus a weak third
    # event, made 16% of quiet stars' inverted light curves report a signal
    # at SDE 15 (experiments.md, stitched hold-out). A real transit is in
    # every event, so the rest still show it. Requires >= 3 events in
    # effect. A vetoed peak has both events masked.
    min_snr_without_two_strongest: float = 0.0
    max_vetoes: int = 5
    # Running-median width as a fraction of the period grid, so the detrended
    # SDE does not depend on how many trials the grid has (1001 of 24000).
    power_detrend_frac: float = 0.042
    durations_hours: tuple[float, ...] = (0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0)
    # Each period only tries durations a transit could physically have. For a
    # circular orbit the central-transit duration is ~13 h (P / 1 yr)^(1/3)
    # (rho / rho_sun)^(-1/3), so a floor on stellar density caps the duration.
    # 0.1 rho_sun admits subgiants (cap 3.9 h at 1 d, 9.2 h at 13 d) and
    # excludes giants. Measured reason: on the 40 no-TCE pilot stars, 7 had
    # "signals", nearly all an 8 h box every 0.5-1.1 d (in transit 30-65% of
    # the time) -- stellar variability, not transits.
    min_stellar_density_solar: float = 0.1
    snr_threshold: float = 7.1  # mirrors SPOC's MES threshold; not the same statistic
    # SNR alone is not enough: the best of ~10^7 trial boxes on PURE white
    # noise reached SNR 6.4-7.8 (12 synthetic sectors), so 7.1 alone accepts
    # about half of noise-only stars. SDE is measured on the DETRENDED
    # periodogram (see detrended_power): noise 6.1-10.0, injected 300-1500 ppm
    # transits at 3.7-23 d 41.8-159.6. ⚠️ PROVISIONAL: real TESS noise is red,
    # so calibrate on real non-TCE stars.
    sde_threshold: float = 11.0
    # Red-noise-aware power (experimental, off by default). BLS weighs every
    # point by the white scatter, so on a pulsator the residual variability
    # makes any transit-length box look significant and the periodogram is
    # full of it. With this on, each trial duration's power and SNR are
    # divided by beta^2 and beta, where beta = the scatter of the flux binned
    # to that duration over what white noise predicts (Pont et al. 2006;
    # floored at 1), measured on the points still searched.
    red_noise_power: bool = False
    max_signals: int = 3
    # False-alarm test (Kepler DR25's inverted light curves): search 2 - flux
    # instead of flux. A transit becomes a brightening, which a dip search
    # cannot find, so everything an inverted search reports is a false alarm
    # from noise, systematics or variability -- measured on real light
    # curves rather than synthetic noise. ⚠️ Not perfectly symmetric: a flare
    # inverts into a dip that the normal search never sees, and pre-whitening's
    # dip mask then protects what were brightenings.
    invert_flux: bool = False
    mask_width_factor: float = 1.5  # mask this many durations around each transit

    def validate(self) -> None:
        if self.detrend_method not in ("median", "biweight", "spline", "none"):
            raise ValueError(f"unknown detrend_method {self.detrend_method!r}")
        if self.detrend_per_duration and self.detrend_window_factor < 3:
            raise ValueError("detrend_window_factor under 3 would flatten the transits")
        longest = max(self.durations_hours) / 24.0
        if self.prewhiten_max_period_days and self.detrend_method != "spline":
            raise ValueError("pre-whitening is only implemented for spline detrending")
        if self.prewhiten_max_period_days and not (
            0 < self.prewhiten_min_period_days < self.prewhiten_max_period_days
        ):
            raise ValueError(
                "need 0 < prewhiten_min_period_days < prewhiten_max_period_days"
            )
        if self.detrend_method == "spline" and self.detrend_per_duration:
            raise ValueError("spline detrending has its own layers; no per_duration")
        if self.detrend_method == "spline":
            # A cubic spline bends within ~one knot interval; it has to be
            # longer than the longest transit or it fits the transit too.
            if min(self.spline_knot_days) < 1.5 * longest:
                raise ValueError(
                    f"spline knot spacing {min(self.spline_knot_days)} d is under 1.5x "
                    f"the longest duration searched ({longest:.3f} d)"
                )
        elif (
            self.detrend_method != "none"
            and not self.detrend_per_duration
            and self.detrend_window_days < 3 * longest
        ):
            raise ValueError(
                f"detrend window {self.detrend_window_days} d is under 3x the longest "
                f"duration searched ({longest:.3f} d); it would flatten those transits"
            )
        if self.min_stellar_density_solar <= 0:
            raise ValueError("min_stellar_density_solar must be positive")
        if not 0 < self.min_period_days < self.max_period_days:
            raise ValueError("need 0 < min_period_days < max_period_days")
        if self.min_transits < 1:
            raise ValueError("min_transits must be >= 1")
        if self.min_snr_without_two_strongest < 0:
            raise ValueError("min_snr_without_two_strongest must be >= 0 (0 = off)")
        if self.max_signals < 1:
            raise ValueError("max_signals must be >= 1")

    def override(self, pairs: list[str]) -> SearchConfig:
        """A copy with ``FIELD=VALUE`` settings applied (the CLIs' ``--set``), validated.

        Values are parsed by the field's current type; a tuple takes a comma list.
        """
        types = {f.name: type(getattr(self, f.name)) for f in fields(self)}
        kw: dict[str, Any] = {}
        for pair in pairs:
            name, sep, value = pair.partition("=")
            if not sep or name not in types:
                raise ValueError(
                    f"--set {pair!r}: expected FIELD=VALUE, FIELD a SearchConfig field"
                )
            kind = types[name]
            if kind is tuple:
                kw[name] = tuple(float(x) for x in value.split(","))
            elif kind is bool:
                if value not in ("true", "false"):
                    raise ValueError(f"--set {pair!r}: a boolean is true or false")
                kw[name] = value == "true"
            else:
                kw[name] = kind(value)
        out = replace(self, **kw)
        out.validate()
        return out

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["durations_hours"] = list(self.durations_hours)
        return d


@dataclass
class Signal:
    """One periodic transit-like signal. Times BTJD, period/epoch in days."""

    rank: int  # 1 = found first (strongest)
    period_days: float
    epoch_btjd: float
    duration_hours: float
    depth_ppm: float
    snr: float  # depth / depth uncertainty, from BLS with the measured scatter
    sde: float  # (peak - mean) / std of the periodogram
    n_transits: int  # transits with in-transit data
    snr_red: float | None = None  # combined per-event SNR, red-noise uncertainties
    snr_without_strongest: float | None = None  # the same, strongest event left out
    snr_without_two_strongest: float | None = None  # ... two strongest left out


@dataclass
class Veto:
    """A peak that cleared the thresholds but failed per-transit consistency."""

    period_days: float
    epoch_btjd: float
    duration_hours: float
    snr: float
    sde: float
    snr_red: float
    snr_without_strongest: float
    strongest_event_btjd: float  # masked before searching again
    # "without_strongest" | "without_two_strongest"; the second event is
    # masked too in the second case.
    reason: str = "without_strongest"
    snr_without_two_strongest: float | None = None
    second_event_btjd: float | None = None


@dataclass
class SearchResult:
    tic_id: int
    sector: int
    config: dict[str, Any]
    n_points: int
    scatter_ppm: float  # robust (MAD) scatter of the detrended flux
    signals: list[Signal] = field(default_factory=list)
    # "below_threshold" | "max_signals" | "too_few_points" | "no_multi_transit_peak"
    # | "max_vetoes"
    stop_reason: str = ""
    # The strongest peak that did NOT clear the thresholds: how close the search
    # came, and WHERE -- so a missed TCE can be told apart as "found at its
    # period but rejected" from "never the strongest peak".
    final_snr: float | None = None
    final_sde: float | None = None
    final_period_days: float | None = None
    final_epoch_btjd: float | None = None
    vetoes: list[Veto] = field(default_factory=list)
    # Spline detrending only: the knot spacing BIC chose for each layer, and
    # the fast periodicities pre-whitening divided out before it.
    knot_days: list[float] = field(default_factory=list)
    prewhitened: list[Prewhitened] = field(default_factory=list)


def running_median(
    time: np.ndarray, flux: np.ndarray, window_days: float
) -> np.ndarray:
    """Median of ``flux`` within ±window/2 of each point, by TIME, not index.

    By time so that a data gap does not pull distant points into the window.
    """
    half = window_days / 2.0
    lo = np.searchsorted(time, time - half, side="left")
    hi = np.searchsorted(time, time + half, side="right")
    return np.array([np.median(flux[a:b]) for a, b in zip(lo, hi)])


def _segments(time: np.ndarray, gap_days: float) -> list[slice]:
    """Contiguous stretches of data, split where the cadence stops for > gap."""
    cuts = np.flatnonzero(np.diff(time) > gap_days) + 1
    edges = [0, *cuts.tolist(), len(time)]
    return [slice(a, b) for a, b in pairwise(edges) if b > a]


def _spline_fit(
    t: np.ndarray, f: np.ndarray, knot_days: float, keep: np.ndarray
) -> np.ndarray:
    """Least-squares cubic spline through the ``keep`` points, evaluated at ``t``."""
    from scipy.interpolate import LSQUnivariateSpline

    tk, fk = t[keep], f[keep]
    n_int = int((tk[-1] - tk[0]) // knot_days)
    if n_int < 1 or len(tk) < 20:
        # Too short for a spline at this spacing: a quadratic, as a spline
        # with no interior knots but robust to a 3-point segment.
        deg = min(2, max(0, len(tk) - 1))
        return np.polyval(np.polyfit(tk - tk[0], fk, deg), t - tk[0])
    knots = tk[0] + knot_days * np.arange(1, n_int + 1)
    knots = knots[knots < tk[-1] - 0.5 * knot_days]
    # Schoenberg-Whitney: every knot interval must hold data, else the fit is
    # singular. Drop knots in gaps rather than fail.
    counts = np.histogram(tk, bins=[tk[0], *knots, tk[-1]])[0]
    while len(knots) and (counts < 4).any():
        bad = int(np.argmin(counts))
        knots = np.delete(knots, min(bad, len(knots) - 1))
        counts = np.histogram(tk, bins=[tk[0], *knots, tk[-1]])[0]
    return LSQUnivariateSpline(tk, fk, knots, k=3)(t)


def spline_trend(
    time: np.ndarray,
    flux: np.ndarray,
    knot_days: tuple[float, ...],
    clip_sigma: float = 3.0,
    max_iter: int = 5,
    gap_days: float = 0.5,
) -> tuple[np.ndarray, float]:
    """``(trend, knot spacing chosen)``: the clipped spline with the lowest BIC.

    Each candidate spacing is fitted per data segment, clipping points more
    than ``clip_sigma`` robust sigmas from the fit and refitting, so a deep
    transit does not drag the trend down. BIC = chi^2 + n_params ln N with
    the WHITE-noise sigma from point-to-point differences, so a spacing only
    wins by explaining variability the noise cannot (Vanderburg+ 2016).
    """
    segs = _segments(time, gap_days)
    sigma_w = 1.4826 * np.median(np.abs(np.diff(flux))) / np.sqrt(2.0)
    best: tuple[float, np.ndarray, float] | None = None
    for kd in knot_days:
        trend = np.empty_like(flux)
        keep_all = np.ones(len(flux), dtype=bool)
        n_params = 0
        for sl in segs:
            t, f = time[sl], flux[sl]
            keep = np.ones(len(t), dtype=bool)
            for _ in range(max_iter):
                fit = _spline_fit(t, f, kd, keep)
                r = f - fit
                sig = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep])))
                new = np.abs(r) < clip_sigma * sig
                if np.array_equal(new, keep) or new.sum() < 20:
                    break
                keep = new
            trend[sl], keep_all[sl] = fit, keep
            n_params += int((t[-1] - t[0]) // kd) + 4
        resid = (flux - trend)[keep_all]
        bic = float(np.sum(resid**2) / sigma_w**2 + n_params * np.log(keep_all.sum()))
        if best is None or bic < best[0]:
            best = (bic, trend, kd)
    assert best is not None
    return best[1], best[2]


@dataclass
class Prewhitened:
    period_days: float
    semi_amplitude_ppm: float  # half the peak-to-peak of the folded, binned data
    symmetry: float  # brightening / dimming of the fitted shape (1 = sinusoid)


def _harmonic_fit(
    t: np.ndarray, r: np.ndarray, keep: np.ndarray, freq: float, n_harm: int
) -> np.ndarray:
    """Least-squares Fourier series at ``freq`` through ``keep`` points, at all ``t``."""
    cols = [np.ones_like(t)]
    for k in range(1, n_harm + 1):
        w = 2 * np.pi * k * freq * t
        cols += [np.cos(w), np.sin(w)]
    a = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(a[keep], r[keep], rcond=None)
    return a @ coef - coef[0]


def prewhiten(
    time: np.ndarray, flux: np.ndarray, cfg: SearchConfig
) -> tuple[np.ndarray, list[Prewhitened]]:
    """``flux`` with significant, two-sided fast periodicities divided out."""
    from astropy.timeseries import LombScargle

    f = flux.copy()
    removed: list[Prewhitened] = []
    cadence = float(np.median(np.diff(time)))
    hour = 1.0 / 24.0
    sigma_w = 1.4826 * np.median(np.abs(np.diff(flux))) / np.sqrt(2.0) / np.median(flux)
    white_hour = sigma_w / np.sqrt(max(1.0, hour / cadence))
    for _ in range(cfg.prewhiten_max_terms):
        r = f / spline_trend(time, f, cfg.spline_knot_days)[0] - 1.0
        if cfg.prewhiten_min_red_excess > 0:
            excess = binned_scatter(time, r, hour) / white_hour
            if excess < cfg.prewhiten_min_red_excess:
                break  # what is left is close to white noise: nothing to model
        ls = LombScargle(time, r)
        freq, power = ls.autopower(
            minimum_frequency=1.0 / cfg.prewhiten_max_period_days,
            maximum_frequency=1.0 / cfg.prewhiten_min_period_days,
            samples_per_peak=10,
        )
        i = int(np.argmax(power))
        fap = float(
            ls.false_alarm_probability(
                power[i],
                minimum_frequency=1.0 / cfg.prewhiten_max_period_days,
                maximum_frequency=1.0 / cfg.prewhiten_min_period_days,
                method="baluev",
            )
        )
        if fap > cfg.prewhiten_fap:
            break
        # Two-sidedness is measured on the DATA folded at the period: the
        # brightest and faintest of 50 phase bins against their median. A
        # transit train dips without brightening; variability does both. (A
        # 3-harmonic fitted shape rings: a pure 0.75 d, 1 h transit read as
        # 0.43 two-sided and was removed -- measured, synthetic.)
        period = 1.0 / float(freq[i])
        if any(abs(period / p.period_days - 1) < 1e-3 for p in removed):
            break  # the same peak again: nothing more this term can remove
        phase = ((time - time[0]) / period) % 1.0
        idx = np.minimum((phase * 50).astype(int), 49)
        counts = np.bincount(idx, minlength=50)
        binned = (
            np.bincount(idx, weights=r, minlength=50)[counts > 0] / counts[counts > 0]
        )
        mid = float(np.median(binned))
        up, down = float(binned.max()) - mid, mid - float(binned.min())
        symmetry = up / down if down > 0 else 0.0
        if (
            min(symmetry, 1 / symmetry if symmetry else 0.0)
            < cfg.prewhiten_min_symmetry
        ):
            break  # one-sided: transit- or eclipse-like, not stellar variability
        # Subtract a fit without points far BELOW the rest, so a transit
        # inside the variability does not pull the removed shape down.
        sig = 1.4826 * np.median(np.abs(r - np.median(r)))
        keep = r > np.median(r) - 3 * sig
        model = _harmonic_fit(time, r, keep, float(freq[i]), cfg.prewhiten_harmonics)
        if cfg.prewhiten_mask_dips:
            # A transit is a STRETCH of low points, each within the noise when
            # the star is noisy: find dips in 1-hour running means of what the
            # first fit leaves, and refit without them (twice).
            half = max(1, round(0.5 * hour / cadence))
            kernel = np.ones(2 * half + 1) / (2 * half + 1)
            for _ in range(2):
                smooth = np.convolve(r - model, kernel, mode="same")
                s_sig = 1.4826 * np.median(np.abs(smooth - np.median(smooth)))
                dip = smooth < np.median(smooth) - 3 * s_sig
                dip = np.convolve(dip, np.ones(2 * half + 1), mode="same") > 0  # widen
                model = _harmonic_fit(
                    time, r, keep & ~dip, float(freq[i]), cfg.prewhiten_harmonics
                )
        f = f / (1.0 + model)
        removed.append(
            Prewhitened(1.0 / float(freq[i]), 0.5 * (up + down) * 1e6, symmetry)
        )
    return f, removed


def flatten(
    time: np.ndarray, flux: np.ndarray, cfg: SearchConfig
) -> tuple[np.ndarray, float | None, list[Prewhitened]]:
    """``(flat, knot spacing or None, removed periodicities)`` with ``cfg``'s method."""
    if cfg.detrend_method != "spline":
        return (
            detrend(time, flux, cfg.detrend_window_days, cfg.detrend_method),
            None,
            [],
        )
    removed: list[Prewhitened] = []
    if cfg.prewhiten_max_period_days > 0:
        flux, removed = prewhiten(time, flux, cfg)
    trend, knot = spline_trend(time, flux, cfg.spline_knot_days)
    return flux / trend, knot, removed


def detrend(
    time: np.ndarray,
    flux: np.ndarray,
    window_days: float,
    method: str = "median",
    knot_days: tuple[float, ...] = (1.0,),
) -> np.ndarray:
    """``flux`` divided by its trend. May contain NaN (biweight, sparse edges)."""
    if np.any(np.diff(time) < 0):
        raise ValueError("time must be sorted ascending")
    if method == "none":  # already flattened (multisector: per sector, then stitched)
        return np.array(flux, dtype=float)
    if method == "spline":
        return flux / spline_trend(time, flux, knot_days)[0]
    if method == "median":
        return flux / running_median(time, flux, window_days)
    if method == "biweight":
        from wotan import flatten

        return np.asarray(
            flatten(
                time, flux, method="biweight", window_length=window_days, edge_cutoff=0
            )
        )
    raise ValueError(f"unknown detrend method {method!r}")


@dataclass
class _Layer:
    """One flattening of the light curve, and the trial durations searched on it."""

    durations: np.ndarray  # days
    flat: np.ndarray
    valid: np.ndarray  # finite after flattening
    scatter: float
    window_days: float
    prewhitened: list[Prewhitened] = field(default_factory=list)


def _layers(time: np.ndarray, flux: np.ndarray, cfg: SearchConfig) -> list[_Layer]:
    durations = np.sort(np.array(cfg.durations_hours) / 24.0)
    if not cfg.detrend_per_duration:
        groups = [(cfg.detrend_window_days, durations)]
    else:
        groups = [
            (
                max(cfg.detrend_min_window_days, cfg.detrend_window_factor * d),
                np.array([d]),
            )
            for d in durations
        ]
    out = []
    removed: list[Prewhitened] = []
    for window, durs in groups:
        if cfg.detrend_method == "spline":
            flat, window, removed = flatten(time, flux, cfg)
        else:
            flat = detrend(time, flux, window, cfg.detrend_method)
        valid = np.isfinite(flat)
        out.append(_Layer(durs, flat, valid, _mad_sigma(flat[valid]), window, removed))
    return out


def _mad_sigma(x: np.ndarray) -> float:
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def max_duration_days(period_days: np.ndarray | float, min_density_solar: float) -> Any:
    """Longest physically plausible transit at ``period_days`` (circular, b=0)."""
    hours = (
        13.0
        * (np.asarray(period_days) / 365.25) ** (1 / 3)
        * min_density_solar ** (-1 / 3)
    )
    return hours / 24.0


def _period_chunks(
    periods: np.ndarray, durations: np.ndarray, min_density: float
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Split the period grid so each chunk tries only the durations allowed there.

    Returns ``(mask over periods, durations allowed)`` per chunk.

    The allowed set only grows with period, so the chunks are contiguous and
    every chunk's durations are shorter than its periods (BLS requires it).
    """
    cap = max_duration_days(periods, min_density)
    n_allowed = np.searchsorted(np.sort(durations), cap, side="right")
    chunks = []
    for n in np.unique(n_allowed):
        if n == 0:
            continue  # no searched duration is short enough at these periods
        chunks.append((n_allowed == n, np.sort(durations)[:n]))
    return chunks


def _in_transit(
    time: np.ndarray, period: float, epoch: float, half_width: float
) -> np.ndarray:
    phase = (time - epoch + 0.5 * period) % period - 0.5 * period
    return np.abs(phase) < half_width


def detrended_power(power: np.ndarray, window: int) -> np.ndarray:
    """BLS power minus its running median over the (log-spaced) grid, ÷ robust std.

    Raw BLS power rises with period: few transits fit noise better, and one
    deep dip fits EVERY long period that puts the next transit in a gap. Left
    in, that plateau inflated the mean and std a plain SDE is measured
    against, and a real 2-transit planet at SNR 49 scored SDE 6.5 (measured,
    synthetic 18 d). The running median is computed on every 10th trial and
    interpolated — the trend is smooth, the peaks are not.
    """
    step = 10
    coarse_x = np.arange(0, len(power), step)
    half = max(1, window // (2 * step))
    coarse = power[::step]
    lo = np.clip(np.arange(len(coarse)) - half, 0, None)
    hi = np.clip(np.arange(len(coarse)) + half + 1, None, len(coarse))
    med = np.array([np.median(coarse[a:b]) for a, b in zip(lo, hi)])
    resid = power - np.interp(np.arange(len(power)), coarse_x, med)
    return resid / (1.4826 * np.median(np.abs(resid - np.median(resid))))


def _n_transits(t: np.ndarray, period: float, epoch: float, duration: float) -> int:
    """Distinct transits with at least one in-transit cadence."""
    in_tr = _in_transit(t, period, epoch, duration / 2.0)
    return len(np.unique(np.round((t[in_tr] - epoch) / period)))


@dataclass
class EventStats:
    depths: np.ndarray  # per-event depth (fraction), one per transit with data
    sigmas: np.ndarray  # per-event uncertainty (red-noise aware)
    centres: np.ndarray  # per-event mid-time, BTJD
    snr_red: float
    snr_without_strongest: float
    strongest: int  # index into the arrays above
    snr_without_two_strongest: float = 0.0
    second: int | None = None  # index of the second-strongest event


def binned_scatter(t: np.ndarray, f: np.ndarray, width: float) -> float:
    """Robust std of ``f`` averaged in time bins of ``width`` (half-full bins kept).

    For white noise this is sigma/sqrt(n per bin); red noise makes it larger,
    which is exactly the uncertainty a transit-length average really has.
    """
    idx = np.floor((t - t[0]) / width).astype(int)
    counts = np.bincount(idx)
    sums = np.bincount(idx, weights=f)
    full = counts >= 0.5 * counts.max()
    means = sums[full] / counts[full]
    if len(means) < 5:
        raise ValueError(f"only {len(means)} usable bins of {width:.3f} d")
    return float(1.4826 * np.median(np.abs(means - np.median(means))))


def red_noise_beta(
    t: np.ndarray, f: np.ndarray, durations: np.ndarray, scatter: float
) -> dict[float, float]:
    """Per duration: binned scatter over the white-noise prediction, floored at 1."""
    cadence = float(np.median(np.diff(t)))
    out = {}
    for d in durations:
        try:
            red = binned_scatter(t, f, float(d))
        except ValueError:
            red = 0.0
        white = scatter / np.sqrt(max(1.0, float(d) / cadence))
        out[float(d)] = max(1.0, red / white)
    return out


def event_stats(
    t: np.ndarray, f: np.ndarray, period: float, epoch: float, duration: float
) -> EventStats:
    """Per-transit depths and the consistency statistics built from them."""
    phase = (t - epoch + 0.5 * period) % period - 0.5 * period
    in_tr = np.abs(phase) < duration / 2.0
    sigma_dur = binned_scatter(t[~in_tr], f[~in_tr], duration)
    k = np.round((t - epoch) / period).astype(int)
    events = np.unique(k[in_tr])
    cadence = float(np.median(np.diff(t)))
    n_full = max(1.0, duration / cadence)
    # Depth against the out-of-transit MEAN, as BLS measures it -- not against
    # 1.0. After a signal is masked, a variable star's remaining flux does not
    # sit at 1, and measuring from 1 scored a strong EB peak (TIC 149308175)
    # at per-event SNR -16 and vetoed it.
    baseline = float(np.mean(f[~in_tr]))
    depths, sigmas, centres = [], [], []
    for e in events:
        sel = in_tr & (k == e)
        n = sel.sum()
        depths.append(baseline - float(np.mean(f[sel])))
        sigmas.append(sigma_dur * np.sqrt(max(1.0, n_full / n)))
        centres.append(epoch + e * period)
    d, sg = np.array(depths), np.array(sigmas)

    def combined(mask: np.ndarray) -> float:
        if not mask.any():
            return 0.0
        w = 1.0 / sg[mask] ** 2
        return float(np.sum(d[mask] * w) / np.sum(w) * np.sqrt(np.sum(w)))

    order = np.argsort(d / sg)[::-1]
    strongest = int(order[0])
    rest = np.ones(len(d), dtype=bool)
    rest[strongest] = False
    rest2 = rest.copy()
    second = int(order[1]) if len(order) > 1 else None
    if second is not None:
        rest2[second] = False
    return EventStats(
        depths=d,
        sigmas=sg,
        centres=np.array(centres),
        snr_red=combined(np.ones(len(d), dtype=bool)),
        snr_without_strongest=combined(rest),
        strongest=strongest,
        snr_without_two_strongest=combined(rest2),
        second=second,
    )


def _best_with_transits(
    t: np.ndarray,
    power: np.ndarray,
    periods: np.ndarray,
    epochs: np.ndarray,
    durations: np.ndarray,
    min_transits: int,
    max_tries: int,
) -> int | None:
    """Index of the strongest peak showing ``min_transits`` transits with data."""
    for i in np.argsort(power)[::-1][:max_tries]:
        if _n_transits(t, periods[i], epochs[i], durations[i]) >= min_transits:
            return int(i)
    return None


def search(
    time: np.ndarray,
    flux: np.ndarray,
    tic_id: int,
    sector: int,
    config: SearchConfig | None = None,
) -> SearchResult:
    """Detrend ``flux`` (normalised to ~1) and run the iterative BLS search."""
    from astropy.timeseries import BoxLeastSquares

    cfg = config or SearchConfig()
    cfg.validate()
    if cfg.invert_flux:
        flux = 2.0 - np.asarray(flux)
    layers = _layers(time, flux, cfg)
    longest = layers[-1]  # the widest window: the closest to the raw light curve
    result = SearchResult(
        tic_id=tic_id,
        sector=sector,
        config=cfg.to_dict(),
        n_points=len(time),
        scatter_ppm=longest.scatter * 1e6,
        prewhitened=list(longest.prewhitened),
        knot_days=[ly.window_days for ly in layers]
        if cfg.detrend_method == "spline"
        else [],
    )
    periods = np.exp(
        np.linspace(
            np.log(cfg.min_period_days), np.log(cfg.max_period_days), cfg.n_periods
        )
    )
    durations = np.array(cfg.durations_hours) / 24.0
    chunks = _period_chunks(periods, durations, cfg.min_stellar_density_solar)
    keep = np.ones(len(time), dtype=bool)
    fields = ("power", "transit_time", "duration", "depth", "depth_snr")

    while True:
        if keep.sum() < 100:
            result.stop_reason = "too_few_points"
            return result
        # One periodogram per layer, on the period grid; -inf where a layer
        # has no duration allowed at that period.
        grid = {k: np.full((len(layers), len(periods)), -np.inf) for k in fields}
        for li, layer in enumerate(layers):
            use = keep & layer.valid
            bls = BoxLeastSquares(
                time[use], layer.flat[use], dy=np.full(use.sum(), layer.scatter)
            )
            beta = (
                red_noise_beta(time[use], layer.flat[use], layer.durations, layer.scatter)
                if cfg.red_noise_power
                else {}
            )  # fmt: skip
            for sel, allowed in chunks:
                d = np.intersect1d(allowed, layer.durations)
                if len(d) == 0:
                    continue
                if not beta:
                    r = bls.power(periods[sel], d, objective="likelihood")
                    for k in fields:
                        grid[k][li, sel] = np.asarray(getattr(r, k))
                    continue
                # One duration at a time, so each is scaled by its own beta
                # before the best duration at a period is chosen.
                idx = np.flatnonzero(sel) if sel.dtype == bool else np.asarray(sel)
                for dd in d:
                    r = bls.power(periods[sel], [dd], objective="likelihood")
                    bt = beta[float(dd)]
                    pw = np.asarray(r.power) / bt**2
                    better = pw > grid["power"][li, idx]
                    grid["power"][li, idx[better]] = pw[better]
                    grid["depth_snr"][li, idx[better]] = (
                        np.asarray(r.depth_snr)[better] / bt
                    )
                    for k in ("transit_time", "duration", "depth"):
                        grid[k][li, idx[better]] = np.asarray(getattr(r, k))[better]
        # Likelihood power ~ SNR^2 / 2 with each layer's own scatter, so the
        # layers are comparable: keep the best layer at each period.
        best_layer = np.argmax(grid["power"], axis=0)
        cols = np.arange(len(periods))
        res = {k: grid[k][best_layer, cols] for k in fields}
        raw = res["power"]
        raw[~np.isfinite(raw)] = np.min(raw[np.isfinite(raw)])
        power = detrended_power(
            raw, max(3, round(cfg.power_detrend_frac * cfg.n_periods))
        )
        i = _best_with_transits(
            time[keep], power, periods, res["transit_time"], res["duration"],
            cfg.min_transits, cfg.max_peaks_tried,
        )  # fmt: skip
        if i is None:
            result.stop_reason = "no_multi_transit_peak"
            return result
        layer = layers[best_layer[i]]
        snr = float(res["depth_snr"][i])
        sde = float(power[i])
        p, t0, dur = (
            float(periods[i]),
            float(res["transit_time"][i]),
            float(res["duration"][i]),
        )
        if (
            snr < cfg.snr_threshold
            or sde < cfg.sde_threshold
            or not res["depth"][i] > 0
        ):
            result.stop_reason = "below_threshold"
            result.final_snr, result.final_sde = snr, sde
            result.final_period_days, result.final_epoch_btjd = p, t0
            return result
        use = keep & layer.valid
        ev = event_stats(time[use], layer.flat[use], p, t0, dur)
        fail1 = ev.snr_without_strongest < cfg.min_snr_without_strongest
        fail2 = (
            cfg.min_snr_without_two_strongest > 0
            and ev.snr_without_two_strongest < cfg.min_snr_without_two_strongest
        )
        if fail1 or fail2:
            centre = float(ev.centres[ev.strongest])
            second = (
                float(ev.centres[ev.second])
                if fail2 and not fail1 and ev.second is not None
                else None
            )
            result.vetoes.append(
                Veto(
                    p,
                    t0,
                    dur * 24.0,
                    snr,
                    sde,
                    ev.snr_red,
                    ev.snr_without_strongest,
                    centre,
                    reason="without_strongest" if fail1 else "without_two_strongest",
                    snr_without_two_strongest=ev.snr_without_two_strongest,
                    second_event_btjd=second,
                )
            )
            if len(result.vetoes) > cfg.max_vetoes:
                result.stop_reason = "max_vetoes"
                return result
            for c in (centre, second):
                if c is not None:
                    keep &= np.abs(time - c) >= cfg.mask_width_factor * dur / 2.0
            continue
        result.signals.append(
            Signal(
                rank=len(result.signals) + 1,
                period_days=p,
                epoch_btjd=t0,
                duration_hours=dur * 24.0,
                depth_ppm=float(res["depth"][i]) * 1e6,
                snr=snr,
                sde=sde,
                n_transits=len(ev.depths),
                snr_red=ev.snr_red,
                snr_without_strongest=ev.snr_without_strongest,
                snr_without_two_strongest=ev.snr_without_two_strongest,
            )
        )
        if len(result.signals) >= cfg.max_signals:
            result.stop_reason = "max_signals"
            return result
        keep &= ~_in_transit(time, p, t0, cfg.mask_width_factor * dur / 2.0)
