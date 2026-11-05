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

from dataclasses import asdict, dataclass, field
from itertools import pairwise
from typing import Any

import numpy as np


@dataclass(frozen=True)
class SearchConfig:
    # How the light curve is flattened before the search.
    #   detrend_method: "median" (running median) or "biweight" (Tukey's
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
    max_signals: int = 3
    mask_width_factor: float = 1.5  # mask this many durations around each transit

    def validate(self) -> None:
        if self.detrend_method not in ("median", "biweight", "spline"):
            raise ValueError(f"unknown detrend_method {self.detrend_method!r}")
        if self.detrend_per_duration and self.detrend_window_factor < 3:
            raise ValueError("detrend_window_factor under 3 would flatten the transits")
        longest = max(self.durations_hours) / 24.0
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
        elif not self.detrend_per_duration and self.detrend_window_days < 3 * longest:
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
        if self.max_signals < 1:
            raise ValueError("max_signals must be >= 1")

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
    # Spline detrending only: the knot spacing BIC chose for each layer.
    knot_days: list[float] = field(default_factory=list)


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
    for window, durs in groups:
        if cfg.detrend_method == "spline":
            trend, window = spline_trend(time, flux, cfg.spline_knot_days)
            flat = flux / trend
        else:
            flat = detrend(time, flux, window, cfg.detrend_method)
        valid = np.isfinite(flat)
        out.append(_Layer(durs, flat, valid, _mad_sigma(flat[valid]), window))
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

    strongest = int(np.argmax(d / sg))
    rest = np.ones(len(d), dtype=bool)
    rest[strongest] = False
    return EventStats(
        depths=d,
        sigmas=sg,
        centres=np.array(centres),
        snr_red=combined(np.ones(len(d), dtype=bool)),
        snr_without_strongest=combined(rest),
        strongest=strongest,
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
    layers = _layers(time, flux, cfg)
    longest = layers[-1]  # the widest window: the closest to the raw light curve
    result = SearchResult(
        tic_id=tic_id,
        sector=sector,
        config=cfg.to_dict(),
        n_points=len(time),
        scatter_ppm=longest.scatter * 1e6,
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
            for sel, allowed in chunks:
                d = np.intersect1d(allowed, layer.durations)
                if len(d) == 0:
                    continue
                r = bls.power(periods[sel], d, objective="likelihood")
                for k in fields:
                    grid[k][li, sel] = np.asarray(getattr(r, k))
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
        if ev.snr_without_strongest < cfg.min_snr_without_strongest:
            centre = float(ev.centres[ev.strongest])
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
                )
            )
            if len(result.vetoes) > cfg.max_vetoes:
                result.stop_reason = "max_vetoes"
                return result
            keep &= np.abs(time - centre) >= cfg.mask_width_factor * dur / 2.0
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
            )
        )
        if len(result.signals) >= cfg.max_signals:
            result.stop_reason = "max_signals"
            return result
        keep &= ~_in_transit(time, p, t0, cfg.mask_width_factor * dur / 2.0)
