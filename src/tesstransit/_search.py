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
from typing import Any

import numpy as np


@dataclass(frozen=True)
class SearchConfig:
    detrend_window_days: float = 1.0
    min_period_days: float = 0.25
    max_period_days: float = 13.0
    n_periods: int = 20000
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
    # SNR alone is not enough: the peak of ~10^7 trial boxes on PURE white noise
    # reached SNR 6.4-7.8 (12 synthetic sectors), so 7.1 alone would accept
    # about half of all noise-only stars. Their SDE was 3.6-5.1; injected
    # 300-1500 ppm transits gave 12.5-14.3. 7.0 sits in that gap. ⚠️
    # PROVISIONAL: real TESS noise is red, so calibrate on real non-TCE stars.
    sde_threshold: float = 7.0
    max_signals: int = 3
    mask_width_factor: float = 1.5  # mask this many durations around each transit

    def validate(self) -> None:
        longest = max(self.durations_hours) / 24.0
        if self.detrend_window_days < 3 * longest:
            raise ValueError(
                f"detrend window {self.detrend_window_days} d is under 3x the longest "
                f"duration searched ({longest:.3f} d); it would flatten those transits"
            )
        if self.min_stellar_density_solar <= 0:
            raise ValueError("min_stellar_density_solar must be positive")
        if not 0 < self.min_period_days < self.max_period_days:
            raise ValueError("need 0 < min_period_days < max_period_days")
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


@dataclass
class SearchResult:
    tic_id: int
    sector: int
    config: dict[str, Any]
    n_points: int
    scatter_ppm: float  # robust (MAD) scatter of the detrended flux
    signals: list[Signal] = field(default_factory=list)
    stop_reason: str = ""  # "below_threshold" | "max_signals" | "too_few_points"
    final_snr: float | None = None  # the peak that did NOT clear the thresholds
    final_sde: float | None = None


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


def detrend(time: np.ndarray, flux: np.ndarray, window_days: float) -> np.ndarray:
    if np.any(np.diff(time) < 0):
        raise ValueError("time must be sorted ascending")
    return flux / running_median(time, flux, window_days)


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

    The allowed set only grows with period, so the chunks are contiguous and
    every chunk's durations are shorter than its periods (BLS requires it).
    """
    cap = max_duration_days(periods, min_density)
    n_allowed = np.searchsorted(np.sort(durations), cap, side="right")
    chunks = []
    for n in np.unique(n_allowed):
        if n == 0:
            continue  # no searched duration is short enough at these periods
        sel = n_allowed == n
        chunks.append((periods[sel], np.sort(durations)[:n]))
    return chunks


def _in_transit(
    time: np.ndarray, period: float, epoch: float, half_width: float
) -> np.ndarray:
    phase = (time - epoch + 0.5 * period) % period - 0.5 * period
    return np.abs(phase) < half_width


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
    flat = detrend(time, flux, cfg.detrend_window_days)
    scatter = _mad_sigma(flat)
    result = SearchResult(
        tic_id=tic_id,
        sector=sector,
        config=cfg.to_dict(),
        n_points=len(time),
        scatter_ppm=scatter * 1e6,
    )
    periods = np.exp(
        np.linspace(
            np.log(cfg.min_period_days), np.log(cfg.max_period_days), cfg.n_periods
        )
    )
    durations = np.array(cfg.durations_hours) / 24.0
    chunks = _period_chunks(periods, durations, cfg.min_stellar_density_solar)
    keep = np.ones(len(time), dtype=bool)

    while True:
        if keep.sum() < 100:
            result.stop_reason = "too_few_points"
            return result
        t, f = time[keep], flat[keep]
        dy = np.full(len(t), scatter)
        bls = BoxLeastSquares(t, f, dy=dy)
        parts = [bls.power(p, d, objective="likelihood") for p, d in chunks]
        res = {
            k: np.concatenate([np.asarray(getattr(r, k)) for r in parts])
            for k in (
                "power",
                "period",
                "transit_time",
                "duration",
                "depth",
                "depth_snr",
            )
        }
        power = res["power"]
        i = int(np.argmax(power))
        snr = float(res["depth_snr"][i])
        sde = float((power[i] - power.mean()) / power.std())
        if (
            snr < cfg.snr_threshold
            or sde < cfg.sde_threshold
            or not res["depth"][i] > 0
        ):
            result.stop_reason = "below_threshold"
            result.final_snr = snr
            result.final_sde = sde
            return result
        p, t0, dur = (
            float(res["period"][i]),
            float(res["transit_time"][i]),
            float(res["duration"][i]),
        )
        in_tr = _in_transit(t, p, t0, dur / 2.0)
        epochs_seen = np.unique(np.round((t[in_tr] - t0) / p))
        result.signals.append(
            Signal(
                rank=len(result.signals) + 1,
                period_days=p,
                epoch_btjd=t0,
                duration_hours=dur * 24.0,
                depth_ppm=float(res["depth"][i]) * 1e6,
                snr=snr,
                sde=sde,
                n_transits=len(epochs_seen),
            )
        )
        if len(result.signals) >= cfg.max_signals:
            result.stop_reason = "max_signals"
            return result
        keep &= ~_in_transit(time, p, t0, cfg.mask_width_factor * dur / 2.0)
