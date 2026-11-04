"""SPOC's TCE table, and whether one of our signals is the same thing as one of its TCEs.

The reference is the per-sector DV statistics table MAST publishes, e.g.::

    https://archive.stsci.edu/missions/tess/catalogs/tce/
        tess2018206190142-s0001-s0001_dvr-tcestats.csv

(sector 1: 1,457 TCEs on 971 stars). It opens with ``#`` comment lines. Units
are mixed and easy to get wrong: period in DAYS, ``tce_time0bt`` in BTJD,
duration in **HOURS**, depth in ppm.

Every TCE and every one of our signals lands in exactly ONE bucket:

=====================  ======================================================
``recovered``          a signal at the TCE's period and epoch
``recovered_alias``    a signal at 2x, 1/2x, 3x or 1/3x the TCE's period whose
                       transits coincide with the TCE's (``alias`` records k)
``missed``             a TCE with no matching signal; ``cause`` says what we
                       KNOW, never a guess (see below)
``unmatched``          our signal with no TCE. UNVETTED: most will be
                       systematics or eclipsing binaries, never "new planets"
=====================  ======================================================

Miss causes, deliberately limited to what the search can establish:

- ``not_searched`` — the star has no search result at all
- ``outside_search_grid`` — the TCE's period lies outside our period grid
- ``search_exhausted`` — we stopped at ``max_signals`` before reaching it
- ``not_detected`` — searched, inside the grid, nothing matched. This does NOT
  say why (below threshold, removed by detrending, masked by quality flags);
  telling those apart needs injection tests, not a label.
"""

from __future__ import annotations

import csv
import io
from collections.abc import Iterable
from dataclasses import dataclass

from ._search import SearchResult, Signal

DEFAULT_ALIASES: tuple[float, ...] = (1.0, 2.0, 0.5, 3.0, 1.0 / 3.0)


@dataclass(frozen=True)
class TCE:
    tce_id: str
    tic_id: int
    planet_num: int
    period_days: float
    epoch_btjd: float
    duration_hours: float
    depth_ppm: float | None
    snr: float | None  # SPOC's model-fit SNR (tce_model_snr)
    n_transits: int | None


def _f(v: str | None) -> float | None:
    v = (v or "").strip()
    return float(v) if v else None


def _i(v: str | None) -> int | None:
    x = _f(v)
    return int(x) if x is not None else None


def parse_tce_csv(text: str) -> list[TCE]:
    """Parse a ``*_dvr-tcestats.csv`` body. Rows missing period/epoch/duration raise."""
    lines = (ln for ln in io.StringIO(text) if not ln.startswith("#"))
    out: list[TCE] = []
    for row in csv.DictReader(lines):
        period = _f(row.get("tce_period"))
        epoch = _f(row.get("tce_time0bt"))
        dur = _f(row.get("tce_duration"))
        if period is None or epoch is None or dur is None:
            raise ValueError(f"TCE {row.get('tceid')!r} lacks period/epoch/duration")
        out.append(
            TCE(
                tce_id=row["tceid"],
                tic_id=int(row["ticid"]),
                planet_num=int(row["tce_plnt_num"]),
                period_days=period,
                epoch_btjd=epoch,
                duration_hours=dur,
                depth_ppm=_f(row.get("tce_depth")),
                snr=_f(row.get("tce_model_snr")),
                n_transits=_i(row.get("tce_num_transits")),
            )
        )
    return out


@dataclass(frozen=True)
class MatchConfig:
    period_rel_tol: float = 0.005
    epoch_tol_durations: float = 0.5  # in units of the longer of the two durations
    aliases: tuple[float, ...] = DEFAULT_ALIASES


@dataclass
class MatchRow:
    tic_id: int
    bucket: str
    tce: TCE | None = None
    signal: Signal | None = None
    alias: float | None = None  # our period / TCE period, when matched
    cause: str | None = None  # for "missed"


def _alias_of(sig: Signal, tce: TCE, cfg: MatchConfig) -> float | None:
    """The alias factor k with P_ours ≈ k·P_tce and coincident transits, or None."""
    hits: list[tuple[float, float]] = []
    for k in cfg.aliases:
        rel = abs(sig.period_days / (k * tce.period_days) - 1.0)
        if rel > cfg.period_rel_tol:
            continue
        # Wrap the epoch offset by the SHORTER period: at 2x every one of our
        # transits is a TCE transit; at 1/2x every TCE transit is one of ours.
        wrap = min(sig.period_days, tce.period_days)
        dt = (sig.epoch_btjd - tce.epoch_btjd + 0.5 * wrap) % wrap - 0.5 * wrap
        tol = (
            cfg.epoch_tol_durations * max(sig.duration_hours, tce.duration_hours) / 24.0
        )
        if abs(dt) > tol:
            continue
        hits.append((k, rel))
    if not hits:
        return None
    return min(hits, key=lambda h: (h[0] != 1.0, h[1]))[0]


def match_star(
    tic_id: int,
    result: SearchResult | None,
    tces: Iterable[TCE],
    cfg: MatchConfig | None = None,
) -> list[MatchRow]:
    """Bucket one star's TCEs and signals. Each signal matches at most one TCE."""
    cfg = cfg or MatchConfig()
    tces = sorted(tces, key=lambda t: t.planet_num)
    if result is None:
        return [MatchRow(tic_id, "missed", tce=t, cause="not_searched") for t in tces]

    search_cfg = result.config
    free = list(result.signals)
    rows: list[MatchRow] = []
    for t in tces:
        hits = [(s, k) for s in free if (k := _alias_of(s, t, cfg)) is not None]
        if hits:
            # Prefer an un-aliased match, then the strongest signal.
            s, k = min(hits, key=lambda h: (h[1] != 1.0, h[0].rank))
            free.remove(s)
            rows.append(
                MatchRow(
                    tic_id,
                    "recovered" if k == 1.0 else "recovered_alias",
                    tce=t,
                    signal=s,
                    alias=k,
                )
            )
            continue
        if (
            not search_cfg["min_period_days"]
            <= t.period_days
            <= search_cfg["max_period_days"]
        ):
            cause = "outside_search_grid"
        elif result.stop_reason == "max_signals":
            cause = "search_exhausted"
        else:
            cause = "not_detected"
        rows.append(MatchRow(tic_id, "missed", tce=t, cause=cause))
    rows.extend(MatchRow(tic_id, "unmatched", signal=s) for s in free)
    return rows


def match_sector(
    results: dict[int, SearchResult],
    tces: Iterable[TCE],
    cfg: MatchConfig | None = None,
) -> list[MatchRow]:
    """Bucket every searched star and every TCE in the sector.

    A TCE on a star we never searched is reported (``not_searched``), never
    dropped: a comparison that silently covers only the stars it reached would
    overstate the recovery rate.
    """
    by_tic: dict[int, list[TCE]] = {}
    for t in tces:
        by_tic.setdefault(t.tic_id, []).append(t)
    rows: list[MatchRow] = []
    for tic in sorted(set(results) | set(by_tic)):
        rows.extend(match_star(tic, results.get(tic), by_tic.get(tic, []), cfg))
    return rows


def summarise(rows: Iterable[MatchRow]) -> dict[str, int]:
    """Counts per bucket, and per miss cause as ``missed:<cause>``."""
    out: dict[str, int] = {}
    for r in rows:
        out[r.bucket] = out.get(r.bucket, 0) + 1
        if r.bucket == "missed":
            key = f"missed:{r.cause}"
            out[key] = out.get(key, 0) + 1
    return out
