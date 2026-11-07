"""Shape tests: is a periodic dip a planet or an eclipsing binary?

    python -m tesstransit.vet --runs runs --planets 1=s1-planets-pwB ... \\
        --workdir lc --out runs/vet

Three standard tests on the flattened light curve at a signal's own
ephemeris (Kepler's Robovetter and SPOC DV make the same three):

- **odd/even**: an eclipsing binary of two similar stars, found at half its
  period, alternates primary and secondary eclipses. ``odd_even_sigma`` is
  the difference of the mean odd and even depths over its uncertainty.
- **secondary**: a dip at phase 0.5, with the same duration. A planet's
  secondary eclipse is tens of ppm at most, so ``secondary_snr`` well above
  noise means a self-luminous companion.
- **shape**: a planet crossing a star gives a flat-bottomed transit (short
  ingress); grazing eclipses are V-shaped. ``shape_ratio`` is the mean depth
  over the whole box divided by the mean depth in its central half: ~0.9 for
  a planet, 2/3 for a pure V.

Every uncertainty is from the flux binned to the transit duration, so red
noise is included (``_search.event_stats``). The scoring CLI recomputes the
tests for every TOI a planet run recovered whose TFOPWG disposition is a
planet (CP, KP) or an eclipsing binary (EB), so each test's separating power
is measured on labelled TESS data, not assumed.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from . import _lightcurve as lc
from . import _manifest as mf
from . import _search as sr


@dataclass
class Vet:
    odd_even_sigma: float  # |odd - even| / combined sigma
    depth_odd_ppm: float
    depth_even_ppm: float
    secondary_snr: float  # dip at phase 0.5, same duration
    secondary_ppm: float
    shape_ratio: float  # full-box depth / central-half depth
    n_odd: int
    n_even: int


def _depth(
    t: np.ndarray, f: np.ndarray, period: float, epoch: float, half: float
) -> tuple[np.ndarray, np.ndarray]:
    """``(inside, phase)``: points within ``half`` days of a transit centre."""
    phase = (t - epoch + 0.5 * period) % period - 0.5 * period
    return np.abs(phase) < half, phase


def vet(
    t: np.ndarray, f: np.ndarray, period: float, epoch: float, duration: float
) -> Vet:
    """The three tests at ``period``/``epoch``/``duration`` (days) on flat ``f``."""
    ev = sr.event_stats(t, f, period, epoch, duration)
    k = np.round((ev.centres - epoch) / period).astype(int)
    odd = k % 2 == 1

    def mean(mask: np.ndarray) -> tuple[float, float]:
        if not mask.any():
            return float("nan"), float("inf")
        w = 1.0 / ev.sigmas[mask] ** 2
        return float(np.sum(ev.depths[mask] * w) / np.sum(w)), float(
            1.0 / np.sqrt(np.sum(w))
        )

    (d_o, s_o), (d_e, s_e) = mean(odd), mean(~odd)
    oe = abs(d_o - d_e) / np.hypot(s_o, s_e) if odd.any() and (~odd).any() else 0.0
    in_full, _ = _depth(t, f, period, epoch, duration / 2.0)
    # The primary is masked first: left in, it drags the out-of-eclipse
    # baseline down and a planet reads as a negative secondary (measured).
    prim, _ = _depth(t, f, period, epoch, 0.75 * duration)
    try:
        sec = sr.event_stats(t[~prim], f[~prim], period, epoch + 0.5 * period, duration)
        w = 1.0 / sec.sigmas**2
        sec_snr, sec_depth = sec.snr_red, float(np.sum(sec.depths * w) / np.sum(w))
    except ValueError:  # no data at phase 0.5 (it fell in the gaps)
        sec_snr, sec_depth = float("nan"), float("nan")
    in_core, _ = _depth(t, f, period, epoch, duration / 4.0)
    base = float(np.median(f[~in_full]))
    full = base - float(np.mean(f[in_full]))
    core = base - float(np.mean(f[in_core])) if in_core.any() else float("nan")
    return Vet(
        odd_even_sigma=float(oe),
        depth_odd_ppm=d_o * 1e6,
        depth_even_ppm=d_e * 1e6,
        secondary_snr=float(sec_snr),
        secondary_ppm=sec_depth * 1e6,
        shape_ratio=float(full / core) if core and core > 0 else float("nan"),
        n_odd=int(odd.sum()),
        n_even=int((~odd).sum()),
    )


# Thresholds measured on the TOIs sectors 1-13 recovered (403 CP/KP planets,
# 265 EBs), chosen to flag at most 3% of planets, and the same when tuned on
# odd sectors and tested on even ones or the reverse: 10 of 403 planets (2.5%)
# and 76 of 265 binaries (29%) flagged. Each test needs a SIZE as well as a
# significance -- on a deep hot Jupiter a 2% odd/even difference is already
# 3 sigma. Most EB TOIs are not flagged because they are diluted background
# binaries, which look like planets in shape; that needs a centroid test.
ODD_EVEN_MIN_SIGMA = 3.0
ODD_EVEN_MIN_FRACTION = 0.3  # |odd - even| / mean depth
SECONDARY_MIN_SNR = 3.0
SECONDARY_MIN_FRACTION = 0.05  # secondary / primary depth
SHAPE_MAX_RATIO = 0.7


def eb_like(v: Vet, depth_ppm: float) -> list[str]:
    """The tests ``v`` fails, by name; empty = consistent with a planet."""
    out = []
    mean = max(abs(v.depth_odd_ppm + v.depth_even_ppm) / 2.0, 1.0)
    if (
        v.odd_even_sigma > ODD_EVEN_MIN_SIGMA
        and abs(v.depth_odd_ppm - v.depth_even_ppm) / mean > ODD_EVEN_MIN_FRACTION
    ):
        out.append("odd_even")
    if (
        v.secondary_snr > SECONDARY_MIN_SNR
        and v.secondary_ppm / max(depth_ppm, 1.0) > SECONDARY_MIN_FRACTION
    ):
        out.append("secondary")
    if v.shape_ratio < SHAPE_MAX_RATIO:
        out.append("v_shape")
    return out


_client: Any = None


def _vet_star(
    tic: int,
    sector: int,
    key: str,
    workdir: str,
    cfg: dict[str, Any],
    signals: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    global _client
    if _client is None:
        _client = lc.s3_client()
    out = []
    try:
        path = lc.fetch_lc(_client, tic, sector, Path(workdir), key=key)
        c = lc.load_lc(path)
        scfg = sr.SearchConfig(
            **{
                k: tuple(v) if isinstance(v, list) else v
                for k, v in cfg.items()
                if k in sr.SearchConfig.__dataclass_fields__
            }
        )
        flat, _, _ = sr.flatten(c.time, c.flux, scfg)
        ok = np.isfinite(flat)
        for s in signals:
            v = vet(
                c.time[ok], flat[ok], s["period_days"], s["epoch_btjd"],
                s["duration_hours"] / 24.0,
            )  # fmt: skip
            fails = eb_like(v, s["row"]["sig_depth_ppm"])
            out.append({**s["row"], **asdict(v), "fails": "+".join(fails), "ok": True})
        path.unlink()
    except Exception as exc:  # noqa: BLE001 — recorded per star, never swallowed
        out = [{**s["row"], "ok": False, "error": f"{type(exc).__name__}: {exc}"}
               for s in signals]  # fmt: skip
    return out


def _jobs(runs: Path, run: str, sector: int, workdir: Path) -> list[tuple]:
    """Per star: the recovered CP/KP/EB TOIs, at our matched signal's ephemeris."""
    results = {
        r["tic_id"]: r
        for r in map(json.loads, (runs / run / "results.jsonl").read_text().splitlines())
        if r["ok"]
    }  # fmt: skip
    targets = mf.parse((workdir / f"manifest-s{sector:04d}.json").read_text(), sector)
    by_star: dict[int, list[dict[str, Any]]] = {}
    for row in csv.DictReader((runs / run / "tois.csv").open()):
        if not row["outcome"].startswith("recovered"):
            continue
        if row["disposition"] not in ("CP", "KP", "EB"):
            continue
        tic = int(row["tic_id"])
        sigs = results[tic]["result"]["signals"]
        sp = float(row["sig_period"])
        s = min(sigs, key=lambda x: abs(x["period_days"] - sp))
        by_star.setdefault(tic, []).append(
            {
                "period_days": s["period_days"],
                "epoch_btjd": s["epoch_btjd"],
                "duration_hours": s["duration_hours"],
                "row": {
                    "sector": sector,
                    "toi": row["toi"],
                    "tic_id": tic,
                    "disposition": row["disposition"],
                    "alias": row["alias"],
                    "sig_period": s["period_days"],
                    "sig_depth_ppm": s["depth_ppm"],
                    "sig_snr": s["snr"],
                },
            }
        )
    return [
        (tic, sector, targets[tic], str(workdir), results[tic]["result"]["config"], sg)
        for tic, sg in sorted(by_star.items())
    ]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--planets", action="append", required=True, metavar="SECTOR=RUN")
    ap.add_argument("--workdir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args(argv)
    jobs = []
    for spec in args.planets:
        s, _, run = spec.partition("=")
        jobs += _jobs(args.runs, run, int(s), args.workdir)
    rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for fut in as_completed([pool.submit(_vet_star, *j) for j in jobs]):
            rows.extend(fut.result())
    rows.sort(key=lambda r: (r["sector"], r["toi"]))
    args.out.mkdir(parents=True, exist_ok=True)
    fields = sorted({k for r in rows for k in r}, key=lambda k: (k not in rows[0], k))
    with (args.out / "vet.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    bad = [r for r in rows if not r["ok"]]
    print(f"{len(rows)} signals vetted on {len(jobs)} stars; {len(bad)} failed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
