"""Injection-recovery: what fraction of transits the search finds, by SNR and period.

Usage::

    python -m tesstransit.inject --sector 1 --tce-csv s0001_tcestats.csv \\
        --workdir /path/to/scratch --sample-run runs/s1-pilot-pw2 \\
        --out runs/s1-inject

Synthetic transits go into the REAL PDCSAP flux of a sample run's stars with no
TCE, before detrending, so the measurement includes everything the pipeline
does to a transit: the spline, pre-whitening, the quality mask, the gaps. Each
light curve carries ONE injection and is searched with the unchanged search; a
star gets ``--per-star`` of them, drawn independently.

What is drawn, per injection:

- period log-uniform in ``[0.5, 13]`` d, the range a 27-day sector can show
  twice; an epoch uniform over one period, redrawn until at least two transits
  fall on data
- impact parameter uniform in ``[0, 0.9]``, and the duration a solar-density
  star gives that period and impact parameter
- an **expected SNR** log-uniform in ``[4, 40]``. The depth is set from it with
  the star's own per-cadence scatter: ``depth = snr * sigma / sqrt(n_in)``,
  ``sigma`` the MAD scatter of the detrended, un-injected flux and ``n_in``
  the cadences inside the flat bottom. This is the white-noise SNR; a star's
  red noise is what makes the measured curve fall short of it.

⚠️ The transit is a trapezoid (ingress from the radius ratio and impact
parameter), not a limb-darkened model: a box search loses a little more on a
rounder transit, so the curve is slightly optimistic for grazing geometry.

⚠️ Stars with no TCE only. A TCE star already has a strong signal that would
compete with the injection, and the sample's no-TCE stars are a random draw, so
they stand for the bulk of the sector.

Recovered = a signal at the injected period and epoch (``recovered``), or at
2x, 1/2x, 3x or 1/3x with coincident transits (``recovered_alias``), with the
same tolerances used to match SPOC's TCEs.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from . import _lightcurve as lc
from . import _manifest as mf
from . import _match as mt
from . import _search as sr

MIN_PERIOD_DAYS = 0.5
MAX_PERIOD_DAYS = 13.0
MIN_SNR = 4.0
MAX_SNR = 40.0
MAX_IMPACT = 0.9
_client: Any = None


@dataclass(frozen=True)
class Injection:
    index: int
    period_days: float
    epoch_btjd: float
    duration_hours: float  # first to fourth contact
    impact: float
    radius_ratio: float
    depth_ppm: float  # flat-bottom depth
    snr_expected: float
    n_transits: int


def solar_duration_hours(period_days: float, impact: float, k: float = 0.0) -> float:
    """First-to-fourth-contact duration for a circular orbit around a Sun-like star."""
    a_rs = 215.03 * (period_days / 365.25) ** (2.0 / 3.0)  # a / R_sun, Kepler III
    chord = max((1.0 + k) ** 2 - impact**2, 0.0)
    return (
        24.0 * period_days / np.pi * float(np.arcsin(min(np.sqrt(chord) / a_rs, 1.0)))
    )


def trapezoid(
    t: np.ndarray,
    period: float,
    epoch: float,
    dur_h: float,
    ingress_h: float,
    depth: float,
) -> np.ndarray:
    """Fractional flux loss of a trapezoidal transit (``depth`` at the flat bottom)."""
    half = dur_h / 48.0
    ing = min(ingress_h / 24.0, half)
    x = np.abs((t - epoch + 0.5 * period) % period - 0.5 * period)
    out = np.zeros_like(t)
    flat = x <= half - ing
    out[flat] = depth
    edge = (x > half - ing) & (x < half)
    if ing > 0:
        out[edge] = depth * (half - x[edge]) / ing
    return out


def _ingress_hours(dur_h: float, impact: float, k: float) -> float:
    """Ingress duration from the ratio of the inner to the outer chord."""
    outer = (1.0 + k) ** 2 - impact**2
    inner = max((1.0 - k) ** 2 - impact**2, 0.0)
    if outer <= 0:
        return dur_h / 2.0
    return 0.5 * dur_h * (1.0 - np.sqrt(inner / outer))


def draw(
    rng: np.random.Generator, t: np.ndarray, sigma: float, index: int
) -> Injection:
    """One injection whose transits fall on ``t`` at least twice."""
    for _ in range(1000):
        period = float(
            np.exp(rng.uniform(np.log(MIN_PERIOD_DAYS), np.log(MAX_PERIOD_DAYS)))
        )
        epoch = float(t[0] + rng.uniform(0.0, period))
        impact = float(rng.uniform(0.0, MAX_IMPACT))
        snr = float(np.exp(rng.uniform(np.log(MIN_SNR), np.log(MAX_SNR))))
        dur_h = solar_duration_hours(period, impact)
        n_tr = sr._n_transits(t, period, epoch, dur_h / 24.0)
        if n_tr < 2:
            continue
        # The radius ratio and the depth depend on each other only through the
        # ingress, so two fixed-point passes settle it.
        k = 0.05
        for _ in range(3):
            ing = _ingress_hours(dur_h, impact, k)
            flat = trapezoid(t, period, epoch, dur_h, ing, 1.0) >= 1.0
            n_in = int(flat.sum())
            if n_in == 0:
                break
            depth = snr * sigma / np.sqrt(n_in)
            k = float(np.sqrt(depth))
        if n_in == 0 or depth >= 0.5:
            continue
        dur_h = solar_duration_hours(period, impact, k)
        return Injection(
            index, period, epoch, dur_h, impact, k, float(depth * 1e6), snr, n_tr
        )
    raise ValueError("no injection with two transits on data in 1000 draws")


def recovered(
    result: sr.SearchResult, inj: Injection
) -> tuple[str, sr.Signal | None, float | None]:
    """``(bucket, signal, alias)`` for one injection, with the TCE-match tolerances."""
    ref = mt.TCE(
        tce_id=f"inj-{inj.index}",
        tic_id=result.tic_id,
        planet_num=1,
        period_days=inj.period_days,
        epoch_btjd=inj.epoch_btjd,
        duration_hours=inj.duration_hours,
        depth_ppm=inj.depth_ppm,
        snr=inj.snr_expected,
        n_transits=inj.n_transits,
    )
    row = next(
        r for r in mt.match_star(result.tic_id, result, [ref]) if r.tce is not None
    )
    return row.bucket, row.signal, row.alias


def _inject_star(
    tic: int,
    sector: int,
    key: str,
    workdir: str,
    keep: bool,
    cfg: sr.SearchConfig,
    per_star: int,
    seed: int,
) -> list[dict[str, Any]]:
    global _client
    if _client is None:
        _client = lc.s3_client()
    try:
        path = lc.fetch_lc(_client, tic, sector, Path(workdir), key=key)
        curve = lc.load_lc(path)
        if not keep:
            path.unlink()
        flat = sr.flatten(curve.time, curve.flux, cfg)[0]
        sigma = sr._mad_sigma(flat[np.isfinite(flat)] - 1.0)
    except Exception as exc:  # noqa: BLE001 — recorded per star, never swallowed
        return [{"tic_id": tic, "ok": False, "error": f"{type(exc).__name__}: {exc}"}]
    rng = np.random.default_rng([seed, sector, tic])
    out: list[dict[str, Any]] = []
    for i in range(per_star):
        t0 = time.time()
        try:
            inj = draw(rng, curve.time, sigma, i)
            ing = _ingress_hours(inj.duration_hours, inj.impact, inj.radius_ratio)
            loss = trapezoid(
                curve.time,
                inj.period_days,
                inj.epoch_btjd,
                inj.duration_hours,
                ing,
                inj.depth_ppm * 1e-6,
            )
            res = sr.search(curve.time, curve.flux * (1.0 - loss), tic, sector, cfg)
            bucket, sig, alias = recovered(res, inj)
            out.append({
                "tic_id": tic, "ok": True, "sigma_ppm": sigma * 1e6,
                "injection": asdict(inj), "bucket": bucket, "alias": alias,
                "signal": asdict(sig) if sig else None,
                "n_signals": len(res.signals), "stop_reason": res.stop_reason,
                "seconds": time.time() - t0,
            })  # fmt: skip
        except Exception as exc:  # noqa: BLE001
            out.append(
                {
                    "tic_id": tic,
                    "ok": False,
                    "index": i,
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
    return out


def no_tce_stars(sample_run: Path, tces: list[mt.TCE]) -> list[int]:
    """The sample run's searched stars that carry no TCE in this sector."""
    with_tce = {t.tic_id for t in tces}
    stars = {
        json.loads(ln)["tic_id"]
        for ln in (sample_run / "results.jsonl").read_text().splitlines()
    }
    return sorted(stars - with_tce)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument("--tce-csv", type=Path, required=True)
    ap.add_argument("--workdir", type=Path, required=True, help="light-curve scratch")
    ap.add_argument("--sample-run", type=Path, required=True, help="a pilot run dir")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--per-star", type=int, default=6)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--keep-lc", action="store_true")
    ap.add_argument(
        "--prewhiten",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="first remove two-sided stellar periodicities faster than 1.5 d",
    )
    args = ap.parse_args(argv)

    cfg = sr.SearchConfig(prewhiten_max_period_days=1.5 if args.prewhiten else 0.0)
    cfg.validate()
    args.out.mkdir(parents=True, exist_ok=True)
    client = lc.s3_client()
    targets = mf.load_or_list(
        client, args.sector, args.workdir / f"manifest-s{args.sector:04d}.json"
    )
    tces = mt.parse_tce_csv(args.tce_csv.read_text())
    stars = no_tce_stars(args.sample_run, tces)
    print(f"{len(stars)} no-TCE stars x {args.per_star} injections", file=sys.stderr)

    records: list[dict[str, Any]] = []
    t1 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = [
            pool.submit(
                _inject_star, tic, args.sector, targets[tic], str(args.workdir),
                args.keep_lc, cfg, args.per_star, args.seed,
            )
            for tic in stars
        ]  # fmt: skip
        for i, fut in enumerate(as_completed(futs), 1):
            records.extend(fut.result())
            if i % 10 == 0:
                print(f"  {i}/{len(stars)} ({time.time() - t1:.0f}s)", file=sys.stderr)
    records.sort(
        key=lambda r: (
            r["tic_id"],
            r.get("injection", {}).get("index", r.get("index", -1)),
        )
    )
    with (args.out / "injections.jsonl").open("w") as fh:
        for r in records:
            fh.write(json.dumps(r) + "\n")

    with (args.out / "injections.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["tic_id", "index", "period_days", "epoch_btjd", "duration_hours", "impact",
             "radius_ratio", "depth_ppm", "snr_expected", "n_transits", "sigma_ppm",
             "bucket", "alias", "sig_rank", "sig_snr", "sig_sde"]
        )  # fmt: skip
        for r in records:
            if not r["ok"]:
                continue
            inj, sig = r["injection"], r["signal"] or {}
            w.writerow(
                [r["tic_id"], inj["index"], inj["period_days"], inj["epoch_btjd"],
                 inj["duration_hours"], inj["impact"], inj["radius_ratio"], inj["depth_ppm"],
                 inj["snr_expected"], inj["n_transits"], r["sigma_ppm"], r["bucket"], r["alias"],
                 sig.get("rank", ""), sig.get("snr", ""), sig.get("sde", "")]
            )  # fmt: skip

    ok = [r for r in records if r["ok"]]
    found = [r for r in ok if r["bucket"] != "missed"]
    summary = {
        "sector": args.sector,
        "sample_run": args.sample_run.name,
        "stars": len(stars),
        "per_star": args.per_star,
        "seed": args.seed,
        "injections": len(ok),
        "failed": [r for r in records if not r["ok"]],
        "recovered": sum(r["bucket"] == "recovered" for r in ok),
        "recovered_alias": sum(r["bucket"] == "recovered_alias" for r in ok),
        "fraction": len(found) / len(ok) if ok else None,
        "search_config": cfg.to_dict(),
        "draws": {
            "period_days": [MIN_PERIOD_DAYS, MAX_PERIOD_DAYS],
            "snr_expected": [MIN_SNR, MAX_SNR],
            "impact": [0.0, MAX_IMPACT],
        },
        "seconds": {"wall": time.time() - t1, "workers": args.workers},
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(
        json.dumps({k: v for k, v in summary.items() if k != "search_config"}, indent=2)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
