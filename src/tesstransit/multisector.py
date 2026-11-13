"""Search stars observed in many sectors as ONE stitched light curve.

    python -m tesstransit.multisector --sectors 1-13 \\
        --tce-csv s0001-s0013_tcestats.csv --toi-csv tois.csv \\
        --workdir lc --out runs/ms-pilot --min-sectors 9

Each sector is flattened on its own, exactly as the single-sector search does
(spline with knots by BIC, gated pre-whitening): a star's variability and the
spacecraft systematics change from sector to sector, so one trend through a
year of data would fit neither. The flattened sectors are stitched, binned to
``--bin-minutes`` and searched with ``detrend_method="none"`` from 1 d to half
the baseline, so a planet needs two transits anywhere in the year, not in one
sector.

The period grid is set by the baseline: a trial period wrong by dP moves the
last transit by ~T dP / P, so the log spacing is the shortest duration over
``4 T`` (the single-sector grid's own ratio). Over sectors 1-13 that is
~350,000 periods against 24,000 for one sector.

The reference is SPOC's multi-sector search over the same sectors
(``<stamp>-s0001-s0013_dvr-tcestats.csv``), scored with the same matcher, so
the numbers are comparable with the single-sector pilot's.

``--sample-of RUN`` re-searches an earlier run's stars (e.g. sectors 1-13
stars over 1-13 + 27-39), a paired comparison rather than a new draw.

⚠️ Declared limits: durations stop at 8 h (the 0.5 d spline knots), so
long-period planets around large stars are searched with too short a box;
and a planet whose transits fall in only one sector is a single event here,
which this search does not report.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

import numpy as np

from . import _lightcurve as lc
from . import _manifest as mf
from . import _match as mt
from . import _search as sr
from .pilot import _result_from_dict, _toi_for, parse_tois, select_sample, stars_of

_client: Any = None


def parse_sectors(text: str) -> list[int]:
    """``"1-13"`` or ``"1,2,5"`` -> sorted sector numbers."""
    out: set[int] = set()
    for part in text.split(","):
        a, _, b = part.partition("-")
        out.update(range(int(a), int(b or a) + 1))
    return sorted(out)


def bin_curve(
    t: np.ndarray, f: np.ndarray, minutes: float
) -> tuple[np.ndarray, np.ndarray]:
    """Mean time and flux in bins of ``minutes``; empty bins dropped."""
    if minutes <= 0:
        return t, f
    width = minutes / 1440.0
    idx = np.floor((t - t[0]) / width).astype(np.int64)
    _, inv, counts = np.unique(idx, return_inverse=True, return_counts=True)
    return (
        np.bincount(inv, weights=t) / counts,
        np.bincount(inv, weights=f) / counts,
    )


def stitch(
    curves: list[lc.LightCurve], cfg: sr.SearchConfig, bin_minutes: float
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    """Flatten each sector with ``cfg``, stitch, bin: ``(time, flat, per_sector)``."""
    ts, fs, per = [], [], []
    for c in sorted(curves, key=lambda c: c.time[0]):
        flat, knot, removed = sr.flatten(c.time, c.flux, cfg)
        ok = np.isfinite(flat)
        ts.append(c.time[ok])
        fs.append(flat[ok])
        per.append(
            {
                "sector": c.sector,
                "n_points": int(ok.sum()),
                "knot_days": knot,
                "prewhitened": [asdict(p) for p in removed],
            }
        )
    t, f = bin_curve(np.concatenate(ts), np.concatenate(fs), bin_minutes)
    return t, f, per


def stitched_config(
    cfg: sr.SearchConfig, baseline_days: float, min_period_days: float
) -> sr.SearchConfig:
    """The search over the stitched curve: no detrending, a baseline-sized grid."""
    p_max = baseline_days / 2.0
    d_min = min(cfg.durations_hours) / 24.0
    n = math.ceil(math.log(p_max / min_period_days) * 4.0 * baseline_days / d_min)
    return replace(
        cfg,
        detrend_method="none",
        prewhiten_max_period_days=0.0,
        min_period_days=min_period_days,
        max_period_days=p_max,
        n_periods=n,
    )


# Thresholds where the stitched search differs from a single sector's. It
# tries ~15x the periods, so coincident noise dips clear the single-sector
# cuts: at SNR 7.1 / SDE 11 the inverted light curves of 22 of 50 quiet stars
# reported a signal (pilot), at 3 transits / SDE 15 still 8 of 50 (hold-out).
# The two-event veto plus SDE 13 was fixed BEFORE a third, fresh sample
# (seed 3, experiments.md) and gave 3 of 50, all 9 in-grid planet TOIs kept.
STITCHED_DEFAULTS: dict[str, Any] = {
    "sde_threshold": 13.0,
    "min_snr_without_two_strongest": 3.0,
}


# A stitch spanning more than a year (sectors 1-13 + 27-39, ~3 years) tries
# ~3x the periods again. At SDE 13 the inverted curves of 25 of 104 quiet
# stars gave a signal (ms26-pair-inv); SDE 18, fixed before a fresh sample
# (seed 26), gave 7 of 100, and all 9 planets of the paired 26-sector run
# have SDE >= 18. Validated only on that 1-13 + 27-39 layout.
LONG_STITCH_SECTOR_SPAN = 14
LONG_STITCH_SDE = 18.0


def base_config(
    prewhiten: bool = True, invert: bool = False, sectors: list[int] | None = None
) -> sr.SearchConfig:
    """The stitched search's SearchConfig before any ``--set`` overrides.
    ``sectors`` spanning more than a year raise the SDE threshold."""
    defaults = dict(STITCHED_DEFAULTS)
    if sectors and max(sectors) - min(sectors) >= LONG_STITCH_SECTOR_SPAN:
        defaults["sde_threshold"] = LONG_STITCH_SDE
    cfg = replace(
        sr.SearchConfig(),
        prewhiten_max_period_days=1.5 if prewhiten else 0.0,
        invert_flux=invert,
        **defaults,
    )
    cfg.validate()
    return cfg


def _search_star(
    tic: int,
    keys: dict[int, str],
    workdir: str,
    keep: bool,
    cfg: sr.SearchConfig,
    bin_minutes: float,
    min_period_days: float,
) -> dict[str, Any]:
    global _client
    if _client is None:
        _client = lc.s3_client()
    t0 = time.time()
    try:
        curves, paths = [], []
        for sector, key in sorted(keys.items()):
            path = lc.fetch_lc(_client, tic, sector, Path(workdir), key=key)
            paths.append(path)
            curves.append(lc.load_lc(path))
        t_fetch = time.time() - t0
        flux_cfg = cfg
        if cfg.invert_flux:  # invert before flattening, as search() does
            for c in curves:
                c.flux = 2.0 - c.flux
            flux_cfg = replace(cfg, invert_flux=False)
        t, f, per = stitch(curves, flux_cfg, bin_minutes)
        baseline = float(t[-1] - t[0])
        scfg = stitched_config(flux_cfg, baseline, min_period_days)
        res = sr.search(t, f, tic, 0, scfg)
        if not keep:
            for p in paths:
                p.unlink()
        return {
            "tic_id": tic,
            "ok": True,
            "result": asdict(res),
            "sectors": per,
            "baseline_days": baseline,
            "n_binned": len(t),
            "seconds": {"fetch": t_fetch, "total": time.time() - t0},
        }
    except Exception as exc:  # noqa: BLE001 — recorded per star, never swallowed
        return {"tic_id": tic, "ok": False, "error": f"{type(exc).__name__}: {exc}"}


def eligible(
    manifests: dict[int, dict[int, str]], min_sectors: int
) -> dict[int, dict[int, str]]:
    """``{tic: {sector: key}}`` for stars with 2-min data in >= ``min_sectors``."""
    by_star: dict[int, dict[int, str]] = {}
    for sector, targets in manifests.items():
        for tic, key in targets.items():
            by_star.setdefault(tic, {})[sector] = key
    return {k: v for k, v in by_star.items() if len(v) >= min_sectors}


def paired_sample(
    prior: frozenset[int], tces: list[mt.TCE], stars: dict[int, dict[int, str]]
) -> tuple[list[int], list[int], list[int]]:
    """The stars of earlier runs, re-grouped against THIS run's TCE table.

    ``(tce_stars, quiet_stars, dropped)``: a star with a TCE in the new table
    is a TCE star even if it was quiet before (a longer baseline can give
    SPOC a TCE it did not have), and a star no longer eligible is dropped and
    counted, never silently replaced.
    """
    with_tce = {t.tic_id for t in tces}
    kept = sorted(t for t in prior if t in stars)
    dropped = sorted(t for t in prior if t not in stars)
    return (
        [t for t in kept if t in with_tce],
        [t for t in kept if t not in with_tce],
        dropped,
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sectors", required=True, help='e.g. "1-13"')
    ap.add_argument(
        "--tce-csv", type=Path, required=True, help="SPOC multi-sector TCEs"
    )
    ap.add_argument("--toi-csv", type=Path)
    ap.add_argument("--workdir", type=Path, required=True, help="light-curve scratch")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--min-sectors", type=int, default=9)
    ap.add_argument("--n-tce", type=int, default=150)
    ap.add_argument("--n-quiet", type=int, default=50)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--bin-minutes", type=float, default=10.0)
    ap.add_argument("--min-period", type=float, default=1.0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--keep-lc", action="store_true")
    ap.add_argument(
        "--holdout-of", type=Path, action="append", default=[],
        help="earlier run dir(s): draw a sample sharing none of their stars",
    )  # fmt: skip
    ap.add_argument(
        "--sample-of", type=Path, action="append", default=[],
        help="earlier run dir(s): search the SAME stars (a paired comparison)",
    )  # fmt: skip
    ap.add_argument("--invert", action="store_true", help="false-alarm run (2 - flux)")
    ap.add_argument(
        "--no-prewhiten", dest="prewhiten", action="store_false",
        help="flatten each sector with the spline alone",
    )  # fmt: skip
    ap.add_argument(
        "--set", action="append", default=[], metavar="FIELD=VALUE",
        help="override one SearchConfig field (repeatable), e.g. snr_threshold=8",
    )  # fmt: skip
    args = ap.parse_args(argv)

    sectors = parse_sectors(args.sectors)
    cfg = base_config(
        prewhiten=args.prewhiten, invert=args.invert, sectors=sectors
    ).override(args.set)
    args.out.mkdir(parents=True, exist_ok=True)
    client = lc.s3_client()
    manifests = {
        s: mf.load_or_list(client, s, args.workdir / f"manifest-s{s:04d}.json")
        for s in sectors
    }
    stars = eligible(manifests, args.min_sectors)
    tces = mt.parse_tce_csv(args.tce_csv.read_text())
    dropped: list[int] = []
    if args.sample_of:
        if args.holdout_of:
            raise SystemExit("--sample-of and --holdout-of are opposites")
        tce_stars, quiet, dropped = paired_sample(stars_of(args.sample_of), tces, stars)
        unobserved = sorted({t.tic_id for t in tces} - set(stars))
    else:
        tce_stars, quiet, unobserved = select_sample(
            tces, stars, args.n_tce, args.n_quiet, args.seed, stars_of(args.holdout_of)
        )
    print(
        f"{len(stars)} stars with 2-min data in >= {args.min_sectors} of sectors "
        f"{args.sectors}; sample {len(tce_stars)} TCE + {len(quiet)} quiet",
        file=sys.stderr,
    )
    sample = tce_stars + quiet
    records: list[dict[str, Any]] = []
    t1 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = [
            pool.submit(
                _search_star,
                tic,
                stars[tic],
                str(args.workdir),
                args.keep_lc,
                cfg,
                args.bin_minutes,
                args.min_period,
            )
            for tic in sample
        ]
        for i, fut in enumerate(as_completed(futs), 1):
            records.append(fut.result())
            if i % 10 == 0:
                print(f"  {i}/{len(sample)} ({time.time() - t1:.0f}s)", file=sys.stderr)
    wall = time.time() - t1
    records.sort(key=lambda r: r["tic_id"])
    with (args.out / "results.jsonl").open("w") as fh:
        for r in records:
            fh.write(json.dumps(r) + "\n")

    results = {r["tic_id"]: _result_from_dict(r["result"]) for r in records if r["ok"]}
    failed = {r["tic_id"]: r["error"] for r in records if not r["ok"]}
    tce_set, quiet_set = set(tce_stars), set(quiet)
    sampled = [t for t in tces if t.tic_id in tce_set]
    rows = mt.match_sector(results, sampled)
    tois = parse_tois(args.toi_csv.read_text()) if args.toi_csv else {}
    with (args.out / "matches.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["tic_id", "group", "bucket", "cause", "alias", "tce_id", "tce_period", "tce_snr",
             "tce_depth_ppm", "sig_rank", "sig_period", "sig_snr", "sig_sde", "sig_depth_ppm",
             "sig_duration_h", "toi", "toi_disposition"]
        )  # fmt: skip
        for r in rows:
            period = (
                r.tce.period_days
                if r.tce
                else (r.signal.period_days if r.signal else None)
            )
            toi = _toi_for(r.tic_id, period, tois)
            w.writerow(
                [r.tic_id, "quiet" if r.tic_id in quiet_set else "tce", r.bucket, r.cause, r.alias,
                 r.tce.tce_id if r.tce else "", r.tce.period_days if r.tce else "",
                 r.tce.snr if r.tce else "", r.tce.depth_ppm if r.tce else "",
                 r.signal.rank if r.signal else "", r.signal.period_days if r.signal else "",
                 r.signal.snr if r.signal else "", r.signal.sde if r.signal else "",
                 r.signal.depth_ppm if r.signal else "", r.signal.duration_hours if r.signal else "",
                 toi["toi"] if toi else "", toi["disposition"] if toi else ""]
            )  # fmt: skip
    totals = [r["seconds"]["total"] for r in records if r["ok"]]
    summary = {
        "sectors": sectors,
        "min_sectors": args.min_sectors,
        "eligible_stars": len(stars),
        "bin_minutes": args.bin_minutes,
        "sample": {"tce_stars": len(tce_stars), "quiet_stars": len(quiet), "failed": failed},
        "holdout_of": [str(p) for p in args.holdout_of],
        "sample_of": [str(p) for p in args.sample_of],
        "sample_of_not_eligible": dropped,
        "tce_stars_not_eligible": len(unobserved),
        "tces_in_sample": len(sampled),
        "tce_group": mt.summarise([r for r in rows if r.tic_id in tce_set]),
        "quiet_group": {
            "stars_with_any_signal": len({r.tic_id for r in rows if r.tic_id in quiet_set}),
            "signals": sum(1 for r in rows if r.tic_id in quiet_set),
        },
        "search_config": cfg.to_dict(),
        "match_config": asdict(mt.MatchConfig()),
        "seconds": {
            "wall": wall,
            "per_star_median": float(np.median(totals)) if totals else None,
            "workers": args.workers,
        },
    }  # fmt: skip
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
