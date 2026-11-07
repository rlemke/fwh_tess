"""Phase-0 pilot: re-search a stratified sample of one sector and match to SPOC.

Usage::

    python -m tesstransit.pilot --sector 1 --tce-csv s0001_tcestats.csv \\
        --toi-csv tois.csv --workdir /path/to/scratch --out /path/to/run

The sample has two parts, chosen for what each can measure:

- **TCE stars** (``--n-tce``), spread evenly over SPOC's SNR range, so the
  recovery rate is measured from the faintest TCEs to the strongest rather
  than dominated by easy hot Jupiters.
- **Quiet stars** (``--n-quiet``), a seeded random draw of observed stars with
  NO TCE. Any signal we report on one of these is either something SPOC
  rejected or a false alarm of ours, so this group is what calibrates the
  thresholds.

Light curves are streamed and deleted after the search unless ``--keep-lc``.
A star that fails is recorded with its error, never dropped from the counts.

⚠️ The comparison covers the SAMPLED stars only: TCEs on stars not in the
sample are excluded rather than reported as ``not_searched``, and the summary
says so. This is a pilot, not a sector result.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

from . import _lightcurve as lc
from . import _manifest as mf
from . import _match as mt
from . import _search as sr

_client: Any = None


def forced_stats(
    time_: np.ndarray,
    flux: np.ndarray,
    ephemerides: list[tuple[str, float, float, float]],
    cfg: sr.SearchConfig,
) -> dict[str, dict[str, Any]]:
    """Each known ephemeris measured on the light curve exactly as the search sees it.

    ``ephemerides`` are ``(name, period_d, epoch_btjd, duration_h)``. The flux
    is detrended with the search's own settings, so a transit the detrend
    erases shows up here as a low SNR: that separates "not in the data we
    search" from "in the data, and the search did not find it".
    """
    f = sr.flatten(time_, flux, cfg)[0]
    ok = np.isfinite(f)
    t, f = time_[ok], f[ok]
    out: dict[str, dict[str, Any]] = {}
    for name, period, epoch, dur_h in ephemerides:
        dur = dur_h / 24.0
        n = sr._n_transits(t, period, epoch, dur)
        rec: dict[str, Any] = {"n_transits": n, "snr_red": None, "depth_ppm": None}
        if n:
            try:
                ev = sr.event_stats(t, f, period, epoch, dur)
                rec["snr_red"] = float(ev.snr_red)
                rec["depth_ppm"] = float(np.average(ev.depths) * 1e6)
            except ValueError as exc:  # too few bins for a very long duration
                rec["error"] = str(exc)
        out[name] = rec
    return out


def _search_one(
    tic: int,
    sector: int,
    key: str,
    workdir: str,
    keep: bool,
    cfg: sr.SearchConfig,
    forced: list[tuple[str, float, float, float]] | None = None,
) -> dict[str, Any]:
    global _client
    if _client is None:
        _client = lc.s3_client()
    t0 = time.time()
    try:
        path = lc.fetch_lc(_client, tic, sector, Path(workdir), key=key)
        t_fetch = time.time() - t0
        curve = lc.load_lc(path)
        res = sr.search(curve.time, curve.flux, tic, sector, cfg)
        extra = (
            {"forced": forced_stats(curve.time, curve.flux, forced, cfg)}
            if forced
            else {}
        )
        if not keep:
            path.unlink()
        return {
            "tic_id": tic,
            "ok": True,
            **extra,
            "result": asdict(res),
            "mask": {
                "n_cadences": curve.n_cadences,
                "n_flagged": curve.n_flagged,
                "n_nonfinite": curve.n_nonfinite,
            },
            "t_span": [float(curve.time[0]), float(curve.time[-1])],
            "seconds": {"fetch": t_fetch, "total": time.time() - t0},
        }
    except Exception as exc:  # noqa: BLE001 — recorded per star, never swallowed
        return {"tic_id": tic, "ok": False, "error": f"{type(exc).__name__}: {exc}"}


def _result_from_dict(d: dict[str, Any]) -> sr.SearchResult:
    d = dict(d)
    d["signals"] = [sr.Signal(**s) for s in d["signals"]]
    d["vetoes"] = [sr.Veto(**v) for v in d.get("vetoes", [])]
    d["prewhitened"] = [sr.Prewhitened(**v) for v in d.get("prewhitened", [])]
    return sr.SearchResult(**d)


def select_sample(
    tces: list[mt.TCE],
    targets: dict[int, str],
    n_tce: int,
    n_quiet: int,
    seed: int,
    exclude: frozenset[int] = frozenset(),
) -> tuple[list[int], list[int], list[int]]:
    """``(tce_stars, quiet_stars, tce_stars_not_observed)``.

    With no ``exclude``: TCE stars at evenly spaced SNR ranks (deterministic;
    the development sample). With ``exclude`` (a HOLD-OUT): stars in it are
    removed, the rest split into ``n_tce`` contiguous SNR bins, and one star is
    drawn per bin with ``seed`` -- the same SNR coverage, none of the stars.
    """
    best: dict[int, float] = {}
    for t in tces:
        best[t.tic_id] = max(
            best.get(t.tic_id, -1.0), t.snr if t.snr is not None else -1.0
        )
    unobserved = sorted(tic for tic in best if tic not in targets)
    ranked = sorted(
        (tic for tic in best if tic in targets and tic not in exclude),
        key=lambda k: best[k],
    )
    rng = np.random.default_rng(seed)
    if n_tce >= len(ranked):
        tce_stars = ranked
    elif exclude:
        tce_stars = sorted(int(rng.choice(b)) for b in np.array_split(ranked, n_tce))
    else:
        idx = np.unique(np.linspace(0, len(ranked) - 1, n_tce).round().astype(int))
        tce_stars = [ranked[i] for i in idx]
    pool = sorted(set(targets) - set(best) - set(exclude))
    quiet = sorted(
        int(x) for x in rng.choice(pool, size=min(n_quiet, len(pool)), replace=False)
    )
    return tce_stars, quiet, unobserved


def parse_tois(text: str) -> dict[int, list[dict[str, Any]]]:
    lines = (ln for ln in io.StringIO(text) if not ln.startswith("#"))
    out: dict[int, list[dict[str, Any]]] = {}
    for row in csv.DictReader(lines):
        period = row.get("Orbital Period (days) Value", "").strip()
        out.setdefault(int(row["TIC"]), []).append(
            {
                "toi": row["Full TOI ID"],
                "disposition": row["TOI Disposition"],
                "period_days": float(period) if period else None,
                "comment": row.get("Public Comment", ""),
            }
        )
    return out


def _toi_for(tic: int, period: float | None, tois: dict[int, list[dict[str, Any]]]):
    """The TOI on this star at (an alias of) this period, if any."""
    if period is None:
        return None
    for t in tois.get(tic, []):
        p = t["period_days"]
        if p and any(abs(period / (k * p) - 1) < 0.01 for k in mt.DEFAULT_ALIASES):
            return t
    return None


def stars_of(runs: list[Path]) -> frozenset[int]:
    """Every star an earlier run searched, so a hold-out can avoid them."""
    return frozenset(
        json.loads(line)["tic_id"]
        for run in runs
        for line in (run / "results.jsonl").read_text().splitlines()
    )


def _search_sample(
    args: argparse.Namespace,
    cfg: sr.SearchConfig,
    targets: dict[int, str],
    sample: list[int],
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    t1 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = [
            pool.submit(
                _search_one,
                tic,
                args.sector,
                targets[tic],
                str(args.workdir),
                args.keep_lc,
                cfg,
            )
            for tic in sample
        ]
        for i, fut in enumerate(as_completed(futs), 1):
            records.append(fut.result())
            if i % 10 == 0:
                print(f"  {i}/{len(sample)} ({time.time() - t1:.0f}s)", file=sys.stderr)
    return records


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument("--tce-csv", type=Path, required=True)
    ap.add_argument("--toi-csv", type=Path)
    ap.add_argument("--workdir", type=Path, required=True, help="light-curve scratch")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--n-tce", type=int, default=60)
    ap.add_argument("--n-quiet", type=int, default=40)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--keep-lc", action="store_true")
    ap.add_argument(
        "--holdout-of",
        type=Path,
        action="append",
        default=[],
        help="earlier run dir(s): draw a sample sharing none of their stars",
    )
    ap.add_argument(
        "--detrend-method", choices=["median", "biweight", "spline"], default="spline"
    )
    ap.add_argument(
        "--prewhiten",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="spline: first remove two-sided stellar periodicities faster than 1.5 d "
        "(default on; ignored by the other detrenders)",
    )
    ap.add_argument(
        "--per-duration",
        action="store_true",
        help="one detrend window per trial duration",
    )
    ap.add_argument(
        "--reuse",
        action="store_true",
        help="re-match the run's existing results.jsonl instead of searching again",
    )
    ap.add_argument(
        "--invert",
        action="store_true",
        help="false-alarm test: search 2 - flux, so every signal is a false alarm",
    )
    ap.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="FIELD=VALUE",
        help="override one SearchConfig field (repeatable), e.g. prewhiten_max_terms=5",
    )
    args = ap.parse_args(argv)

    cfg = sr.SearchConfig(
        detrend_method=args.detrend_method,
        detrend_per_duration=args.per_duration,
        invert_flux=args.invert,
        prewhiten_max_period_days=(
            1.5 if args.prewhiten and args.detrend_method == "spline" else 0.0
        ),
    )
    cfg = cfg.override(args.set)
    args.out.mkdir(parents=True, exist_ok=True)
    client = lc.s3_client()
    t0 = time.time()
    targets = mf.load_or_list(
        client, args.sector, args.workdir / f"manifest-s{args.sector:04d}.json"
    )
    print(
        f"manifest: {len(targets)} targets ({time.time() - t0:.0f}s)", file=sys.stderr
    )
    tces = mt.parse_tce_csv(args.tce_csv.read_text())
    exclude = stars_of(args.holdout_of)
    tce_stars, quiet, unobserved = select_sample(
        tces, targets, args.n_tce, args.n_quiet, args.seed, exclude
    )
    print(
        f"sample: {len(tce_stars)} TCE stars + {len(quiet)} quiet; "
        f"{len(unobserved)} TCE stars have no 2-min light curve in the bucket",
        file=sys.stderr,
    )

    sample = tce_stars + quiet
    records: list[dict[str, Any]] = []
    t1 = time.time()
    if args.reuse:
        records = [
            json.loads(ln)
            for ln in (args.out / "results.jsonl").read_text().splitlines()
        ]
        missing = set(sample) - {r["tic_id"] for r in records}
        if missing:
            raise SystemExit(
                f"--reuse: {len(missing)} sample stars not in results.jsonl"
            )
        records = [r for r in records if r["tic_id"] in set(sample)]
    else:
        records = _search_sample(args, cfg, targets, sample)
    wall = time.time() - t1
    records.sort(key=lambda r: r["tic_id"])
    if not args.reuse:
        with (args.out / "results.jsonl").open("w") as fh:
            for r in records:
                fh.write(json.dumps(r) + "\n")

    results = {r["tic_id"]: _result_from_dict(r["result"]) for r in records if r["ok"]}
    failed = {r["tic_id"]: r["error"] for r in records if not r["ok"]}
    tce_set = set(tce_stars)
    sampled_tces = [t for t in tces if t.tic_id in tce_set]
    rows = mt.match_sector(results, sampled_tces)
    tois = parse_tois(args.toi_csv.read_text()) if args.toi_csv else {}
    quiet_set = set(quiet)

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

    tce_rows = [r for r in rows if r.tic_id in tce_set]
    quiet_rows = [r for r in rows if r.tic_id in quiet_set]
    totals = [r["seconds"]["total"] for r in records if r["ok"]]
    summary = {
        "sector": args.sector,
        "scope": "sampled stars only; TCEs on unsampled stars excluded",
        "holdout_of": [str(p) for p in args.holdout_of],
        "seed": args.seed,
        "sample": {
            "tce_stars": len(tce_stars),
            "quiet_stars": len(quiet),
            "failed": failed,
        },
        "tce_stars_without_2min_lc": len(unobserved),
        "tces_in_sample": len(sampled_tces),
        "tce_group": mt.summarise(tce_rows),
        "quiet_group": {
            "stars_with_any_signal": len({r.tic_id for r in quiet_rows}),
            "signals": len(quiet_rows),
        },
        "search_config": cfg.to_dict(),
        "match_config": asdict(mt.MatchConfig()),
        # Per-star times are recorded with each result, so they survive
        # --reuse; the wall time does not -- a re-match did no searching.
        "seconds": {
            "wall": None if args.reuse else wall,
            "per_star_median": float(np.median(totals)) if totals else None,
            "workers": None if args.reuse else args.workers,
        },
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
