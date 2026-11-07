"""Transit-timing variations: does each transit arrive when a fixed period says?

    python -m tesstransit.ttv --sectors 1-13 --toi-csv tois.csv \\
        --workdir lc_ms --out runs/ttv --toi 216.01 --toi 216.02
    python -m tesstransit.ttv --sectors 1-13 --toi-csv tois.csv \\
        --workdir lc_ms --out runs/ttv --from-run runs/ms-pilot

Each star's sectors are flattened and stitched as in ``multisector``. Every
transit with data is timed on its own: a trapezoid with the signal's depth and
duration (ingress 15% of the duration) is slid across the predicted time, and
the mid-time is the chi-square minimum, its uncertainty where chi-square rises
by 1 (a parabola through the minimum). A straight line through the times gives
the best fixed ephemeris; what is left is O - C.

``chi2_dof`` far above 1 means the times wander more than their errors allow:
a second planet tugging (a resonant pair swaps energy and the times swing over
months), or an eclipsing binary's light-travel time. One sector holds too few
transits to show it, which is what a stitched year adds. ⚠️ A timing error
from a spot crossing or a gap at ingress inflates chi-square too, so a large
value is a candidate, not a detection, until the O - C curve is looked at.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from . import _lightcurve as lc
from . import _manifest as mf
from . import _search as sr
from . import multisector as ms


def _trapezoid(dt: np.ndarray, duration: float, depth: float) -> np.ndarray:
    ingress = 0.15 * duration
    x = np.abs(dt)
    edge = duration / 2.0
    return -depth * np.clip((edge - x) / ingress, 0.0, 1.0)


def refine_ephemeris(
    t: np.ndarray, f: np.ndarray, period: float, duration: float, width: float = 0.03
) -> tuple[float, float, float]:
    """``(period, epoch, depth)`` from a narrow BLS search on THIS data.

    A catalogue ephemeris can be years from the data, and for a planet whose
    times swing by hours the mean period does not carry back: TOI 216.02's,
    propagated 58 orbits to sectors 1-13, puts every trial time where there
    is no transit (folded depth ~0 ppm, measured).
    """
    from astropy.timeseries import BoxLeastSquares

    periods = np.linspace(period * (1 - width), period * (1 + width), 12000)
    r = BoxLeastSquares(t, f).power(periods, [duration], objective="snr")
    i = int(np.argmax(r.power))
    return float(r.period[i]), float(r.transit_time[i]), float(r.depth[i])


@dataclass
class TTV:
    period_days: float  # the fixed ephemeris fitted to the times
    epoch_btjd: float
    n_transits: int
    epochs: list[int] = field(default_factory=list)  # transit numbers
    times: list[float] = field(default_factory=list)  # BTJD
    sigmas_min: list[float] = field(default_factory=list)
    o_minus_c_min: list[float] = field(default_factory=list)
    chi2_dof: float = float("nan")
    max_abs_oc_min: float = float("nan")


def time_transits(
    t: np.ndarray,
    f: np.ndarray,
    period: float,
    epoch: float,
    duration: float,
    depth: float,
    min_coverage: float = 0.6,
) -> list[tuple[int, float, float]]:
    """``(n, mid_time, sigma)`` for each transit with ``min_coverage`` of data."""
    cadence = float(np.median(np.diff(t)))
    noise = sr._mad_sigma(f)
    # A resonant pair's times swing by hours (TOI 216: up to ~2 h on a 3.8 h
    # transit), so the trial shifts reach a full duration either way; a
    # minimum at the edge is still rejected rather than trusted.
    half_win = 2.0 * duration
    out = []
    n_all = np.round((t - epoch) / period).astype(int)
    for n in np.unique(n_all):
        tc = epoch + n * period
        sel = np.abs(t - tc) < half_win
        in_tr = np.abs(t - tc) < duration / 2.0
        if (
            in_tr.sum() < min_coverage * duration / cadence
            or sel.sum() < 3 * in_tr.sum()
        ):
            continue
        ts, fs = t[sel], f[sel] - 1.0
        shifts = np.arange(-duration, duration + 1e-9, cadence / 4.0)
        chi2 = (
            np.array(
                [
                    np.sum((fs - _trapezoid(ts - tc - s, duration, depth)) ** 2)
                    for s in shifts
                ]
            )
            / noise**2
        )
        i = int(np.argmin(chi2))
        if i in (0, len(shifts) - 1):
            continue  # the minimum is at the window edge: not timed
        y0, y1, y2 = chi2[i - 1], chi2[i], chi2[i + 1]
        h = shifts[1] - shifts[0]
        curv = (y0 - 2 * y1 + y2) / h**2
        if curv <= 0:
            continue
        best = shifts[i] + 0.5 * h * (y0 - y2) / (y0 - 2 * y1 + y2)
        out.append((int(n), float(tc + best), float(np.sqrt(2.0 / curv))))
    return out


def fit_ttv(timed: list[tuple[int, float, float]]) -> TTV:
    """The weighted linear ephemeris through the times, and what is left."""
    n = np.array([x[0] for x in timed], dtype=float)
    tm = np.array([x[1] for x in timed])
    sg = np.array([x[2] for x in timed])
    w = 1.0 / sg**2
    a = np.vstack([np.ones_like(n), n]).T
    coef = np.linalg.solve(a.T @ (a * w[:, None]), a.T @ (w * tm))
    oc = tm - a @ coef
    dof = max(1, len(n) - 2)
    return TTV(
        period_days=float(coef[1]),
        epoch_btjd=float(coef[0]),
        n_transits=len(n),
        epochs=[int(x) for x in n],
        times=tm.tolist(),
        sigmas_min=(sg * 1440).tolist(),
        o_minus_c_min=(oc * 1440).tolist(),
        chi2_dof=float(np.sum(oc**2 * w) / dof),
        max_abs_oc_min=float(np.max(np.abs(oc)) * 1440),
    )


_client: Any = None


def _ttv_star(
    tic: int,
    keys: dict[int, str],
    workdir: str,
    cfg: sr.SearchConfig,
    targets: list[dict[str, Any]],
    bin_minutes: float,
) -> list[dict[str, Any]]:
    global _client
    if _client is None:
        _client = lc.s3_client()
    try:
        curves = [
            lc.load_lc(lc.fetch_lc(_client, tic, s, Path(workdir), key=k))
            for s, k in sorted(keys.items())
        ]
        t, f, _ = ms.stitch(curves, cfg, bin_minutes)
        out = []
        keep = np.ones(len(t), dtype=bool)
        # Deepest first, each one's transits masked before the next is refined:
        # refined alone, TOI 216.02 locked onto 216.01 at half its period
        # (17.278 vs 34.555 d, measured).
        for x in sorted(targets, key=lambda x: -x["depth_ppm"]):
            dur = x["duration_hours"] / 24.0
            tk, fk = t[keep], f[keep]
            if x.get("refine"):
                p, e, d = refine_ephemeris(tk, fk, x["period_days"], dur)
                x = {**x, "period_days": p, "epoch_btjd": e, "depth_ppm": d * 1e6}
            timed = time_transits(
                tk, fk, x["period_days"], x["epoch_btjd"], dur, x["depth_ppm"] * 1e-6
            )
            phase = (t - x["epoch_btjd"]) % x["period_days"]
            keep &= np.minimum(phase, x["period_days"] - phase) > dur
            for _, tm, _ in timed:
                keep &= np.abs(t - tm) > dur
            if len(timed) < 4:
                out.append(
                    {
                        **x,
                        "tic_id": tic,
                        "ok": True,
                        "n_transits": len(timed),
                        "note": "under 4 transits timed",
                    }
                )
                continue
            out.append({**x, "tic_id": tic, "ok": True, **asdict(fit_ttv(timed))})
        return out
    except Exception as exc:  # noqa: BLE001 — recorded per star, never swallowed
        return [
            {**x, "tic_id": tic, "ok": False, "error": f"{type(exc).__name__}: {exc}"}
            for x in targets
        ]


def _toi_targets(toi_csv: Path, wanted: set[str]) -> dict[int, list[dict[str, Any]]]:
    lines = (ln for ln in toi_csv.open() if not ln.startswith("#"))
    out: dict[int, list[dict[str, Any]]] = {}
    for row in csv.DictReader(lines):
        if row["Full TOI ID"] not in wanted:
            continue
        out.setdefault(int(row["TIC"]), []).append(
            {
                "name": row["Full TOI ID"],
                "period_days": float(row["Orbital Period (days) Value"]),
                "epoch_btjd": float(row["Orbital Epoch Value"]),  # already BTJD
                "duration_hours": float(row["Transit Duration (hours) Value"]),
                "depth_ppm": float(row["Transit Depth Value"]),
                "refine": True,  # the catalogue ephemeris may be years away
            }
        )
    return out


def _run_targets(run: Path, min_transits: int) -> dict[int, list[dict[str, Any]]]:
    out: dict[int, list[dict[str, Any]]] = {}
    for r in map(json.loads, (run / "results.jsonl").read_text().splitlines()):
        if not r["ok"]:
            continue
        for s in r["result"]["signals"]:
            if s["n_transits"] >= min_transits:
                out.setdefault(r["tic_id"], []).append(
                    {
                        "name": f"{r['tic_id']}-{s['rank']}",
                        **{
                            k: s[k]
                            for k in (
                                "period_days",
                                "epoch_btjd",
                                "duration_hours",
                                "depth_ppm",
                            )
                        },
                    }
                )
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sectors", required=True)
    ap.add_argument("--toi-csv", type=Path)
    ap.add_argument(
        "--toi", action="append", default=[], help="a TOI to time, e.g. 216.01"
    )
    ap.add_argument(
        "--from-run", type=Path, help="time every signal of a multisector run"
    )
    ap.add_argument("--min-transits", type=int, default=4)
    ap.add_argument("--workdir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--bin-minutes", type=float, default=10.0)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args(argv)

    targets: dict[int, list[dict[str, Any]]] = {}
    if args.toi:
        targets.update(_toi_targets(args.toi_csv, set(args.toi)))
    if args.from_run:
        for tic, x in _run_targets(args.from_run, args.min_transits).items():
            targets.setdefault(tic, []).extend(x)
    client = lc.s3_client()
    manifests = {
        s: mf.load_or_list(client, s, args.workdir / f"manifest-s{s:04d}.json")
        for s in ms.parse_sectors(args.sectors)
    }
    keys = ms.eligible(manifests, 1)
    cfg = sr.SearchConfig(prewhiten_max_period_days=1.5)
    rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = [
            pool.submit(
                _ttv_star, tic, keys[tic], str(args.workdir), cfg, x, args.bin_minutes
            )
            for tic, x in sorted(targets.items())
            if tic in keys
        ]
        for fut in as_completed(futs):
            rows.extend(fut.result())
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "ttv.jsonl").open("w") as fh:
        for r in sorted(rows, key=lambda r: r["name"]):
            fh.write(json.dumps(r) + "\n")
    for r in sorted(rows, key=lambda r: -r.get("chi2_dof", 0) if r["ok"] else 0)[:20]:
        print(
            r["name"],
            r.get("n_transits"),
            f"chi2/dof {r.get('chi2_dof', float('nan')):.1f}",
            f"max|O-C| {r.get('max_abs_oc_min', float('nan')):.1f} min",
            r.get("error", r.get("note", "")),
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
