"""Planet sample: re-search every star carrying a TOI observed in one sector.

Usage::

    python -m tesstransit.planets --sector 1 --tce-csv s0001_tcestats.csv \\
        --toi-csv tois.csv --workdir /path/to/scratch --out /path/to/run

TCE recovery measures agreement with SPOC, including on SPOC's false alarms.
This measures what bears on planets: the reference is the TOI catalogue, which
human vetters built (often from many sectors), not one sector's TCE list.

Every TOI listing the sector, on a star with a 2-minute light curve in it, is
in the sample. Each is labelled on two independent axes, never collapsed:

- **tce_backed**: SPOC made a TCE for it in THIS sector (period within 1%,
  alias-aware). The rest were found later, in more data or by another
  pipeline, and may simply not be detectable in one sector.
- **forced**: the known ephemeris measured on OUR detrended light curve --
  transits with data and their red-noise SNR. Fewer than 2 transits, or a
  forced SNR below the search threshold, is ``weak_in_sector``: a miss there
  says little about the search. EXCEPT when SPOC detected it in this sector
  at SNR >= the threshold: then the signal is in the data and our detrending
  lost it, which is ``missed:weak_after_detrend`` -- our failure, not the
  data's (measured: DS Tuc Ab, SPOC SNR 25, forced 1.0 on an active star).

The ephemeris used is the TOI's period with, when the TOI is TCE-backed,
SPOC's in-sector epoch. TOI epochs are often fitted years later, and
propagating them back is only as good as the period: TOI 121.01's catalogue
ephemeris lands 12.4 h from its sector-1 transits, 50 orbits back, where both
SPOC and this search found them. A TOI with no sector TCE keeps its catalogue
epoch, and its propagated uncertainty is reported (``epoch_sigma_h``) and
flagged when it exceeds the match tolerance.

The headline is recovery of planet-like TOIs (CP, KP, PC) that are in the
data. EB / V / IS / FP TOIs are searched and reported, but kept out of it.

⚠️ The forced SNR is measured at the known ephemeris, so it is an upper bound
on what a blind search can reach.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

from . import _lightcurve as lc
from . import _manifest as mf
from . import _match as mt
from . import _search as sr
from .pilot import _result_from_dict, _search_one

PLANET_LIKE = ("CP", "KP", "PC")


@dataclass
class TOI:
    toi: str
    tic_id: int
    disposition: str
    period_days: float | None
    epoch_btjd: float
    duration_hours: float
    depth_ppm: float
    sectors: tuple[int, ...]
    pipelines: str
    period_err: float | None = None
    epoch_err: float | None = None
    epoch_source: str = "toi"  # "tce" once replaced by SPOC's in-sector epoch

    def as_tce(self) -> mt.TCE:
        """The TOI in the matcher's shape (planet_num orders multi-TOI stars)."""
        return mt.TCE(
            tce_id=self.toi,
            tic_id=self.tic_id,
            planet_num=int(self.toi.split(".")[1]),
            period_days=self.period_days or 0.0,
            epoch_btjd=self.epoch_btjd,
            duration_hours=self.duration_hours,
            depth_ppm=self.depth_ppm,
            snr=None,
            n_transits=None,
        )


def _float(v: str | None) -> float | None:
    return float(v) if v and v.strip() else None


def epoch_sigma_hours(toi: TOI, t_mid: float) -> float | None:
    """The catalogue epoch's uncertainty propagated to ``t_mid``, in hours."""
    if toi.period_days is None or toi.epoch_err is None or toi.period_err is None:
        return None
    n = abs(round((t_mid - toi.epoch_btjd) / toi.period_days))
    return 24.0 * (toi.epoch_err**2 + (n * toi.period_err) ** 2) ** 0.5


def local_ephemeris(toi: TOI, tce: mt.TCE | None) -> TOI:
    """The TOI with SPOC's in-sector epoch when a sector TCE backs it.

    The TCE's period may be an alias of the TOI's; its epoch is still a real
    transit, so only the epoch is taken and the TOI's period is kept.
    """
    if tce is None:
        return toi
    return replace(toi, epoch_btjd=tce.epoch_btjd, epoch_source="tce")


def parse_tois(text: str) -> list[TOI]:
    """MAST's ``tois.csv``. Epochs are BTJD there (measured: 925-3988, none > 2e6)."""
    lines = (ln for ln in io.StringIO(text) if not ln.startswith("#"))
    out: list[TOI] = []
    for row in csv.DictReader(lines):
        period = row["Orbital Period (days) Value"].strip()
        epoch = float(row["Orbital Epoch Value"])
        if epoch > 2_000_000:  # BJD, not BTJD: convert rather than mis-phase
            epoch -= 2_457_000
        out.append(
            TOI(
                toi=row["Full TOI ID"],
                tic_id=int(row["TIC"]),
                disposition=row["TOI Disposition"],
                # A single-transit TOI has period 0 or blank: not searchable
                # by a >= 2-transit search, reported as such.
                period_days=float(period) if period and float(period) > 0 else None,
                epoch_btjd=epoch,
                duration_hours=float(row["Transit Duration (hours) Value"]),
                depth_ppm=float(row["Transit Depth Value"]),
                sectors=tuple(int(s) for s in re.findall(r"\d+", row["Sectors"])),
                pipelines=row["Detection Pipeline(s)"],
                period_err=_float(row.get("Orbital Period Uncertainty")),
                epoch_err=_float(row.get("Epoch Uncertainty")),
            )
        )
    return out


def _tce_backed(toi: TOI, tces: list[mt.TCE]) -> mt.TCE | None:
    """The sector TCE at (an alias of) this TOI's period, if SPOC made one."""
    if toi.period_days is None:
        return None
    for t in tces:
        if any(
            abs(t.period_days / (k * toi.period_days) - 1) < 0.01
            for k in mt.DEFAULT_ALIASES
        ):
            return t
    return None


def classify(
    toi: TOI,
    row: mt.MatchRow | None,
    forced: dict[str, Any] | None,
    cfg: sr.SearchConfig,
    spoc_snr: float | None = None,
) -> str:
    """One outcome per TOI. ``row`` is None when the star failed or, for a
    single-transit TOI, when there was no period to match against."""
    if toi.period_days is None:
        return "single_transit"
    if row is None:
        return "not_searched"
    if row.bucket.startswith("recovered"):
        return row.bucket
    if not cfg.min_period_days <= toi.period_days <= cfg.max_period_days:
        return "outside_search_grid"
    if forced is None or forced["n_transits"] < cfg.min_transits:
        return "weak_in_sector:few_transits"
    if forced["snr_red"] is None or forced["snr_red"] < cfg.snr_threshold:
        if spoc_snr is not None and spoc_snr >= cfg.snr_threshold:
            return "missed:weak_after_detrend"
        return "weak_in_sector:low_forced_snr"
    if row.cause == "search_exhausted":
        return "missed:search_exhausted"
    return "missed:not_detected"


def _search_all(
    args: argparse.Namespace,
    cfg: sr.SearchConfig,
    targets: dict[int, str],
    by_star: dict[int, list[TOI]],
) -> list[dict[str, Any]]:
    """Search every sample star, forcing each of its TOI ephemerides."""
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
                [
                    (t.toi, t.period_days, t.epoch_btjd, t.duration_hours)
                    for t in ts
                    if t.period_days
                ],
            )
            for tic, ts in sorted(by_star.items())
        ]
        for i, fut in enumerate(as_completed(futs), 1):
            records.append(fut.result())
            if i % 20 == 0:
                print(f"  {i}/{len(futs)} ({time.time() - t1:.0f}s)", file=sys.stderr)
    return records


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sector", type=int, required=True)
    ap.add_argument("--tce-csv", type=Path, required=True)
    ap.add_argument("--toi-csv", type=Path, required=True)
    ap.add_argument("--workdir", type=Path, required=True, help="light-curve scratch")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--keep-lc", action="store_true")
    ap.add_argument(
        "--reuse",
        action="store_true",
        help="re-match the run's existing results.jsonl instead of searching again",
    )
    args = ap.parse_args(argv)

    cfg = sr.SearchConfig()
    cfg.validate()
    args.out.mkdir(parents=True, exist_ok=True)
    client = lc.s3_client()
    targets = mf.load_or_list(
        client, args.sector, args.workdir / f"manifest-s{args.sector:04d}.json"
    )
    tces = mt.parse_tce_csv(args.tce_csv.read_text())
    tces_by: dict[int, list[mt.TCE]] = {}
    for t in tces:
        tces_by.setdefault(t.tic_id, []).append(t)
    all_tois = [
        t for t in parse_tois(args.toi_csv.read_text()) if args.sector in t.sectors
    ]
    catalogue = {t.toi: t for t in all_tois}
    backing = {t.toi: _tce_backed(t, tces_by.get(t.tic_id, [])) for t in all_tois}
    tois = [local_ephemeris(t, backing[t.toi]) for t in all_tois if t.tic_id in targets]
    by_star: dict[int, list[TOI]] = {}
    for t in tois:
        by_star.setdefault(t.tic_id, []).append(t)
    print(
        f"{len(all_tois)} TOIs list sector {args.sector}; {len(tois)} on {len(by_star)} "
        f"stars with a 2-min light curve",
        file=sys.stderr,
    )

    records: list[dict[str, Any]] = []
    t1 = time.time()
    if args.reuse:
        records = [
            json.loads(ln)
            for ln in (args.out / "results.jsonl").read_text().splitlines()
        ]
        missing = set(by_star) - {r["tic_id"] for r in records}
        if missing:
            raise SystemExit(
                f"--reuse: {len(missing)} sample stars not in results.jsonl"
            )
    else:
        records = _search_all(args, cfg, targets, by_star)
    wall = time.time() - t1
    records.sort(key=lambda r: r["tic_id"])
    if not args.reuse:
        with (args.out / "results.jsonl").open("w") as fh:
            for r in records:
                fh.write(json.dumps(r) + "\n")

    rec_by = {r["tic_id"]: r for r in records}
    out_rows: list[dict[str, Any]] = []
    unmatched = 0
    for tic, ts in sorted(by_star.items()):
        r = rec_by[tic]
        result = _result_from_dict(r["result"]) if r["ok"] else None
        searchable = [t for t in ts if t.period_days]
        rows = (
            mt.match_star(tic, result, [t.as_tce() for t in searchable])
            if result
            else []
        )
        by_toi = {m.tce.tce_id: m for m in rows if m.tce}
        unmatched += sum(1 for m in rows if m.bucket == "unmatched")
        for t in ts:
            m = by_toi.get(t.toi)
            forced = (r.get("forced") or {}).get(t.toi)
            tce = backing[t.toi]
            sigma = (
                epoch_sigma_hours(catalogue[t.toi], sum(r["t_span"]) / 2)
                if r["ok"] and t.epoch_source == "toi"
                else None
            )
            dt = None
            if m and m.signal:
                wrap = min(m.signal.period_days, t.period_days)
                dt = (
                    (m.signal.epoch_btjd - t.epoch_btjd + 0.5 * wrap) % wrap
                    - 0.5 * wrap
                ) * 24
            out_rows.append({
                "toi": t.toi, "tic_id": tic, "disposition": t.disposition,
                "planet_like": t.disposition in PLANET_LIKE,
                "tce_backed": tce.tce_id if tce else "",
                "spoc_snr": tce.snr if tce else None,
                "epoch_source": t.epoch_source,
                "epoch_sigma_h": sigma,
                "epoch_uncertain": sigma is not None
                and sigma > mt.MatchConfig().epoch_tol_durations * t.duration_hours,
                "pipelines": t.pipelines,
                "period_days": t.period_days, "duration_hours": t.duration_hours,
                "depth_ppm": t.depth_ppm,
                "forced_n_transits": forced["n_transits"] if forced else None,
                "forced_snr": forced["snr_red"] if forced else None,
                "outcome": classify(
                    t, m if result else None, forced, cfg, tce.snr if tce else None
                ),
                "alias": m.alias if m else None,
                "sig_period": m.signal.period_days if m and m.signal else None,
                "sig_snr": m.signal.snr if m and m.signal else None,
                "epoch_offset_h": dt,
                "error": "" if r["ok"] else r["error"],
            })  # fmt: skip
    with (args.out / "tois.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)

    def table(rows: list[dict[str, Any]]) -> dict[str, int]:
        return dict(sorted(Counter(r["outcome"] for r in rows).items()))

    planet = [r for r in out_rows if r["planet_like"]]
    groups = {
        "planet_like": {
            "tce_backed": table([r for r in planet if r["tce_backed"]]),
            "no_sector_tce": table([r for r in planet if not r["tce_backed"]]),
        },
        "other_dispositions": table([r for r in out_rows if not r["planet_like"]]),
    }
    summary = {
        "sector": args.sector,
        "scope": "every TOI listing the sector, on a star with a 2-min light curve",
        "tois_listing_sector": len(all_tois),
        "tois_in_sample": len(tois),
        "stars": len(by_star),
        "failed": {r["tic_id"]: r["error"] for r in records if not r["ok"]},
        "outcomes": groups,
        "signals_matching_no_toi": unmatched,
        "search_config": cfg.to_dict(),
        "match_config": asdict(mt.MatchConfig()),
        # A --reuse run did no searching: its wall time is not a search time.
        "seconds": None if args.reuse else {"wall": wall, "workers": args.workers},
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary["outcomes"], indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
