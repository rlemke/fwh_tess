"""Build the stand-alone results page from finished runs.

Usage::

    python -m tesstransit.site --runs /path/to/runs --workdir /path/to/scratch \\
        --planets 1=planets-spline7a-v2 --planets 2=s2-planets-v2 \\
        --sample 1=Development=dev-spline7a --sample 1=Hold-out=holdout-spline7a \\
        --sample 2=Sample=s2-pilot \\
        --variant "Pilot 5|1-day running median|planets-s1|pilot5-veto-fixed|holdout-seed2" \\
        --variant "7a|Cubic spline, knots by BIC|planets-spline7a|dev-spline7a|holdout-spline7a" \\
        --default 7a --out /path/to/site/tess-sector1

Every row carries its ``sector``; the ``--variant`` detrending scorecards are
sector-1 experiments.

``--stitched "LABEL|SETTINGS|RUN|INVERTED-RUN|yes"`` adds a sample of the
13-sector stitched search (``tesstransit.multisector``): the run as observed,
the same stars with their light curves turned upside down (no real transit
survives that, so every signal there is a false alarm), and whether the
settings are the adopted ones. It adds three tables, ``stitched_samples``,
``stitched_planets`` and ``stitched_signals``; no sky positions are needed.

Writes one self-contained ``index.html``: the data is embedded, and only d3
and the fonts load from a CDN. It is the same page as the results dashboard,
over a fixed snapshot of the named runs.

Sky positions are not in the run outputs: each star's RA/Dec comes from the
primary header of its light curve (``RA_OBJ``/``DEC_OBJ``), read with one
ranged GET of the first 14 KB, and cached in ``<workdir>/coords-s000N.json``.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from importlib.resources import files
from pathlib import Path
from typing import Any

from . import _lightcurve as lc
from . import _manifest as mf

HEADER_BYTES = 14400  # five FITS blocks: the whole primary header of a SPOC _lc.fits


def _group(outcome: str) -> str:
    if outcome.startswith("recovered"):
        return "Recovered"
    if outcome.startswith("missed"):
        return "Missed"
    if outcome.startswith("weak_in_sector"):
        return "Too weak in one sector"
    if outcome == "outside_search_grid":
        return "Period outside search"
    if outcome == "single_transit":
        return "Single transit"
    return "Not searched"


def _r(x: str | None, nd: int) -> str:
    """A CSV value rounded for display, kept as text like the dashboard's files."""
    if x in ("", None, "None"):
        return ""
    v = round(float(x), nd)
    return str(int(v)) if nd == 0 else str(v)


def star_coords(
    client: Any, targets: dict[int, str], tics: set[int], cache: Path
) -> dict[int, tuple[float, float]]:
    """RA/Dec (degrees) per TIC, from cache or the light curves' headers."""
    from astropy.io import fits

    have: dict[int, tuple[float, float]] = {}
    if cache.exists():
        have = {int(k): tuple(v[:2]) for k, v in json.loads(cache.read_text()).items()}

    def one(tic: int) -> tuple[int, tuple[float, float]]:
        body = client.get_object(
            Bucket=lc.BUCKET, Key=targets[tic], Range=f"bytes=0-{HEADER_BYTES - 1}"
        )["Body"].read()
        end = body.index(b"END" + b" " * 77)
        h = fits.Header.fromstring(body[: end + 80].decode("ascii"))
        return tic, (float(h["RA_OBJ"]), float(h["DEC_OBJ"]))

    todo = sorted(tics - set(have))
    if todo:
        with ThreadPoolExecutor(16) as ex:
            have.update(ex.map(one, todo))
        cache.write_text(json.dumps({str(k): v for k, v in have.items()}))
    return have


def toi_rows(
    sector: int, run: Path, pos: dict[tuple[int, int], tuple[float, float]]
) -> list[dict[str, str]]:
    out = []
    for r in csv.DictReader((run / "tois.csv").open()):
        ra, dec = pos[(sector, int(r["tic_id"]))]
        out.append({
            "sector": str(sector),
            "toi": r["toi"], "tic_id": r["tic_id"], "ra": _r(ra, 4), "dec": _r(dec, 4),
            "disposition": r["disposition"],
            "planet_like": "yes" if r["planet_like"] == "True" else "no",
            "sector_tce": "yes" if r["tce_backed"] else "no",
            "period_days": _r(r["period_days"], 4), "depth_ppm": _r(r["depth_ppm"], 0),
            "duration_hours": _r(r["duration_hours"], 2),
            "forced_snr": _r(r["forced_snr"], 1), "our_snr": _r(r["sig_snr"], 1),
            "outcome": r["outcome"], "group": _group(r["outcome"]),
        })  # fmt: skip
    return out


def sample_rows(
    runs: list[tuple[int, str, Path]], pos: dict[tuple[int, int], tuple[float, float]]
) -> list[dict[str, str]]:
    """Every TCE, every no-TCE signal, and every no-TCE star with no signal."""
    out = []
    for sector, sample, run in runs:
        matches = list(csv.DictReader((run / "matches.csv").open()))
        tce_stars = {int(r["tic_id"]) for r in matches if r["group"] == "tce"}
        with_signal = set()
        for r in matches:
            tic = int(r["tic_id"])
            ra, dec = pos[(sector, tic)]
            base = {
                "sector": str(sector),
                "sample": sample,
                "tic_id": str(tic),
                "ra": _r(ra, 4),
                "dec": _r(dec, 4),
            }
            if r["group"] == "quiet":
                with_signal.add(tic)
                out.append({**base, "kind": "No-TCE star", "tce_id": "",
                            "period_days": _r(r["sig_period"], 4),
                            "depth_ppm": _r(r["sig_depth_ppm"], 0), "spoc_snr": "",
                            "our_snr": _r(r["sig_snr"], 1), "bucket": "unmatched",
                            "toi_disposition": "", "status": "Signal on a no-TCE star"})  # fmt: skip
            elif r["tce_id"]:
                b = r["bucket"]
                status = (
                    "Recovered" if b.startswith("recovered")
                    else "Outside search" if r["cause"] == "outside_search_grid"
                    else "Missed"
                )  # fmt: skip
                out.append({**base, "kind": "TCE", "tce_id": r["tce_id"],
                            "period_days": _r(r["tce_period"], 4),
                            "depth_ppm": _r(r["tce_depth_ppm"], 0),
                            "spoc_snr": _r(r["tce_snr"], 1), "our_snr": _r(r["sig_snr"], 1),
                            "bucket": b + (":" + r["cause"] if r["cause"] else ""),
                            "toi_disposition": r["toi_disposition"], "status": status})  # fmt: skip
        for line in (run / "results.jsonl").read_text().splitlines():
            tic = json.loads(line)["tic_id"]
            if tic in tce_stars or tic in with_signal:
                continue
            ra, dec = pos[(sector, tic)]
            out.append({"sector": str(sector), "sample": sample, "tic_id": str(tic), "ra": _r(ra, 4),
                        "dec": _r(dec, 4), "kind": "No-TCE star", "tce_id": "",
                        "period_days": "", "depth_ppm": "", "spoc_snr": "", "our_snr": "",
                        "bucket": "", "toi_disposition": "", "status": "No signal"})  # fmt: skip
    return out


def variant_rows(
    runs: Path, variants: list[list[str]], default: str
) -> list[dict[str, str]]:
    """One scorecard row per ``name|detrending|planets|dev|holdout`` variant."""

    def planets(d: Path) -> tuple[int, int]:
        rs = [
            r for r in csv.DictReader((d / "tois.csv").open())
            if r["planet_like"] == "True" and r["tce_backed"]
        ]  # fmt: skip
        return sum(r["outcome"].startswith("recovered") for r in rs), len(rs)

    def pilot(d: Path) -> tuple[int, int, int, float]:
        s = json.loads((d / "summary.json").read_text())
        g = s["tce_group"]
        in_grid = sum(v for k, v in g.items() if k in ("recovered", "recovered_alias"))
        in_grid += sum(
            v for k, v in g.items()
            if k.startswith("missed:") and k != "missed:outside_search_grid"
        )  # fmt: skip
        found = g.get("recovered", 0) + g.get("recovered_alias", 0)
        quiet = s["quiet_group"]["stars_with_any_signal"]
        return found, in_grid, quiet, s["seconds"]["per_star_median"]

    out = []
    for i, (name, det, p, d, h) in enumerate(variants, 1):
        pr, pt = planets(runs / p)
        dr, dt, dq, sec = pilot(runs / d)
        hr, ht, hq, _ = pilot(runs / h)
        out.append({
            "order": str(i), "variant": name, "detrending": det,
            "planets_recovered": str(pr), "planets_total": str(pt),
            "dev_tces_recovered": str(dr), "dev_tces_total": str(dt),
            "dev_quiet_with_signal": str(dq),
            "holdout_tces_recovered": str(hr), "holdout_tces_total": str(ht),
            "holdout_quiet_with_signal": str(hq),
            "seconds_per_star": _r(sec, 1), "default": "yes" if name == default else "no",
        })  # fmt: skip
    return out


PLANET_DISPOSITIONS = ("CP", "KP", "PC")


def _wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """The Wilson score interval for ``k`` of ``n``, as fractions."""
    if n == 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / d
    return max(0.0, c - h), min(1.0, c + h)


def _stitched_run(run: Path) -> tuple[list[dict[str, str]], set[int], set[int]]:
    """``(matches rows, stars with no TCE, of those the stars with a signal)``."""
    matches = list(csv.DictReader((run / "matches.csv").open()))
    tce_stars = {int(r["tic_id"]) for r in matches if r["group"] == "tce"}
    searched = {
        json.loads(line)["tic_id"]
        for line in (run / "results.jsonl").read_text().splitlines()
        if json.loads(line)["ok"]
    }
    quiet = searched - tce_stars
    with_signal = {int(r["tic_id"]) for r in matches if r["group"] == "quiet"}
    return matches, quiet, with_signal & quiet


def stitched_rows(runs: Path, specs: list[list[str]]) -> dict[str, list[dict[str, str]]]:
    """The three stitched tables from ``label|settings|run|inverted-run|adopted`` specs.

    A planet is a CP/KP/PC TOI on a sampled star that a TCE or a signal
    matched (``matches.csv``), counted once per sample: found if any of its
    rows was recovered (an alias counts) or is a signal SPOC made no TCE for,
    else missed. TCEs whose period is outside the search grid are left out.
    """
    samples, planets, signals = [], [], []
    for order, (label, settings, run, inv, adopted) in enumerate(specs, 1):
        matches, quiet, sig = _stitched_run(runs / run)
        _, inv_quiet, inv_sig = _stitched_run(runs / inv)
        by_toi: dict[str, dict[str, str]] = {}
        for r in matches:
            if r["toi_disposition"] not in PLANET_DISPOSITIONS or r["group"] != "tce":
                continue
            if r["cause"] == "outside_search_grid":
                continue
            b = r["bucket"]
            found = b.startswith("recovered") or b == "unmatched"
            row = {
                "order": str(order), "sample": label, "toi": r["toi"],
                "disposition": r["toi_disposition"], "tic_id": r["tic_id"],
                "period_days": _r(r["tce_period"] or r["sig_period"], 4),
                "our_period_days": _r(r["sig_period"], 4),
                "our_snr": _r(r["sig_snr"], 1), "our_sde": _r(r["sig_sde"], 1),
                "status": ("Found at an alias" if b == "recovered_alias"
                           else "Found, no SPOC TCE" if b == "unmatched"
                           else "Found" if found else "Missed"),
            }  # fmt: skip
            old = by_toi.get(r["toi"])
            if old is None or (old["status"] == "Missed" and found):
                by_toi[r["toi"]] = row
        planets += by_toi.values()
        n_found = sum(p["status"] != "Missed" for p in by_toi.values())
        lo, hi = _wilson(len(inv_sig), len(inv_quiet))
        samples.append({
            "order": str(order), "sample": label, "settings": settings,
            "adopted": "yes" if adopted == "yes" else "no",
            "quiet_stars": str(len(quiet)), "quiet_signal": str(len(sig)),
            "inv_stars": str(len(inv_quiet)), "inv_signal": str(len(inv_sig)),
            "inv_pct": f"{100 * len(inv_sig) / len(inv_quiet):.1f}" if inv_quiet else "",
            "inv_lo": f"{100 * lo:.1f}", "inv_hi": f"{100 * hi:.1f}",
            "planets_found": str(n_found) if by_toi else "",
            "planets_total": str(len(by_toi)) if by_toi else "",
        })  # fmt: skip
        if adopted != "yes":
            continue  # candidates are listed only from the adopted settings
        for r in matches:
            if r["group"] != "quiet":
                continue
            signals.append({
                "order": str(order), "sample": label, "tic_id": r["tic_id"],
                "period_days": _r(r["sig_period"], 4), "depth_ppm": _r(r["sig_depth_ppm"], 0),
                "duration_hours": _r(r["sig_duration_h"], 1),
                "our_snr": _r(r["sig_snr"], 1), "our_sde": _r(r["sig_sde"], 1),
                "toi": r["toi"], "toi_disposition": r["toi_disposition"],
            })  # fmt: skip
    return {"stitched_samples": samples, "stitched_planets": planets, "stitched_signals": signals}


def render(data: dict[str, Any], built: str) -> str:
    page = files("tesstransit").joinpath("site_template.html").read_text()
    # "</" inside embedded JSON would end the <script> element early.
    blob = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    return page.replace("/*__DATA__*/null", blob).replace("/*__BUILT__*/", built)


def _split(values: list[str], parts: int, flag: str) -> list[list[str]]:
    out = [v.split("=") for v in values]
    if any(len(v) != parts for v in out):
        raise SystemExit(f"{flag} expects {parts} '='-separated parts, got {values}")
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument(
        "--runs", type=Path, default=Path("."), help="base for relative run names"
    )
    ap.add_argument(
        "--workdir", type=Path, required=True, help="manifests + coords caches"
    )
    ap.add_argument("--planets", action="append", required=True, help="SECTOR=RUN")
    ap.add_argument("--sample", action="append", required=True, help="SECTOR=LABEL=RUN")
    ap.add_argument(
        "--variant", action="append", default=[], help="sector-1 scorecard row"
    )
    ap.add_argument("--default", default="")
    ap.add_argument(
        "--stitched", action="append", default=[],
        help="stitched sample: label|settings|run|inverted-run|adopted (yes/no)",
    )  # fmt: skip
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    stitched = [v.split("|") for v in args.stitched]
    if any(len(v) != 5 for v in stitched):
        raise SystemExit("--stitched is label|settings|run|inverted-run|adopted")
    variants = [v.split("|") for v in args.variant]
    if any(len(v) != 5 for v in variants):
        raise SystemExit("--variant is name|detrending|planets-run|dev-run|holdout-run")
    planets = [(int(s), args.runs / r) for s, r in _split(args.planets, 2, "--planets")]
    samples = [
        (int(s), lbl, args.runs / r) for s, lbl, r in _split(args.sample, 3, "--sample")
    ]
    client = lc.s3_client()
    pos: dict[tuple[int, int], tuple[float, float]] = {}
    for sector in sorted({s for s, _ in planets} | {s for s, _, _ in samples}):
        runs = [r for s, r in planets if s == sector]
        runs += [r for s, _, r in samples if s == sector]
        tics = {
            json.loads(line)["tic_id"]
            for run in runs
            for line in (run / "results.jsonl").read_text().splitlines()
        }
        targets = mf.load_or_list(
            client, sector, args.workdir / f"manifest-s{sector:04d}.json"
        )
        cache = args.workdir / f"coords-s{sector:04d}.json"
        pos.update(
            {
                (sector, t): v
                for t, v in star_coords(client, targets, tics, cache).items()
            }
        )
    data = {
        "tois": [row for sector, run in planets for row in toi_rows(sector, run, pos)],
        "tce_sample": sample_rows(samples, pos),
        "variants": variant_rows(args.runs, variants, args.default),
        **stitched_rows(args.runs, stitched),
    }
    built = datetime.now(UTC).strftime("%Y-%m-%d")
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "index.html").write_text(render(data, built))
    print(
        f"{args.out / 'index.html'}: {len(data['tois'])} TOIs, "
        f"{len(data['tce_sample'])} sample rows, {len(data['variants'])} variants, "
        f"{len(data['stitched_samples'])} stitched samples",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
