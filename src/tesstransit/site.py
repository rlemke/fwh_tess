"""Build the stand-alone results page from finished runs.

Usage::

    python -m tesstransit.site --runs /path/to/runs --workdir /path/to/scratch \\
        --planets planets-spline7a --dev dev-spline7a --holdout holdout-spline7a \\
        --variant "Pilot 5|1-day running median|planets-s1|pilot5-veto-fixed|holdout-seed2" \\
        --variant "7a|Cubic spline, knots by BIC|planets-spline7a|dev-spline7a|holdout-spline7a" \\
        --default 7a --out /path/to/site/tess-sector1

Writes one self-contained ``index.html``: the data is embedded, and only d3
and the fonts load from a CDN. It is the same page as the results dashboard,
over a fixed snapshot of the named runs.

Sky positions are not in the run outputs: each star's RA/Dec comes from the
primary header of its light curve (``RA_OBJ``/``DEC_OBJ``), read with one
ranged GET of the first 14 KB, and cached in ``<workdir>/coords-s0001.json``.
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
        return "Too weak in sector 1"
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


def toi_rows(run: Path, pos: dict[int, tuple[float, float]]) -> list[dict[str, str]]:
    out = []
    for r in csv.DictReader((run / "tois.csv").open()):
        ra, dec = pos[int(r["tic_id"])]
        out.append({
            "toi": r["toi"], "tic_id": r["tic_id"], "ra": _r(ra, 4), "dec": _r(dec, 4),
            "disposition": r["disposition"],
            "planet_like": "yes" if r["planet_like"] == "True" else "no",
            "sector1_tce": "yes" if r["tce_backed"] else "no",
            "period_days": _r(r["period_days"], 4), "depth_ppm": _r(r["depth_ppm"], 0),
            "duration_hours": _r(r["duration_hours"], 2),
            "forced_snr": _r(r["forced_snr"], 1), "our_snr": _r(r["sig_snr"], 1),
            "outcome": r["outcome"], "group": _group(r["outcome"]),
        })  # fmt: skip
    return out


def sample_rows(
    runs: list[tuple[str, Path]], pos: dict[int, tuple[float, float]]
) -> list[dict[str, str]]:
    """Every TCE, every no-TCE signal, and every no-TCE star with no signal."""
    out = []
    for sample, run in runs:
        matches = list(csv.DictReader((run / "matches.csv").open()))
        tce_stars = {int(r["tic_id"]) for r in matches if r["group"] == "tce"}
        with_signal = set()
        for r in matches:
            tic = int(r["tic_id"])
            ra, dec = pos[tic]
            base = {
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
            ra, dec = pos[tic]
            out.append({"sample": sample, "tic_id": str(tic), "ra": _r(ra, 4),
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


def render(data: dict[str, Any], built: str) -> str:
    page = files("tesstransit").joinpath("site_template.html").read_text()
    # "</" inside embedded JSON would end the <script> element early.
    blob = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    return page.replace("/*__DATA__*/null", blob).replace("/*__BUILT__*/", built)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument(
        "--workdir", type=Path, required=True, help="manifest + coords cache"
    )
    ap.add_argument("--sector", type=int, default=1)
    ap.add_argument("--planets", required=True)
    ap.add_argument("--dev", required=True)
    ap.add_argument("--holdout", required=True)
    ap.add_argument("--variant", action="append", default=[], required=True)
    ap.add_argument("--default", required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    variants = [v.split("|") for v in args.variant]
    if any(len(v) != 5 for v in variants):
        raise SystemExit("--variant is name|detrending|planets-run|dev-run|holdout-run")
    client = lc.s3_client()
    targets = mf.load_or_list(
        client, args.sector, args.workdir / f"manifest-s{args.sector:04d}.json"
    )
    runs = {k: args.runs / v for k, v in
            (("planets", args.planets), ("dev", args.dev), ("holdout", args.holdout))}  # fmt: skip
    tics = {
        json.loads(line)["tic_id"]
        for run in runs.values()
        for line in (run / "results.jsonl").read_text().splitlines()
    }
    pos = star_coords(
        client, targets, tics, args.workdir / f"coords-s{args.sector:04d}.json"
    )
    data = {
        "tois": toi_rows(runs["planets"], pos),
        "tce_sample": sample_rows(
            [("Development", runs["dev"]), ("Hold-out", runs["holdout"])], pos
        ),
        "variants": variant_rows(args.runs, variants, args.default),
    }
    built = datetime.now(UTC).strftime("%Y-%m-%d")
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "index.html").write_text(render(data, built))
    print(
        f"{args.out / 'index.html'}: {len(data['tois'])} TOIs, "
        f"{len(data['tce_sample'])} sample rows, {len(data['variants'])} variants",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
