"""Signals our search found that no catalogue claims.

Every signal in the given runs is checked against everything published about
its star: every SPOC TCE (single-sector AND multi-sector tables), every TOI
(any disposition, false positives included), every community TOI (CTOI) and
every confirmed planet in the NASA Exoplanet Archive. A signal is CLAIMED when
any of them lists the same star at the same period or an alias of it (x1, x2,
x1/2, x3, x1/3, within 0.5%: the matcher's own tolerance). What is left is
UNCLAIMED: found here and, as far as these catalogues show, nowhere else.

⚠️ Unclaimed is not a discovery. Measured on the same stars searched upside
down, the no-TCE signal rate equals the false-alarm rate (experiments.md), so
AS A GROUP these signals are what noise and stellar variability produce. The
report therefore says, for each one, the two things that separate a signal
worth a look from the expected noise:

- **upside down**: whether the inverted search of the same star found a
  signal at the same period. If it did, the dips are as strong as the bumps:
  variability, not a transit.
- **independent sectors**: in how many separately searched sectors the same
  period came back. Noise does not repeat at one period in independent data;
  a planet does.

Nothing here is vetted (no centroid test, no follow-up), and none of it is
reported as a planet.

    python -m tesstransit.unclaimed --runs RUNS --catalogs DIR \
        --sample 45=s45-pilot-pwC:s45-pilot-inv-pwC --planets 45=s45-planets-pwC \
        --stitched "South A=ms26-fresh:ms26-fresh-inv" --out docs/unclaimed.md --csv docs/data/unclaimed.csv
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import io
import json
import math
import os
import sys
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

# Grouping one star's finds across runs uses the matcher's own aliases and
# tolerance. Deciding whether a catalogue CLAIMS a signal errs broad, since
# calling a known signal new is the worse mistake: integer ratios to 6, the
# common rational ones, and a tolerance no tighter than the transits allow.
ALIASES: tuple[float, ...] = (1.0, 2.0, 0.5, 3.0, 1.0 / 3.0)
CLAIM_ALIASES: tuple[float, ...] = tuple(
    sorted(
        {float(k) for k in range(1, 7)}
        | {1.0 / k for k in range(2, 7)}
        | {1.5, 2.0 / 3.0, 2.5, 0.4}
    )
)
REL_TOL = 0.005


@dataclass
class Claim:
    source: str  # "TCE s0045", "TCE s0001-s0096", "TOI 123.01 (PC)", "CTOI ...", "planet ..."
    period_days: float


@dataclass
class Found:
    """One signal in one run."""

    run: str
    kind: str  # sample | planets | stitched
    sector: str  # "45" or the stitched label
    tic_id: int
    period_days: float
    epoch_btjd: float
    duration_hours: float
    depth_ppm: float
    snr: float
    sde: float
    n_transits: int
    claimed_by: list[str] = field(default_factory=list)
    upside_down: str = (
        "not run"  # "same period" | "other period" | "nothing" | "not run"
    )


def alias_of(
    p: float, ref: float, aliases: tuple[float, ...] = ALIASES, tol: float = REL_TOL
) -> float | None:
    for k in aliases:
        if ref > 0 and abs(p / (k * ref) - 1.0) <= tol:
            return k
    return None


def claim_tol(f: Found) -> float:
    """How precisely f's own transits fix its period: each is timed to about a
    duration, over (n - 1) periods. Two transits of an 8 h box 10.6 d apart pin
    the period to ~3%, not to the matcher's 0.5%."""
    n = max(2, f.n_transits)
    return max(REL_TOL, (f.duration_hours / 24.0) / (f.period_days * (n - 1)))


def _num(v: str | None) -> float | None:
    try:
        x = float((v or "").strip())
    except ValueError:
        return None
    return x if math.isfinite(x) and x > 0 else None


def read_tce_claims(paths: Iterable[Path]) -> dict[int, list[Claim]]:
    """(tic -> claims) from SPOC DV TCE tables; tolerant of odd rows."""
    out: dict[int, list[Claim]] = defaultdict(list)
    for p in paths:
        label = "TCE " + p.name.split("_tcestats")[0]
        lines = (ln for ln in io.StringIO(p.read_text()) if not ln.startswith("#"))
        for row in csv.DictReader(lines):
            per = _num(row.get("tce_period"))
            if per is None:
                continue
            out[int(row["ticid"])].append(Claim(label, per))
    return out


def read_toi_claims(path: Path) -> dict[int, list[Claim]]:
    out: dict[int, list[Claim]] = defaultdict(list)
    lines = (ln for ln in io.StringIO(path.read_text()) if not ln.startswith("#"))
    for row in csv.DictReader(lines):
        per = _num(row.get("Orbital Period (days) Value"))
        if per is not None:
            out[int(row["TIC"])].append(
                Claim(
                    f"TOI {row['Full TOI ID']} ({row.get('TOI Disposition', '').strip() or '?'})",
                    per,
                )
            )
    return out


def read_ctoi_claims(path: Path) -> dict[int, list[Claim]]:
    out: dict[int, list[Claim]] = defaultdict(list)
    for row in csv.DictReader(io.StringIO(path.read_text()), delimiter="|"):
        per = _num(row.get("Period (days)"))
        if per is not None and (row.get("TIC ID") or "").strip().isdigit():
            out[int(row["TIC ID"])].append(Claim(f"CTOI {row['CTOI']}", per))
    return out


def read_planet_claims(path: Path) -> dict[int, list[Claim]]:
    out: dict[int, list[Claim]] = defaultdict(list)
    for row in csv.DictReader(io.StringIO(path.read_text())):
        tic = (row.get("tic_id") or "").replace("TIC", "").strip()
        per = _num(row.get("pl_orbper"))
        if tic.isdigit() and per is not None:
            out[int(tic)].append(Claim(f"planet {row['pl_name']}", per))
    return out


def signals_of(run_dir: Path) -> dict[int, list[dict[str, Any]]]:
    out: dict[int, list[dict[str, Any]]] = {}
    for line in (run_dir / "results.jsonl").read_text().splitlines():
        rec = json.loads(line)
        if rec.get("ok"):
            out[int(rec["tic_id"])] = rec["result"]["signals"]
    return out


def _inverted_one(args: tuple[int, int, dict[str, Any], str]) -> dict[str, Any]:
    """Search one star's sector upside down with the config its result recorded."""
    from dataclasses import fields, replace

    from . import _search as sr
    from . import pilot

    tic, sector, config, workdir = args
    known = {f.name for f in fields(sr.SearchConfig)}
    cfg = sr.SearchConfig(
        **{
            k: (tuple(v) if isinstance(v, list) else v)
            for k, v in config.items()
            if k in known
        }
    )
    cfg = replace(cfg, invert_flux=True)
    rec = pilot._search_one(tic, sector, None, workdir, False, cfg)
    return {
        "tic_id": tic,
        "sector": sector,
        "ok": rec["ok"],
        "error": rec.get("error"),
        "signals": rec["result"]["signals"] if rec["ok"] else [],
    }


def fill_inverted(
    found: list[Found], runs_dir: Path, cache: Path, workdir: Path, workers: int
) -> int:
    """Give every unclaimed single-sector find with no upside-down run one.

    The planet runs have no inverted control, and variability repeats across
    sectors exactly as a planet does, so without this a repeat proves nothing.
    Results are cached by (tic, sector) in `cache`, so a re-run costs nothing.
    """
    from concurrent.futures import ProcessPoolExecutor

    done: dict[tuple[int, int], dict[str, Any]] = {}
    if cache.exists():
        for line in cache.read_text().splitlines():
            r = json.loads(line)
            done[(r["tic_id"], r["sector"])] = r
    need = {
        (f.tic_id, int(f.sector), f.run)
        for f in found
        if not f.claimed_by and f.upside_down == "not run" and f.kind != "stitched"
    }
    todo = []
    configs: dict[str, dict[int, dict[str, Any]]] = {}
    for tic, sector, run in sorted(need):
        if (tic, sector) in done:
            continue
        if run not in configs:
            configs[run] = {}
            for line in (runs_dir / run / "results.jsonl").read_text().splitlines():
                rec = json.loads(line)
                if rec.get("ok"):
                    configs[run][int(rec["tic_id"])] = rec["result"]["config"]
        todo.append((tic, sector, configs[run][tic], str(workdir)))
    if todo:
        print(
            f"searching {len(todo)} star-sectors upside down ({workers} workers)",
            file=sys.stderr,
        )
        workdir.mkdir(parents=True, exist_ok=True)
        with ProcessPoolExecutor(max_workers=workers) as pool, cache.open("a") as fh:
            for i, r in enumerate(pool.map(_inverted_one, todo, chunksize=4), 1):
                fh.write(json.dumps(r) + "\n")
                done[(r["tic_id"], r["sector"])] = r
                if i % 200 == 0:
                    print(f"  {i}/{len(todo)}", file=sys.stderr)
    for f in found:
        if f.claimed_by or f.upside_down != "not run" or f.kind == "stitched":
            continue
        r = done.get((f.tic_id, int(f.sector)))
        if r is None or not r["ok"]:
            continue
        other = r["signals"]
        if any(alias_of(f.period_days, o["period_days"]) for o in other):
            f.upside_down = "same period"
        elif other:
            f.upside_down = "other period"
        else:
            f.upside_down = "nothing"
    return len(todo)


def collect(
    runs_dir: Path,
    samples: list[tuple[str, str, str | None]],
    planets: list[tuple[str, str]],
    stitched: list[tuple[str, str, str | None]],
    claims: dict[int, list[Claim]],
) -> tuple[list[Found], dict[str, dict[str, int]]]:
    """Every signal of every run, claimed or not, plus each run's star counts."""
    found: list[Found] = []
    counts: dict[str, dict[str, int]] = {}
    jobs = (
        [("sample", s, o, i) for s, o, i in samples]
        + [("planets", s, o, None) for s, o in planets]
        + [("stitched", s, o, i) for s, o, i in stitched]
    )
    for kind, sector, obs, inv in jobs:
        sig = signals_of(runs_dir / obs)
        inv_sig = signals_of(runs_dir / inv) if inv else None
        counts[obs] = {
            "stars": len(sig),
            "stars_with_signal": sum(1 for v in sig.values() if v),
        }
        if inv_sig is not None:
            counts[obs]["inverted_stars_with_signal"] = sum(
                1 for v in inv_sig.values() if v
            )
            # The measured false-alarm comparison is the stars WITHOUT a TCE:
            # upside down, a deep binary on a TCE star lights up at P/2, P/3.
            so = json.loads((runs_dir / obs / "summary.json").read_text())
            si = json.loads((runs_dir / inv / "summary.json").read_text())
            n_quiet = so.get("sample", {}).get("quiet_stars")
            if n_quiet is not None:
                counts[obs]["quiet_stars"] = n_quiet
                counts[obs]["quiet_with_signal"] = so["quiet_group"][
                    "stars_with_any_signal"
                ]
                counts[obs]["quiet_inverted_with_signal"] = si["quiet_group"][
                    "stars_with_any_signal"
                ]
        for tic, sigs in sig.items():
            for s in sigs:
                f = Found(
                    run=obs,
                    kind=kind,
                    sector=sector,
                    tic_id=tic,
                    period_days=s["period_days"],
                    epoch_btjd=s["epoch_btjd"],
                    duration_hours=s["duration_hours"],
                    depth_ppm=s["depth_ppm"],
                    snr=s["snr"],
                    sde=s["sde"],
                    n_transits=s["n_transits"],
                )
                tol = claim_tol(f)
                f.claimed_by = sorted(
                    {
                        c.source
                        for c in claims.get(tic, [])
                        if alias_of(f.period_days, c.period_days, CLAIM_ALIASES, tol)
                    }
                )
                if inv_sig is not None:
                    other = inv_sig.get(tic, [])
                    if any(alias_of(f.period_days, o["period_days"]) for o in other):
                        f.upside_down = "same period"
                    elif other:
                        f.upside_down = "other period"
                    else:
                        f.upside_down = "nothing"
                found.append(f)
    return found, counts


@dataclass
class Candidate:
    """Unclaimed signals of one star at one period, over every run that found it."""

    tic_id: int
    period_days: float
    finds: list[Found]
    known_on_star: list[str]
    pileup: bool = False

    @property
    def sectors(self) -> list[str]:
        return sorted(
            {f.sector for f in self.finds if f.kind != "stitched"}, key=lambda s: int(s)
        )

    @property
    def stitched(self) -> list[str]:
        return sorted({f.sector for f in self.finds if f.kind == "stitched"})

    @property
    def best(self) -> Found:
        return max(self.finds, key=lambda f: f.sde)

    @property
    def upside_down(self) -> str:
        seen = {f.upside_down for f in self.finds}
        for v in ("same period", "other period", "nothing"):
            if v in seen:
                return v
        return "not run"

    @property
    def repeated(self) -> bool:
        return len(self.sectors) >= 2 or bool(self.sectors and self.stitched)

    @property
    def tier(self) -> str:
        """A/B: the upside-down search found NOTHING on the star; C: it found
        another period (a variable star); D: the same period (variability);
        X: a depth of 100% or more, which is a broken light curve, not a dip.
        "not run" (a fetch failed) counts as C: untested is not clean."""
        if self.best.depth_ppm >= 1e6:
            return "X"
        if self.upside_down == "same period":
            return "D"
        if self.upside_down != "nothing":
            return "C"
        if self.pileup:
            return "S"
        return "A" if self.repeated else "B"


def group(found: list[Found], claims: dict[int, list[Claim]]) -> list[Candidate]:
    """Unclaimed finds, one candidate per star and period (aliases folded)."""
    by_star: dict[int, list[Candidate]] = defaultdict(list)
    for f in sorted((f for f in found if not f.claimed_by), key=lambda f: -f.sde):
        for c in by_star[f.tic_id]:
            if alias_of(f.period_days, c.period_days):
                c.finds.append(f)
                break
        else:
            by_star[f.tic_id].append(
                Candidate(
                    f.tic_id,
                    f.period_days,
                    [f],
                    sorted({c.source for c in claims.get(f.tic_id, [])}),
                )
            )
    return [c for cs in by_star.values() for c in cs]


PILEUP_BIN = 0.01  # in ln(period): bins 1% wide
PILEUP_MIN = 5
PILEUP_FACTOR = 4.0
PILEUP_NEIGHBOURS = 20  # bins either side


def mark_pileups(cands: list[Candidate]) -> list[dict[str, Any]]:
    """Set aside periods that many unrelated stars share.

    Stars do not choose a period, so a 1%-wide period bin holding far more
    clean candidates than its neighbours is an artefact (measured: 33 at
    0.99-1.01 d, every one with the same 3.0 h box). A bin is a pile-up when
    it holds at least PILEUP_MIN and PILEUP_FACTOR x the median of the
    PILEUP_NEIGHBOURS bins on each side; its A/B candidates become tier S.
    """
    clean = [c for c in cands if c.tier in ("A", "B")]
    bins: dict[int, list[Candidate]] = defaultdict(list)
    for c in clean:
        bins[math.floor(math.log(c.period_days) / PILEUP_BIN)].append(c)
    out = []
    for b, members in sorted(bins.items()):
        around = [
            len(bins.get(b + d, []))
            for d in range(-PILEUP_NEIGHBOURS, PILEUP_NEIGHBOURS + 1)
            if d
        ]
        base = max(1.0, sorted(around)[len(around) // 2])
        if len(members) >= PILEUP_MIN and len(members) >= PILEUP_FACTOR * base:
            for c in members:
                c.pileup = True
            durs = defaultdict(int)
            for c in members:
                durs[c.best.duration_hours] += 1
            out.append(
                {
                    "period_from": round(math.exp(b * PILEUP_BIN), 4),
                    "period_to": round(math.exp((b + 1) * PILEUP_BIN), 4),
                    "candidates": len(members),
                    "neighbour_median": base,
                    "durations_h": dict(sorted(durs.items())),
                }
            )
    return out


def _row(c: Candidate) -> dict[str, Any]:
    b = c.best
    return {
        "tier": c.tier,
        "tic_id": c.tic_id,
        "period_days": round(b.period_days, 5),
        "epoch_btjd": round(b.epoch_btjd, 4),
        "duration_hours": b.duration_hours,
        "depth_ppm": round(b.depth_ppm, 1),
        "snr": round(b.snr, 1),
        "sde": round(b.sde, 1),
        "n_transits": b.n_transits,
        "sectors": " ".join(c.sectors),
        "stitched": "; ".join(c.stitched),
        "upside_down": c.upside_down,
        "best_run": b.run,
        "known_on_star": "; ".join(c.known_on_star),
    }


def write_csv(cands: list[Candidate], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [_row(c) for c in sorted(cands, key=lambda c: (c.tier, -c.best.sde))]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=list(_row(cands[0]).keys()) if cands else ["tier"]
        )
        w.writeheader()
        w.writerows(rows)


STAR_FLAGS = ("(EB)", "(FP)", "(V)", "(FA)", "(IS)", "(O)")


def _star_note(known: str) -> str:
    """'binary/false positive on star' when a TOI on the star is EB/FP/V/FA."""
    return "EB/FP/V TOI on star" if any(f in known for f in STAR_FLAGS) else ""


def _fmt_row(r: dict[str, Any]) -> str:
    where = "; ".join(
        x
        for x in (r["sectors"], f"stitched {r['stitched']}" if r["stitched"] else "")
        if x
    )
    return (
        f"| {r['tic_id']} | {r['period_days']:.5f} | {r['duration_hours']:g} | "
        f"{r['depth_ppm']:.0f} | {r['snr']:.1f} | {r['sde']:.1f} | {r['n_transits']} | "
        f"{where} | {r['upside_down']} | {_star_note(r['known_on_star']) or '-'} |"
    )


HEADER = (
    "| TIC | period (d) | duration (h) | depth (ppm) | SNR | SDE | transits | "
    "sectors | upside down | star |\n|---|---|---|---|---|---|---|---|---|---|"
)


def write_md(data: dict[str, Any], path: Path, csv_name: str, max_b: int = 100) -> None:
    """The write-up: what was checked, the measured false-alarm context, the lists."""
    cands = data["candidates"]
    tiers = data["tiers"]
    runs = data["runs"]
    q = [v for v in runs.values() if "quiet_stars" in v]
    q_n = sum(v["quiet_stars"] for v in q)
    q_obs = sum(v["quiet_with_signal"] for v in q)
    q_inv = sum(v["quiet_inverted_with_signal"] for v in q)
    cat = data["catalogues"]
    tier_a = [c for c in cands if c["tier"] == "A"]
    tier_b = [c for c in cands if c["tier"] == "B"]
    non_toi = data["claimed_only_by_non_toi"]
    by_src: dict[str, int] = defaultdict(int)
    for f in non_toi:
        kinds = {s.split()[0] for s in f["claimed_by"]}
        by_src[" + ".join(sorted(kinds))] += 1
    lines = [
        "# Signals no catalogue claims",
        "",
        (
            "Generated by `python -m tesstransit.unclaimed` "
            f"({data['generated']}); every row is also in [`{csv_name}`]({csv_name})."
        ),
        "",
        "⚠️ **None of these is a planet, and none is claimed as one.** They are",
        "signals our search found on a star where no published catalogue lists",
        "anything at that period. As a group they are what noise and stellar",
        "variability produce. Measured on the stars WITHOUT a TCE in the sample",
        f"runs ({q_n:,} star-searches): **{q_obs:,}** gave a signal the right way up",
        f"and **{q_inv:,}** upside down, where no real transit can survive. The",
        "lists below sort the signals by the two tests that can set one apart",
        "from that background. Neither is vetting: no centroid test, no",
        "follow-up, no check that the star is the one dimming.",
        "",
        '## What "no catalogue claims it" means',
        "",
        "A signal is **claimed** when any of these lists the same star at the same",
        "period or an alias of it (x1, x2, x1/2, x3, x1/3, within 0.5%):",
        "",
        f"- every SPOC TCE: {len(data['tce_tables'])} DV tables, single-sector and",
        f"  multi-sector ({cat['tce']:,} TCEs);",
        f"- every TOI, false positives included ({cat['toi']:,});",
        f"- every community TOI, CTOI ({cat['ctoi']:,});",
        "- every confirmed planet in the NASA Exoplanet Archive with a TIC ID",
        f"  ({cat['confirmed']:,}).",
        "",
        f"Of **{data['signals']:,}** signals in {len(runs)} runs, **{data['unclaimed_signals']:,}**",
        "are unclaimed. One star at one period, over every run that found it, is",
        f"one candidate: **{len(cands):,}** candidates.",
        "",
        "## The two tests",
        "",
        "- **Upside down.** The star's light curve searched as 2 − flux with the",
        "  settings its own result recorded. The sample runs have this control",
        "  already; for every other unclaimed signal it was run for this report.",
        "  If the inverted search finds the same period, the bumps are as strong",
        "  as the dips: variability. If it finds another period, the star varies",
        "  enough to make signals of its own, so a dip on it is not clean either.",
        "- **Independent sectors.** Each sector is searched on its own data, so",
        "  the same period coming back in two sectors is not a coincidence of one",
        "  light curve. ⚠️ Variability repeats too, which is why this counts only",
        "  for stars that pass the first test.",
        "",
        "| tier | rule | candidates |",
        "|---|---|---|",
        f"| **A** | upside down finds nothing on the star; found in 2+ independent sectors (or a sector and a stitch) | {tiers['A']} |",
        f"| **B** | upside down finds nothing on the star; found once | {tiers['B']} |",
        f"| S | as A or B, but at a period many unrelated stars share (a pile-up, below) | {tiers['S']} |",
        f"| C | upside down finds another period: a variable star (or the inverted search could not run) | {tiers['C']} |",
        f"| D | upside down finds the same period: variability | {tiers['D']} |",
        f"| X | depth of 100% or more: a broken light curve, not a dip | {tiers['X']} |",
        "",
        '"star" flags a star with an eclipsing-binary, false-positive or variable',
        "TOI: a new period there is more likely the binary than a planet. Depths",
        "of several percent with two transits are most likely eclipsing binaries.",
        "",
        "## Tier A: repeated, nothing upside down",
        "",
        HEADER,
        *[_fmt_row(r) for r in tier_a],
        "",
        f"## Tier B: found once, nothing upside down (strongest {min(max_b, len(tier_b))} of {len(tier_b)})",
        "",
        "By SDE; the rest are in the CSV, with tiers C, D and X.",
        "",
        HEADER,
        *[_fmt_row(r) for r in tier_b[:max_b]],
        "",
        "## Period pile-ups (tier S)",
        "",
        "Stars do not choose a period, so a 1%-wide period bin holding far more",
        f"clean candidates than its neighbours (at least {PILEUP_MIN}, and {PILEUP_FACTOR:g}x the median",
        f"of the {PILEUP_NEIGHBOURS} bins either side) is an artefact of the data or the search,",
        "not a population. Their candidates are set aside, not deleted. Cause not",
        "yet investigated.",
        "",
        "| period (d) | candidates | neighbouring median | box durations (h: count) |",
        "|---|---|---|---|",
        *[
            f"| {p['period_from']:.4f}-{p['period_to']:.4f} | {p['candidates']} | {p['neighbour_median']:g} | "
            + ", ".join(f"{k:g}: {v}" for k, v in p["durations_h"].items())
            + " |"
            for p in data["pileups"]
        ],
        "",
        "## Found by others, but not a TOI",
        "",
        "Signals a catalogue does list, but only a SPOC TCE table or a community",
        "TOI: someone else's search found them, and nobody has made them a TOI or",
        "a planet. GJ 774 A (TIC 351339274, 4.225 d), from the Nature news story of",
        "2026-11-20, is this kind: a SPOC multi-sector TCE written up by an amateur,",
        "not a TOI. Our stitched search of its four sectors (13, 27, 67, 94) finds",
        "it at 4.225140 d (SDE 46.7, 23 transits, nothing upside down); that search",
        "was run by hand and is not one of the runs counted here.",
        "",
        "| claimed only by | signals |",
        "|---|---|",
        *[f"| {k} | {v:,} |" for k, v in sorted(by_src.items(), key=lambda kv: -kv[1])],
        "",
        "## Runs",
        "",
        "| run | stars | with a signal | no-TCE stars | of them with a signal | ... upside down |",
        "|---|---|---|---|---|---|",
        *[
            f"| {k} | {v['stars']} | {v['stars_with_signal']} | {v.get('quiet_stars', '-')} | "
            f"{v.get('quiet_with_signal', '-')} | {v.get('quiet_inverted_with_signal', '-')} |"
            for k, v in runs.items()
        ],
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument(
        "--catalogs", type=Path, required=True, help="dir with the TCE tables"
    )
    ap.add_argument("--toi", type=Path, required=True)
    ap.add_argument("--ctoi", type=Path, required=True)
    ap.add_argument("--confirmed", type=Path, required=True)
    ap.add_argument("--sample", action="append", default=[], metavar="SECTOR=OBS[:INV]")
    ap.add_argument("--planets", action="append", default=[], metavar="SECTOR=RUN")
    ap.add_argument(
        "--stitched", action="append", default=[], metavar="LABEL=OBS[:INV]"
    )
    ap.add_argument("--csv", type=Path, required=True)
    ap.add_argument("--json", type=Path, help="counts + candidates, for the write-up")
    ap.add_argument(
        "--inverted-cache",
        type=Path,
        help="search upside down every unclaimed find that has no inverted run; cache here",
    )
    ap.add_argument(
        "--workdir", type=Path, help="light-curve scratch for --inverted-cache"
    )
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument(
        "--md", type=Path, help="write the Markdown report here (needs --json)"
    )
    ap.add_argument("--md-max-b", type=int, default=100)
    a = ap.parse_args(argv)

    def pair(s: str) -> tuple[str, str, str | None]:
        k, v = s.split("=", 1)
        o, _, i = v.partition(":")
        return k, o, (i or None)

    tce_paths = sorted(a.catalogs.glob("s*_tcestats.csv"))
    claims: dict[int, list[Claim]] = defaultdict(list)
    sources = {
        "tce": read_tce_claims(tce_paths),
        "toi": read_toi_claims(a.toi),
        "ctoi": read_ctoi_claims(a.ctoi),
        "confirmed": read_planet_claims(a.confirmed),
    }
    for d in sources.values():
        for tic, cs in d.items():
            claims[tic].extend(cs)
    found, counts = collect(
        a.runs,
        [pair(s) for s in a.sample],
        [tuple(s.split("=", 1)) for s in a.planets],  # type: ignore[misc]
        [pair(s) for s in a.stitched],
        claims,
    )
    if a.inverted_cache:
        fill_inverted(
            found,
            a.runs,
            a.inverted_cache,
            a.workdir or a.inverted_cache.parent / "lc_unclaimed",
            a.workers,
        )
    cands = group(found, claims)
    write_csv(cands, a.csv)
    pileups = mark_pileups(cands)
    tiers = {t: sum(1 for c in cands if c.tier == t) for t in "ABSCDX"}
    print(
        f"{len(found)} signals, {sum(1 for f in found if not f.claimed_by)} unclaimed -> "
        f"{len(cands)} candidates ({', '.join(f'{k} {v}' for k, v in tiers.items())}); "
        f"{len(tce_paths)} TCE tables",
        file=sys.stderr,
    )
    data = {
        "generated": dt.datetime.now(dt.UTC).date().isoformat(),
        "tce_tables": [p.name for p in tce_paths],
        "catalogues": {k: sum(len(v) for v in d.values()) for k, d in sources.items()},
        "runs": counts,
        "signals": len(found),
        "unclaimed_signals": sum(1 for f in found if not f.claimed_by),
        "tiers": tiers,
        "pileups": pileups,
        "candidates": [
            _row(c) for c in sorted(cands, key=lambda c: (c.tier, -c.best.sde))
        ],
        "claimed_only_by_non_toi": [
            asdict(f)
            for f in found
            if f.claimed_by
            and not any(s.startswith(("TOI", "planet")) for s in f.claimed_by)
        ],
    }
    if a.json:
        a.json.write_text(json.dumps(data, indent=1))
    if a.md:
        write_md(
            data,
            a.md,
            # relative, never absolute: the report is committed and published,
            # and an absolute path would carry this machine's home directory
            os.path.relpath(a.csv, a.md.parent),
            a.md_max_b,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
