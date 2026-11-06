"""Facetwork handlers for the ``tess`` facets.

Each handler is a thin step over the same code the command-line tools run, so
a workflow run and a CLI run over the same inputs produce the same files. The
split is the one the CLIs already had internally:

1. ``tess.source.TCETable`` / ``tess.source.Targets`` -- a sector's SPOC TCE
   table and the stars it observed at 2-minute cadence.
2. ``tess.Plan`` -- choose the stars (every TOI star, or a TCE/no-TCE sample)
   and cut them into chunks.
3. ``tess.SearchChunk`` -- search one chunk; the fan-out runs these.
4. ``tess.MatchSector`` -- gather the chunks into ``results.jsonl`` and run the
   CLI's matcher over it (``--reuse``), which writes the run directory.
5. ``tess.BuildSite`` -- the results page over finished run directories.

⚠️ Paths may be local or ``s3://``. Work happens on local scratch; inputs are
localized and outputs uploaded, so no step depends on a disk another host
cannot see. A chunk's results are rewritten whole by every attempt, so a
re-run (at-least-once delivery) replaces rather than duplicates them.
"""

from __future__ import annotations

import contextlib
import io
import json
import os
import re
import shutil
import tempfile
import time
import urllib.request
from collections.abc import Callable, Iterator
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from . import _lightcurve as lc
from . import _manifest as mf
from . import _match as mt
from . import _search as sr

TCE_INDEX = "https://archive.stsci.edu/missions/tess/catalogs/tce/"
RUN_FILES = ("results.jsonl", "matches.csv", "tois.csv", "summary.json")


def _say(params: dict[str, Any]) -> Callable[..., None]:
    cb = params.get("_step_log")
    return cb if callable(cb) else (lambda *a, **k: None)


def _store() -> Any:
    from facetwork.domains.storage import domain_storage

    return domain_storage("tess")


def _remote(path: str) -> bool:
    from facetwork.domains.storage import is_remote

    return is_remote(path)


def _join(base: str, *parts: str) -> str:
    return "/".join([base.rstrip("/"), *parts])


def _read(path: str) -> str:
    return _store().read_text(path)


def _write(path: str, body: str) -> None:
    _store().write_text(path, body)


def _exists(path: str) -> bool:
    return _store().exists(path)


@contextlib.contextmanager
def _scratch() -> Iterator[Path]:
    root = os.environ.get("FW_LOCAL_SCRATCH") or None
    d = Path(tempfile.mkdtemp(prefix="fw-tess-", dir=root))
    try:
        yield d
    finally:
        shutil.rmtree(d, ignore_errors=True)


def _pull(path: str, into: Path) -> Path:
    """A local copy of one file: the path itself when local."""
    if not _remote(path):
        return Path(path)
    local = into / Path(path).name
    local.write_bytes(_store().read_bytes(path))
    return local


def _pull_run(run_dir: str, into: Path) -> Path:
    """A local copy of a run directory (just the files the next step reads)."""
    if not _remote(run_dir):
        return Path(run_dir)
    local = into / Path(run_dir.rstrip("/")).name
    local.mkdir(parents=True, exist_ok=True)
    for name in RUN_FILES:
        src = _join(run_dir, name)
        if _exists(src):
            (local / name).write_bytes(_store().read_bytes(src))
    return local


def _push_dir(local: Path, dest: str) -> None:
    if not _remote(dest) and Path(dest).resolve() == local.resolve():
        return
    for f in local.iterdir():
        if f.is_file():
            _store().write_bytes(_join(dest, f.name), f.read_bytes())


def _cfg(prewhiten: bool) -> sr.SearchConfig:
    # Same rule as the CLIs: pre-whitening is on unless switched off.
    cfg = sr.SearchConfig(prewhiten_max_period_days=1.5 if prewhiten else 0.0)
    cfg.validate()
    return cfg


def _run_dir(work_dir: str, sector: int, kind: str, label: str) -> str:
    return _join(work_dir, f"s{sector}-{kind}" + (f"-{label}" if label else ""))


def _manifest_name(sector: int) -> str:
    return f"manifest-s{sector:04d}.json"


def _fetch(url: str, timeout: int = 300) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "facetwork-tess/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


# ---------------------------------------------------------------------------
# sources
# ---------------------------------------------------------------------------


def tce_url(sector: int, index_html: str) -> str:
    """The single-sector DV TCE table for ``sector`` from MAST's directory page.

    ⚠️ The file name carries a processing stamp (``tess<stamp>-s0001-s0001``)
    that cannot be derived, so it is read from the listing. If SPOC ever
    reprocesses a sector there are two; the newest stamp wins, and which one
    was used is returned with the table.
    """
    pat = re.compile(
        rf'href="(tess(\d+)-s{sector:04d}-s{sector:04d}_dvr-tcestats\.csv)"'
    )
    found = sorted(pat.findall(index_html), key=lambda m: m[1])
    if not found:
        raise ValueError(
            f"no single-sector TCE table for sector {sector} at {TCE_INDEX}"
        )
    return TCE_INDEX + found[-1][0]


def handle_tce_table(params: dict[str, Any]) -> dict[str, Any]:
    say = _say(params)
    sector = int(params["sector"])
    dest = _join(str(params["dest_dir"]), f"s{sector:04d}_tcestats.csv")
    url = tce_url(sector, _fetch(TCE_INDEX).decode())
    marker = dest + ".url"
    # A stamped TCE table never changes; a new stamp is a new URL.
    if _exists(dest) and _exists(marker) and _read(marker).strip() == url:
        say(f"sector {sector}: TCE table current ({url})")
    else:
        say(f"sector {sector}: fetching {url}")
        _write(dest, _fetch(url).decode())
        _write(marker, url + "\n")
    n = len(mt.parse_tce_csv(_read(dest)))
    return {"path": dest, "url": url, "tces": n}


def handle_targets(params: dict[str, Any]) -> dict[str, Any]:
    say = _say(params)
    sector = int(params["sector"])
    dest = _join(str(params["dest_dir"]), _manifest_name(sector))
    if _exists(dest):
        targets = mf.parse(_read(dest), sector, dest)
    else:
        t0 = time.time()
        targets = mf.list_sector(lc.s3_client(), sector)
        _write(dest, mf.dumps(sector, targets))
        say(
            f"sector {sector}: listed {len(targets):,} light curves ({time.time() - t0:.0f}s)"
        )
    return {"path": dest, "stars": len(targets)}


# ---------------------------------------------------------------------------
# plan -> search -> match
# ---------------------------------------------------------------------------


def handle_plan(params: dict[str, Any]) -> dict[str, Any]:
    from . import pilot, planets

    say = _say(params)
    sector = int(params["sector"])
    kind = str(params["kind"])
    size = int(params.get("chunk_size") or 50)
    if size < 1:
        raise ValueError("chunk_size must be at least 1")
    run_dir = _run_dir(
        str(params["work_dir"]), sector, kind, str(params.get("label") or "")
    )
    targets = mf.parse(
        _read(str(params["targets_path"])), sector, str(params["targets_path"])
    )
    tces = mt.parse_tce_csv(_read(str(params["tce_path"])))

    jobs: list[dict[str, Any]]
    if kind == "planets":
        sel = planets.select(sector, targets, tces, _read(str(params["toi_path"])))
        jobs = [
            {
                "tic_id": tic,
                "key": targets[tic],
                "forced": planets.forced_ephemerides(ts),
            }
            for tic, ts in sorted(sel.by_star.items())
        ]
        say(
            f"sector {sector}: {len(sel.all_tois)} TOIs list it, {len(jobs)} stars observed"
        )
    elif kind == "sample":
        with _scratch() as tmp:
            held = [
                _pull_run(str(r), tmp / f"h{i}")
                for i, r in enumerate(params.get("holdout_of") or [])
            ]
            exclude = pilot.stars_of(held)
        tce_stars, quiet, unobserved = pilot.select_sample(
            tces, targets, int(params.get("n_tce") or 60), int(params.get("n_quiet") or 40),
            int(params.get("seed") or 1), exclude,
        )  # fmt: skip
        jobs = [
            {"tic_id": t, "key": targets[t], "forced": []} for t in tce_stars + quiet
        ]
        say(f"sector {sector}: {len(tce_stars)} TCE stars + {len(quiet)} no-TCE; "
            f"{len(unobserved)} TCE stars without a 2-min light curve")  # fmt: skip
    else:
        raise ValueError(f"kind must be 'planets' or 'sample', got {kind!r}")

    chunks = []
    for i in range(0, len(jobs), size):
        path = _join(run_dir, "chunks", f"c{i // size:04d}.json")
        _write(path, json.dumps({"sector": sector, "jobs": jobs[i : i + size]}))
        chunks.append(path)
    return {"chunks": chunks, "stars": len(jobs), "run_dir": run_dir}


def search_jobs(
    sector: int,
    jobs: list[dict[str, Any]],
    cfg: sr.SearchConfig,
    workdir: Path,
    workers: int,
) -> list[dict[str, Any]]:
    """Search each job's star (forcing its ephemerides), as the CLIs do."""
    from .pilot import _search_one

    def args(j: dict[str, Any]) -> tuple[Any, ...]:
        forced = [tuple(e) for e in j["forced"]] or None
        return (j["tic_id"], sector, j["key"], str(workdir), False, cfg, forced)

    if workers <= 1:
        return [_search_one(*args(j)) for j in jobs]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(_search_one, *args(j)) for j in jobs]
        return [f.result() for f in as_completed(futs)]


def handle_search_chunk(params: dict[str, Any]) -> dict[str, Any]:
    say = _say(params)
    chunk_path = str(params["chunk_path"])
    chunk = json.loads(_read(chunk_path))
    cfg = _cfg(bool(params.get("prewhiten", True)))
    t0 = time.time()
    with _scratch() as tmp:
        records = search_jobs(
            int(chunk["sector"]),
            chunk["jobs"],
            cfg,
            tmp,
            int(params.get("workers") or 4),
        )
    records.sort(key=lambda r: r["tic_id"])
    dest = chunk_path.replace("/chunks/", "/parts/").removesuffix(".json") + ".jsonl"
    _write(dest, "".join(json.dumps(r) + "\n" for r in records))
    failed = [r for r in records if not r["ok"]]
    for r in failed:
        # ⚠️ Not a debug line: a failed star is missing from every result built on it.
        say(f"TIC {r['tic_id']} failed: {r['error']}", level="warning")
    say(f"{len(records)} stars in {time.time() - t0:.0f}s, {len(failed)} failed")
    return {"results_path": dest, "searched": len(records), "failed": len(failed)}


def handle_match_sector(params: dict[str, Any]) -> dict[str, Any]:
    from . import pilot, planets

    say = _say(params)
    sector = int(params["sector"])
    kind = str(params["kind"])
    run_dir = str(params["run_dir"])
    parts = [str(p) for p in params["results_paths"]]
    records = [json.loads(ln) for p in parts for ln in _read(p).splitlines() if ln]
    records.sort(key=lambda r: r["tic_id"])
    if len({r["tic_id"] for r in records}) != len(records):
        raise ValueError(f"sector {sector} {kind}: a star appears in two chunks")

    with _scratch() as tmp:
        out = tmp / "run"
        out.mkdir()
        (out / "results.jsonl").write_text(
            "".join(json.dumps(r) + "\n" for r in records)
        )
        lcdir = tmp / "lc"
        lcdir.mkdir()
        shutil.copy(
            _pull(str(params["targets_path"]), tmp), lcdir / _manifest_name(sector)
        )
        argv = [
            "--sector", str(sector),
            "--tce-csv", str(_pull(str(params["tce_path"]), tmp)),
            "--toi-csv", str(_pull(str(params["toi_path"]), tmp)),
            "--workdir", str(lcdir), "--out", str(out), "--reuse",
        ]  # fmt: skip
        if not params.get("prewhiten", True):
            argv.append("--no-prewhiten")
        if kind == "planets":
            main = planets.main
        elif kind == "sample":
            main = pilot.main
            argv += ["--n-tce", str(int(params.get("n_tce") or 60)),
                     "--n-quiet", str(int(params.get("n_quiet") or 40)),
                     "--seed", str(int(params.get("seed") or 1))]  # fmt: skip
            for i, r in enumerate(params.get("holdout_of") or []):
                argv += ["--holdout-of", str(_pull_run(str(r), tmp / f"h{i}"))]
        else:
            raise ValueError(f"kind must be 'planets' or 'sample', got {kind!r}")
        with contextlib.redirect_stdout(io.StringIO()):
            rc = main(argv)
        if rc:
            raise RuntimeError(f"sector {sector} {kind}: matcher exited {rc}")
        summary = json.loads((out / "summary.json").read_text())
        _push_dir(out, run_dir)

    failed = summary.get("failed") or summary.get("sample", {}).get("failed") or {}
    if kind == "planets":
        g = summary["outcomes"]["planet_like"]["tce_backed"]
        found = sum(v for k, v in g.items() if k.startswith("recovered"))
        say(
            f"sector {sector}: {found} of {sum(g.values())} planets SPOC detected found again"
        )
        site_arg = f"{sector}={run_dir}"
    else:
        g = summary["tce_group"]
        found = g.get("recovered", 0) + g.get("recovered_alias", 0)
        say(f"sector {sector}: {found} sample TCEs found again")
        site_arg = f"{sector}={params.get('sample_label') or 'Sample'}={run_dir}"
    if failed:
        say(f"sector {sector} {kind}: {len(failed)} star(s) failed", level="warning")
    return {
        "run_dir": run_dir,
        "summary_path": _join(run_dir, "summary.json"),
        "stars": len(records),
        "failed": len(failed),
        "found": found,
        "site_arg": site_arg,
    }


# ---------------------------------------------------------------------------
# results page
# ---------------------------------------------------------------------------


def handle_build_site(params: dict[str, Any]) -> dict[str, Any]:
    from . import site

    say = _say(params)
    out_dir = str(params["out_dir"])
    with _scratch() as tmp:

        def local(arg: str) -> str:
            *head, run = arg.split("=")
            return "=".join(
                [*head, str(_pull_run(run, tmp / f"r{len(os.listdir(tmp))}"))]
            )

        argv = [
            "--workdir",
            str(_pull_lc(str(params["lc_dir"]), tmp)),
            "--out",
            str(tmp / "site"),
        ]
        for a in params["planets"]:
            argv += ["--planets", local(str(a))]
        for a in params["samples"]:
            argv += ["--sample", local(str(a))]
        for v in params.get("variants") or []:
            argv += ["--variant", str(v)]
        if params.get("default"):
            argv += ["--default", str(params["default"])]
        rc = site.main(argv)
        if rc:
            raise RuntimeError(f"site build exited {rc}")
        page = (tmp / "site" / "index.html").read_text()
    dest = _join(out_dir, "index.html")
    _write(dest, page)
    say(
        f"results page: {len(params['planets'])} planet runs, {len(params['samples'])} samples"
    )
    return {"output_dir": out_dir, "page_path": dest, "bytes": len(page.encode())}


def _pull_lc(lc_dir: str, tmp: Path) -> Path:
    """Manifests and sky-position caches; local when the lc dir is local."""
    if not _remote(lc_dir):
        Path(lc_dir).mkdir(parents=True, exist_ok=True)
        return Path(lc_dir)
    local = tmp / "lc"
    local.mkdir()
    for name in _store().list_dir(lc_dir):
        base = Path(name).name
        if base.startswith(("manifest-", "coords-")):
            (local / base).write_bytes(_store().read_bytes(_join(lc_dir, base)))
    return local


# ---------------------------------------------------------------------------
# registration
# ---------------------------------------------------------------------------

_DISPATCH: dict[str, Callable[[dict[str, Any]], dict[str, Any]]] = {
    "tess.source.TCETable": handle_tce_table,
    "tess.source.Targets": handle_targets,
    "tess.Plan": handle_plan,
    "tess.SearchChunk": handle_search_chunk,
    "tess.MatchSector": handle_match_sector,
    "tess.BuildSite": handle_build_site,
}


def handle(payload: dict[str, Any]) -> dict[str, Any]:
    """RegistryRunner entrypoint: route on the facet name the runner supplies."""
    facet = payload.get("_facet_name", "")
    fn = _DISPATCH.get(facet)
    if fn is None:
        raise ValueError(f"no tess handler for facet {facet!r}")
    return fn(payload)


def register_handlers(runner: Any) -> None:
    """RegistryRunner registration -- one entrypoint per facet."""
    for facet_name in _DISPATCH:
        runner.register_handler(
            facet_name=facet_name,
            module_uri=f"file://{os.path.abspath(__file__)}",
            entrypoint="handle",
        )


def register_poller(poller: Any) -> None:
    for facet_name, fn in _DISPATCH.items():
        poller.register(facet_name, fn)
