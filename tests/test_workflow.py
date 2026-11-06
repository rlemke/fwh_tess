"""The FFL and its handlers, run end to end in memory with the search faked.

What is under test is the plumbing: the plan is chunked, every chunk is
searched exactly once, the chunks are gathered, and the CLI's own matcher
writes the run directory. The search itself is ``test_tess.py``'s business,
so here it is replaced by a stub that "finds" one planet.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

import pytest

from tesstransit import _search as sr
from tesstransit import handlers as h
from tesstransit import pilot

from ._drive import program, run

TCE_CSV = """# TESS DV TCE statistics
tceid,ticid,tce_plnt_num,tce_period,tce_time0bt,tce_duration,tce_depth,tce_model_snr,tce_num_transits
00000000042-01,42,1,3.7,1326.1,3.0,1500,40,7
00000000077-01,77,1,2.0,1325.9,2.0,800,12,13
"""
TOI_CSV = (
    "TIC,Full TOI ID,TOI Disposition,Orbital Epoch Value,Orbital Period (days) Value,"
    "Transit Duration (hours) Value,Transit Depth Value,Sectors,Detection Pipeline(s)\n"
    '42,101.01,CP,1326.1,3.7,3.0,1500,"[1, 27]",SPOC\n'
    '77,102.01,PC,1325.9,2.0,2.0,800,"[1]",SPOC\n'
    '99,103.01,PC,1325.0,5.0,2.0,800,"[2]",SPOC\n'
)
STARS = [42, 77, 5, 6, 7, 8, 9]


def _fake_search(tic, sector, key, workdir, keep, cfg, forced=None):
    signals = (
        [sr.Signal(1, 3.7, 1326.1, 3.0, 1500.0, 30.0, 20.0, 7)] if tic == 42 else []
    )
    res = sr.SearchResult(
        tic_id=tic, sector=sector, config=cfg.to_dict(), n_points=1000,
        scatter_ppm=500.0, signals=signals, stop_reason="below_threshold",
    )  # fmt: skip
    rec = {
        "tic_id": tic, "ok": True, "result": asdict(res),
        "mask": {"n_cadences": 1000, "n_flagged": 0, "n_nonfinite": 0},
        "t_span": [1325.3, 1352.9], "seconds": {"fetch": 0.1, "total": 1.0},
    }  # fmt: skip
    if forced:
        rec["forced"] = {n: {"n_transits": 7, "snr_red": 30.0, "depth_ppm": 1500.0}
                         for n, *_ in forced}  # fmt: skip
    return rec


@pytest.fixture
def inputs(tmp_path, monkeypatch):
    monkeypatch.setattr(pilot, "_search_one", _fake_search)
    (tmp_path / "tce.csv").write_text(TCE_CSV)
    (tmp_path / "tois.csv").write_text(TOI_CSV)
    (tmp_path / "targets.json").write_text(
        json.dumps({"sector": 1, "targets": {str(t): f"key-{t}" for t in STARS}})
    )
    return {
        "sector": 1, "work_dir": str(tmp_path / "work"),
        "targets_path": str(tmp_path / "targets.json"),
        "tce_path": str(tmp_path / "tce.csv"), "toi_path": str(tmp_path / "tois.csv"),
        "chunk_size": 1, "workers": 1,
    }  # fmt: skip


HANDLERS = {name: fn for name, fn in h._DISPATCH.items()}


def test_ffl_declares_exactly_the_handled_facets():
    events = {
        f"{ns['name']}.{d['name']}"
        for ns in program()["declarations"]
        if ns.get("type") == "Namespace" and ns["name"].startswith("tess")
        for d in ns["declarations"]
        if d.get("type") == "EventFacetDecl"
    }
    assert events == set(h._DISPATCH)


def test_planet_run_chunks_searches_once_and_writes_the_run(inputs):
    out, ran = run("tess.workflows.SectorRun", inputs, HANDLERS)
    # TOIs 101 and 102 list sector 1 on observed stars; 103 is sector 2 only.
    assert ran.count("tess.SearchChunk") == 2
    run_dir = Path(out["run_dir"])
    rows = list(csv.DictReader((run_dir / "tois.csv").open()))
    assert {r["toi"]: r["outcome"] for r in rows} == {
        "101.01": "recovered", "102.01": "missed:not_detected",
    }  # fmt: skip
    assert out["stars"] == 2 and out["failed"] == 0 and out["found"] == 1
    assert out["site_arg"] == f"1={run_dir}"
    assert (run_dir / "results.jsonl").read_text().count("\n") == 2


def test_sample_run_matches_what_the_cli_would(inputs, tmp_path):
    out, ran = run(
        "tess.workflows.SectorRun",
        {**inputs, "kind": "sample", "n_tce": 2, "n_quiet": 3, "seed": 1},
        HANDLERS,
    )
    assert ran.count("tess.SearchChunk") == 5
    summary = json.loads((Path(out["run_dir"]) / "summary.json").read_text())
    assert summary["sample"] == {"tce_stars": 2, "quiet_stars": 3, "failed": {}}
    assert summary["tce_group"]["recovered"] == 1
    assert out["site_arg"] == f"1=Sample={out['run_dir']}"

    # The same sample, searched by the CLI in one go, gives the same matches.
    lcdir = tmp_path / "lc"
    lcdir.mkdir()
    (lcdir / "manifest-s0001.json").write_text((tmp_path / "targets.json").read_text())
    cli = tmp_path / "cli"
    rc = pilot.main([
        "--sector", "1", "--tce-csv", inputs["tce_path"], "--toi-csv", inputs["toi_path"],
        "--workdir", str(lcdir), "--out", str(cli), "--n-tce", "2", "--n-quiet", "3",
        "--workers", "1",
    ])  # fmt: skip
    assert rc == 0
    assert (cli / "matches.csv").read_text() == (
        Path(out["run_dir"]) / "matches.csv"
    ).read_text()


def test_a_star_in_two_chunks_is_refused(inputs, tmp_path):
    part = tmp_path / "p.jsonl"
    part.write_text(
        json.dumps(_fake_search(42, 1, "k", "", False, sr.SearchConfig())) + "\n"
    )
    with pytest.raises(ValueError, match="two chunks"):
        h.handle_match_sector({
            "sector": 1, "kind": "planets", "run_dir": str(tmp_path / "r"),
            "results_paths": [str(part), str(part)], "targets_path": inputs["targets_path"],
            "tce_path": inputs["tce_path"], "toi_path": inputs["toi_path"],
        })  # fmt: skip


def test_tce_url_takes_the_newest_stamp_for_the_sector():
    page = (
        '<a href="tess2018206190142-s0001-s0001_dvr-tcestats.csv">'
        '<a href="tess2020101010101-s0001-s0001_dvr-tcestats.csv">'
        '<a href="tess2018234235059-s0002-s0002_dvr-tcestats.csv">'
        '<a href="tess2019170095531-s0001-s0013_dvr-tcestats.csv">'
    )
    assert h.tce_url(1, page).endswith("tess2020101010101-s0001-s0001_dvr-tcestats.csv")
    with pytest.raises(ValueError, match="sector 3"):
        h.tce_url(3, page)


def test_unknown_facet_is_an_error():
    with pytest.raises(ValueError, match="no tess handler"):
        h.handle({"_facet_name": "tess.Nope"})


def test_every_ffl_block_in_the_docs_compiles(tmp_path):
    """The README snippet and the gallery are compiled, so they cannot drift."""
    import re
    import subprocess
    import sys

    import facetwork

    root = Path(__file__).resolve().parents[1]
    lib = Path(facetwork.__file__).parent / "ffl"
    text = (root / "README.md").read_text() + (
        root / "docs" / "ffl-examples.md"
    ).read_text()
    blocks = re.findall(r"```ffl\n(.*?)```", text, re.DOTALL)
    assert len(blocks) >= 7
    for i, block in enumerate(blocks):
        src = tmp_path / f"b{i}.ffl"
        src.write_text(block)
        argv = [
            sys.executable,
            "-m",
            "facetwork.cli",
            "compile",
            "--check",
            "--primary",
            str(src),
        ]
        for library in (root / "src/tesstransit/ffl/tess.ffl", lib / "http.ffl",
                        lib / "compare.ffl", lib / "provenance.ffl"):  # fmt: skip
            argv += ["--library", str(library)]
        done = subprocess.run(argv, capture_output=True, text=True, check=False)
        assert done.returncode == 0, f"block {i}:\n{block}\n{done.stderr}{done.stdout}"
