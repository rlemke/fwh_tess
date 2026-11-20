"""The FFL for the checks on finished runs, end to end in memory.

The upside-down control of a sector sample, injection-recovery, shape tests,
transit timing and the page's stitched sections each became a step over the
CLI that made them. What is under test is the plumbing: inputs localized, the
CLI's own ``main`` run, its files uploaded, and the counts read back. The
light curves are faked; the per-star work is the CLIs' own tests' business.
"""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from tesstransit import handlers as h
from tesstransit import inject, site, ttv, vet

from ._drive import run
from .test_workflow import HANDLERS, STARS, TCE_CSV, TOI_CSV, inputs  # noqa: F401


@pytest.fixture
def lc_layout(inputs, tmp_path):  # noqa: F811
    """TCE table and manifest where TCETable and Targets put them."""
    work = Path(inputs["work_dir"])
    (work / "lc").mkdir(parents=True)
    (work / "s0001_tcestats.csv").write_text(TCE_CSV)
    (work / "lc" / "manifest-s0001.json").write_text(
        (tmp_path / "targets.json").read_text()
    )
    return inputs


def _serial_pools(monkeypatch, *mods):
    # The CLIs fork a process pool; threads see the monkeypatched workers.
    for m in mods:
        monkeypatch.setattr(m, "ProcessPoolExecutor", ThreadPoolExecutor)
        monkeypatch.setattr(m.lc, "s3_client", lambda: None)


def test_the_upside_down_control_is_recorded_as_such(inputs, monkeypatch):  # noqa: F811
    seen = []
    real = h.search_jobs

    def spy(sector, jobs, cfg, workdir, workers):
        seen.append(cfg)
        return real(sector, jobs, cfg, workdir, workers)

    monkeypatch.setattr(h, "search_jobs", spy)
    out, _ = run(
        "tess.workflows.SectorRun",
        {**inputs, "kind": "sample", "label": "inv", "n_tce": 2, "n_quiet": 3,
         "invert": True, "set": ["sde_threshold=12"], "sample_label": "Upside down"},
        HANDLERS,
    )  # fmt: skip
    assert seen and all(c.invert_flux and c.sde_threshold == 12.0 for c in seen)
    summary = json.loads((Path(out["run_dir"]) / "summary.json").read_text())
    assert summary["search_config"]["invert_flux"] is True
    assert summary["search_config"]["sde_threshold"] == 12.0
    assert out["run_dir"].endswith("s1-sample-inv")
    assert out["site_arg"] == f"1=Upside down={out['run_dir']}"


def test_an_upside_down_planet_run_is_refused(inputs):  # noqa: F811
    with pytest.raises(ValueError, match="'sample' only"):
        h.handle_match_sector({
            "sector": 1, "kind": "planets", "run_dir": "unused", "results_paths": [],
            "targets_path": inputs["targets_path"], "tce_path": inputs["tce_path"],
            "toi_path": inputs["toi_path"], "invert": True,
        })  # fmt: skip


def test_completeness_injects_into_each_sample_run(lc_layout, monkeypatch):
    def fake(tic, sector, key, workdir, keep, cfg, per_star, seed, replay=None):
        return [
            {"tic_id": tic, "ok": True, "sigma_ppm": 500.0, "bucket": b, "alias": "",
             "signal": None, "n_signals": 0, "stop_reason": "x", "seconds": 0.1,
             "injection": {"index": i, "period_days": 2.0, "epoch_btjd": 1326.0,
                           "duration_hours": 2.0, "impact": 0.1, "radius_ratio": 0.05,
                           "depth_ppm": 900.0, "snr_expected": 10.0, "n_transits": 9}}
            for i, b in enumerate(["recovered", "missed"][:per_star])
        ]  # fmt: skip

    monkeypatch.setattr(inject, "_inject_star", fake)
    _serial_pools(monkeypatch, inject)
    sample, _ = run(
        "tess.workflows.SectorRun",
        {**lc_layout, "kind": "sample", "n_tce": 2, "n_quiet": 3},
        HANDLERS,
    )
    out, ran = run(
        "tess.workflows.Completeness",
        {"samples": [sample["site_arg"]], "work_dir": lc_layout["work_dir"],
         "per_star": 2, "workers": 1},
        HANDLERS,
    )  # fmt: skip
    assert ran == ["tess.check.Inject"]
    run_dir = Path(out["runs"][0])
    assert run_dir.name == "s1-sample-inject"
    summary = json.loads((run_dir / "summary.json").read_text())
    # The 3 no-TCE stars (42 and 77 carry TCEs), 2 injections each.
    assert summary["stars"] == 3 and summary["injections"] == 6
    assert out["injections"] == [6] and out["found"] == [3]
    assert (run_dir / "injections.csv").exists()


def test_vetting_runs_over_each_planet_run(lc_layout, monkeypatch):
    calls = []

    def fake(tic, sector, key, workdir, cfg, signals):
        calls.append(tic)
        return [{**s["row"], "odd_even_sigma": 0.5, "fails": "", "ok": True}
                for s in signals]  # fmt: skip

    monkeypatch.setattr(vet, "_vet_star", fake)
    _serial_pools(monkeypatch, vet)
    planets, _ = run("tess.workflows.SectorRun", lc_layout, HANDLERS)
    out, ran = run(
        "tess.workflows.VetPlanets",
        {
            "planets": [planets["site_arg"]],
            "work_dir": lc_layout["work_dir"],
            "workers": 1,
        },
        HANDLERS,
    )
    # TOI 101.01 (CP) is the one recovered planet; 102.01 (PC) is not vetted.
    assert ran == ["tess.check.Vet"] and calls == [42]
    assert out == {"csv_path": str(Path(lc_layout["work_dir"]) / "vet" / "vet.csv"),
                   "signals": 1, "flagged": 0, "failed": 0}  # fmt: skip


def test_timing_times_the_named_tois(lc_layout, monkeypatch, tmp_path):
    def fake(tic, keys, workdir, cfg, targets, bin_minutes):
        return [{**x, "tic_id": tic, "ok": True, "n_transits": 7, "chi2_dof": 1.2}
                for x in targets]  # fmt: skip

    monkeypatch.setattr(ttv, "_ttv_star", fake)
    _serial_pools(monkeypatch, ttv)
    out = h.handle_timing({
        "sectors": [1], "lc_dir": str(Path(lc_layout["work_dir"]) / "lc"),
        "out_dir": str(tmp_path / "ttv"), "tois": ["101.01"],
        "toi_path": lc_layout["toi_path"], "workers": 1,
    })  # fmt: skip
    rows = [json.loads(ln) for ln in (tmp_path / "ttv" / "ttv.jsonl").open()]
    assert [r["name"] for r in rows] == ["101.01"]
    assert out["signals"] == 1 and out["timed"] == 1 and out["failed"] == 0
    with pytest.raises(ValueError, match="toi_path"):
        h.handle_timing(
            {"sectors": [1], "lc_dir": "x", "out_dir": "y", "tois": ["1.01"]}
        )
    with pytest.raises(FileNotFoundError, match="sector 2"):
        h.handle_timing({"sectors": [1, 2], "lc_dir": str(Path(lc_layout["work_dir"]) / "lc"),
                         "out_dir": "y", "from_run": "r"})  # fmt: skip


def test_the_page_takes_both_stitched_sections(lc_layout, monkeypatch, tmp_path):
    seen = {}

    def fake_main(argv):
        seen["argv"] = argv
        out = Path(argv[argv.index("--out") + 1])
        out.mkdir(parents=True)
        (out / "index.html").write_text("<html></html>")
        return 0

    monkeypatch.setattr(site, "main", fake_main)
    a, b = tmp_path / "ms", tmp_path / "ms-inv"
    out, _ = run(
        "tess.workflows.BuildPage",
        {"planets": ["1=" + str(tmp_path / "p")], "samples": [],
         "work_dir": lc_layout["work_dir"], "out_dir": str(tmp_path / "site"),
         "stitched": [f"A|SDE 13|{a}|{b}|yes"],
         "long_stitched": [f"South A|SDE 18|{a}|{b}|yes"]},
        HANDLERS,
    )  # fmt: skip
    argv = seen["argv"]
    assert argv[argv.index("--stitched") + 1] == f"A|SDE 13|{a}|{b}|yes"
    assert argv[argv.index("--long-stitched") + 1] == f"South A|SDE 18|{a}|{b}|yes"
    assert Path(out["page_path"]).read_text() == "<html></html>"
    with pytest.raises(ValueError, match="label\\|settings"):
        h.handle_build_site({"planets": [], "samples": [], "lc_dir": str(tmp_path),
                             "out_dir": str(tmp_path / "x"), "stitched": ["A|b|c"]})  # fmt: skip


def test_every_tess_workflow_compiles_with_its_libraries():
    """All workflows, including the site/publish ones, are reachable by name."""
    from ._drive import program

    names = {
        f"{ns['name']}.{d['name']}"
        for ns in program()["declarations"]
        if ns.get("type") == "Namespace"
        for d in ns["declarations"]
        if d.get("type") == "WorkflowDecl"
    }
    assert {
        "tess.workflows.SectorRun", "tess.workflows.SearchSectors",
        "tess.workflows.BuildResults", "tess.workflows.StitchedSample",
        "tess.workflows.FalseAlarms", "tess.workflows.Completeness",
        "tess.workflows.VetPlanets", "tess.workflows.TimeTOIs",
        "tess.workflows.BuildPage",
    } <= names  # fmt: skip
    assert STARS  # the shared fixture's stars
    assert "Full TOI ID" in TOI_CSV
