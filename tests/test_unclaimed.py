"""The rules that decide what counts as found here and nowhere else."""

from __future__ import annotations

import json

from tesstransit import unclaimed as u


def _found(**kw):
    base = {
        "run": "s5-pilot", "kind": "sample", "sector": "5", "tic_id": 1, "period_days": 3.0,
        "epoch_btjd": 1.0, "duration_hours": 2.0, "depth_ppm": 500.0, "snr": 10.0,
        "sde": 20.0, "n_transits": 8,
    }  # fmt: skip
    return u.Found(**{**base, **kw})


def test_harmonics_beyond_three_count_as_claimed():
    """SPOC at 0.27864 d, ours at 1.11456 d: x4, outside the matcher's aliases."""
    assert u.alias_of(1.11456, 0.27864) is None
    assert u.alias_of(1.11456, 0.27864, u.CLAIM_ALIASES, 0.005) == 4.0


def test_claim_tolerance_follows_the_transits():
    """Two 8 h transits 10.6 d apart fix the period to ~3%, not 0.5%."""
    f = _found(period_days=10.586, duration_hours=8.0, n_transits=2)
    assert u.claim_tol(f) > 0.03
    assert u.alias_of(10.586, 10.6912, u.CLAIM_ALIASES, u.claim_tol(f)) == 1.0
    many = _found(period_days=1.0, duration_hours=1.0, n_transits=25)
    assert u.claim_tol(many) == u.REL_TOL


def _cand(**kw):
    f = _found(**kw)
    return u.Candidate(f.tic_id, f.period_days, [f], [])


def test_tiers():
    assert _cand(upside_down="nothing").tier == "B"
    assert _cand(upside_down="other period").tier == "C"
    assert _cand(upside_down="same period").tier == "D"
    assert _cand(upside_down="not run").tier == "C"  # untested is not clean
    assert _cand(upside_down="nothing", depth_ppm=2e6).tier == "X"
    c = _cand(upside_down="nothing")
    c.finds.append(_found(sector="9", upside_down="nothing"))
    assert c.tier == "A"
    # the same sector twice (a sample run and a planet run) is not a repeat
    d = _cand(upside_down="nothing")
    d.finds.append(_found(run="s5-planets", kind="planets", upside_down="nothing"))
    assert d.tier == "B"


def test_a_shared_period_is_a_pileup_not_a_population():
    cands = [
        _cand(tic_id=i, period_days=1.0 + i * 1e-5, upside_down="nothing")
        for i in range(8)
    ]
    spread = [
        _cand(tic_id=100 + i, period_days=0.5 * 1.02**i, upside_down="nothing")
        for i in range(60)
    ]
    found = u.mark_pileups(cands + spread)
    assert len(found) == 1 and found[0]["candidates"] >= 8
    assert all(c.tier == "S" for c in cands)
    assert sum(c.tier == "S" for c in spread) <= 1


def test_end_to_end(tmp_path):
    runs = tmp_path / "runs"
    for name, sigs in {
        "s5-pilot": [3.0, 7.0],
        "s5-pilot-inv": [],
        "s9-pilot": [3.0],
        "s9-pilot-inv": [3.0],
    }.items():
        d = runs / name
        d.mkdir(parents=True)
        res = {"tic_id": 1, "ok": True, "result": {"config": {}, "signals": [
            {"period_days": p, "epoch_btjd": 1.0, "duration_hours": 2.0, "depth_ppm": 500.0,
             "snr": 10.0, "sde": 20.0, "n_transits": 8} for p in sigs]}}  # fmt: skip
        (d / "results.jsonl").write_text(json.dumps(res) + "\n")
        (d / "summary.json").write_text(json.dumps(
            {"sample": {"quiet_stars": 1}, "quiet_group": {"stars_with_any_signal": int(bool(sigs))}}
        ))  # fmt: skip
    claims = {1: [u.Claim("TOI 1.01 (PC)", 14.0)]}  # 7.0 = 14.0 / 2: claimed
    found, counts = u.collect(
        runs,
        [("5", "s5-pilot", "s5-pilot-inv"), ("9", "s9-pilot", "s9-pilot-inv")],
        [],
        [],
        claims,
    )
    assert [f.period_days for f in found if f.claimed_by] == [7.0]
    cands = u.group(found, claims)
    assert len(cands) == 1  # 3.0 d in sectors 5 and 9 is one candidate
    # sector 9's upside-down run finds 3.0 d too, so the candidate is variability
    assert cands[0].tier == "D"
    assert counts["s5-pilot"]["quiet_inverted_with_signal"] == 0


def _run(tmp_path, name, stars):
    """A run dir: stars = {tic: [period, ...]}; 12 transits of a 3 h box."""
    d = tmp_path / "runs" / name
    d.mkdir(parents=True)
    lines = [
        json.dumps(
            {
                "tic_id": tic,
                "ok": True,
                "result": {
                    "config": {},
                    "signals": [
                        {
                            "period_days": p,
                            "epoch_btjd": 1.0,
                            "duration_hours": 3.0,
                            "depth_ppm": 160.0,
                            "snr": 9.0,
                            "sde": 14.0,
                            "n_transits": 12,
                        }
                        for p in ps
                    ],
                },
            }
        )
        for tic, ps in stars.items()
    ]
    (d / "results.jsonl").write_text("".join(ln + "\n" for ln in lines))
    (d / "summary.json").write_text(json.dumps({}))
    return d


def test_a_planet_with_timing_variations_claims_more_widely(tmp_path):
    """TOI-1408 c: 2.1664 d confirmed, 2.2222 d as a TOI; one sector gave 2.1777 d,
    0.522% from the confirmed period against a tolerance of 0.522%. The archive
    flags the planet for timing variations, so the star's claims widen."""
    _run(tmp_path, "s76", {1: [2.1777136], 2: [2.1777136]})
    plain = [u.Claim("planet X c", 2.1664), u.Claim("TOI 9.02 (CP)", 2.2222)]
    ttv = [u.Claim("planet Y c", 2.1664, True), u.Claim("TOI 8.02 (CP)", 2.2222)]
    found, _ = u.collect(tmp_path / "runs", [], [("76", "s76")], [], {1: plain, 2: ttv})
    by = {f.tic_id: f.claimed_by for f in found}
    assert by[1] == []
    assert by[2] == ["TOI 8.02 (CP)", "planet Y c"]


def test_planet_claims_read_the_ttv_flag(tmp_path):
    p = tmp_path / "ps.csv"
    p.write_text(
        "pl_name,tic_id,hostname,pl_orbper,disc_facility,disc_year,ttv_flag\n"
        '"A b","TIC 5","A",2.1664,"TESS",2024,1\n'
        '"B b","TIC 6","B",3.0,"TESS",2024,0\n'
    )
    claims = u.read_planet_claims(p)
    assert claims[5][0].ttv and not claims[6][0].ttv


def test_the_csv_carries_the_pileup_tier(tmp_path):
    """The CSV was written before pile-ups were marked, so a pile-up candidate
    read tier A or B there and tier S in the report."""
    stars = {i: [1.0 + i * 1e-5] for i in range(8)}
    stars.update({100 + i: [0.5 * 1.02**i] for i in range(60)})
    _run(tmp_path, "s5", stars)
    _run(tmp_path, "s5-inv", {})
    cat = tmp_path / "cat"
    cat.mkdir()
    (tmp_path / "toi.csv").write_text(
        "TIC,Full TOI ID,TOI Disposition,Orbital Period (days) Value\n"
    )
    (tmp_path / "ctoi.csv").write_text("CTOI|TIC ID|Period (days)\n")
    (tmp_path / "ps.csv").write_text("pl_name,tic_id,hostname,pl_orbper\n")
    out = tmp_path / "u.csv"
    assert u.main([
        "--runs", str(tmp_path / "runs"), "--catalogs", str(cat),
        "--toi", str(tmp_path / "toi.csv"), "--ctoi", str(tmp_path / "ctoi.csv"),
        "--confirmed", str(tmp_path / "ps.csv"), "--sample", "5=s5:s5-inv",
        "--csv", str(out), "--json", str(tmp_path / "u.json"),
    ]) == 0  # fmt: skip
    tiers = [ln.split(",")[0] for ln in out.read_text().splitlines()[1:]]
    data = json.loads((tmp_path / "u.json").read_text())
    assert tiers.count("S") == data["tiers"]["S"] == 8
