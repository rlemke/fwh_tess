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
