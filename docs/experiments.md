# Experiment log

Every change to the search, and what it did, measured on the SAME sample:
sector 1, 60 TCE stars spread evenly over SPOC's SNR range (87 TCEs) + 40
seeded-random stars with no TCE (`python -m tesstransit.pilot --sector 1`,
`--seed 1`). Terms are defined in [definitions.md](definitions.md).

⚠️ Every threshold below was chosen while looking at this sample. Its numbers
are a development record, not a result: a threshold must be confirmed on a
sample it was not tuned on before it is quoted.

| Pilot | Change | Recovered | Alias | Outside grid | Not detected | Exhausted | No-TCE stars with a signal | s/star (8 workers) |
|---|---|---|---|---|---|---|---|---|
| 1 | Baseline: 1-d running median, BLS 0.5–13 d, durations 1–8 h, SNR ≥ 7.1, SDE ≥ 7 | 14 | 13 | 35 | 23 | 2 | 7 / 40 | 6.8 |
| — | Matcher fix: direct matches before aliases (TIC 350479101) | 14 | 13 | 35 | 23 | 2 | — | — |
| 2 | Duration capped by period (Kepler's law, ρ ≥ 0.1 ρ☉); grid to 0.25 d, + 0.5 h | **26** | 4 | 33 | 18 | 6 | **2 / 40** | 4.8 |
| 3 | Periods to 26 d, ≥ 2 transits with data; SDE on detrended periodogram, ≥ 11 | **32** | 4 | **0** | 34 | 17 | 8 / 40 | 6.4 |
| 4 | Per-transit veto (events other than the strongest ≥ SNR 3, red-noise σ) | 32 | 3 | 0 | 41 | 11 | 8 / 40 | 7.5 |
| 5 | Veto bug fixed: event depth from the out-of-transit mean, not 1.0 | 32 | 4 | 0 | 36 | 15 | 8 / 40 | 7.7 |
| 6a | Biweight (wotan) instead of running median, fixed 1-d window | 31 | 4 | 0 | 38 | 14 | 7 / 40 | 5.8 |
| 6b | Running median, one window per trial duration (3×) | 27 | 5 | 0 | 41 | 14 | 8 / 40 | 16.2 |
| 6c | Biweight, one window per trial duration (3×) | 27 | 6 | 0 | 43 | 11 | 7 / 40 | 15.2 |

**Current default: pilot 5** (running median, fixed 1-day window).

## What each step taught

1. **The grid decided more than the search.** The 13-day cap excluded 40% of
   TCEs before anything was searched: SPOC reports two-transit TCEs out to
   ~24 d.
2. **An 8-hour box every 0.5–1.1 days is not a transit.** Most early
   no-TCE "signals" were in transit 30–65% of the time. A physical duration
   cap removed them and, as a side effect, moved 10 eclipsing binaries from
   an alias to their true period.
3. **SNR alone accepts noise; SDE must be detrended.** On white noise the best
   of ~10⁷ boxes reaches SNR 6.4–7.8. Raw BLS power rises with period, so a
   plain SDE marked a real 2-transit planet (SNR 49) at 6.5.
4. **White-noise calibration does not transfer.** On real light curves the
   SDE gate failed both ways: no-TCE signals at SDE 11–76, missed TCEs whose
   best peak had SNR 10–75 but SDE 5.6–10.9.
5. **The veto is correct but is not what limits false alarms.** The no-TCE
   signals are consistent across events (the dip is present without the
   strongest event, SNR 3.9–9.5). Its first version vetoed real eclipsing
   binaries — a baseline bug, caught because the pilot recorded every veto.
6. **Detrending method barely matters; a per-duration window hurts.** A 3×
   window is safe for a box that matches the transit, but transits are
   longer than the box that finds them: four of five lost TCEs were 4–8 h
   transits at 10–14 d previously found with 3–4 h boxes. And the short-period
   no-TCE signals (TIC 260239272 at 0.356 d, 276936458 at 0.460 d) survived a
   6-hour window — evidence they are real periodic stars that SPOC did not
   turn into TCEs, not leftover variability. Unvetted either way.

## Operational

- Pilots 6b/6c first ran during a network loss: 33 and 100 of 100 fetches
  failed with `EndpointConnectionError`. Each failure was recorded per star,
  so the runs were visibly invalid rather than silently short. The S3 client
  now retries 10× with backoff; both were re-run clean.

## Hold-out check (2026-11-04)

Same configuration as pilot 5, on a sample sharing **no** stars with it
(`--seed 2 --holdout-of runs/pilot5-veto-fixed`): 60 TCE stars, one drawn
per SNR bin (91 TCEs), + 40 new no-TCE stars.

| Sample | In-grid TCEs | Recovered (incl. alias) | Not detected | Exhausted | No-TCE stars with a signal |
|---|---|---|---|---|---|
| Development (pilot 5) | 87 | 36 (41%) | 36 | 15 | 8 / 40 |
| Hold-out | 86 | 34 (40%) | 37 | 15 | 7 / 40 |

The development numbers hold on unseen stars: no sign that the tuning fitted
the sample. The hold-out's no-TCE signals are of the same kind — short
periods (0.29–2.1 d), many transits, some plainly eclipsing binaries (TIC
220507474: 6,939 ppm every 0.963 d).

## Why TCEs are missed

From what each search recorded (strongest rejected peak, vetoes, signals):

| Class | Development | Hold-out |
|---|---|---|
| Never the strongest peak | 29 | 29 |
| Found at the TCE's period, rejected by the SDE gate | 4 | 2 |
| A signal within 3% of the period but outside the match tolerance | 3 | 6 |

Thresholds are not the main limit. Of the 29 development misses that were never
the strongest peak, 20 have periods of 10–24 d (two or three transits), 6 have
durations above 8 h — beyond the longest box searched (**182 of 1,457
sector-1 TCEs are longer than 8 h**, 92 longer than 12 h) — and 3 have a SPOC
model SNR below 7.1.

## What became of the TCEs: the TOI catalogue

TCEs are unvetted detections. Joining each to the TOI catalogue (same star,
period within 1%, alias-aware):

| Outcome | Development | Hold-out |
|---|---|---|
| Recovered | 36 — 2 TOIs (1 CP, 1 KP) | 34 — 4 TOIs (3 CP, 1 EB) |
| Not detected | 36 — 1 TOI (PC: TOI 119.01) | 37 — 0 TOIs |
| Search exhausted | 15 — 0 TOIs | 15 — 0 TOIs |

**Every confirmed or known planet in both samples (5 of 5 CP/KP) was
recovered. Of 73 not-detected TCEs, 1 became a TOI.** Most of what the search
misses are TCEs SPOC's own vetting never promoted. ⚠️ "No TOI" is not proof of
a false alarm (vetting also uses later sectors), and seven TOIs are far too
few to quote a planet recovery rate — the sample was drawn across TCE SNR, not
across planets. Recovery against TCEs measures pipeline agreement; recovery
against TOIs is what bears on planets, and needs a sample built for it.

## Planet sample: every sector-1 TOI (2026-11-04)

`python -m tesstransit.planets --sector 1` — pilot-5 configuration, no tuning.
445 TOIs list sector 1; **230 are on 186 stars with a 2-minute light curve**,
and all 186 were searched (0 failures, 143 s on 8 workers). Each TOI is
matched against its own ephemeris and measured at it ("forced") on our
detrended flux; outcomes are defined in [definitions.md §3.6](definitions.md).

| Planet-like TOIs (CP, KP, PC) | n | Recovered | Missed, in data | Weak in sector | Outside grid / single |
|---|---|---|---|---|---|
| SPOC made a sector-1 TCE | 64 | **59 (92%)** | 5 | 0 | 0 |
| No sector-1 TCE | 138 | 3 | 1 | 84 | 45 + 5 |

**Against planets, the search recovers 92% of what SPOC detected in the same
sector, against 41% of raw TCEs.** CP + KP: 36 of 38; PC: 23 of 26. Of 28
EB / V / IS TOIs, 14 recovered (1 as an alias), 9 weak, 5 outside the grid.

The five TCE-backed misses:

| TOI | SPOC SNR | Forced SNR (ours) | Outcome | Reading |
|---|---|---|---|---|
| 200.01 (DS Tuc Ab, CP) | 25.4 | 0.9 | weak_after_detrend | Young, active star; the 1-day median leaves variability that swamps a 3.3 h transit. Our search reported three rotation residues at SNR 85–144 instead. |
| 140.01 (PC) | 11.5 | 2.8 | weak_after_detrend | A 1.12 d variable on the same star (SNR 55) dominates the noise. |
| 119.01 (PC) | 8.0 | 5.7 | weak_after_detrend | Near threshold in both pipelines. |
| 210.01 (PC) | 7.9 | 6.8 | weak_after_detrend | Near threshold in both pipelines. |
| 139.01 (CP) | 10.5 | 7.2 | not_detected | Two transits; our strongest peak was 14.3 d. |

So the remaining planet losses are **stellar variability, not the grid or the
thresholds**: the two clear ones are active or variable stars where SPOC's
detrending keeps a transit ours erases. That points at detrending that models
the star (e.g. a GP or iterative masked fit) rather than at the search.

Without a sector-1 TCE, almost nothing is in the sector-1 data — 84 of 88
in-grid TOIs are weak in the sector by our own forced measurement, which is
why SPOC needed more sectors for them. Of the 4 that are not, we recovered
3 that SPOC's sector-1 run did not make TCEs for (TOI 209.01, 487.02,
708.01, at our SNR ≥ 7.1 though forced 5.7–7.0 — near threshold) and missed
TOI 174.01 (CP, 2 transits, forced 11.3), which SPOC also missed in sector 1.

⚠️ Caveats. The forced SNR is measured at the known ephemeris, an upper bound
for a blind search. "Weak in sector" uses OUR noise model, which is the thing
that failed on DS Tuc; only where SPOC detected the TOI in this sector can a
low forced SNR be called our loss. 52 signals on these stars match no TOI;
they are unvetted. One sector, 64 TCE-backed planet-like TOIs: a 92% rate has
a binomial 1σ of about ±3.4 points.

### Two bugs the planet sample found

1. **Matcher epoch drift.** `_alias_of` wrapped the raw epoch offset by our
   period; with a TOI epoch ~350 orbits away, a 0.01% period difference
   became hours. Hot Jupiters found at SNR 340 were "not detected" (15 of 64
   TCE-backed TOIs). The reference epoch is now propagated with its own
   period first. Re-matching pilot 5 and the hold-out gives identical counts,
   as expected for in-sector TCE epochs.
2. **Catalogue ephemerides don't always reach back.** TOI 121.01's
   ephemeris, propagated 50 orbits, lands 12.4 h from transits both
   pipelines found. TCE-backed TOIs now use SPOC's in-sector epoch with the
   TOI period; others keep the catalogue epoch with its propagated
   uncertainty reported (1 of 151 exceeds the match tolerance).

## Detrending that follows the star (2026-11-05)

Each variant scored on all three samples: TCE-backed planet-like TOIs
(planet sample), TCEs recovered incl. alias and no-TCE stars with a signal
(development and hold-out samples).

| Variant | Planets (of 64) | Planet-sample signals matching no TOI | Dev: TCEs / no-TCE stars with signal | Hold-out: TCEs / no-TCE | s/star |
|---|---|---|---|---|---|
| Pilot 5: 1-day running median | 59 | 52 | 36 / 8 | 34 / 7 | 7.7 |
| **7a: clipped cubic spline, knots 0.5–3 d by BIC per star** | **61** | **29** | 34 / **4** | 35 / **5** | 5.1 |
| 7b: 7a + a 0.25–0.4 d knot layer for ≤ 2 h boxes on stars at the 0.5 d floor | 61 | 29 | 35 / 4 | 34 / 5 | 4.4 |

**New default: 7a.** Planets +2 net (gained DS Tuc Ab at SNR 59, 139.01,
210.01; lost 4186.01), false alarms on no-TCE stars roughly halved, and
signals on TOI stars that match no TOI down 44%. TCE agreement is flat
(dev −2, neither a TOI; hold-out +3 −2).

- **DS Tuc Ab** shows why: the 1-day median left flat-bottomed residuals ~1%
  deep on a 2.85 d, 2% rotation; the 0.5 d spline leaves ~0.2% wiggles, and
  the forced SNR went 0.9 → 7.1.
- **4186.01** was not eaten by the detrend (forced 13.4 → 11.1); an unrelated
  11 d peak (SNR 14) outranked its two transits. A near-threshold swap.
- **7b was dropped, and its code removed.** On TOI 140.01's star (a 1.12 d,
  0.25% oscillation) the flexible layer does expose the transit (forced 3.1
  → 7.9 at 0.3 d knots), but the stiff layer's residual harmonics at 0.887 d
  took all three signal slots first. And the 0.5 d floor is hit on 28–40% of
  stars, not a few: the layer added trials everywhere for no measured gain.
  Fixing 140.01 means not searching the stiff layer on such stars, which
  re-opens the long-transit loss of 6b/6c.
- On white noise the spline default keeps the noise peaks of the three
  test seeds below SNR 7.1 (5.3–6.6, vs 7.2–7.8 with the median).

Remaining TCE-backed planet misses with 7a: 140.01 (above), 119.01 (SPOC SNR
8.0, ours 6.5), 4186.01.
