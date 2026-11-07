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

## Sector 2: the first sector nothing was tuned on (2026-11-05)

Spline default (7a), unchanged. `python -m tesstransit.planets --sector 2`
and `python -m tesstransit.pilot --sector 2` (seed 1: 60 TCE stars across
SPOC's SNR range + 40 random no-TCE stars). SPOC's sector-2 table: 1,062
TCEs on 726 stars; 15,994 2-minute light curves. 0 failed stars in either run.

| | Sector 1 | Sector 2 |
|---|---|---|
| TCE-backed planet-like TOIs recovered | 61 / 64 (95%) | **61 / 62 (98%)** |
| Confirmed + known planets among them | 38 / 38 | **40 / 40** |
| Planet-like TOIs with no sector TCE, recovered | 5 | **15** |
| TCE sample recovered (incl. alias) | 34 / 87 | 53 / 87 |
| No-TCE stars with a signal | 4 / 40 | 6 / 40 (10 signals) |
| Signals on TOI stars matching no TOI | 29 | 40 |

**Planet recovery holds on unseen data.** The one TCE-backed miss is TOI
214.01 (PC, 18.6 d, two transits, forced SNR 7.7; SPOC 8.4): near threshold
in both pipelines.

**TCE recovery is not comparable between sectors.** Sector 1's sample has a
median TCE period of 10.7 d (48 of 87 above 10 d); sector 2's 5.7 d (29 of
87). Short periods mean many transits. Sector 1, the mission's first, carries
many long-period TCEs that vetting discarded; the rate measures the
reference as much as the search.

**Fifteen planet candidates SPOC made no sector-2 TCE for were found** at
our SNR 7.6–13 (forced 5.0–11.4), five of them confirmed planets (125.02,
178.03, 206.01, 220.01, 237.01). SPOC reached these only by combining
sectors. A chance match is implausible (period within 0.5% AND transits
within half a duration), but these sit near threshold, where the search's
optimistic SNR (best of many trials) is largest.

### A bug the second sector found

`_tce_backed` paired each TOI with the FIRST TCE within 1% of any alias.
On L 98-59 (TOI 175) it gave planet d (7.4507 d) planet c's TCE (3.6905 d;
2x is 0.94% off) while d's own TCE (7.451 d) went unused, so d was scored
against c's transit times: "missed, weak after detrend", forced SNR 1.3,
for a planet SPOC found at SNR 23.7. Backing is now assigned per star like
the matcher (direct first, closest period, each TCE once, 0.5%). It changed
one TOI per sector: 175.02 is now recovered; in sector 1, 184.01 (an EB)
moved from a direct to an alias match. Sector 1's headline is unchanged.

## Sectors 1–13: the first year of the southern sky (2026-11-05)

Spline default, no further tuning; sectors 3–13 run back to back on 8
workers in 37 minutes. Per sector: every TOI listing it on a 2-minute star,
plus one TCE sample (60 TCE stars across SPOC's SNR range + 40 random no-TCE
stars). One star failed of ~4,000 (sector 4, TIC 300038935: median PDCSAP
flux negative, refused at load rather than searched).

| Sector | Planet TOIs SPOC detected, found | CP + KP | No-sector-TCE planets in data, found | TCE sample | No-TCE stars with a signal |
|---|---|---|---|---|---|
| 1 | 61 / 64 | 38 / 38 | 5 / 6 | 34 / 87 | 4 / 40 |
| 2 | 61 / 62 | 40 / 40 | 15 / 16 | 53 / 87 | 6 / 40 |
| 3 | 33 / 36 | 22 / 23 | 5 / 6 | 49 / 82 | 4 / 40 |
| 4 | 43 / 47 | 26 / 28 | 8 / 12 | 47 / 80 | 9 / 40 |
| 5 | 47 / 52 | 31 / 32 | 11 / 12 | 53 / 90 | 7 / 40 |
| 6 | 46 / 50 | 32 / 35 | 7 / 8 | 52 / 83 | 7 / 40 |
| 7 | 47 / 47 | 27 / 27 | 5 / 6 | 54 / 92 | 10 / 40 |
| 8 | 34 / 35 | 18 / 19 | 9 / 9 | 57 / 83 | 7 / 40 |
| 9 | 76 / 78 | 36 / 36 | 4 / 4 | 50 / 83 | 6 / 40 |
| 10 | 64 / 68 | 29 / 33 | 6 / 7 | 48 / 95 | 6 / 40 |
| 11 | 46 / 47 | 25 / 26 | 10 / 11 | 38 / 85 | 7 / 40 |
| 12 | 47 / 47 | 17 / 17 | 12 / 12 | 44 / 75 | 14 / 40 |
| 13 | 71 / 74 | 28 / 29 | 7 / 7 | 51 / 78 | 5 / 40 |

Pooled, with Wilson 95% intervals:

| Measure | Count | Rate | 95% interval |
|---|---|---|---|
| Planet-like TOIs SPOC detected in the same sector | 676 / 707 | **95.6%** | 93.8–96.9% |
| Confirmed + known planets among them | 369 / 383 | **96.3%** | 94.0–97.8% |
| Planet TOIs with no sector TCE, strong enough in the sector | 104 / 116 | 89.7% | 82.8–94.0% |
| TCE sample (in the period range) | 630 / 1,100 | 57.3% | 54.3–60.2% |
| No-TCE stars with a signal (unvetted) | 92 / 520 | 17.7% | 14.7–21.2% |

Sector 1–2's 100% on confirmed planets did not hold: 14 of 383 were missed
over the year. The rate is now measured to about ±2 points.

### The 31 TCE-backed planet misses

| Cause | Misses | Examples |
|---|---|---|
| Fast stellar variability (< 1 d) the 0.5 d knot floor cannot follow; its residual takes the signal slots and swamps the transit | ~10 | 1066.01 (KP, 13,793 ppm, forced 2.7: a 0.63 d variable), 6551.01, 540.01 (twice), 837.01, 140.01 |
| Near threshold in both pipelines (forced SNR ≲ 8.5, SPOC ≲ 9) | ~12 | 214.01, 697.01, 536.01, 1080.01 |
| Reached, then refused by the SDE gate | ≥ 2 | 216.02 (CP, TTV pair in resonance: best peak SNR 42 at a 4.27 d blend of both planets, SDE 8.3 < 11) |
| Multi-planet star: signal slots used up first | 2 | 561.01, 175.03 (L 98-59 b) |
| Other / not yet diagnosed | ~5 | 179.01 (forced 12.6), 1233.01 (forced 11.4; a peak at 14.16 d was rejected) |

The first two dominate, and only the first is a fixable loss: variability
faster than ~1 day needs either a variability model (a periodogram-driven
sinusoid fit, or a GP) or a knot floor tied to the trial duration, which
pilots 6b/6c and 7b showed costs long transits. Counts are from an automatic
classification and need a per-star look before quoting.

## Pre-whitening fast stellar variability, sectors 1–13 (2026-11-05)

`--prewhiten` (off by default): before the spline, divide out up to three
Lomb-Scargle periodicities faster than 1.5 d, each only if significant (FAP <
1e-4) and two-sided in the folded data (brightest-bin rise ≥ 0.3 × faintest-bin
dip), so transits and eclipses, which only dim the star, are left in. Same
stars as the default runs.

| Sectors 1–13 | Default (spline) | Pre-whitened |
|---|---|---|
| Planet-like TOIs SPOC detected, found | **676 / 707 (95.6%, 93.8–96.9)** | 657 / 707 (92.9%, 90.8–94.6) |
| Confirmed + known | 369 / 383 | 368 / 383 |
| No-TCE stars with a signal | 92 / 520 (17.7%) | **44 / 520 (8.5%)** |
| TCE sample | 630 / 1,100 | 487 / 1,093 |
| EB / V / IS TOIs found | 281 | 244 |
| Signals on TOI stars matching no TOI | 491 | 348 |

**Rejected as a default.** It did what it was built for: 10 planets gained, all
on fast variables (1066.01, 837.01, 540.01 twice, 4325.01, 4180.01, 140.01,
179.01, 4364.01, 175.03). But it lost 29, including pi Men c (144.01, sector 8)
and TOI 407.01 in five sectors, for a net −19.

Why, from the losses:

- **It fires on stars that barely vary.** On TOI 396's bright, quiet star it
  removed 40–50 ppm terms that only just pass the false-alarm test, which was
  enough to push 180–233 ppm planets under threshold.
- **The two-sidedness test does not separate.** A real win (1066.01's 0.63 d
  variable) scored 0.43; a loss (4181.01, a 24,380 ppm dip every 0.6 d,
  probably an EB) scored 0.48. No threshold keeps one and drops the other.
- **The fit absorbs part of the transit** on variable stars (407.01: forced
  SNR 14 → 5.8): clipping individual low points does not protect a transit
  whose depth is comparable to the per-point scatter.
- **It removes what SPOC also reports**: contact binaries and short-period
  variables are TCEs and EB/V TOIs, so agreement drops (630 → 487 TCEs).

What it does well, it does very well: no-TCE stars with a signal halve (92 →
44), and stray signals on TOI stars fall 29%. A variant worth testing gates
it on amplitude relative to the noise (so quiet stars are left alone) and
masks the strongest BLS events before fitting; or use it as a second pass only
on stars where the default leaves a sub-day variability peak.

### Version 2: gated, with transit masking (2026-11-05)

Pre-whitening only while the hour-binned spline residual scatters ≥ 2× white
noise, and the variability fit refitted without stretches dipping below 3σ in
1-hour running means (b29d381). Same stars, sectors 1–13.

| Sectors 1–13 | Default (spline) | Pre-whitening v1 | **Pre-whitening v2** |
|---|---|---|---|
| Planet-like TOIs SPOC detected, found | 676 / 707 (95.6%) | 657 (92.9%) | **682 / 707 (96.5%, 94.8–97.6)** |
| Confirmed + known | 369 / 383 | 368 | **373 / 383** |
| No-TCE stars with a signal | 92 / 520 (17.7%) | 44 (8.5%) | 70 / 520 (13.5%) |
| EB / V / IS TOIs found | 281 | 244 | 280 |
| Signals on TOI stars matching no TOI | 491 | 348 | 410 |
| TCE sample | 630 / 1,100 | 487 | 539 / 1,097 |

**Planets: 8 gained, 2 lost (+6).** Gained, all on fast variables: 140.01,
540.01 (twice), 865.01, 4325.01, 4180.01, 837.01, 1066.01. Lost: TOI 252.01 in
sectors 2 and 11, a 1.002 d planet on a star varying at 0.998 d, which no
periodicity filter can separate. The gate restored every v1 loss checked.

**The TCE-sample drop is SPOC's false alarms, not planets.** All 109 TCEs v2
stops recovering never became TOIs, and 107 have periods under 1.5 d (median
0.36 d): stellar variability SPOC flagged and its vetters rejected. 18 TCEs
were gained, one a known planet. EB/V TOIs are unchanged (281 → 280).

The planet gain is small next to the interval (+0.9 points, intervals
overlap), but it is paired, same stars, and every gain and loss has a
mechanism. A sign test on 8 gained vs 2 lost gives p ≈ 0.055 (one-sided).

**Adopted as the default (2026-11-06).** `tesstransit.pilot` and
`tesstransit.planets` now pre-whiten by default with the spline
(`--no-prewhiten` turns it off; the median and biweight detrenders ignore
it). `SearchConfig` itself still defaults to off, so library callers opt in
with `prewhiten_max_period_days=1.5`. The published page and dashboard show
the v2 runs, as scorecard variant 7c. Its sector-1 hold-out (seed 2, same
stars as `holdout-spline7a`) found 30 of 86 TCEs against 7a's 35, and
signals on 2 of 40 no-TCE stars against 5, the same trade as sectors 1-13.

## FFL workflow reproduces the CLI runs (2026-11-06)

`tess.workflows.BuildResults` over sectors 1-13, run in memory with the real
handlers and the real search (8 workers per chunk, chunks of 50), against the
same TOI catalogue copy as the published runs: 85 chunks, each searched once,
46 minutes. Against the CLI's pre-whitened runs:

- 25 of 26 run directories are identical: `tois.csv` / `matches.csv`,
  `results.jsonl` (timings aside) and `summary.json`.
- The 26th, sector 4's planets, differs only in the text of one error: TIC
  300038935 (negative median PDCSAP flux) fails in both, but the message
  carried the local scratch path. Errors now name only the file.
- The page's embedded data is identical: 3,149 TOI rows, the sector 2-13 TCE
  sample rows, and the scorecard. Sector 1 matches the published development
  sample; the published page also shows a hold-out, which `BuildResults` does
  not run.

## False-alarm rate from inverted light curves (2026-11-06)

The thresholds (SDE 11, SNR 7.1) were set on synthetic white noise. To
measure the false-alarm rate on real TESS noise, the 1,300 sample stars of
sectors 1-13 (the same stars, seeds and pre-whitened search as `s{N}-pilot-pw2`)
were searched again with the flux inverted (`2 - flux`, `pilot --invert`, Kepler
DR25's test). A transit becomes a brightening that a dip search cannot find,
so every inverted signal is a false alarm. 1,298 stars searched; ~1 min per
sector of 100 stars.

**No-TCE stars: 73 of 520 (14.0%, 11.3-17.3%) give an inverted signal; the
normal search gives a signal on 70 of the same 520 (13.5%, 10.8-16.7%).** The
rates are the same, so the no-TCE signals are, as a group, what false alarms
alone produce. They also look alike: under 1 day for 87 of 122 inverted and
98 of 118 normal signals, with many transits each. That is periodic stellar
variability surviving pre-whitening, not white noise, which a threshold
cannot separate from a transit.

**Stars with a deep real signal are not neutral under inversion.** Inverted,
75% of TCE stars give a signal. Of those, 81 match a TCE's ephemeris and all
81 are TCEs that never became TOIs (variability SPOC flagged). Recovered
planets show the mechanism: every one whose period reappears inverted is deep
(4,000-37,000 ppm: TOI 418.01, 187.01, 507.01, 656.01, 781.01, 1084.01 ...),
and the inverted signal sits at P, P/2 or P/3, wider and much shallower -- the
flattening around the inverted transit, now a large bump, leaves dips beside
it. Eclipsing binaries do the same. So the false-alarm rate is taken from the
no-TCE stars only.

**A stricter cut costs planets faster than it removes false alarms**
(false alarms per no-TCE star, inverted | planets SPOC detected in the sector,
found again; recovery at each cut from the matched signal and every signal
ranked above it):

| cut | false alarms | planets |
|---|---|---|
| SDE 11, SNR 7.1 (current) | 14.0% | 682 / 707 (96.5%) |
| SDE 12 | 11.9% | 676 (95.6%) |
| SDE 13 | 10.8% | 667 (94.3%) |
| SDE 14 | 9.8% | 657 (92.9%) |
| SDE 16 | 8.1% | 639 (90.4%) |
| SDE 20 | 5.4% | 593 (83.9%) |
| SDE 11, SNR 10 | 10.4% | 634 (89.7%) |

The current cut is kept. Looser cuts cannot be scored from these runs (a
signal the current cut rejected was never matched), but they only add false
alarms: 16.0% at SDE 10, 21.5% at SDE 9.

**A symmetry veto was tried and rejected.** Dropping a normal signal whose
period (or 2x, 1/2, 3x, 1/3) also appears in the same star's inverted search
would remove 46% of no-TCE signals, but also 10 of 57 recovered planet signals in the
sample (18%), all deep hot Jupiters, for the reason above.

**For a full-sector run** this means about 14% of the ~16,000 stars per sector
-- some 2,200 -- will report a signal that is a false alarm. Those signals
stay unvetted; separating them needs a shape test (odd/even depths, secondary
eclipse, transit shape), not a stricter threshold.

## Completeness from injected transits (2026-11-06)

`python -m tesstransit.inject` puts ONE synthetic transit into the real
PDCSAP flux of a star, before detrending, and searches it with the unchanged
pipeline (pre-whitening on). The stars are the 520 no-TCE stars of the
`s{N}-pilot-pw2` samples, 6 independent injections each: 3,120 searches, none
failed, ~70 s per sector on 12 workers. Per injection: period log-uniform
0.5-13 d (at least two transits on data), impact parameter 0-0.9, the duration
a Sun-like star gives, a trapezoidal shape, and an **expected SNR**
log-uniform 4-40 that sets the depth from the star's own scatter
(`depth = snr * sigma / sqrt(n_in)`, `sigma` the MAD of the detrended,
un-injected flux). Recovered = matched at the injected period and epoch, or
an alias, with the TCE-matching tolerances (43 of 1,887 were aliases).

**Recovery by expected SNR** (all stars | stars the normal search did NOT
pre-whiten, 82% of injections):

| expected SNR | all | without fast variability |
|---|---|---|
| 4-6 | 3.7% | 4.3% |
| 6-7 | 18.4% | 22.1% |
| 7-8 | 43.3% | 51.2% |
| 8-10 | 56.0% | 68.8% |
| 10-12 | 75.0% | 92.7% |
| 12-15 | 76.8% | 91.3% |
| 15-20 | 83.3% | 98.7% |
| 20-40 | 90.1% | 98.8% |

A logistic fit puts 50% at SNR 7.9 and 90% at 11.8 for ordinary stars (all
stars: 9.6 and 21.2). The search's measured SNR on a recovered injection is
0.95 of the expected (16-84%: 0.84-1.09), so the 7.1 threshold sits at an
expected SNR of ~7.5, where the curve crosses 50%: on white-ish noise the
pipeline loses almost nothing to detrending. Period, impact parameter, number
of transits and duration barely matter above SNR 15 (84-92% in every bin; 2
transits 84%, 20+ 90%).

**The 10% loss at high SNR is the pulsating stars.** Of 92 injections missed
at expected SNR >= 20, 83 are on the 93 stars (18%) whose normal search
removed a periodicity; on every other star 9 of 775 are missed. Measured at
the injected ephemeris with red-noise uncertainties, a sample of those misses
has SNR 2-13 (median ~5) against an expected 20-40: the white-noise SNR
overstates what can be detected on a pulsator, and the search stops on the
star's own variability first (74 of 92 stop below threshold). In 6 of 25
checked, the search had the right period at SNR > 7.1 but SDE 7-11: the
pulsations widen the periodogram, so the SDE cut is what rejects them.

**Pre-whitening helps on exactly those stars.** The same 558 injections on the
93 pre-whitened stars, searched with it on and off: at expected SNR >= 20,
48% found with it, 36% without. At SNR >= 15 it rescues 35 injections and
loses 7 (3 of them where a removed period is a sub-harmonic of the injected
one); where both find the transit the SNR is the same (median ratio 1.01).

**So**: on the ~80% of stars without fast variability the search is complete
to ~99% above expected SNR 15 and at the level its threshold promises. On
pulsating stars completeness is about half even for strong transits, and a
better treatment of their variability (a fuller pre-whitening, or a red-noise
aware SDE) is the lever, not a lower threshold, which the false-alarm test
shows is already spent on the same stars.

## Pre-whitening v3: more terms, faster periods (2026-11-06)

The injections showed pulsating stars are where completeness is lost, so the
pre-whitening was tuned on them, then checked on everything. `--set
FIELD=VALUE` (on `pilot`, `planets` and `inject`) overrides any `SearchConfig`
field, and `inject --replay` re-searches the identical injected transits, so
each variant is paired with the default star by star.

**Screen** (69 no-TCE stars whose default search used all three terms, 414
injections, the same stars inverted, and the 18 TCE-backed planets on
pre-whitened stars): ten terms instead of three found 68 of 171 injections at
expected SNR >= 15 against 52; adding periods down to 0.01 d (delta Scuti
modes under an hour) gave 73; lowering the gate to 1.5 gave 74 with false
alarms 23 -> 17 of 69. Twenty terms or fifteen gained nothing over ten.
Extending the period range to 3 d lost 2 of the 18 planets and was dropped.

**Full check, sectors 1-13**, both candidates against the default (pw2):

| | Default | pwA: 10 terms, 0.01 d, gate 1.5 | **pwB: 10 terms, 0.01 d** |
|---|---|---|---|
| planets (TCE-backed TOIs) | 682/707 | 681 (+2 −3) | **682 (+1 −1)** |
| false alarms (inverted, 520 no-TCE stars) | 73 (14.0%) | 57 (11.0%) | 70 (13.5%) |
| signals on no-TCE stars (normal search) | 70 | 53 | 66 |
| pilot TCEs recovered | 539/1097 | 519/1098 | 536/1099 |
| injections found, pulsators, SNR >= 15 | 88/229 | 115 (+38 −11) | **112 (+31 −7)** |
| injections found, other stars, SNR >= 15 | 1064/1077 | 1058 (+2 −8) | **1065 (+1 −0)** |
| planet run wall time | 25 min | 25 min | 25 min |

The lower gate removes more false alarms, but it also starts pre-whitening
quieter stars: 8 injections lost on ordinary stars, a planet net lost, and 8 of
its 36 lost pilot TCEs are on stars only it pre-whitens, often at the TCE's own
period (0.604 d TCE, 0.6035 d term removed). pwB keeps the gate at 2, so it
changes only stars already pre-whitened, and is the new default
(`prewhiten_max_terms=10`, `prewhiten_min_period_days=0.01`). Every pilot TCE
either variant loses is a non-TOI signal at 0.25-1.4 d, i.e. variability or a
contact binary that the pre-whitening now removes as intended.

Completeness with the new default (groups as above, by the default run's
pre-whitening):

| | 50% found at | 90% found at | found, expected SNR >= 20 |
|---|---|---|---|
| all stars | SNR 9.4 (was 9.6) | SNR 19.7 (21.2) | 92.0% (90.1%) |
| stars without fast variability | 7.9 (7.9) | 11.7 (11.8) | 99.0% (98.8%) |
| pulsating stars | 24.7 (28.8) | 56.2 (60.6) | 57.9% (47.8%) |

Pulsators remain the weak spot: ~40% of strong transits on them are still
missed. What is left is not periodic enough to fit with sinusoids
(amplitude-modulated modes, beating) and needs a red-noise-aware detection
statistic rather than more terms.
