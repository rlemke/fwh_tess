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

## Red-noise-aware power (2026-11-06; kept off)

`red_noise_power` divides each trial duration's BLS power and SNR by beta^2
and beta, beta being the scatter of the flux binned to that duration over
the white-noise prediction (Pont et al. 2006), so a box the length of a
pulsation residue no longer looks significant. Screened on the 69 no-TCE
stars whose default search used every pre-whitening term (414 injections,
the same stars inverted, 18 TCE-backed planets), with v3 pre-whitening:

| | strong injections found (SNR >= 15, of 171) | false alarms (of 69) | planets (of 18) | s / search |
|---|---|---|---|---|
| default | 73 | 20 | 14 | 2.4 |
| red-noise power | 68 (+2 −7) | 13 | 14 | 8.0 |
| red-noise power, SDE 9 | 68 | 16 | 14 | 8.8 |
| red-noise power, SNR 6 | 76 | 18 | 15 | 8.7 |

Measured honestly, the transits still missed on pulsators are weak: scaling
the significance by the red noise cuts false alarms by a third and finds no
more of them. Lowering the cut to SNR 6 to compensate is a wash inside the
counting noise, at 3.6x the run time, and would loosen the cut on ordinary
stars (beta ~ 1). So the option stays off; it is a false-alarm lever if one is
needed, not a completeness one. The ~40% of strong transits missed on
pulsators are, at their real red-noise SNR, mostly below any threshold that
keeps the false-alarm rate.

## Shape tests: planet or eclipsing binary? (2026-11-06)

`python -m tesstransit.vet` recomputes three tests at our own signal's
ephemeris on the flattened light curve, for every TOI the sectors 1-13
planet runs (`s{N}-planets-pwB`) recovered whose TFOPWG disposition is a
planet (CP, KP) or an eclipsing binary (EB): 672 signals on 639 stars.
Uncertainties come from the flux binned to the transit duration, so red
noise is included.

- **odd/even**: mean depth of odd and of even transits; an EB of two similar
  stars, found at half its period, alternates.
- **secondary**: the dip at phase 0.5, same duration, with the primary
  masked first (left in, it pulled the baseline down and every hot Jupiter
  read as a negative secondary).
- **V-shape**: full-box depth over central-half depth, ~1 for a flat bottom,
  2/3 for a pure V.

Significance alone over-flags planets: at 3 sigma, 45 of 403 planets failed
something, because on a deep hot Jupiter a 2% odd/even difference or a real
thermal secondary is already significant. Each test therefore also needs a
SIZE. A grid over the size cuts, capped at 3% of planets flagged, picked
the same thresholds tuned on odd sectors (tested on even) and on even
(tested on odd): odd/even > 3 sigma AND > 30% of the depth; secondary > 3
sigma AND > 5% of the primary; shape ratio < 0.7.

| held-out sectors | planets flagged | EBs flagged |
|---|---|---|
| even (tuned on odd) | 4/185 | 36/123 |
| odd (tuned on even) | 6/218 | 40/142 |
| all sectors | 10/403 (2.5%) | 76/265 (29%) |

(Four stars first failed with no data at phase 0.5; with that handled, 405
planets and 267 EBs, the same 10 and 76 flagged.) The secondary eclipse
does most of the work (58 of the 76). Most EB TOIs
pass: their median depth is 3,200 ppm against 6,200 ppm for the planets,
the mark of a binary diluted by a brighter neighbour, whose dip has a
planet's shape. Telling those apart needs the target pixels (a centroid
shift), which this search does not read.

Applied to what nobody has labelled: 19 of 385 recovered PC TOIs (4.9%) and
18 of 105 signals on no-TCE pilot stars (17%) fail a test.

## Planets SPOC's single-sector search did not flag (2026-11-06)

The planet runs search every TOI on a star with a 2-minute light curve,
including planet-like TOIs with no SPOC TCE in that sector: 1,896 such
TOI-sectors over sectors 1-13. Most are out of reach in one sector (830 too
weak at their known ephemeris, 261 with too few transits, 631 outside the
period grid, 53 single transits). Of those within reach:

| found by | found | missed | |
|---|---|---|---|
| QLP only (MIT's FFI pipeline) | 8 | 1 | 89% |
| SPOC, other searches (other sectors, multi-sector) | 98 | 10 | 91% |

So on planets SPOC's same-sector search did not report, the search does as
well as on those it did. 200 of the 209 QLP-only TOI-sectors are out of
reach in one sector, which is what stitching sectors addresses.

## Stitched sectors 1-13: a pilot (2026-11-06)

`python -m tesstransit.multisector` flattens each sector with the
single-sector spline (pre-whitening included), stitches them, bins to 10
minutes and searches 1 d to half the baseline, about 355,000 periods per
star against 24,000 for a sector. Sample: 6,911 stars have 2-min data in at
least nine of the thirteen sectors. Of those, 150 were drawn from stars with
a TCE in SPOC's own sectors 1-13 search (265 TCEs) and 50 from stars with
none. Each search took ~250 s per star with 12 in parallel. The same 200
stars were run again inverted.

Of the 265 TCEs, 62 were found. Most of the rest were out of reach by
design:

| TCE period | found | why the rest were missed |
|---|---|---|
| under 1 d | 22 of 110 (at a multiple) | 88 below the grid |
| 1-13 d | 29 of 54 | 14 not detected, 11 search exhausted |
| 13-26 d | 6 of 6 | |
| 26-180 d | 5 of 56 | 38 not detected, 12 search exhausted |
| over 180 d | 0 of 39 | above the grid (half the baseline) |

All 9 TOIs in the sample were found: the TOI 700 system (9.98, 16.05 and
37.43 d), 713.01 (36.0 d) and .02, 187.01, 867.01 (15.4 d) and 705.01
(47.0 d). Most SPOC multi-sector TCEs at 26-180 d are near its own
threshold (SNR 7-10) and not TOIs.

⚠️ **False alarms at the single-sector thresholds (SNR 7.1, SDE 11).** On
the 50 stars with no TCE, the search put a signal on 21 of them, and on 22
of them when inverted. In both runs these are mostly two-transit signals at
40-170 d. With 15 times the trial periods, two coincident dips in a year of
data clear SDE 11 easily. Cuts measured on the same runs (signals from the
full run, kept if they pass):

| cut | stars without a TCE flagged, normal | inverted | TCE matches kept |
|---|---|---|---|
| none (SNR 7.1, SDE 11, 2 transits) | 21 / 50 | 22 / 50 | 62 / 62 |
| 3 transits | 13 / 50 | 8 / 50 | 61 / 62 |
| 3 transits, SNR 10 | 4 / 50 | 3 / 50 | 55 / 62 |
| 3 transits, SDE 15 | 5 / 50 | 3 / 50 | 60 / 62 |
| 3 transits, SDE 15, SNR 12 | 2 / 50 | 0 / 50 | 49 / 62 |

Three transits and SDE 15 gives 6% inverted (Wilson 2-16%), against 13.5%
for a single sector. That was chosen on this pilot, and the search itself
was run at SDE 11, so a raised threshold could change which later signals
are tried. It needs a fresh sample run at those settings before it is a
default; the code still searches at SDE 11. On the TCE stars the inverted
run flags 53-130 of 200 under any of these cuts. Those are mostly
short-period signals on variable stars and binaries, where an inverted
light curve still has periodic dips, so they say nothing about quiet stars.

### Hold-out at 3 transits, SDE 15 (2026-11-07): the cut does not hold

A fresh 200 stars (seed 2, `--holdout-of runs/ms-pilot`, no star shared
with the pilot; 150 with a TCE, 305 TCEs, and 50 without), this time
searched AT the cut (`--set sde_threshold=15 --set min_transits=3`), once
normally and once inverted. Runs `ms-holdout` and `ms-holdout-inv`, ~48 min
each with 12 workers (light curves cached in `lc_ms_ho`, 4.1 GB).

| | pilot (cut applied afterwards) | hold-out (searched at the cut) |
|---|---|---|
| stars without a TCE flagged, normal | 5 / 50 | 4 / 50 |
| stars without a TCE flagged, inverted | 3 / 50 (Wilson 2-16%) | **8 / 50 (Wilson 8-28%)** |
| planet-like TOIs found | 9 / 9 | 5 / 6 |
| TCEs found, 26-180 d | 5 / 56 (at SDE 11) | 1 / 71 |

The inverted rate is 16%, no better than a single sector's 14%. Six of the
eight inverted signals are still long-period, 27-119 d, so a third transit
does not stop coincident dips once there are 355,000 trial periods. The
other two are at 1.01-1.02 d. The 3 / 50 on the pilot was a cut chosen on
the data it was scored on.

The planets are fine: TOI 813.01 at 84 d was found, as were 210.01, 789.01
and .02. The miss is 789.03 (8.04 d), the third planet of that star, which
went to search exhaustion. But at SDE 15 only 1 of 71 TCEs at 26-180 d was
found, against 5 of 56 at SDE 11, because SPOC's own long-period TCEs sit
near its threshold.

**Not adopted.** The stitched search keeps the single-sector defaults, and
its long-period signals stay unvetted: about one quiet star in six gets a
false one. A threshold alone does not fix this. Two candidates are an
individual-transit test (each event should be significant on its own,
which a pair of coincident noise dips is not) and calibrating the SDE
against the trial count per period range instead of using one value for
the whole grid.

## Transit-timing variations: TOI 216 (2026-11-06)

`python -m tesstransit.ttv` times every transit on the stitched sectors
1-13 light curve (a trapezoid slid across the predicted time, the error
from chi-square + 1) and fits one fixed period through the times. TOI 216 is
a pair near 2:1 resonance that the single-sector search missed because its
transit times wander.

Two traps, both measured on it, are now handled in the code. The
catalogue ephemeris does not carry back: TOI 216.02's, 58 orbits from
sectors 1-13, puts every trial time where there is no transit (folded depth
~0 ppm). So the ephemeris is first refined by a narrow BLS search (±3% in
period) on this data. And refined alone, 216.02 locked onto 216.01 at half
its period, so on a multi-planet star the deepest is timed first and its
transits masked.

| | period (this data) | transits timed | O − C range | timing error | chi² / dof |
|---|---|---|---|---|---|
| TOI 216.01 | 34.554 d | 8 | −23 to +22 min | ~1 min | 362 |
| TOI 216.02 | 17.073 d | 15 | −90 to +56 min | 3-8 min | 100 |

Both curves are smooth over the year and opposite in phase, 216.02 early
while 216.01 is late: the exchange of energy between two resonant planets,
the lighter inner one swinging more. No fixed period fits either, which is
why each looked weak to a fixed-period search in single sectors.

### Null test: chi² / dof does not single TOI 216 out

The same timing was run on every signal of the stitched pilot (above) with
at least four transits: 201 signals, most of them not TTV systems. Median
chi² / dof is 33, the 90th percentile 87, and 16 of 201 exceed 100. TOI
216.02's 100 is ordinary among them; only 216.01's 362 sits in the tail, and
deep eclipsing binaries go higher still (11,500 and 2,400). Two causes,
split by the signal's SNR per transit (its SNR over √transits):

| SNR per transit | signals | median chi² / dof | median max \|O − C\| |
|---|---|---|---|
| under 2 | 140 | 33 | 197 min |
| 2-5 | 37 | 39 | 161 min |
| 5-10 | 8 | 62 | 51 min |
| 10-30 | 5 | 10 | 20 min |
| over 30 | 11 | 24 | 9 min |

- Below ~5 a single transit is not detected on its own. The "time" is the
  lowest noise point within ±1 duration (O − C near the window's limit),
  and chi² + 1 on a noisy chi² curve understates the error.
- Above ~10 the times are real, but chi² / dof still sits at 4-25. The
  confirmed single planets here (TOI 700.01 and .02, 705.01, 713.01) give
  3-7, so the errors are about 2× too small even when timing works. The
  likely causes are 10-minute binning, the fixed 15% ingress and red noise.

So chi² / dof is not a detection statistic as it stands. TOI 216 was
recognised by its O − C curves: smooth over a year and opposite in phase,
which noise does not produce. Before timing could screen a sample, it would
need errors calibrated on these planets (or bootstrapped), a per-transit SNR
floor of ~10, and a test for a smooth O − C rather than a large chi².

### The two-event veto (2026-11-07): adopted for the stitched search

The hold-out's eight inverted false alarms were nearly all two strong dips
plus a weak third event. The existing veto drops ONE event (the strongest)
and asks the rest to show the dip at SNR 3, which the second dip alone
satisfies. `min_snr_without_two_strongest` drops the TWO strongest. A real
transit is in every event, so the rest still show it, and a dip pair does
not. A vetoed peak has both events masked and the search goes on.

Scored afterwards on the stored signals (each event's SNR recomputed from the
cached light curves, the star's other signals masked), stars without a TCE
with an inverted signal, of 50:

| | without | with the two-event test at SNR 3 |
|---|---|---|
| pilot, SDE 11 | 22 | 7 |
| pilot, SDE 13 | 7 | 2 |
| hold-out, SDE 15 (searched at 3 transits) | 8 | 1 |

Planet TOIs kept: 9 of 9 and 4 of 4. TCE matches kept: 48 of 50 (one is a
two-transit TCE, which the test rejects by design) and 63 of 63. A ratio
form (the weaker events' mean SNR over all events') did not separate: many
of the pilot's inverted signals with 3+ events are as even as a planet.

Both samples had now been looked at, so the decision went to a third.
**Fixed before running:** SDE 13 with the two-event test at SNR 3, on 200
fresh stars (seed 3, `--holdout-of` both earlier runs); adopt for the
stitched search if at most 4 of 50 stars without a TCE give an inverted
signal and no planet TOI is lost except by design. Runs `ms-v3`,
`ms-v3-inv` (~60 min each; cache `lc_ms_v3`).

| | result |
|---|---|
| stars without a TCE, inverted signal | **3 / 50** (6%, Wilson 2-16%) |
| stars without a TCE, normal signal | 5 / 50 |
| planet-like TOIs in the grid | 9 / 9, incl. 1209.01 at 40.7 d |
| TCE matches (incl. aliases) | 66 / 289 (hold-out 68 / 305 at SDE 15) |

The missed planet TOI is 199.01 at 314.6 d, above the grid (half the
baseline). Of the three inverted false alarms, two are 1.04 and 1.15 d
signals that appear in both directions (variability, which the single-sector
search also reports), and one is a 37.9 d, 5,000 ppm event at SNR 30. The
long-period dip pairs are gone. **Adopted** as `multisector.STITCHED_DEFAULTS`
(the single-sector `SearchConfig` keeps it off). 50 stars cannot pin the rate
down; its 95% range still reaches the single-sector 14%.

A targeted check, not part of any sample: TOI 214 (TIC 167415965, 12
sectors), whose two planets were in the archive's November 2026 additions,
was missed in 11 of 12 single sectors. Stitched, both are found at either
setting: 18.5523 d at SNR 28.5 and 9.6958 d at SNR 11.1, with the rest at
SNR 19.0 and 7.6 after the two strongest transits are dropped.

## The false-alarm rate on 250 fresh stars (2026-11-07): 10%, not 6%

The adopted default (SDE 13, two-event test at SNR 3) run unchanged on 250
stars without a TCE, none in an earlier sample (seed 4, `--n-tce 0`,
`--holdout-of` all three earlier runs), searched inverted and as observed.
Runs `ms-fa250-inv`, `ms-fa250` (~55 min each; cache `lc_ms_fa`).

| | stars with a signal, of 250 | Wilson 95% |
|---|---|---|
| inverted | **26** (10.4%) | 7.2-14.8% |
| inverted, pooled with `ms-v3` (29 / 300) | 9.7% | 6.8-13.5% |
| as observed | 17 (6.8%) | 4.3-10.6% |

So the 3 / 50 of `ms-v3` was the low end of its range. The stitched search's
rate is lower than a single sector's 14% (10.8-16.7%) but the two ranges
overlap. The false alarms are two different populations:

| period | inverted, stars | as observed, stars |
|---|---|---|
| 1.0-2.5 d | 14 | 14 |
| 4.7 d | 1 | 0 |
| 33-150 d | 12 | 3 |

**Short periods: stellar variability.** Most of these are the same stars at the
same period in both directions (e.g. TIC 309657720 at 1.286 d, 150437899 at
1.319 d, 220411074 at 1.60 d, 142018697 at 1.85 d, 150391740 at 1.87 d), with
depths of 30-1,100 ppm and 3-4 h boxes. A signal found in both directions
is symmetric, and a transit is not. These stars are periodic, either just above
the 1.5 d pre-whitening limit or below it but not removed by its gates.

**Long periods: data edges.** The long-period signals pass the two-event
test only barely (11 of 14 inverted have SNR 3.1-4.8 with the two strongest
transits left out). Their events sit at the start or end of a data segment
(a sector start or the mid-sector downlink gap): **17 of 59 events lie within
0.5 d of an edge**, where only 9.7% of the data is, and 12 of the 17
long-period signals have at least one such event. Several share a time
across unrelated stars (four at BTJD 1491.6-1491.8, the start of sector 7).
The ramp at a segment edge, which the spline does not follow, is mostly a
brightening: inverted, it becomes a dip, which explains 12 stars inverted
against 3 as observed.

Not yet tested, each to be fixed in advance and scored on a fresh sample:
(1) leave out the first and last ~0.5 d of each segment, or veto an event
that lies within 0.5 d of an edge; (2) a both-directions test, a signal whose
period the inverted search also finds is variability, not a transit.
Scored afterwards on these stars, only as a guide: the edge rule would leave
3 of the 12 inverted long-period stars, and inverted 26 → about 17.

### Edge trim and both-directions test: checked afterwards, then fixed in advance

Two options in `SearchConfig`, both off by default: `edge_trim_days` (drop
data within that much of a segment edge) and `opposite_snr_ratio` (the
both-directions test; see the field's comment). Re-searched afterwards on
stars already seen (the 26 + 17 flagged stars without a TCE from the
250-star run, and every star where `ms-pilot`, `ms-holdout` or `ms-v3`
recovered a planet TOI or TCE):

| | inverted, of 26 | as observed, of 17 |
|---|---|---|
| edge trim 0.5 d | 14 | 13 |
| both-directions, ratio 0.5 | 16 | 9 |
| both | 9 | 7 |

What this check found:
- **The both-directions test needs many events.** TOI 813.01 (3 events at
  83.8 d) was vetoed: the brightening is the best of ~1,000 trial phases, so
  noise reached SNR 8 against the dip's 9.7. Now only signals with at least
  10 events (`opposite_min_events`) are tested; all the variability false
  alarms were under 5 d.
- **A ratio of 0.5 is too low.** TOI 787.01 (2.13 d) had a brightening at
  0.52 of its dip. Scored afterwards: inverted stars left 9/11/13 of 26 at
  ratios 0.6/0.7/0.8, and known signals kept 54/65/66 of 77. Most of the
  known signals lost are SPOC TCEs at 0.3-1.7 d with no TOI, the same
  variability. 0.7 was chosen for its margin above 787.01.
- **A bug in peak selection, exposed by the trim.** On TOI 802.01, one strong
  single dip made a peak at every period that put its other transits in
  gaps, and these 500 neighbouring grid points used up `max_peaks_tried`
  before the planet's peak was reached ("no_multi_transit_peak"). The budget
  now counts distinct failed epochs, so 802.01 is found again (SNR 10.1). This
  changes the default search as well.

**Fixed before running:** `--set edge_trim_days=0.5 --set
opposite_snr_ratio=0.7` (`opposite_min_events` 10, `opposite_min_repeat_snr` 3)
on the stitched defaults, on 150 TCE stars + 250 stars without a TCE, none in
an earlier sample (seed 5, `--holdout-of` all four earlier runs), inverted and
as observed. Runs `ms-v4-inv`, `ms-v4`; cache `lc_ms_v4`. Adopt as the
stitched default if (a) at most 15 of 250 stars without a TCE give an inverted
signal (6%, against 26 of 250 before), and (b) every planet TOI (CP/KP/PC) in
the grid with at least 3 transits is found.

**Result (2026-11-08): fails (a), not adopted.** Runs `ms-v4-inv`, `ms-v4`
(400 stars; one failed, TIC 55497082, with negative PDCSAP flux in sector 13).

| | result | condition |
|---|---|---|
| stars without a TCE, inverted signal | **21 / 250** (8.4%, Wilson 5.6-12.5%) | at most 15: **fails** |
| stars without a TCE, normal signal | 23 / 250 | |
| planet TOIs in the grid, 3+ transits | 8 / 8 (216.01, 220.01, 723.01, 798.01 at 95.2 d, 794.01, 283.01, 708.01; 801.01 at its 2x alias) | all: passes |

Of the 21 inverted false alarms, 13 are 1.0-1.9 d variability again, most
of them on stars that also give a signal as observed. The both-directions
test misses them because their brightening is weaker than their dip: the
ratio is 0.43-0.70 (median 0.63), just below the 0.7 cut, and TOI 787.01's
0.52 lies inside that range. A box brightening at the dip's duration does not
separate this variability from planets. The other 8 are 24-162 d. Whether
the edge trim helped cannot be told from 8 against the 12 of a different
250 stars. So both options stay off. The peak-budget fix stays in, since it
is a bug fix.

The variability left is periodic at 1-2.5 d, just above the 1.5 d where each
sector's pre-whitening stops. Untested next step: pre-whiten to a longer
period for the stitched search, judged on planets lost as well as false
alarms removed.

## Short-period variability: why pre-whitening misses it, and a variability veto (2026-11-08)

**Longer per-sector pre-whitening does not help.** Re-searched on stars already
seen, with `prewhiten_max_period_days` 2.5 or 3.5 instead of 1.5:

| stars without a TCE | default | 2.5 d | 3.5 d |
|---|---|---|---|
| 250-star run, inverted | 26 | 24 | 24 |
| `ms-v4`, inverted | 21 | 19 | 20 |
| `ms-v4`, as observed | 23 | 17 | 18 |

(2.5 d kept all 30 planet TOIs and 212 of 216 recovered TCEs.) Most of the
surviving short-period signals, 32 of 77, are under 1.5 d, inside the range
pre-whitening already covers. On TIC 167247077 (inverted signal at 1.269 d),
10 of 12 flattened sectors still hold a 0.634 d wave with Lomb-Scargle
false-alarm probability 1e-8 to 1e-18, and pre-whitening removed nothing in
them: the wave is ~60 ppm, so the hour-binned scatter stays under twice white
noise and the red-excess gate stops it. That gate is there to protect
shallow planets (TOI 396). Stitched, the wave's false-alarm probability is
1e-59, and BLS fits a box at twice its period.

**Pre-whitening the stitched curve loses planets.** The same tests without the
red-excess gate, on the stitched curve up to 3 d (not kept): inverted false
alarms 26 → 11 and 21 → 8, but TOI 787.01 and 801.01 (PC) were lost, with the
eclipsing binaries TOI 284.01 and 2246.01, and 117 of 216 SPOC TCEs with no
TOI, nearly all at 1-2 d. What it removed on the planets were harmonics of the
transit period (787.01: P/4, P/3; 801.01: P/4, P/6). Folded at P/k a shallow
transit is diluted k times to the noise, and the one-sided test cannot see it.
Median time per star 210 s (was 120 s).

**A variability veto, measured afterwards.** For each signal, mask its
transits (± one duration) and find the strongest Lomb-Scargle peak within
0.5% of j/P, j = 1-4, on the stitched curve. A transit train's own power at
j/P goes with the masked points; variability stays. Scored on every recorded
signal of the false-alarm stars and of the recovered TCE stars of all four
earlier runs (`coh.py`; 311 star runs):

| veto when | false-alarm signals | planet + EB signals | other TCE signals |
|---|---|---|---|
| FAP < 1e-4 | 97 / 137 | 17 / 46 | 349 / 374 |
| FAP < 1e-10, two-sided ≥ 0.3 | 82 / 137 | 6 / 46 | 253 / 374 |
| … and semi-amplitude ≥ 0.5 × box depth | ~77 / 137 | 1 / 46 | |

The amplitude condition is what separates them: in a false alarm the box IS
the wave (ratio median ~0.95), while an eclipse is much deeper than any wave
on its star. The one planet signal lost is TOI 708.01 (PC, 1.996 d), whose
host varies at P/4 with 279 ppm against a 330 ppm transit: a planet host can
vary at a harmonic of the orbit, as 787.01's does at P/2. Per star, counting
a star flagged if any signal is not vetoed (an approximation, since a vetoed
signal is masked and the search continues): inverted 26 → ~16 (250-star run)
and 21 → ~12 (`ms-v4`).

Implemented as `coherent_veto_fap` (off by default) with
`coherent_harmonics` 4, `coherent_rel_window` 0.005, `coherent_min_symmetry`
0.3, `coherent_min_amp_ratio` 0.5; vetoed signals have reason `"coherent"`
and are masked like found ones. The wave is measured on every valid point,
not only those left unmasked: earlier vetoes at the same period mask the
wave's troughs, and without them the fold looks one-sided (a synthetic
sinusoid was vetoed three times, then accepted).

Checked on three stars with the implementation: TIC 167247077 inverted gives
no signal (1.269 d vetoed twice, FAP 1e-29, 1e-34); TOI 787.01 is kept (SNR
11.3); TOI 708.01 is lost as predicted, and its star then reports a 73.5 d
signal. There the veto fired six times (1.996 d, then 1.49 d five times),
exhausting `max_vetoes`: on a strongly varying star, vetoes at the
variability's period can use up the budget before a planet is reached. The
fresh sample measures whether that costs planets.

**Fixed before running:** `--set coherent_veto_fap=1e-10` on the stitched
defaults (edge trim and both-directions test off), on 150 TCE stars + 250
stars without a TCE, none in an earlier sample (seed 6, `--holdout-of` all
five earlier runs), inverted and as observed. Runs `ms-v5-inv`, `ms-v5`; cache
`lc_ms_v5`. Adopt as the stitched default if (a) at most 15 of 250 stars
without a TCE give an inverted signal (6%), and (b) every planet TOI
(CP/KP/PC) in the grid with at least 3 transits is found. If it fails, the
stitched false-alarm rate stays at ~8-10% and false-alarm work stops here.

**Result (2026-11-09): fails (a) and (b), not adopted.** Runs `ms-v5-inv`,
`ms-v5` (400 stars, none failed; median 145-158 s per star).

| | result | condition |
|---|---|---|
| stars without a TCE, inverted signal | **17 / 250** (6.8%, Wilson 4.3-10.6%) | at most 15: **fails** |
| stars without a TCE, normal signal | 20 / 250 | |
| planet TOIs in the grid, 3+ transits | **11 / 13** | all: **fails** |

The veto fired 200 times on the inverted run. Of the 17 inverted false alarms,
9 are still at 1.1-2.8 d, and 8 are 7.9-95 d; 5 of the 17 stars also give a
signal as observed. Pooled with the two earlier fresh 250-star samples, the
stitched inverted rate is 7-10% whatever is switched on, so false-alarm work
on the stitched search stops here: a stitched signal on a star with no TCE
is unvetted.

The two planets missed are **not** the veto's doing (no veto fired on either
star): TOI 6902.01 (PC, 36.4 d, 9.3 h, 5 transits in the data) and TOI 810.02
(PC, 90.8 d, 11.1 h, 3 transits). Both transits are longer than the longest
box tried (8 h), and the first reading of this result blamed that. It is
wrong; the next section has the check.

## Long transits: the 8 h box is not why TOI 6902.01 and 810.02 are missed (2026-11-09)

The longest box is tied to the detrending: `validate()` requires the spline
knot spacing to be at least 1.5x the longest duration, so the 0.5 d knots
allow at most 8 h. Longer boxes need coarser knots, which flatten stellar
variability less. Both stars were searched again (stitched, sectors 1-13) with
knots of 1-3 d, first with boxes up to 8 h and then up to 14 h:

| star | knots 1-3 d, boxes to 8 h | knots 1-3 d, boxes to 14 h |
|---|---|---|
| TOI 6902.01 as observed | no signal | no signal |
| TOI 6902.01 inverted | 1.671 d, SNR 24 (false alarm) | no signal |
| TOI 810.02 as observed | 28.31 d only (TOI 810.01) | 28.31 d only (TOI 810.01) |
| TOI 810.02 inverted | no signal | no signal |

Longer boxes find neither planet. Folding the stitched curve at the TOI
ephemerides shows why: neither miss has anything to do with the box length.

- **TOI 810.02** is in the data at full depth: the three transits measure
  365, 373 and 444 ppm (TOI: 314 ppm), each about SNR 5 against the binned
  point scatter, about 9.6 combined. The knot spacing makes no difference.
  The original search's best peak on this star was already at 90.82 d
  (SNR 7.2, SDE 10.4), below the stitched SDE 13 cut. It is a threshold miss:
  three transits at SNR 5 each do not clear the look-elsewhere cost of a
  13-sector period grid.
- **TOI 6902.01** is not cleanly in the data: the star is very bright
  (Tmag 5.7) and its five transits measure 31, 118, 14, 94 and 181 ppm with
  0.5-3 d knots, but -110, 113, 278, 173 and 81 ppm with 1-3 d knots
  (TOI: 384 ppm; point scatter 96 and 136 ppm). The trend left by detrending
  is as large as the transit, so no box length recovers it. SPOC found it
  in a much longer data set (s01-s69).

So longer boxes are not tested on a fresh sample: no planet in hand needs
them, and coarser knots let more stellar variability through. Whether the
8 h limit costs planets with much longer transits is not measured.

## Lowering the stitched SDE cut: checked on the earlier samples, not adopted (2026-11-09)

TOI 810.02 misses only the SDE 13 cut, so the question is what a lower cut
buys. Each stitched result records the peak that stopped the search
(`final_snr`, `final_sde`, `final_period_days`). A star with no signal whose
stopping peak has SNR >= 7.1 and SDE between the lower cut and 13 is an upper
bound on what that cut adds: the peak still has to pass the two-event veto,
and a cut lower than 13 could also add a second or third signal on a star.
The planet side is a bound too, since only the stopping peak is recorded.
Across the four samples run with the adopted stitched settings (`ms-fa250`,
`ms-v3`, `ms-v4`, `ms-v5`; 799 quiet stars, plus the TCE stars for TOIs):

| cut | quiet stars added, inverted | quiet stars added, as observed | planet TOIs (CP/KP/PC) gained |
|---|---|---|---|
| SDE 11 | 59, 10/50, 56, 63 per 250 | 54, 12/50, 43, 54 | 0 |
| SDE 10 | 96, 19/50, 89, 108 per 250 | 90, 20/50, 84, 101 | 1 (TOI 810.02, SDE 10.4) |
| SDE 9 | 128, 25/50, 121, 143 per 250 | 123, 25/50, 112, 137 | 1 |

The noise sits just below the cut: inverted, about a quarter of quiet stars
stop at SDE 11-13 against 7-10% that clear 13, and the as-observed counts are
the same size, so these are noise peaks, not hidden planets. One planet in
four samples is not worth even a fraction of that. SDE 13 stays, and no
fresh sample is run: the bound already decides it.

## Sectors 14-26: the second year, untuned (2026-11-09)

The second year (northern ecliptic hemisphere) was searched with the
single-sector default unchanged: nothing was tuned on it. Per sector, the
same three runs as for sectors 1-13: every TOI on a 2-minute star
(`s{N}-planets-pwC`), the 60 TCE + 40 quiet sample (`s{N}-pilot-pwC`), and
that sample upside down (`s{N}-pilot-inv-pwC`). SPOC TCE tables
`s0014`-`s0026`, 726-1,591 TCEs each; TOI table as of 2026-11-04.

Sectors 1-13 were re-run on the same code first, since the
`_best_with_transits` fix came after their published runs (`pwB`): 4 of
5,215 star searches changed, all on TCE stars (two that used to raise now
finish, two eclipsing binaries gain a signal), and every pooled number
below is identical. So the published sector 1-13 figures stand and the two
years are compared on one code version.

| | sectors 1-13 | sectors 14-26 |
|---|---|---|
| planet TOIs (CP/KP/PC) with a TCE in the sector, found | 682/707 = 96.5% (94.8-97.6) | 670/696 = 96.3% (94.6-97.4) |
| planet TOIs with no TCE in the sector, found | 108 | 165 |
| ... of those detectable in the sector (found + not detected + exhausted) | 108/120 | 165/184 |
| sample TCEs found (incl. aliases) | 536/1,125 | 627/1,136 |
| quiet stars with a signal | 66/520 = 12.7% | 61/520 = 11.7% |
| quiet stars with a signal, upside down | 70/520 = 13.5% (10.8-16.7) | 59/520 = 11.3% (8.9-14.4) |

The second year reproduces the first within the intervals. The larger
no-TCE count is volume: 2,382 planet TOIs list a year-2 sector with no TCE
from it, against 1,896; most are too weak in one sector either way
(1,515 and 1,091 below the forced-SNR or transit-count bar). The
false-alarm conclusion is unchanged: signals on no-TCE stars occur at the
rate the upside-down search produces, so they are not candidates.

Per sector, planets with a TCE missed: 1-5 (sector 20: 5 of 55; sector 21:
4 of 60). Causes are the year-1 ones: weak after detrending (14), not
detected (10), one search that ran out of signal slots, and one TOI too weak
even at its own ephemeris (forced SNR below the bar).

Two stars failed. TIC 232564158 (sector 16) has a negative median PDCSAP
flux and is rightly refused, as TIC 300038935 is in sector 4. TIC 259006185
(sector 26, a 2.5% eclipsing binary with TCEs) failed only upside down: the
inverted eclipses drew 2 h boxes at 0.277 d, a seventh of the 1.94 d orbit,
and after two accepted signals only 4 cadences lay outside the next peak's
transits, so `event_stats` raised. The search now stops with
`too_few_out_of_transit` (fwh_tess c63b205); the sector 26 inverted run was
repeated, and its quiet-star count is unchanged (10/40).

## Sectors 27-39: the third year, untuned (2026-11-10)

The southern hemisphere again, two years after sectors 1-13, searched with
the single-sector default unchanged (`s{N}-{planets,pilot,pilot-inv}-pwC`).
SPOC TCE tables `s0027`-`s0039`, 909-2,287 TCEs each; TOI table as of
2026-11-04.

| | sectors 1-13 | sectors 14-26 | sectors 27-39 |
|---|---|---|---|
| planet TOIs (CP/KP/PC) with a TCE in the sector, found | 682/707 = 96.5% | 670/696 = 96.3% | 1,341/1,383 = 97.0% (95.9-97.7) |
| planet TOIs with no TCE in the sector, found | 108 | 165 | 203 |
| ... of those detectable in the sector | 108/120 | 165/184 | 203/236 |
| sample TCEs found (incl. aliases) | 536/1,125 | 627/1,136 | 606/1,091 |
| quiet stars with a signal | 66/520 = 12.7% | 61/520 = 11.7% | 81/520 = 15.6% |
| quiet stars with a signal, upside down | 70/520 = 13.5% (10.8-16.7) | 59/520 = 11.3% (8.9-14.4) | 93/520 = 17.9% (14.8-21.4) |

Planet recovery holds; the year has twice the planet TOIs per sector of the
first two (the TOI catalogue grew). The 42 misses: weak after detrending
(16), not detected (14), outside the search grid (12, about one a sector).
Two stars are refused for a negative median PDCSAP flux (TIC 264979636 in
sector 34, TIC 447061717 in sector 37).

⚠️ **The upside-down rate rose**, 17.9% against 11-13%, outside both
earlier intervals; per sector 4-12 of 40 (highest: sectors 32, 34, 38).
Observed and inverted still agree (15.6% vs 17.9%), so the conclusion that
no-TCE signals are false alarms as a group stands, but the cause is not yet
known. Unchecked candidates: year-3 2-minute targets chosen more often for
variability, or a processing change. Not investigated before travel.

## Sectors 1-13 + 27-39 stitched: a paired run (2026-11-10)

The same stars as the published stitched samples (`ms-v3`, `ms-fa250`),
re-searched with both southern years, so the comparison is paired, not two
draws. New option `multisector --sample-of RUN` (repeatable) searches an
earlier run's stars, re-grouped against the NEW TCE table (a star quiet in
sectors 1-13 can have a sectors 1-39 TCE) and counts the stars no longer
eligible instead of replacing them.

- Reference: SPOC's sectors 1-39 multi-sector TCEs
  (`tess2018206190142-s0001-s0039_dvr-tcestats.csv`, 9,699 TCEs). It also
  searched sectors 14-26, which these southern stars barely touch.
- Eligibility: 2-minute data in at least 18 of the 26 sectors: 3,311 stars.
  Of the 6,911 stars with 9+ sectors in year 1, 3,593 have **no** 2-minute
  data in year 3 (targets are re-chosen each sector). 216 of the 450 sampled
  stars remain: 81 with a TCE, 135 without.
- Cost: baseline ~1,064 d, ~81,500 binned points, ~3.7x the periods. One
  star took 28.7 min against 2.3 min for sectors 1-13 (12x); the run of 216
  took 12.2 h on 12 workers.
- Runs: `ms26-pair` (observed), then `ms26-pair-inv` (the 135 quiet-sample
  stars upside down; not yet run). The driver skips a run whose
  `summary.json` exists, so re-running it resumes at the inverted run.
  ⚠️ A run writes its results only at the end; an interrupted run starts
  over.

Results (`ms26-pair`, observed; same stars, earlier result from `ms-v3` /
`ms-fa250`):

| same stars | sectors 1-13 | sectors 1-13 + 27-39 |
|---|---|---|
| planet TOIs (CP/KP/PC) found | 8/9 | 9/9 |
| quiet stars (no s1-39 TCE) with a signal | 5/135 = 3.7% | **26/135 = 19.3%** |

- The added planet is TOI 199.01 (SPOC period 314.6 d), found at a third of
  its period (104.9 d, an alias). At 314.6 d it cannot be found in sectors
  1-13, whose search stops at half the ~350 d baseline.
- ⚠️ False alarms rose five-fold on the same stars. The search tries ~3.7x
  the periods of the 1-13 search, and SDE 13 was set on the 1-13 search, so
  the threshold does not carry over to a longer baseline. The inverted run
  (`ms26-pair-inv`, 135 quiet stars) has not been run; it is the test of
  whether these 26 are noise. Until then the 26-sector stitch is NOT
  adopted, and nothing from it goes on the site. Next: the inverted run,
  then an SDE cut set on it and checked on a fresh sample.
