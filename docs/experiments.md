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
