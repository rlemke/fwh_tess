# Definitions

What every term, field and outcome in fwh_tess means. Two vocabularies meet
here, and keeping them apart is the point of this document:

- **Ours** — the search (`SearchConfig`, `Signal`, `SearchResult`) and the
  comparison (`bucket`, `cause`, `alias`). Defined by this repo's code; the
  code is the authority and this file must change with it.
- **SPOC's** — the `tce_*` columns of MAST's per-sector TCE table. Defined by
  the SPOC Data Validation pipeline; this file records what each column means
  and **where that definition comes from**.

Where a definition could not be traced to a source, it says so.

---

## 1. Core terms

| Term | Meaning |
|---|---|
| **TESS** | NASA's Transiting Exoplanet Survey Satellite. Observes the sky in ~27-day **sectors**. |
| **Sector** | One ~27-day pointing of TESS's four cameras. Our unit of comparison: one TIC ID × one sector. |
| **TIC ID** | Identifier of a star in the TESS Input Catalog. |
| **2-minute target** | A star whose pixels were downloaded at 2-minute cadence and processed by SPOC into a light curve. Sector 1: 15,889 of them. |
| **Cadence** | One exposure (one row of a light curve), 2 minutes here. |
| **SPOC** | Science Processing Operations Center, NASA Ames — the pipeline that produces the light curves AND the TCEs we compare against. |
| **PDCSAP flux** | Pre-search Data Conditioning Simple Aperture Photometry: SPOC's light curve with instrumental systematics removed. The flux we search. |
| **QUALITY** | Per-cadence bit flags set by SPOC (momentum dumps, scattered light, …). We drop every flagged cadence by default. |
| **BTJD** | Barycentric TESS Julian Date = BJD − 2457000.0. The clock of light-curve `TIME` and of `tce_time0bt`. |
| **TPS** | Transiting Planet Search — SPOC's search: a wavelet-based adaptive matched filter. **Not** what we use. |
| **DV** | Data Validation — the SPOC stage that fits and tests every TCE; the source of the TCE table. |
| **TCE** | Threshold Crossing Event: a periodic signal in one star's light curve that SPOC's TPS search found above its detection threshold (7.1σ MES) and passed to DV. **A TCE is a detection, not a planet** — many are eclipsing binaries or systematics. |
| **TOI** | TESS Object of Interest: a TCE (or a QLP detection) that human vetters promoted to a candidate. Context only here, never the reference. |
| **MES** | Multiple Event Statistic — TPS's detection statistic for a folded, multi-transit signal. Not computed by us; comparable to our SNR in purpose only. |
| **BLS** | Box Least Squares — fits a flat-bottomed box dip at each trial period, epoch and duration. **Our** search. |
| **EB** | Eclipsing binary — two stars eclipsing each other. Deep dips (often > 1%), frequently a primary + secondary eclipse, which makes period aliases common. |
| **Alias** | A signal at a simple multiple (2×, 3×) or fraction (½, ⅓) of the true period. |
| **Duty cycle** | Duration ÷ period: the fraction of the time a star is "in transit". |

---

## 2. Our search

### 2.1 `SearchConfig` — every parameter that shapes a result

Saved with every `SearchResult` (`config`), so any reported signal traces to
the exact search that produced it.

| Field | Default | Meaning |
|---|---|---|
| `detrend_window_days` | 1.0 | Width of the running median divided out of the flux. Anything slower is removed — **including a transit longer than about a third of it**, which is why it must be ≥ 3× the longest duration searched. |
| `min_period_days` / `max_period_days` | 0.25 / 26.0 | The period grid. One sector spans ~27 d, so 26 d is the longest period that can show two transits. A TCE outside this range is `outside_search_grid`. |
| `n_periods` | 24000 | Trial periods, log-spaced across the grid. |
| `min_transits` | 2 | A peak is only accepted if this many transits have data. One transit fits every long period that puts the next in a gap, so it is not a periodic detection. |
| `max_peaks_tried` | 500 | Peaks checked for `min_transits`, strongest first, before giving up (`no_multi_transit_peak`). |
| `power_detrend_frac` | 0.042 | Width of the running median subtracted from the periodogram before `sde` is measured, as a fraction of the period grid. |
| `durations_hours` | 0.5, 1, 2, 3, 4, 6, 8 | Trial transit durations (box widths). |
| `min_stellar_density_solar` | 0.1 | Caps the duration tried at each period at the longest a transit could physically last around a star at least this dense (circular orbit, central transit): **T_max ≈ 13 h × (P / 1 yr)^⅓ × (ρ / ρ☉)^−⅓**. At 0.1 ρ☉: 3.9 h at 1 d, 9.2 h at 13 d. Admits subgiants, excludes giants. |
| `snr_threshold` | 7.1 | Minimum `snr` to report a signal. Numerically SPOC's MES threshold; **not the same statistic**. |
| `sde_threshold` | 11.0 | Minimum `sde` to report a signal. ⚠️ Provisional (see §2.3). |
| `max_signals` | 3 | Stop after this many signals per star. |
| `mask_width_factor` | 1.5 | After a signal is found, mask this many durations around each of its transits before searching again. |

### 2.2 `Signal` — one thing our search found

| Field | Unit | Meaning |
|---|---|---|
| `rank` | — | Order found in: 1 = the strongest, found first. |
| `period_days` | days | Best-fit period. |
| `epoch_btjd` | BTJD | Mid-time of one transit. |
| `duration_hours` | **hours** | Best-fit box width (one of `durations_hours`). |
| `depth_ppm` | ppm | Box depth. Expected **shallower** than SPOC's limb-darkened fit (π Men c: 259 vs 302 ppm) — a box is the wrong shape for a real transit. |
| `snr` | — | Depth ÷ its uncertainty, with the uncertainty from the measured scatter. |
| `sde` | — | Signal Detection Efficiency, on the **detrended** periodogram: (power − running median of power) ÷ robust std. How far the peak stands above nearby trials. |
| `n_transits` | — | Transits that have in-transit data (gaps can hide some). |

### 2.3 Why both `snr` and `sde` — and why `sde` is detrended

On synthetic white noise the strongest of ~10⁷ trial boxes reached SNR 6.4–7.8
(12 sectors), so an SNR cut at 7.1 alone accepted about half of noise-only
stars.

Raw BLS power rises with period: one deep dip fits every long period that puts
the next transit in a gap, so the long-period end is a high plateau. A plain
SDE measured against that plateau scored a real synthetic 2-transit 18 d
planet (SNR 49) at 6.5. On the detrended periodogram: noise 6.1–10.0, injected
300–1500 ppm transits at 3.7–23 d 41.8–159.6.

⚠️ On REAL light curves the separation does not hold (pilot 3, 2026-11-04):
no-TCE stars produced signals with SDE 11–76, while missed TCEs had their best
peak at SNR 10–75 but SDE 5.6–10.9. Red noise and systematics, not white noise,
set the real false-alarm floor. Both thresholds are open.

### 2.4 `SearchResult`

| Field | Meaning |
|---|---|
| `n_points` | Cadences searched (after the quality mask). |
| `scatter_ppm` | Robust scatter (1.4826 × MAD) of the detrended flux. |
| `signals` | The `Signal`s, in `rank` order. |
| `stop_reason` | Why the search stopped — see below. |
| `final_snr`, `final_sde` | The strongest peak that did **not** clear the thresholds. How close the search came to finding something more. |

| `stop_reason` | Meaning |
|---|---|
| `below_threshold` | The next-strongest peak failed `snr_threshold` or `sde_threshold`. The normal ending. |
| `max_signals` | `max_signals` reached; there may be more to find. |
| `too_few_points` | Masking left fewer than 100 cadences. |
| `no_multi_transit_peak` | None of the `max_peaks_tried` strongest peaks had `min_transits` transits with data. |

### 2.5 Light-curve mask counts

| Field | Meaning |
|---|---|
| `n_cadences` | Rows in the file. |
| `n_flagged` | Dropped because a QUALITY bit was set. |
| `n_nonfinite` | Dropped for NaN time/flux/error, among cadences NOT already flagged (each cadence counted once). |
| `n_kept` | Searched. `n_cadences = n_flagged + n_nonfinite + n_kept`. |

---

## 3. The comparison

### 3.1 When a signal matches a TCE

All of:

1. **Same star** (TIC ID).
2. **Period**: `|P_ours / (k · P_TCE) − 1| ≤ period_rel_tol` (0.5%), for some
   `k` in the alias set {1, 2, ½, 3, ⅓}.
3. **Epoch**: the offset between the two epochs, wrapped by the **shorter** of
   the two periods, is within `epoch_tol_durations` (0.5) × the longer of the
   two durations. Wrapping by the shorter period is what makes aliases work:
   at 2× every one of our transits is a TCE transit; at ½× every TCE transit
   is one of ours.

Assignment, per star: **direct matches (k = 1) are assigned across all TCEs
first, aliases only after**, each to the strongest free signal. Each signal
matches at most one TCE. (A single pass in TCE order let one TCE claim, as an
alias, a signal that was a direct match for a later one — measured on TIC
350479101.)

### 3.2 Buckets — every TCE and every signal lands in exactly one

| `bucket` | Row has | Meaning |
|---|---|---|
| `recovered` | TCE + signal | A signal at the TCE's period (k = 1) and epoch. |
| `recovered_alias` | TCE + signal | A signal at 2×, ½×, 3× or ⅓× the TCE's period with coincident transits. `alias` = k. Mostly EBs. |
| `missed` | TCE only | No signal matched. `cause` says what we know. |
| `unmatched` | signal only | Our signal with no TCE. **Unvetted** — something SPOC rejected or never triggered on, or a false alarm of ours. Never a "new planet". |

### 3.3 `alias`

Our period ÷ the TCE's, for a match: `1.0` (direct), `2.0` (we found twice the
TCE's period), `0.5` (half), `3.0`, `0.333…`.

### 3.4 Miss causes — limited to what the search can establish

| `cause` | Meaning |
|---|---|
| `not_searched` | The star has no search result (not in the run, or it failed). In a full-sector run this keeps an incomplete run from overstating recovery. |
| `outside_search_grid` | The TCE's period is outside `[min_period_days, max_period_days]`. We could not have found it at that period; not a search failure. |
| `search_exhausted` | The search stopped at `max_signals` before reaching it — often an EB whose harmonics used up the signal budget. |
| `not_detected` | Searched, inside the grid, nothing matched. **Does not say why** — below threshold, removed by detrending, masked by quality flags, or a period/epoch outside the match tolerances. Telling those apart needs injection tests, not a label. |

### 3.5 Pilot sample groups

| `group` | Meaning |
|---|---|
| `tce` | Stars with ≥ 1 SPOC TCE, chosen evenly across SPOC's SNR range so recovery is measured from the faintest TCEs to the strongest. |
| `quiet` | A seeded random draw of observed stars with **no** TCE. Any signal here is something SPOC rejected or our false alarm — the group that calibrates the thresholds. |

---

## 4. SPOC's TCE table — the `tce_*` columns

One row per TCE, 123 columns, extracted by SPOC from its Data Validation XML.
Download: `https://archive.stsci.edu/missions/tess/catalogs/tce/<stamp>-s00NN-s00NN_dvr-tcestats.csv`
(6 leading `#` lines). Sector 1: 1,457 TCEs on 971 stars.

**Sources.** [M] = MAST's own column file,
[`TCE_csv_description.txt`](https://archive.stsci.edu/missions/tess/catalogs/tce/TCE_csv_description.txt),
which maps every column to a label, unit and DV XML element — the primary
source. Its labels are copied from Kepler with TESS notes in brackets
("[BTJD FOR TESS]"); **the bracket is the TESS meaning**. [K] = the NASA
Exoplanet Archive's
[Kepler TCE column definitions](https://exoplanetarchive.ipac.caltech.edu/docs/API_tce_columns.html),
whose names the TESS file reuses. [T18] = Twicken et al. 2018, *Kepler Data
Validation I* ([arXiv:1803.04526](https://arxiv.org/abs/1803.04526)).
"Measured" = checked against the sector-1 rows. `±` = a `_err` column with the
same units exists.

### 4.1 The columns this repo reads

| Column | Meaning | Unit | Used as |
|---|---|---|---|
| `tceid` | TIC ID + `-NN` planet-number suffix | — | `TCE.tce_id` |
| `ticid` | TIC ID | — | `TCE.tic_id` |
| `tce_plnt_num` | Which TCE on this star (1 = found first) | — | `TCE.planet_num` |
| `tce_period` ± | Period from the all-transits fit | days | matching |
| `tce_time0bt` ± | Epoch: centre of the first transit | **BTJD** | matching (same clock as light-curve `TIME`) |
| `tce_duration` ± | First to last contact | **hours** | matching |
| `tce_depth` ± | Flux decrement at transit minimum | ppm | reported |
| `tce_model_snr` | Transit SNR: depth ÷ mean in-transit flux uncertainty | — | "SPOC SNR" in the pilot's sample stratification |
| `tce_num_transits` | Observed transit count ([K]: excludes transits wholly in gaps) | count | reported |

### 4.2 Identity and bookkeeping

| Column | Meaning | Unit | Src |
|---|---|---|---|
| `sectors` | Delivery sectors of this search (e.g. `s0001`) | — | M |
| `lastUpdate` | Date of last update | date | M |
| `tce_ntoi` | ⚠️ **Number of TCEs on this star** (XML `planetCandidateCount`). Despite the name, **not** a TOI count | count | M |
| `tce_sectors` | Sectors observed (string, e.g. `01`) | — | M |

### 4.3 Transit model fit (limb-darkened, all transits)

| Column | Meaning | Unit | Src |
|---|---|---|---|
| `tce_time0` ± | The same epoch as full BJD (= `tce_time0bt` + 2457000.0). Measured: `_err` equals `tce_time0bt_err` (the offset was not added, rightly) | BJD | M |
| `tce_ror` ± | Planet radius ÷ stellar radius | — | M,K |
| `tce_dor` ± | Semi-major axis ÷ stellar radius (strictly, star–planet distance at mid-transit ÷ R★) | — | M,K |
| `tce_incl` ± | Orbital inclination | deg | M,K |
| `tce_impact` ± | Impact parameter | R★ (unit inferred) | M,K |
| `tce_ingress` ± | First to second contact | hours | M,K |
| `tce_eccen` ±, `tce_longp` ± | Eccentricity; longitude of periastron. Measured: 0 in sampled rows (circular fit) | —; deg | M,K |
| `tce_limbdark_mod` | Limb-darkening model name | — | M,K |
| `tce_ldm_coeff1..4` | Limb-darkening coefficients (held fixed) | — | M,K |
| `tce_trans_mod` | Transit model name (Mandel–Agol) | — | M,K |
| `tce_full_conv` | Fit fully converged | bool | M,K |
| `tce_model_chisq`, `tce_model_dof` | χ² of the fit and its degrees of freedom (plain, not reduced) | — | M,K |

### 4.4 Derived planet and stellar parameters

| Column | Meaning | Unit | Src |
|---|---|---|---|
| `tce_prad` ± | Planet radius (= ror × R★) | R⊕ | M,K |
| `tce_sma` ± | Semi-major axis, from period + stellar mass | AU | M,K |
| `tce_eqt` ± | Equilibrium temperature (Bond albedo 0.3, full redistribution) | K | M,K |
| `tce_insol` ± | Insolation | Earth flux | M,K |
| `tce_steff` ± | Stellar effective temperature | K | M,K |
| `tce_slogg` ± | Stellar surface gravity | log₁₀(cm s⁻²) | M,K |
| `tce_smet` ± | Stellar [Fe/H] | dex | M,K |
| `tce_sradius` ± | Stellar radius | R☉ | M,K |
| `tce_sdensity` ± | Stellar density (measured: consistent with M/R³ from logg and R★) | ρ☉ | M |
| `tce_*_prov` | Provenance of each stellar parameter (`TIC8.1`, `TIC8.1-Der`, `Solar`) | — | M,K |

### 4.5 Trapezoidal fit — `tcet_*`

A second, simpler model (XML `trapezoidalFit`) with the period fixed to the
search value [K]. Measured: every `tcet_*_err` is **−1** (no uncertainty).

| Column | Unit |
|---|---|
| `tcet_period` | days |
| `tcet_time0bt` / `tcet_time0` | BTJD / BJD |
| `tcet_duration`, `tcet_ingress` | hours |
| `tcet_depth` | ppm |
| `tcet_full_conv`, `tcet_model_chisq`, `tcet_model_dof` | — |

### 4.6 Detection statistics (TPS)

| Column | Meaning | Unit | Src |
|---|---|---|---|
| `tce_max_sngle_ev` | Maximum single event statistic (SES) | σ (inferred) | M,K |
| `tce_max_mult_ev` | Maximum **multiple event statistic (MES)** — the 7.1σ detection statistic | σ (inferred) | M,K |
| `tce_robstat` | Robust statistic: consistency of transit depths across events | — | M,K |
| `tce_chisq2`, `tce_dof2` | Second χ² veto discriminator and its dof | — | M,K |
| `tce_chisqgof`, `tce_chisqgofdof` | χ² goodness of fit between the TPS and DV signal, and its dof | — | M,K |

### 4.7 Weak secondary test — `wst_*`, `tce_ws_*`

Searches the folded light curve, with the primary transit removed, for a
secondary eclipse — the main EB tell. XML `weakSecondary`.

| Column | Meaning | Unit | Src |
|---|---|---|---|
| `wst_robstat` | Robust statistic of the strongest secondary event | — | M,K |
| `wst_depth` | Depth of the strongest secondary | ppm | M,K |
| `tce_ws_mesmedian`, `tce_ws_mesmad` | Median and MAD of the MES over all phases | σ | M,K |
| `tce_ws_maxmes`, `tce_ws_minmes` | MES of the strongest secondary; of the strongest positive excursion | σ | M,K |
| `tce_ws_maxmesd`, `tce_ws_minmesd` | Phase of those two events from the primary | days | M,K |

### 4.8 Binary-discrimination and other vetting tests

| Column | Meaning | Unit | Src |
|---|---|---|---|
| `tce_bin_oedp_stat` | Odd/even transit depth comparison | χ² statistic | M,T18 |
| `tce_bin_spc_stat`, `tce_bin_lpc_stat` | Period comparison with the TCE on the same star with the next **shorter** / **longer** period | χ² statistic | M,T18 |
| `tce_albedo` ± | Geometric albedo implied by the secondary depth; > 1 suggests a self-luminous companion | — | M,K |
| `tce_ptemp` ± | Planet temperature implied by the secondary depth | K | M,K |
| `tce_albedo_stat` | (albedo − 1) | σ | M,K |
| `tce_ptemp_stat` | (implied − equilibrium temperature) | σ | M,K |
| `tce_cap_stat`, `tce_hap_stat` | Ghost diagnostic: core- and halo-aperture correlation; core < halo suggests an optical ghost | ~N(0,1) | M,K,T18 |

⚠️ The three `tce_bin_*_stat` columns are the **test statistic**, not its
significance ([M] reads the XML `value`; [T18] reports the significance
separately). The Kepler table's `tce_bin_oedp_sig` is a significance in
percent — **not interchangeable**.

### 4.9 Bootstrap false-alarm test

| Column | Meaning | Unit |
|---|---|---|
| `boot_fap` | Bootstrap false-alarm probability (XML `significance`) | probability |
| `boot_mesthresh` | MES threshold giving the same false-alarm probability as 7.1σ Gaussian noise | σ |
| `boot_mesmean`, `boot_messtd` | Mean and std of the bootstrap null MES distribution | σ |
| `bootstrap_transit_count` | Transits used (−1 = not run) | count |

### 4.10 Centroid offsets — is the dip on this star at all?

| Column | Meaning | Unit |
|---|---|---|
| `tce_dicco_mra/mdec/msky` ± | Difference-image centroid offset from the **out-of-transit centroid**, weighted mean over sectors | arcsec |
| `tce_ditco_mra/mdec/msky` ± | Offset from the **TIC catalogue position** ([M] labels it "KIC"; TESS note: TIC and sectors) | arcsec |

A large offset means the dip comes from a neighbouring star blended in the
aperture.

### 4.11 Missing values

Sentinels, not blanks: **−1** (bootstrap fields, `tcet_*_err`) and **0.0** in
some centroid fields. Measured: TIC 24695725 has every bootstrap field −1.
Treat both as "not available", never as a measurement.

### 4.12 Not in the TESS file

Kepler columns with no TESS counterpart here: `tce_chisq1`/`tce_dof1`, the
rolling-band columns, the flux-weighted centroid columns.
