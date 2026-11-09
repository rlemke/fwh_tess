# fwh_tess

An independent transit re-search of TESS 2-minute light curves, compared
against the SPOC pipeline's single-sector Threshold Crossing Events (TCEs).

The question is reproduction, not discovery: does a different, fully declared
search over the **same** light curves recover SPOC's detections, and where the
two disagree, is the cause the method, the data, or a defect?

⚠️ SPOC searches with TPS (a wavelet-based matched filter, 7.1σ MES); this
package uses Box Least Squares. Some disagreement is expected from the method
alone. Signals with no TCE are **unvetted** and are never reported as planets.

## Status

Sectors 1–13 searched. Of the planets and candidates SPOC detected in the same
sector, **682 of 707 (96.5%)** are found again. Results page:
<https://rlemke.github.io/facetwork-maps/space/tess-sector1/>.

A Facetwork domain (`tess`): the command-line tools and the FFL workflows run
the same code, so either way gives the same files.

## FFL at a glance

```ffl
namespace demo {
    workflow OneSector(work_dir: String = "/data/tess") => (found: Long, run_dir: String) andThen {
        tois    = fw.http.Fetch(url = "https://archive.stsci.edu/missions/tess/catalogs/toi/tois.csv",
                                dest = $.work_dir ++ "/tois.csv", max_age_hours = 24.0)
        tces    = tess.source.TCETable(sector = 1, dest_dir = $.work_dir)
        stars   = tess.source.Targets(sector = 1, dest_dir = $.work_dir ++ "/lc")
        planets = tess.workflows.SectorRun(sector = 1, work_dir = $.work_dir,
                                           targets_path = stars.path, tce_path = tces.path,
                                           toi_path = tois.path)
        yield OneSector(found = planets.found, run_dir = planets.run_dir)
    }
}
```

Every sector end to end: `tess.workflows.SearchSectors` (runs) or
`tess.workflows.BuildResults` (runs + page); `tess.publish.PublishResults` also
publishes. Full gallery: [`docs/ffl-examples.md`](docs/ffl-examples.md).

## Data (measured 2026-11-04)

| What | Where |
|---|---|
| Light curves | `s3://stpubdata/tess/public/tid/s00NN/<TIC split 4/4/4/4>/` — MAST's public copy, anonymous reads. Sector 1: 15,889 `_lc.fits`, 32.4 GB |
| TCE reference | `https://archive.stsci.edu/missions/tess/catalogs/tce/<stamp>-s00NN-s00NN_dvr-tcestats.csv` — sector 1: 1,457 TCEs on 971 stars. ⚠️ duration is in **hours** |

Light curves are streamed and can be discarded: only per-star results are kept.

Each experiment and what it changed is in [docs/experiments.md](docs/experiments.md).
Every term, field and outcome — ours and SPOC's `tce_*` columns — is defined in
[docs/definitions.md](docs/definitions.md).

## Layout

| Module | |
|---|---|
| `_lightcurve.py` | locate / fetch / load one light curve; the quality mask is counted, not hidden |
| `_search.py` | spline detrend (knots by BIC, pre-whitening below 1.5 d) + iterative BLS, durations capped by period (Kepler's law); every parameter in `SearchConfig`, which travels with the result |
| `_manifest.py` | which stars the sector observed, from the bucket itself |
| `pilot.py` | `python -m tesstransit.pilot` — stratified sample → search → match |
| `inject.py` | `python -m tesstransit.inject` — synthetic transits into the real light curves of a sample run's no-TCE stars → search → recovered or not |
| `multisector.py` | `python -m tesstransit.multisector` — each sector flattened, then stitched and searched as one light curve (1 d to half the baseline), scored against SPOC's multi-sector TCEs |
| `vet.py` | `python -m tesstransit.vet` — odd/even, secondary-eclipse and V-shape tests at a signal's ephemeris; thresholds measured on recovered CP/KP planets vs EB TOIs |
| `ttv.py` | `python -m tesstransit.ttv` — each transit timed on the stitched light curve, O − C against the best fixed period |
| `planets.py` | `python -m tesstransit.planets` — every TOI in the sector → search → match against the TOI, with each ephemeris forced on our flux |
| `site.py` | `python -m tesstransit.site` — the stand-alone results page (per-sector figures, period–depth chart, sky map, experiments) from finished runs of any number of sectors; published at rlemke.github.io/facetwork-maps/space/tess-sector1/ |
| `_match.py` | the TCE table and the matching rules: `recovered` / `recovered_alias` / `missed` (with a cause) / `unmatched` |
| `handlers.py` | the Facetwork steps: `TCETable`, `Targets`, `Plan`, `SearchChunk`, `MatchSector`, `BuildSite` — each a thin layer over the modules above |
| `ffl/tess.ffl` | facets and workflows (`SectorRun`, `SearchSectors`, `BuildResults`); `ffl/tess_site.ffl` adds `PublishResults` (needs fwh_census_us) |
| `domain.py` | the `facetwork.domains` entry point, kept apart so the search runs without Facetwork |

## The false-alarm rate is measured

Searching the same light curves turned upside down (`pilot --invert`) turns
every transit into a brightening, so whatever the inverted search still
reports is a false alarm. Over the 520 no-TCE sample stars of sectors 1–13, at
the current cut (SDE 11, SNR 7.1): **13.5% (10.8–16.7%) give a false alarm**.
The normal search gives a signal on 12.7% of the same stars, so signals on
stars without a TCE are, as a group, what false alarms alone produce. Most are
periodic variability under a day that survives pre-whitening, not white
noise, and a stricter cut costs planets faster than it removes them (SDE 14,
measured with the earlier three-term pre-whitening: 9.8% false alarms, 92.9%
of planets). Details: `docs/experiments.md`.

## Completeness is measured

Synthetic transits injected into the real light curves of the same 520 no-TCE
stars (`python -m tesstransit.inject`, 6 per star, 3,120 in all) and searched
with the unchanged pipeline. Depth is set from an **expected SNR** (the
white-noise SNR from the star's own scatter), drawn from 4 to 40.

| | 50% found at | 90% found at | found, expected SNR ≥ 20 |
|---|---|---|---|
| all stars | SNR 9.4 | SNR 19.7 | 92.0% |
| stars without fast variability (82%) | SNR 7.9 | SNR 11.7 | 99.0% |
| pulsating stars, pre-whitened (18%) | SNR 24.7 | SNR 56.2 | 57.9% |

On ordinary stars the search does what its SNR 7.1 threshold promises. On
pulsating stars the residual red noise buries the transit. Pre-whitening up
to ten periodicities down to 0.01 d (tuned on these injections; it was three
down to 0.05 d) raised strong transits found there from 48% to 58% with the
same planets found and no change on other stars. Details: `docs/experiments.md`.

## Planet or eclipsing binary?

`vet.py` applies three standard shape tests to a signal: alternating
odd/even depths, an eclipse at phase 0.5, and a V-shaped dip. Each needs a
size as well as a significance. The thresholds were measured on the TOIs
sectors 1-13 recovered, and held up tuned on odd sectors and tested on
even ones:

| | flagged as binary-like |
|---|---|
| confirmed or known planets (CP, KP) | 10 of 405 (2.5%) |
| TOIs dispositioned EB | 76 of 267 (28%) |
| planet candidates (PC) | 19 of 385 (4.9%) |
| signals on stars with no TCE | 18 of 105 (17%) |

Most EB TOIs pass because they are diluted background binaries, which look
like planets in shape and need a centroid test. Details: `docs/experiments.md`.

## Searching all thirteen sectors at once

`multisector.py` flattens each sector, stitches the year of data and
searches 1 d to half the baseline. In a pilot of 200 stars (150 with a TCE
from SPOC's sectors 1-13 search, 50 without), with at least nine sectors
each, it found all 9 TOIs in the sample, including TOI 700.02 (37.4 d),
713.01 (36.0 d) and 705.01 (47.0 d). A single sector holds at most one
transit of any of these.

At the single-sector thresholds it is not usable. Inverted, it put a
signal on 22 of the 50 stars without a TCE, because it tries 15 times as many
periods and a pair of noise dips is easy to line up. Three transits and SDE
15 looked like a fix on that pilot (3 of 50) but failed a fresh sample (8 of
50). What works is a per-transit test: with the two strongest transits left
out, the rest must still show the dip at SNR 3. Together with SDE 13, which
is now the stitched search's default, a third sample fixed in advance put an
inverted signal on 3 of its 50 stars without a TCE, and kept all 9 in-grid
planet TOIs (including TOI 1209.01 at 40.7 d). ⚠️ On 250 more stars without a
TCE, chosen fresh, the inverted rate is 10.4% (7.2-14.8%), lower than a
single sector's 14% but not clearly so. Half of these false alarms are 1-2.5 d
stellar variability, found at the same period in both directions. The other
half are long-period signals whose events sit at data-segment edges. A trim of
the segment edges plus a both-directions test, fixed in advance, kept every
planet TOI on a fresh sample but still flagged 21 of 250 (8.4%), so neither is
on by default. Nor is a variability veto (a signal is rejected when, with its
transits masked, a two-sided wave at a harmonic of its period is at least
half its depth): fixed in advance, it flagged 17 of 250 (6.8%) against the
15 allowed. So the stitched false-alarm rate stands at about 7-10%: treat a
stitched signal on a star with no TCE as unvetted. The two planets that
sample missed (TOI 810.02 at 90.8 d, TOI 6902.01 at 36.4 d) are not missed
for their 9-11 h transits: 810.02 is in the data but falls below the SDE 13
cut, and 6902.01's transits on a very bright star are as large as the
detrending residuals. Details:
`docs/experiments.md`.

## Tests

```bash
pip install -e '.[dev]'
pytest
```

No test touches the network.
