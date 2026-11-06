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
the current cut (SDE 11, SNR 7.1): **14.0% (11.3–17.3%) give a false alarm**.
The normal search gives a signal on 13.5% of the same stars, so signals on
stars without a TCE are, as a group, what false alarms alone produce. Most are
periodic variability under a day that survives pre-whitening, not white
noise, and a stricter cut costs planets faster than it removes them (SDE 14:
9.8% false alarms, 92.9% of planets). Details: `docs/experiments.md`.

## Tests

```bash
pip install -e '.[dev]'
pytest
```

No test touches the network.
