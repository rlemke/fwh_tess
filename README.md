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

Phase 0 (core + offline tests). Not yet a Facetwork domain: no FFL or
handlers yet.

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
| `_search.py` | running-median detrend + iterative BLS, durations capped by period (Kepler's law); every parameter in `SearchConfig`, which travels with the result |
| `_manifest.py` | which stars the sector observed, from the bucket itself |
| `pilot.py` | `python -m tesstransit.pilot` — stratified sample → search → match |
| `planets.py` | `python -m tesstransit.planets` — every TOI in the sector → search → match against the TOI, with each ephemeris forced on our flux |
| `site.py` | `python -m tesstransit.site` — the stand-alone results page (per-sector figures, period–depth chart, sky map, experiments) from finished runs of any number of sectors; published at rlemke.github.io/facetwork-maps/space/tess-sector1/ |
| `_match.py` | the TCE table and the matching rules: `recovered` / `recovered_alias` / `missed` (with a cause) / `unmatched` |

## Thresholds are provisional

On synthetic white noise the strongest BLS peak reached SNR 6.4–7.8, so an
SNR cut at 7.1 alone accepts about half of noise-only stars. An SDE ≥ 7.0 gate
separates them (noise 3.6–5.1, injected 300–1500 ppm transits 12.5–14.3). Real
TESS noise is red, so both thresholds are to be calibrated on real stars that
have no TCE.

## Tests

```bash
pip install -e '.[dev]'
pytest
```

No test touches the network.
