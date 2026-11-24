# Every tool as FFL

Each command-line tool in this package has an FFL workflow that runs the same
code, so a run from either writes the same files. The handlers are thin: a
step localizes its inputs, calls the CLI's own functions (or its `main`), and
uploads what it wrote. Paths may be local or `s3://`.

## Which workflow stands for which tool

| Tool | What it does | Facets | Workflow |
|---|---|---|---|
| `pilot` | a sector's TCE/no-TCE sample: search, match against the TCEs | `tess.Plan` → `tess.SearchChunk` → `tess.MatchSector` | `tess.workflows.SectorRun(kind = "sample")` |
| `pilot --invert` | the false-alarm control: the sample searched upside down | the same, `invert = true` | `SectorRun(kind = "sample", invert = true)`; every sector: `FalseAlarms` |
| `planets` | every TOI star of a sector, each ephemeris forced on our flux | the same, `kind = "planets"` | `SectorRun(kind = "planets")` |
| `pilot`/`planets` over many sectors | both kinds for each sector | `tess.source.TCETable`, `tess.source.Targets` + the above | `SearchSectors` |
| `--set FIELD=VALUE`, `--no-prewhiten` | SearchConfig overrides | `set`, `prewhiten` on `SectorRun` / `SearchChunk` / `MatchSector` | (any of the above) |
| `inject` | injection-recovery in the no-TCE stars of a sample run | `tess.check.Inject` | `Completeness` |
| `vet` | odd/even, secondary eclipse and V-shape tests | `tess.check.Vet` | `VetPlanets` |
| `ttv` | each transit timed on the stitched light curve, O − C | `tess.check.Timing` | `TimeTOIs` |
| `multisector` (and `--invert`, `--sample-of`, `--holdout-of`) | a stitched sample, as observed and upside down | `tess.source.StitchTCETable`, `tess.stitch.Plan` → `tess.stitch.SearchChunk` → 2 × `tess.stitch.Match` | `StitchedSample` |
| `site` | the results page over finished runs | `tess.BuildSite` | `BuildPage` (finished runs); `BuildResults` (searches first) |
| publishing the page | the page on the maps site | `census.Publish.PublishWebBundle` | `tess.publish.PublishPage` / `PublishResults` (`tess_site.ffl`) |

The run directories are named as the CLI runs were: `<work_dir>/s<N>-<kind>[-<label>]`
for a sector run, `<run>-inject` for injection-recovery, `<label>` and
`<label>-inv` for a stitched sample.

Not separate workflows, on purpose:

- **`--reuse`** (re-match an existing `results.jsonl`): that is what
  `MatchSector` is. It gathers the chunks and runs the CLI's matcher with
  `--reuse`.
- **The pilot 1-6 detrending variants** (`--detrend-method median|biweight`,
  `--per-duration`): they are SearchConfig fields, reached through `set`
  (`set = ["detrend_method=median"]`). They are the experiment log's history,
  not defaults.
- **`--keep-lc`**: a step's light curves live in its own scratch directory
  and are gone when it ends.
- **One stitched pass alone**: `StitchedSample` always searches both, because
  a stitched result without its upside-down control has no false-alarm rate.
- The paired re-searches in the experiment log (the same stars searched in
  two years) were one-off comparisons over `pilot._search_one`, not tools.

## The CLI and the FFL give the same files (measured 2026-11-20)

Each pair below was run on the same day with the same code: the CLI, then the
FFL in memory with the real handlers (`tests/_drive.py`), sector 1 or the
stitched year.

| Check | CLI | FFL | Result |
|---|---|---|---|
| upside-down control, sector 1 (100 stars) | `pilot --invert` | `SectorRun(kind = "sample", invert = true)` | `matches.csv` identical; every search result identical; `summary.json` same config and counts |
| shape tests, sector 1 planet run (51 signals) | `vet` | `VetPlanets` | `vet.csv` identical |
| TOI 216 timing, sectors 1-13 | `ttv --toi 216.01 --toi 216.02` | `tess.check.Timing` | `ttv.jsonl` identical (2 planets, both timed) |
| injection-recovery, sector 1 sample (2 per star) | `inject --per-star 2` | `Completeness(per_star = 2)` | `injections.jsonl` identical except timings, `injections.csv` identical (80 injections, 52 found) |
| the published page (sectors 1-55, 4 stitched, 3 two-year samples) | `site` (the 2026-11-16 build) | `BuildPage` | 6,804,839 bytes both; identical except the snapshot date |

Earlier: `BuildResults` over sectors 1-13 reproduced 25 of 26 CLI runs byte
for byte (2026-11-06); the 26th differed only in an error message carrying a
scratch path, since fixed. The stitched search has run only as FFL on the
fleet since 2026-11-14, and its results match the CLI run on one host (below).

## Rebuilding the published page

The page at `space/tess-sector1/` on the maps site is `BuildPage` (or, to
publish it too, `tess.publish.PublishPage`) over these runs, all under one
data directory `<data>` with the light-curve manifests and sky-position caches
in `<data>/lc`:

- `planets`: `"<N>=<data>/runs/s<N>-planets-pwC"` for N = 1-55
- `samples`: `"1=Development=<data>/runs/s1-pilot-pwC"`,
  `"1=Hold-out=<data>/runs/holdout-pwC"`, and `"<N>=Sample=<data>/runs/s<N>-pilot-pwC"`
  for N = 2-55
- `variants` (the sector-1 detrending scorecard), `default = "7d"`:
  - `Pilot 5|1-day running median|planets-s1|pilot5-veto-fixed|holdout-seed2`
  - `7a|Cubic spline, knots by BIC|planets-spline7a|dev-spline7a|holdout-spline7a`
  - `7b|7a + flexible short-box layer|planets-spline7b|dev-spline7b|holdout-spline7b`
  - `7c|7a + gated pre-whitening under 1.5 d|s1-planets-pw2|s1-pilot-pw2|holdout-pw2`
  - `7d|7c with up to 10 terms, down to 0.01 d|s1-planets-pwB|s1-pilot-pwB|holdout-pwB`

  each run as `<data>/runs/<name>`
- `stitched` (one year, sectors 1-13), each `label|settings|run|inverted run|adopted`:
  - `A|Adopted: two-event test, SDE 13|ms-v3|ms-v3-inv|yes`
  - `B|Adopted: two-event test, SDE 13|ms-fa250|ms-fa250-inv|yes`
  - `C|Adopted + gap-edge trim + upside-down test (not adopted)|ms-v4|ms-v4-inv|no`
  - `D|Adopted + variability veto (not adopted)|ms-v5|ms-v5-inv|no`
- `long_stitched` (two years), all `SDE 18 (adopted for two years)` and `yes`:
  - `South A`: `ms26-fresh` / `ms26-fresh-inv`
  - `North A`: `n4-fresh` / `n4-fresh-inv`
  - `South B`: `south-y13-fleet` / `south-y13-fleet-inv`

Every run named here was itself made by a tool in the first table, with the
settings the experiment log records for it.

## Running one on the fleet

```bash
fw ffl run --primary src/tesstransit/ffl/tess.ffl \
    --library <facetwork>/facetwork/ffl/http.ffl \
    --workflow tess.workflows.StitchedSample \
    --inputs '{"sectors": [1, 2, …, 13, 27, …, 39], "work_dir": "s3://…/stitched",
               "label": "south-c", "min_sectors": 18, "n_tce": 150, "n_quiet": 150,
               "seed": 41, "chunk_size": 3, "workers": 3}'
```

The run appears in the dashboard like any other. What the fleet does with the
stitched chunks (CPU booking, fast-host placement, the reserve a desktop host
keeps) and how three re-runs of one sample went is in
[experiments.md](experiments.md#the-stitched-search-on-the-fleet-2026-11-17-to-11-20).

## Keeping up with new sectors

`fw svc tess-update` (in the facetwork repo) checks MAST every day and, for each
sector released since its last run, submits `tess.workflows.SearchSectors` and
`tess.workflows.FalseAlarms` for it, with the defaults the CLI runs use. A
sector counts as released only when its SPOC single-sector TCE table is listed
AND its 2-minute light curves are in the bucket. Results go to
`s3://afl-cache/cache/tess/sectors`, one run directory per sector and kind.

- MAST gives no release calendar. Measured in 2026, a new sector arrived every
  3-6 weeks (s104 07-13, s105 08-17, s106 09-10, s107 10-01, s108 11-12, dates
  the first file was written).
- **The TOI catalogue lags the data by months.** On 2026-11-23 no TOI listed a
  sector after 96 and few listed one after 70, so a new sector's planet sample
  is usually empty and only the TCE/no-TCE sample says anything. Re-run the
  planet sample of a sector later, against a newer catalogue, with
  `fw svc tess-update --sectors N`.
- `--check` alarms when no new sector has arrived for 70 days, and MAST
  reprocessing an old sector is reported, never re-run unattended.

⚠️ An empty planet sample first ran on 2026-11-23 (sector 108) and failed:
`SearchEach` returned `results: Json` only through `+=` in its `foreach`, so with
zero chunks nothing set it, and `MatchSector`'s `searched.results` was "not found".
The fan-out returns are now array-typed (`[String]`, `[Long]`), which the runtime
sets to `[]` for a zero-iteration `foreach` (that runtime step also had a bug for
facets called as steps, fixed in facetwork), and `planets.py` writes `tois.csv`
with a fixed header instead of the first row's keys.
