# FFL examples — `tess`

Complete, compile-checked scenarios against this domain's own facets. Every
fenced `ffl` block below compiles:

```bash
fw ffl compile --check --primary <block>.ffl \
  --library ~/fw_handlers/fwh_tess/src/tesstransit/ffl/tess.ffl \
  --library facetwork/ffl/http.ffl --library facetwork/ffl/compare.ffl \
  --library facetwork/ffl/provenance.ffl
```

The facets in play:

| Namespace | Facet | Returns |
|---|---|---|
| `tess.source` | `TCETable(sector, dest_dir)` | `path, url, tces` |
| `tess.source` | `Targets(sector, dest_dir)` | `path, stars` |
| `tess` | `Plan(sector, kind, work_dir, targets_path, tce_path, toi_path, …)` | `chunks, stars, run_dir` |
| `tess` | `SearchChunk(chunk_path, prewhiten, workers)` | `results_path, searched, failed` |
| `tess` | `MatchSector(sector, kind, run_dir, results_paths, …)` | `run_dir, summary_path, stars, failed, found, site_arg` |
| `tess` | `BuildSite(planets, samples, lc_dir, out_dir, variants, default)` | `output_dir, page_path, bytes` |

And the workflows built from them: `tess.workflows.SectorRun` (one sector, one
kind), `SearchSectors` (many sectors, both kinds), `BuildResults` (those plus
the page), and `tess.publish.PublishResults` (plus publishing; in
`tess_site.ffl`, which needs fwh_census_us).

`kind` is `"planets"` (every TOI star in the sector, with each TOI's ephemeris
measured on our flux) or `"sample"` (`n_tce` TCE stars and `n_quiet` stars
without a TCE, drawn with `seed`). Each `(sector, kind, label)` gets its own
run directory, `<work_dir>/s<N>-<kind>[-<label>]`.

---

## 1. The minimal workflow

One sector's planet sample. The three inputs have no dependency on each other,
so they run at once; `SectorRun` waits for all three because it reads them.

```ffl
namespace ex.minimal {
    workflow Sector1Planets(work_dir: String = "/data/tess") => (found: Long, stars: Long) andThen {
        tois  = fw.http.Fetch(url = "https://archive.stsci.edu/missions/tess/catalogs/toi/tois.csv",
                              dest = $.work_dir ++ "/tois.csv", max_age_hours = 24.0)
        tces  = tess.source.TCETable(sector = 1, dest_dir = $.work_dir)
        stars = tess.source.Targets(sector = 1, dest_dir = $.work_dir ++ "/lc")
        run   = tess.workflows.SectorRun(sector = 1, work_dir = $.work_dir,
                                         targets_path = stars.path, tce_path = tces.path,
                                         toi_path = tois.path)
        yield Sector1Planets(found = run.found, stars = run.stars)
    }
}
```

## 2. A hold-out sample

The second sample may share no star with the first, so a change tuned on one
is scored on stars it never saw. `holdout_of` takes the first run's directory;
the reference is a value, so the order is already implied.

```ffl
namespace ex.holdout {
    workflow DevAndHoldout(work_dir: String, targets_path: String, tce_path: String, toi_path: String) => (dev_found: Long, holdout_found: Long) andThen {
        dev = tess.workflows.SectorRun(sector = 1, kind = "sample", label = "dev", seed = 1,
                                       work_dir = $.work_dir, targets_path = $.targets_path,
                                       tce_path = $.tce_path, toi_path = $.toi_path)
        held = tess.workflows.SectorRun(sector = 1, kind = "sample", label = "holdout", seed = 2,
                                        holdout_of = [dev.run_dir], sample_label = "Hold-out",
                                        work_dir = $.work_dir, targets_path = $.targets_path,
                                        tce_path = $.tce_path, toi_path = $.toi_path)
        yield DevAndHoldout(dev_found = dev.found, holdout_found = held.found)
    }
}
```

## 3. Many sectors, bounded

`SearchSectors` is a `foreach` over sectors, and inside each a `foreach` over
chunks. `sector_width` bounds the outer fan-out and `width` the inner one; the
total in flight is at most their product. A light curve is ~2 MB and deleted
after its search, so the bound is CPU, not disk.

```ffl
namespace ex.many {
    workflow FirstYear(work_dir: String) => (planets: Json, toi_sha256: String) andThen {
        runs = tess.workflows.SearchSectors(sectors = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
                                            work_dir = $.work_dir, sector_width = 2, width = 8,
                                            chunk_size = 50, workers = 4)
        yield FirstYear(planets = runs.planets, toi_sha256 = runs.toi_sha256)
    }
}
```

⚠️ The TOI catalogue moves: new TOIs appear and dispositions change. A re-run
is against the catalogue of its day, and `toi_sha256` records which one.

## 4. Does pre-whitening change the answer?

The same sector with pre-whitening on and off, each in its own run directory,
then the two `tois.csv` compared row by row on the TOI. `fw.compare.Tabular`
reports the highest level at which they agree, so `keys` means the same TOIs
were searched and some outcomes differ, which is the question being asked.

```ffl
namespace ex.ablation {
    workflow PrewhitenAblation(sector: Int, work_dir: String, targets_path: String,
                               tce_path: String, toi_path: String) => (level: String, agree: Boolean, on_found: Long, off_found: Long) andThen {
        on = tess.workflows.SectorRun(sector = $.sector, label = "pw", prewhiten = true,
                                      work_dir = $.work_dir, targets_path = $.targets_path,
                                      tce_path = $.tce_path, toi_path = $.toi_path)
        off = tess.workflows.SectorRun(sector = $.sector, label = "nopw", prewhiten = false,
                                       work_dir = $.work_dir, targets_path = $.targets_path,
                                       tce_path = $.tce_path, toi_path = $.toi_path)
        cmp = fw.compare.Tabular(expected = off.run_dir ++ "/tois.csv",
                                 actual = on.run_dir ++ "/tois.csv",
                                 key_columns = ["toi"],
                                 report = $.work_dir ++ "/prewhiten-ablation.json")
        yield PrewhitenAblation(level = cmp.level, agree = cmp.agree,
                                on_found = on.found, off_found = off.found)
    }
}
```

## 5. Gate the page on a clean run

A star that fails (a corrupt light curve, a negative median flux) is recorded,
not dropped, and counted in `failed`. Here the page is built only when nothing
failed; otherwise the workflow returns the count and stops.

```ffl
namespace ex.gate {
    workflow PageIfClean(work_dir: String, targets_path: String, tce_path: String, toi_path: String) => (page: String, failed: Long) andThen {
        run = tess.workflows.SectorRun(sector = 1, work_dir = $.work_dir,
                                       targets_path = $.targets_path, tce_path = $.tce_path,
                                       toi_path = $.toi_path) andThen when {
            case $.failed == 0 => {
                site = tess.BuildSite(planets = [$.site_arg], samples = [],
                                      lc_dir = $$.work_dir ++ "/lc", out_dir = $$.work_dir ++ "/site")
                yield PageIfClean(page = site.page_path, failed = 0)
            }
            case _ => {
                yield PageIfClean(page = "", failed = $.failed)
            }
        }
    }
}
```

## 6. Record what produced the page

`fw.provenance.Manifest` writes the software versions, the inputs' recorded
provenance and a digest of every output. The TOI catalogue was fetched with
`fw.http.Fetch`, so it carries a sidecar; the TCE table and star list do not,
and the manifest says so in `incomplete` rather than leaving them out.

```ffl
namespace ex.provenance {
    workflow PageWithManifest(work_dir: String) => (page: String, complete: Boolean, incomplete: [String]) andThen {
        built = tess.workflows.BuildResults(sectors = [1, 2], work_dir = $.work_dir)
        man = fw.provenance.Manifest(outputs = built.output_dir,
                                     inputs = [$.work_dir ++ "/tois.csv",
                                               $.work_dir ++ "/s0001_tcestats.csv",
                                               $.work_dir ++ "/s0002_tcestats.csv"],
                                     dest = built.output_dir ++ "/provenance.json")
        yield PageWithManifest(page = built.page_path, complete = man.complete,
                               incomplete = man.incomplete)
    }
}
```

---

## Syntax cheat sheet

| | |
|---|---|
| step | `name = ns.Facet(arg = value)` |
| reference | `step.field`, never a bare step name |
| container attributes | `$.param` (this block), `$$.param` (one up) |
| string concatenation | `$.work_dir ++ "/tois.csv"` |
| fan-out | `andThen foreach x in $.list limit $.width { … yield F(out += step.field) }` |
| branch | `step = F(…) andThen when { case $.n == 0 => { … } case _ => { … } }` |
| ordering without a value | `b = F(…) after a` |
