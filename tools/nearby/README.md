# Nearby M dwarfs with no known planet (sample for the stitched search)

How `docs/data/nearby-mdwarfs-300.json` was built (2026-11-24):

1. `gaia_near.py OUT.csv`: Gaia DR3 stars with parallax > 25 mas (closer than 40 pc) and parallax_over_error > 10 (23,642 rows).
2. `tl_fetch.sh`: the TESS 2-minute target lists for sectors 1-96 (tess.mit.edu), which give each target's TIC position.
3. `xmatch.py GAIA.csv OUT.json`: Gaia positions are moved back to epoch 2000 with their proper motions and matched to the 2-minute targets within 4 arcsec (7,151 matches). Sectors per star come from the sector manifests. Stars hosting a TOI, a CTOI, or a confirmed planet (archive pscomppars) are flagged. `TESS_DATA` is the directory holding the manifests and catalogue files.
4. `pick.py`: of those with data in at least 3 of sectors 1-96, no TOI, CTOI or confirmed planet, and Gaia BP-RP > 1.8 (M dwarfs), the 300 nearest (all within 22.6 pc).

Of all stars within 40 pc with 3+ sectors, 3,333 have no TOI, CTOI or known planet; 2,258 of them are M dwarfs. This list is the first 300 of those M dwarfs, nearest first.

SPOC already searched these stars over sectors 1-96 (its multi-sector TCE table). A signal is new only if it has no SPOC TCE, single- or multi-sector, as well as no TOI, CTOI or planet, and it must not appear in the same star's upside-down search.
