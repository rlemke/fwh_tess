"""fwh_tess — an independent transit re-search of TESS, compared against SPOC.

One question: does a different, fully declared search over the SAME light
curves reproduce the SPOC pipeline's single-sector detections, and where the
two disagree, why?

- :mod:`._lightcurve` — locate, fetch and load one SPOC 2-minute light curve
  from MAST's public S3 copy, with the quality mask reported, not hidden.
- :mod:`._search` — detrend, then iterative Box Least Squares.
- :mod:`._match` — the SPOC TCE table, and the rules that decide whether one
  of our signals is the same thing as one of their TCEs.

⚠️ SPOC searches with TPS (a wavelet-based matched filter, 7.1σ MES); this
package searches with BLS. Disagreement is expected from the method alone, so
every comparison must say which side of that line a difference falls on.

Not yet a Facetwork domain: the entry point is added with the FFL.
"""
