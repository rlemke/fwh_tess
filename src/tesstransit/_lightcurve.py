"""One SPOC 2-minute light curve: where it lives, how to fetch it, how to load it.

MAST mirrors TESS to a public, anonymously readable bucket. A target's files
for one sector sit in one directory keyed by its zero-padded TIC ID, split
4/4/4/4 (measured 2026-11-04, sector 1)::

    s3://stpubdata/tess/public/tid/s0001/0000/0002/6113/6679/
        tess2018206045859-s0001-0000000261136679-0120-s_lc.fits
        tess2018206045859-s0001-0000000261136679-0120-s_tp.fits
        ...-00106_dvr.pdf / _dvr.xml / _dvm.pdf      (only where SPOC found a TCE)

The light-curve file name carries a per-sector timestamp we cannot derive, so
the key is found by LISTING the target's directory, not constructed.

⚠️ Quality masking is reported, never silent: :class:`LightCurve` carries how
many cadences were dropped and why, because a mask that removes too much can
erase a transit, and that has to be visible when a TCE is missed.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

BUCKET = "stpubdata"
REGION = "us-east-1"

# Every QUALITY bit set means SPOC itself flagged the cadence. We drop any
# flagged cadence by default; the mask is a parameter so it can be relaxed and
# the effect measured.
DEFAULT_QUALITY_BITMASK = 0xFFFFFFFF


def target_prefix(tic_id: int, sector: int) -> str:
    """The S3 directory holding one target's files for one sector."""
    if tic_id <= 0:
        raise ValueError(f"TIC ID must be positive, got {tic_id}")
    if sector <= 0:
        raise ValueError(f"sector must be positive, got {sector}")
    t = f"{tic_id:016d}"
    return f"tess/public/tid/s{sector:04d}/{t[0:4]}/{t[4:8]}/{t[8:12]}/{t[12:16]}/"


def s3_client() -> Any:
    """An anonymous client for MAST's public bucket (no credentials needed)."""
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config

    # Retries with backoff: a pilot lost its network mid-run and 133 of 200
    # fetches failed with EndpointConnectionError; the default retry budget
    # gave up within seconds.
    return boto3.client(
        "s3",
        region_name=REGION,
        config=Config(
            signature_version=UNSIGNED,
            retries={"max_attempts": 10, "mode": "standard"},
            connect_timeout=20,
            read_timeout=60,
        ),
    )


def find_lc_key(client: Any, tic_id: int, sector: int) -> str:
    """The key of the 2-minute light curve, or raise if the target has none.

    Raises ``FileNotFoundError`` rather than returning ``None``: a target in
    the sector's manifest with no light curve is a finding, not a skip.
    """
    prefix = target_prefix(tic_id, sector)
    resp = client.list_objects_v2(Bucket=BUCKET, Prefix=prefix)
    keys = [o["Key"] for o in resp.get("Contents", []) if o["Key"].endswith("_lc.fits")]
    if not keys:
        raise FileNotFoundError(
            f"no 2-minute light curve for TIC {tic_id} sector {sector} "
            f"under s3://{BUCKET}/{prefix}"
        )
    if len(keys) > 1:
        raise ValueError(
            f"TIC {tic_id} sector {sector}: {len(keys)} light curves, expected 1: {keys}"
        )
    return keys[0]


def fetch_lc(
    client: Any, tic_id: int, sector: int, dest_dir: Path, key: str | None = None
) -> Path:
    """Download one light curve into ``dest_dir``; reuse it if already there.

    ``key`` comes from the sector manifest when there is one, saving a listing
    per star. Reuse is by exact size against the object's listing, so a
    truncated earlier download is fetched again rather than trusted.
    """
    key = key or find_lc_key(client, tic_id, sector)
    dest = Path(dest_dir) / key.rsplit("/", 1)[1]
    size = client.head_object(Bucket=BUCKET, Key=key)["ContentLength"]
    if dest.exists() and dest.stat().st_size == size:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    client.download_file(BUCKET, key, str(tmp))
    tmp.replace(dest)
    return dest


@dataclass
class LightCurve:
    """Normalised PDCSAP flux with what was removed to get it."""

    tic_id: int
    sector: int
    time: np.ndarray  # BTJD (BJD - 2457000), the same clock as the TCE table
    flux: np.ndarray  # PDCSAP, divided by its median
    flux_err: np.ndarray
    n_cadences: int  # rows in the file
    n_flagged: int  # dropped by the quality bitmask
    n_nonfinite: int  # dropped for NaN time or flux (after the quality mask)

    @property
    def n_kept(self) -> int:
        return len(self.time)


def load_lc(
    path: Path,
    quality_bitmask: int = DEFAULT_QUALITY_BITMASK,
) -> LightCurve:
    """Read a SPOC ``_lc.fits`` into a :class:`LightCurve`."""
    from astropy.io import fits

    with fits.open(path) as hdul:
        hdr = hdul[0].header
        data = hdul[1].data
        time = np.asarray(data["TIME"], dtype=float)
        flux = np.asarray(data["PDCSAP_FLUX"], dtype=float)
        err = np.asarray(data["PDCSAP_FLUX_ERR"], dtype=float)
        quality = np.asarray(data["QUALITY"], dtype=np.int64)
        tic_id = int(hdr["TICID"])
        sector = int(hdr["SECTOR"])

    n = len(time)
    flagged = (quality & quality_bitmask) != 0
    finite = np.isfinite(time) & np.isfinite(flux) & np.isfinite(err)
    keep = ~flagged & finite
    if not keep.any():
        raise ValueError(
            f"{path}: no usable cadences ({n} rows, {flagged.sum()} flagged)"
        )
    med = np.median(flux[keep])
    if not med > 0:
        raise ValueError(f"{path}: median PDCSAP flux {med} is not positive")
    return LightCurve(
        tic_id=tic_id,
        sector=sector,
        time=time[keep],
        flux=flux[keep] / med,
        flux_err=err[keep] / med,
        n_cadences=n,
        n_flagged=int(flagged.sum()),
        n_nonfinite=int((~flagged & ~finite).sum()),
    )
