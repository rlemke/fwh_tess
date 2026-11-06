"""Which stars a sector actually observed at 2-minute cadence.

Taken from what is in MAST's bucket, not from a target catalogue: the CTL
lists CANDIDATE targets, and only a light curve on disk proves a star was
observed and processed. One listing of the sector prefix (sector 1: 15,889
light curves) gives every TIC ID and the exact key to fetch it by.

Light-curve names carry the TIC ID zero-padded to 16 digits::

    tess2018206045859-s0001-0000000261136679-0120-s_lc.fits
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from ._lightcurve import BUCKET

_LC_NAME = re.compile(r"-s(\d{4})-(\d{16})-\d{4}-s_lc\.fits$")


def tic_from_key(key: str) -> tuple[int, int]:
    """``(sector, tic_id)`` from a light-curve key; raise if it is not one."""
    m = _LC_NAME.search(key)
    if not m:
        raise ValueError(f"not a 2-minute light-curve key: {key}")
    return int(m.group(1)), int(m.group(2))


def list_sector(client: Any, sector: int) -> dict[int, str]:
    """Every 2-minute light curve in the sector, as ``{tic_id: key}``.

    Raises if a TIC ID appears twice: one star, one light curve per sector is
    what the rest of the pipeline assumes, so a duplicate must be seen.
    """
    prefix = f"tess/public/tid/s{sector:04d}/"
    out: dict[int, str] = {}
    for page in client.get_paginator("list_objects_v2").paginate(
        Bucket=BUCKET, Prefix=prefix
    ):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if not key.endswith("_lc.fits"):
                continue
            s, tic = tic_from_key(key)
            if s != sector:
                raise ValueError(
                    f"key {key} is under sector {sector} but names sector {s}"
                )
            if tic in out:
                raise ValueError(f"TIC {tic} has two light curves in sector {sector}")
            out[tic] = key
    return out


def load_or_list(client: Any, sector: int, path: Path) -> dict[int, str]:
    """The manifest from ``path`` if present, else list the bucket and save it."""
    path = Path(path)
    if path.exists():
        return parse(path.read_text(), sector, str(path))
    targets = list_sector(client, sector)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".part")
    tmp.write_text(dumps(sector, targets))
    tmp.replace(path)
    return targets


def dumps(sector: int, targets: dict[int, str]) -> str:
    return json.dumps({"sector": sector, "targets": targets})


def parse(text: str, sector: int, where: str = "manifest") -> dict[int, str]:
    """``{tic_id: key}`` from a saved manifest; raise if it is another sector's."""
    data = json.loads(text)
    if data["sector"] != sector:
        raise ValueError(
            f"{where} is the manifest for sector {data['sector']}, not {sector}"
        )
    return {int(k): v for k, v in data["targets"].items()}
