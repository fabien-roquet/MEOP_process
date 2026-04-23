"""Load CORA reference dataset tiles for a given bounding box.

CORA tiles are stored as individual NetCDF files following the naming convention::

    CORA_lon{LON}E_lat{LAT}N.nc   (East/North)
    CORA_lon{LON}W_lat{LAT}S.nc   (West/South)

where LON and LAT are rounded to the nearest 10°.  Each tile covers a 10°×10° cell
and contains variables TEMP, PSAL on a (N_PROF, PRES_GRID) grid together with
per-profile LATITUDE, LONGITUDE, and JULD.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np


def _tile_lon_name(lon_deg10: int) -> str:
    """Return the longitude component of a tile filename (e.g. 'lon00E', 'lon10W')."""
    if lon_deg10 >= 0:
        return f"lon{lon_deg10:03d}E"
    return f"lon{abs(lon_deg10):03d}W"


def _tile_lat_name(lat_deg10: int) -> str:
    """Return the latitude component of a tile filename (e.g. 'lat60S', 'lat00N')."""
    if lat_deg10 >= 0:
        return f"lat{lat_deg10:02d}N"
    return f"lat{abs(lat_deg10):02d}S"


def _tile_filename(lon_deg10: int, lat_deg10: int) -> str:
    return f"CORA_{_tile_lon_name(lon_deg10)}_{_tile_lat_name(lat_deg10)}.nc"


def _tiles_for_bbox(
    lon_min: float, lon_max: float, lat_min: float, lat_max: float
) -> list[str]:
    """Return the list of tile filenames that intersect the bounding box.

    Tile centres are on multiples of 10°.  A tile (lon10, lat10) covers
    [lon10, lon10+10) × [lat10, lat10+10).
    """
    import math

    lon_lo = int(math.floor(lon_min / 10.0)) * 10
    lon_hi = int(math.floor(lon_max / 10.0)) * 10
    lat_lo = int(math.floor(lat_min / 10.0)) * 10
    lat_hi = int(math.floor(lat_max / 10.0)) * 10

    filenames: list[str] = []
    lon = lon_lo
    while lon <= lon_hi:
        lat = lat_lo
        while lat <= lat_hi:
            filenames.append(_tile_filename(lon, lat))
            lat += 10
        lon += 10
    return filenames


def load_cora_tiles(
    cora_dir: Path | str,
    *,
    lon_min: float,
    lon_max: float,
    lat_min: float,
    lat_max: float,
) -> dict[str, np.ndarray]:
    """Return merged CORA profiles that intersect *[lon_min, lon_max] × [lat_min, lat_max]*.

    Returns a dict with keys:
        ``lat``   — 1-D array, shape (N,)
        ``lon``   — 1-D array, shape (N,)
        ``juld``  — 1-D array, shape (N,)
        ``temp``  — 2-D array, shape (N, PRES_GRID)
        ``psal``  — 2-D array, shape (N, PRES_GRID)
        ``pres``  — 1-D array, shape (PRES_GRID,)  — shared pressure grid

    Only profiles whose *LATITUDE* and *LONGITUDE* lie within the bounding box
    are returned.  Returns arrays with N=0 if no tiles exist or no profiles qualify.
    """
    try:
        import xarray as xr  # type: ignore
    except ImportError as exc:
        raise ImportError("xarray is required to load CORA tiles") from exc

    cora_dir = Path(cora_dir)
    tile_files = _tiles_for_bbox(lon_min, lon_max, lat_min, lat_max)

    lat_parts: list[np.ndarray] = []
    lon_parts: list[np.ndarray] = []
    juld_parts: list[np.ndarray] = []
    temp_parts: list[np.ndarray] = []
    psal_parts: list[np.ndarray] = []
    pres: np.ndarray | None = None

    for fname in tile_files:
        path = cora_dir / fname
        if not path.exists():
            continue
        with xr.open_dataset(path, decode_times=False) as ds:
            tile_lat = np.asarray(ds["LATITUDE"].values, dtype=np.float64)
            tile_lon = np.asarray(ds["LONGITUDE"].values, dtype=np.float64)
            mask = (
                (tile_lat >= lat_min)
                & (tile_lat <= lat_max)
                & (tile_lon >= lon_min)
                & (tile_lon <= lon_max)
            )
            if not np.any(mask):
                continue
            lat_parts.append(tile_lat[mask])
            lon_parts.append(tile_lon[mask])
            juld_parts.append(np.asarray(ds["JULD"].values, dtype=np.float64)[mask])
            tile_temp = np.asarray(ds["TEMP"].values, dtype=np.float64)[mask]
            tile_psal = np.asarray(ds["PSAL"].values, dtype=np.float64)[mask]
            temp_parts.append(tile_temp)
            psal_parts.append(tile_psal)
            if pres is None:
                pres = np.asarray(ds["PRES_GRID"].values, dtype=np.float64)

    if not lat_parts:
        n_pres = int(pres.size) if pres is not None else 0
        return {
            "lat": np.empty(0, dtype=np.float64),
            "lon": np.empty(0, dtype=np.float64),
            "juld": np.empty(0, dtype=np.float64),
            "temp": np.empty((0, n_pres), dtype=np.float64),
            "psal": np.empty((0, n_pres), dtype=np.float64),
            "pres": pres if pres is not None else np.empty(0, dtype=np.float64),
        }

    return {
        "lat": np.concatenate(lat_parts),
        "lon": np.concatenate(lon_parts),
        "juld": np.concatenate(juld_parts),
        "temp": np.concatenate(temp_parts, axis=0),
        "psal": np.concatenate(psal_parts, axis=0),
        "pres": pres,
    }
