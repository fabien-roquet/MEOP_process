"""Load CORA reference dataset tiles for a given bounding box.

CORA tiles are stored as individual NetCDF files following the naming convention::

    CORA_lon{LON}E_lat{LAT}N.nc   (East/North)
    CORA_lon{LON}W_lat{LAT}S.nc   (West/South)

where LON and LAT are rounded to the nearest 10°.  Archives in use by MEOP contain
both zero-padded (``lon040W``) and minimally padded (``lon40W``) longitude fields.
The loader therefore inventories and parses the files that actually exist rather than
constructing one assumed spelling.  Each tile covers a 10°×10° cell and contains
variables TEMP, PSAL on a (N_PROF, PRES_GRID) grid together with per-profile
LATITUDE, LONGITUDE, and JULD.
"""

from __future__ import annotations

import logging
from pathlib import Path
import re
from typing import Callable, Iterable

import numpy as np

logger = logging.getLogger(__name__)


_TILE_FILENAME_PATTERN = re.compile(
    r"^CORA_lon(?P<lon>\d+)(?P<east_west>[EW])_lat(?P<lat>\d+)(?P<north_south>[NS])\.nc$"
)


def _tile_lon_name(lon_deg10: int) -> str:
    """Return the legacy zero-padded longitude component of a tile filename."""
    if lon_deg10 >= 0:
        return f"lon{lon_deg10:03d}E"
    return f"lon{abs(lon_deg10):03d}W"


def _tile_lat_name(lat_deg10: int) -> str:
    """Return the latitude component of a tile filename (e.g. 'lat60S', 'lat00N')."""
    if lat_deg10 >= 0:
        return f"lat{lat_deg10:02d}N"
    return f"lat{abs(lat_deg10):02d}S"


def _tile_filename(lon_deg10: int, lat_deg10: int) -> str:
    """Return the legacy canonical tile name used by older tests and archives."""
    return f"CORA_{_tile_lon_name(lon_deg10)}_{_tile_lat_name(lat_deg10)}.nc"


def _parse_tile_filename(filename: str) -> tuple[int, int] | None:
    """Return the signed 10-degree cell encoded in *filename*.

    The numeric fields may be minimally or zero padded.  ``None`` is returned for
    filenames outside the CORA tile convention.
    """
    match = _TILE_FILENAME_PATTERN.match(filename)
    if match is None:
        return None
    lon = int(match.group("lon"))
    lat = int(match.group("lat"))
    if match.group("east_west") == "W":
        lon = -lon
    if match.group("north_south") == "S":
        lat = -lat
    return lon, lat


def discover_cora_tiles(cora_dir: Path | str) -> dict[tuple[int, int], Path]:
    """Inventory CORA tiles by their signed ``(longitude, latitude)`` cell.

    Discovery deliberately uses actual filenames so that ``lon40W`` and
    ``lon040W`` resolve to the same cell.  If an archive contains two spellings for
    one cell, the shortest filename is selected deterministically and a warning is
    logged; a tile is never counted twice.
    """
    root = Path(cora_dir)
    candidates: dict[tuple[int, int], list[Path]] = {}
    if not root.exists():
        return {}
    for path in sorted(root.glob("CORA_lon*_lat*.nc")):
        cell = _parse_tile_filename(path.name)
        if cell is None:
            logger.warning("Ignoring unrecognised CORA tile filename %s", path.name)
            continue
        candidates.setdefault(cell, []).append(path)

    inventory: dict[tuple[int, int], Path] = {}
    for cell, paths in sorted(candidates.items()):
        ordered = sorted(paths, key=lambda item: (len(item.name), item.name))
        inventory[cell] = ordered[0]
        if len(ordered) > 1:
            logger.warning(
                "Multiple CORA tiles resolve to cell %s; using %s and ignoring %s",
                cell,
                ordered[0].name,
                ", ".join(path.name for path in ordered[1:]),
            )
    return inventory


def _tile_cells_for_bbox(
    lon_min: float, lon_max: float, lat_min: float, lat_max: float
) -> list[tuple[int, int]]:
    """Return the 10-degree cells intersecting the bounding box."""
    import math

    lon_lo = int(math.floor(lon_min / 10.0)) * 10
    lon_hi = int(math.floor(lon_max / 10.0)) * 10
    lat_lo = int(math.floor(lat_min / 10.0)) * 10
    lat_hi = int(math.floor(lat_max / 10.0)) * 10

    cells: list[tuple[int, int]] = []
    lon = lon_lo
    while lon <= lon_hi:
        lat = lat_lo
        while lat <= lat_hi:
            cells.append((lon, lat))
            lat += 10
        lon += 10
    return cells


def _tiles_for_bbox(
    lon_min: float, lon_max: float, lat_min: float, lat_max: float
) -> list[str]:
    """Return the list of tile filenames that intersect the bounding box.

    Tile centres are on multiples of 10°.  A tile (lon10, lat10) covers
    [lon10, lon10+10) × [lat10, lat10+10).
    """
    cells = _tile_cells_for_bbox(lon_min, lon_max, lat_min, lat_max)
    return [_tile_filename(lon, lat) for lon, lat in cells]


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
    tile_cells = _tile_cells_for_bbox(lon_min, lon_max, lat_min, lat_max)

    def bbox_mask(
        _cell: tuple[int, int], tile_lat: np.ndarray, tile_lon: np.ndarray
    ) -> np.ndarray:
        return (
            (tile_lat >= lat_min)
            & (tile_lat <= lat_max)
            & (tile_lon >= lon_min)
            & (tile_lon <= lon_max)
        )

    return _load_selected_cora_tiles(cora_dir, tile_cells, mask_builder=bbox_mask)


def _empty_cora_data(pres: np.ndarray | None = None) -> dict[str, np.ndarray]:
    pressure = np.asarray(pres, dtype=np.float64) if pres is not None else np.empty(0)
    return {
        "lat": np.empty(0, dtype=np.float64),
        "lon": np.empty(0, dtype=np.float64),
        "juld": np.empty(0, dtype=np.float64),
        "temp": np.empty((0, pressure.size), dtype=np.float64),
        "psal": np.empty((0, pressure.size), dtype=np.float64),
        "pres": pressure,
    }


def _load_selected_cora_tiles(
    cora_dir: Path | str,
    tile_cells: Iterable[tuple[int, int]],
    *,
    mask_builder: Callable[[tuple[int, int], np.ndarray, np.ndarray], np.ndarray],
) -> dict[str, np.ndarray]:
    """Load selected profile rows from explicit CORA cells.

    Latitude and longitude are intentionally read first.  TEMP and PSAL are indexed
    on ``N_PROF`` before they are materialised, preserving the archive's 128-profile
    NetCDF chunking instead of reading an entire tile and masking it afterwards.
    """
    try:
        import xarray as xr  # type: ignore
    except ImportError as exc:
        raise ImportError("xarray is required to load CORA tiles") from exc

    tile_inventory = discover_cora_tiles(cora_dir)

    lat_parts: list[np.ndarray] = []
    lon_parts: list[np.ndarray] = []
    juld_parts: list[np.ndarray] = []
    temp_parts: list[np.ndarray] = []
    psal_parts: list[np.ndarray] = []
    pres: np.ndarray | None = None

    for cell in sorted(set(tile_cells)):
        path = tile_inventory.get(cell)
        if path is None:
            continue
        try:
            ds_context = xr.open_dataset(path, decode_times=False)
        except Exception as exc:
            logger.warning("Skipping unreadable CORA tile %s: %s", path.name, exc)
            continue
        with ds_context as ds:
            required = {"LATITUDE", "LONGITUDE", "JULD", "TEMP", "PSAL", "PRES_GRID"}
            missing = sorted(required.difference(ds.variables))
            if missing:
                logger.warning("Skipping malformed CORA tile %s: missing %s", path.name, ", ".join(missing))
                continue
            tile_lat = np.asarray(ds["LATITUDE"].values, dtype=np.float64)
            tile_lon = np.asarray(ds["LONGITUDE"].values, dtype=np.float64)
            mask = np.asarray(mask_builder(cell, tile_lat, tile_lon), dtype=bool)
            if mask.shape != tile_lat.shape:
                raise ValueError(f"CORA profile selector returned the wrong shape for {path.name}")
            profile_indices = np.flatnonzero(mask)
            if profile_indices.size == 0:
                continue
            selected = ds.isel(N_PROF=profile_indices)
            lat_parts.append(tile_lat[profile_indices])
            lon_parts.append(tile_lon[profile_indices])
            juld_parts.append(np.asarray(selected["JULD"].values, dtype=np.float64))
            tile_temp = np.asarray(selected["TEMP"].values, dtype=np.float64)
            tile_psal = np.asarray(selected["PSAL"].values, dtype=np.float64)
            temp_parts.append(tile_temp)
            psal_parts.append(tile_psal)
            if pres is None:
                pres = np.asarray(ds["PRES_GRID"].values, dtype=np.float64)

    if not lat_parts:
        return _empty_cora_data(pres)

    return {
        "lat": np.concatenate(lat_parts),
        "lon": np.concatenate(lon_parts),
        "juld": np.concatenate(juld_parts),
        "temp": np.concatenate(temp_parts, axis=0),
        "psal": np.concatenate(psal_parts, axis=0),
        "pres": pres,
    }


def _point_longitude_intervals(longitude: float, margin: float) -> tuple[tuple[float, float], ...]:
    centre = ((float(longitude) + 180.0) % 360.0) - 180.0
    if margin >= 180.0:
        return ((-180.0, 180.0),)
    lower = centre - margin
    upper = centre + margin
    if lower < -180.0:
        return ((lower + 360.0, 180.0), (-180.0, upper))
    if upper > 180.0:
        return ((lower, 180.0), (-180.0, upper - 360.0))
    return ((lower, upper),)


def _track_cell_targets(
    latitudes: np.ndarray,
    longitudes: np.ndarray,
    *,
    margin: float,
) -> dict[tuple[int, int], np.ndarray]:
    """Map buffered track cells to the target positions that require each cell."""
    if margin < 0.0:
        raise ValueError("track margin must be non-negative")
    lat = np.asarray(latitudes, dtype=np.float64).reshape(-1)
    lon = np.asarray(longitudes, dtype=np.float64).reshape(-1)
    if lat.shape != lon.shape:
        raise ValueError("target latitude and longitude arrays must have the same shape")
    valid = np.isfinite(lat) & np.isfinite(lon) & (lat >= -90.0) & (lat <= 90.0)
    lat = lat[valid]
    lon = ((lon[valid] + 180.0) % 360.0) - 180.0
    mapping: dict[tuple[int, int], list[int]] = {}
    for index, (latitude, longitude) in enumerate(zip(lat, lon, strict=True)):
        lat_min = max(-90.0, float(latitude) - margin)
        lat_max = min(90.0, float(latitude) + margin)
        for lon_min, lon_max in _point_longitude_intervals(float(longitude), margin):
            for cell in _tile_cells_for_bbox(lon_min, lon_max, lat_min, lat_max):
                mapping.setdefault(cell, []).append(index)
    return {
        cell: np.asarray(sorted(set(indices)), dtype=np.int64)
        for cell, indices in mapping.items()
    }


def cora_cells_for_track(
    latitudes: np.ndarray,
    longitudes: np.ndarray,
    *,
    margin: float = 5.0,
) -> tuple[tuple[int, int], ...]:
    """Return only cells intersecting the buffered target positions, not their envelope."""
    return tuple(sorted(_track_cell_targets(latitudes, longitudes, margin=margin)))


def _profiles_near_track(
    reference_latitudes: np.ndarray,
    reference_longitudes: np.ndarray,
    target_latitudes: np.ndarray,
    target_longitudes: np.ndarray,
    *,
    margin: float,
) -> np.ndarray:
    """Select reference positions inside any target profile's local degree window."""
    ref_lat = np.asarray(reference_latitudes, dtype=np.float64)
    ref_lon = np.asarray(reference_longitudes, dtype=np.float64)
    target_lat = np.asarray(target_latitudes, dtype=np.float64)
    target_lon = np.asarray(target_longitudes, dtype=np.float64)
    selected = np.zeros(ref_lat.shape, dtype=bool)
    reference_chunk = 4096
    target_chunk = 128
    for ref_start in range(0, ref_lat.size, reference_chunk):
        ref_stop = min(ref_start + reference_chunk, ref_lat.size)
        chunk_selected = np.zeros(ref_stop - ref_start, dtype=bool)
        chunk_lat = ref_lat[ref_start:ref_stop, None]
        chunk_lon = ref_lon[ref_start:ref_stop, None]
        for target_start in range(0, target_lat.size, target_chunk):
            target_stop = min(target_start + target_chunk, target_lat.size)
            latitude_close = np.abs(
                chunk_lat - target_lat[None, target_start:target_stop]
            ) <= margin
            longitude_delta = np.abs(
                (
                    chunk_lon
                    - target_lon[None, target_start:target_stop]
                    + 180.0
                )
                % 360.0
                - 180.0
            )
            chunk_selected |= np.any(
                latitude_close & (longitude_delta <= margin), axis=1
            )
            if np.all(chunk_selected):
                break
        selected[ref_start:ref_stop] = chunk_selected
    return selected


def load_cora_track(
    cora_dir: Path | str,
    *,
    latitudes: np.ndarray,
    longitudes: np.ndarray,
    margin: float = 5.0,
) -> dict[str, np.ndarray]:
    """Load CORA profiles within local windows along an accepted target track."""
    lat = np.asarray(latitudes, dtype=np.float64).reshape(-1)
    lon = np.asarray(longitudes, dtype=np.float64).reshape(-1)
    valid = np.isfinite(lat) & np.isfinite(lon) & (lat >= -90.0) & (lat <= 90.0)
    lat = lat[valid]
    lon = ((lon[valid] + 180.0) % 360.0) - 180.0
    cell_targets = _track_cell_targets(lat, lon, margin=margin)
    if not cell_targets:
        return _empty_cora_data()

    def track_mask(
        cell: tuple[int, int], tile_lat: np.ndarray, tile_lon: np.ndarray
    ) -> np.ndarray:
        indices = cell_targets[cell]
        return _profiles_near_track(
            tile_lat,
            tile_lon,
            lat[indices],
            lon[indices],
            margin=margin,
        )

    return _load_selected_cora_tiles(cora_dir, cell_targets, mask_builder=track_mask)
