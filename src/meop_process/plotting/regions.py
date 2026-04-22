"""Geographic region labelling using regionmask AR6 ocean basins.

Provides :func:`label_region` (single position) and :func:`label_regions`
(arrays) that map (lon, lat) coordinates to simplified ocean basin names.
Falls back gracefully when ``regionmask`` or ``scipy`` are not installed.
"""
from __future__ import annotations

import numpy as np

try:
    import regionmask  # type: ignore
    from scipy.interpolate import RegularGridInterpolator  # type: ignore

    _AVAILABLE = True
except Exception:
    _AVAILABLE = False


# Mapping from AR6 region names to simplified ocean basin labels.
# Land/coastal AR6 regions are mapped to the nearest ocean basin, which is
# useful when a tag median position falls slightly inland.
_AR6_TO_OCEAN: dict[str, str] = {
    # Arctic / high latitude
    "Greenland/Iceland": "North Atlantic",
    "Arctic-Ocean": "Arctic Ocean",
    "Russian-Arctic": "Arctic Ocean",
    # North Atlantic
    "N.W.North-America": "North Atlantic",
    "N.E.North-America": "North Atlantic",
    "E.North-America": "North Atlantic",
    "N.Europe": "North Atlantic",
    "West&Central-Europe": "North Atlantic",
    "E.Europe": "North Atlantic",
    "N.Atlantic-Ocean": "North Atlantic",
    # Tropical Atlantic
    "Caribbean": "Tropical Atlantic",
    "N.W.South-America": "Tropical Atlantic",
    "N.South-America": "Tropical Atlantic",
    "N.E.South-America": "Tropical Atlantic",
    "South-American-Monsoon": "Tropical Atlantic",
    "Equatorial.Atlantic-Ocean": "Tropical Atlantic",
    # South Atlantic
    "S.W.South-America": "South Atlantic",
    "S.E.South-America": "South Atlantic",
    "S.South-America": "South Atlantic",
    "S.Atlantic-Ocean": "South Atlantic",
    # North Pacific
    "W.North-America": "North Pacific",
    "C.North-America": "North Pacific",
    "N.Central-America": "North Pacific",
    "S.Central-America": "North Pacific",
    "Russian-Far-East": "North Pacific",
    "E.Siberia": "North Pacific",
    "W.Siberia": "North Pacific",
    "N.Pacific-Ocean": "North Pacific",
    # Tropical Pacific
    "Equatorial.Pacific-Ocean": "Tropical Pacific",
    # South Pacific
    "E.Australia": "South Pacific",
    "S.Australia": "South Pacific",
    "New-Zealand": "South Pacific",
    "S.Pacific-Ocean": "South Pacific",
    # Indian Ocean
    "Tibetan-Plateau": "Indian Ocean",
    "Arabian-Peninsula": "Indian Ocean",
    "S.Asia": "Indian Ocean",
    "S.E.Asia": "Indian Ocean",
    "N.Australia": "Indian Ocean",
    "C.Australia": "Indian Ocean",
    "Madagascar": "Indian Ocean",
    "S.Eastern-Africa": "Indian Ocean",
    "N.Eastern-Africa": "Indian Ocean",
    "Arabian-Sea": "Indian Ocean",
    "Bay-of-Bengal": "Indian Ocean",
    "Equatorial.Indic-Ocean": "Indian Ocean",
    "S.Indic-Ocean": "Indian Ocean",
    # Africa (map to nearest ocean basin)
    "Sahara": "North Atlantic",
    "Western-Africa": "Tropical Atlantic",
    "Central-Africa": "Tropical Atlantic",
    "W.Southern-Africa": "South Atlantic",
    "E.Southern-Africa": "Indian Ocean",
    # Central/East Asia
    "W.C.Asia": "Indian Ocean",
    "E.C.Asia": "North Pacific",
    "E.Asia": "North Pacific",
    "Mediterranean": "North Atlantic",
    # Southern Ocean and Antarctica
    "E.Antarctica": "Southern Ocean",
    "W.Antarctica": "Southern Ocean",
    "Southern-Ocean": "Southern Ocean",
}

_UNKNOWN = "Unknown"
_interpolator: object | None = None


def _build_interpolator():
    """Build and cache the RegularGridInterpolator from the AR6 mask grid."""
    global _interpolator  # noqa: PLW0603
    if _interpolator is not None:
        return _interpolator
    if not _AVAILABLE:
        return None
    basins = regionmask.defined_regions.ar6.all
    lons = np.arange(-179.5, 180.0)
    lats = np.arange(-89.5, 90.0)
    # mask() returns (nlat, nlev) — transpose so the grid is (nlon, nlat)
    mask_grid = basins.mask(lons, lats).transpose().values
    _interpolator = RegularGridInterpolator(
        (lons, lats),
        mask_grid,
        method="nearest",
        bounds_error=False,
        fill_value=np.nan,
    )
    return _interpolator


def label_regions(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """Return simplified ocean region names for arrays of (lon, lat) positions.

    Parameters
    ----------
    lon:
        Longitude values in degrees East (−180 … 180).
    lat:
        Latitude values in degrees North (−90 … 90).

    Returns
    -------
    numpy.ndarray of str with the same shape as *lon*/*lat*.
    Entries are one of the simplified ocean basin names (e.g. ``"Southern
    Ocean"``, ``"North Atlantic"``) or ``"Unknown"`` when the position cannot
    be classified.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    original_shape = lon.shape
    lon_flat = lon.ravel()
    lat_flat = lat.ravel()

    interp = _build_interpolator()
    if interp is None:
        return np.full(original_shape, _UNKNOWN, dtype=object)

    points = np.stack([lon_flat, lat_flat], axis=-1)
    indices = interp(points)

    basins = regionmask.defined_regions.ar6.all
    names = basins.names

    results: list[str] = []
    for idx in indices:
        if np.isnan(idx):
            results.append(_UNKNOWN)
        else:
            ar6_name = names[int(round(float(idx)))] if 0 <= int(round(float(idx))) < len(names) else _UNKNOWN
            results.append(_AR6_TO_OCEAN.get(ar6_name, ar6_name))

    return np.array(results, dtype=object).reshape(original_shape)


def label_region(lon: float, lat: float) -> str:
    """Return the simplified ocean region name for a single (lon, lat) position.

    Returns ``"Unknown"`` when the position cannot be classified or when
    ``regionmask`` is not installed.
    """
    result = label_regions(np.asarray([lon]), np.asarray([lat]))
    return str(result[0])
