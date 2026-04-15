from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import xarray as xr


DEFAULT_REFERENCE_DATE = datetime(1950, 1, 1, tzinfo=timezone.utc)


def open_meop_netcdf(path: str | Path, *, decode_times: bool = False) -> xr.Dataset:
    """Open a MEOP netCDF product with a tolerant engine fallback sequence.

    The historical toolbox writes a mix of netCDF3 and netCDF4-compatible outputs,
    so the helper tries the most common backends in order. The returned dataset is
    fully loaded into memory, which keeps plotting/test code simple and avoids open
    file handles leaking across workflow steps.
    """

    dataset_path = Path(path)
    errors: list[Exception] = []
    for engine in ("h5netcdf", "scipy", None):
        try:
            if engine is None:
                return xr.load_dataset(dataset_path, decode_times=decode_times)
            return xr.load_dataset(dataset_path, engine=engine, decode_times=decode_times)
        except Exception as exc:  # pragma: no cover - backend availability varies by environment
            errors.append(exc)
    raise OSError(f"Unable to open MEOP netCDF file: {dataset_path}") from errors[-1]


def decode_text(value: object) -> str:
    """Decode MEOP text values stored as bytes/object scalars or arrays."""

    if isinstance(value, bytes):
        return value.decode("utf-8", errors="ignore").strip()
    if isinstance(value, str):
        return value.strip()
    if np.isscalar(value):
        return str(value).strip()

    array = np.asarray(value)
    if array.ndim == 0:
        return decode_text(array.item())
    if array.dtype.kind in {"S", "U", "O"}:
        flattened = [decode_text(item) for item in array.ravel()]
        return "".join(flattened).strip()
    return str(value).strip()


def decode_text_vector(values: object) -> list[str]:
    array = np.asarray(values)
    if array.ndim == 0:
        return [decode_text(array.item())]
    return [decode_text(item) for item in array]


def juld_to_datetime(values: Iterable[float] | np.ndarray, *, reference: datetime = DEFAULT_REFERENCE_DATE) -> np.ndarray:
    """Convert JULD-style day offsets to timezone-aware UTC datetimes."""

    out: list[datetime | None] = []
    for value in np.asarray(values, dtype=float).ravel():
        if np.isfinite(value):
            out.append(reference + timedelta(days=float(value)))
        else:
            out.append(None)
    return np.asarray(out, dtype=object).reshape(np.asarray(values).shape)
