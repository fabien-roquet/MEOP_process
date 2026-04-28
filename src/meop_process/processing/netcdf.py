from __future__ import annotations

import logging
from pathlib import Path

import xarray as xr

logger = logging.getLogger(__name__)

DEFAULT_FORMAT = "NETCDF4_CLASSIC"


def _sanitize_string(value: object) -> object:
    """Remove surrogate characters from strings that cannot be encoded as UTF-8.
    
    Surrogate characters (U+D800–U+DFFF) result from misinterpreting Latin-1 bytes as UTF-8.
    They are invalid in Python strings and cause UnicodeEncodeError when writing to NetCDF.
    """
    if not isinstance(value, str):
        return value
    
    # Check for surrogate characters
    try:
        value.encode('utf-8')
        return value
    except UnicodeEncodeError:
        # Replace surrogates with replacement character
        sanitized = value.encode('utf-8', errors='replace').decode('utf-8', errors='replace')
        logger.warning(f"Removed invalid UTF-8 surrogate characters from attribute: {value!r} -> {sanitized!r}")
        return sanitized


def _sanitize_dataset_strings(dataset: xr.Dataset) -> xr.Dataset:
    """Sanitize all string-valued attributes and data in a dataset to remove surrogates."""
    # Sanitize global attributes
    if dataset.attrs:
        dataset.attrs = {k: _sanitize_string(v) for k, v in dataset.attrs.items()}
    
    # Sanitize variable attributes
    for var_name in dataset.data_vars:
        var = dataset[var_name]
        if var.attrs:
            var.attrs = {k: _sanitize_string(v) for k, v in var.attrs.items()}
    
    for coord_name in dataset.coords:
        coord = dataset[coord_name]
        if coord.attrs:
            coord.attrs = {k: _sanitize_string(v) for k, v in coord.attrs.items()}
    
    return dataset


def save_dataset_with_compression(
    dataset: xr.Dataset,
    path: Path,
    *,
    format: str = DEFAULT_FORMAT,
) -> None:
    # Sanitize strings before writing to avoid UTF-8 encoding errors
    dataset = _sanitize_dataset_strings(dataset)
    
    encoding = {
        name: {"zlib": True, "complevel": 5}
        for name, variable in dataset.data_vars.items()
        if getattr(variable.dtype, "kind", "") not in {"O", "S", "U"}
    }
    try:
        import netCDF4  # noqa: F401

        dataset.to_netcdf(path, engine="netcdf4", format=format, encoding=encoding)
    except (ImportError, RuntimeError, ValueError):
        dataset.to_netcdf(path, engine="h5netcdf", encoding=encoding)