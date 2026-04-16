from __future__ import annotations

from pathlib import Path

import xarray as xr


DEFAULT_FORMAT = "NETCDF4_CLASSIC"


def save_dataset_with_compression(
    dataset: xr.Dataset,
    path: Path,
    *,
    format: str = DEFAULT_FORMAT,
) -> None:
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