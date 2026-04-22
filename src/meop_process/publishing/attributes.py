"""Update global attributes on published NC files."""
from __future__ import annotations

import re
from datetime import datetime, timezone
from pathlib import Path

import xarray as xr

from ..processing.netcdf import save_dataset_with_compression


_DATE_FORMAT = "%Y%m%dT%H%M%SZ"
_DATA_TYPE = "Argo profile"


def _today_str() -> str:
    return datetime.now(tz=timezone.utc).strftime(_DATE_FORMAT)


def _smru_name_from_path(path: Path) -> str:
    name = path.name
    # <smru_name>_all_prof.nc  or  <smru_name>_<qf>_prof.nc
    return name.split("_")[0]


def update_global_attributes(
    output_dir: Path,
    *,
    version: str = "",
    verbose: bool = False,
) -> list[Path]:
    """Patch standard global attributes on all *_prof.nc files in output_dir.

    Sets data_type, date_update, platform_code, and optionally version_database.
    Returns list of patched file paths.
    """
    patched: list[Path] = []
    date_str = _today_str()
    nc_files = sorted(output_dir.glob("*_prof.nc"))
    if not nc_files:
        nc_files = sorted(output_dir.rglob("*_prof.nc"))

    for path in nc_files:
        try:
            with xr.open_dataset(path) as ds:
                result = ds.load().copy(deep=True)

            new_attrs = dict(result.attrs)
            new_attrs["data_type"] = _DATA_TYPE
            new_attrs["date_update"] = date_str
            smru = _smru_name_from_path(path)
            if "platform_code" not in new_attrs or not str(new_attrs.get("platform_code", "")).strip():
                new_attrs["platform_code"] = smru
            if version:
                new_attrs["version_database"] = version
            result.attrs = new_attrs

            save_dataset_with_compression(result, path)
            patched.append(path)
            if verbose:
                print(f"  updated attributes: {path.name}")
        except Exception as exc:
            if verbose:
                print(f"  WARNING: could not update {path.name}: {exc}")

    return patched
