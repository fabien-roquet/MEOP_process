from __future__ import annotations

from pathlib import Path
from typing import Iterable

from ..catalog.filenames import fname_prof
from ..catalog.tables import non_empty_attributes, read_csv_rows
from ..data.layout import resolve_table_path
from ..models import MeopConfig


def _write_attrs_in_place(path: Path, attributes: dict[str, str]) -> None:
    try:
        import netCDF4 as nc4
        opener = nc4.Dataset
        write_attr = lambda handle, key, value: handle.setncattr(key, value)
    except ImportError:
        import h5netcdf
        opener = h5netcdf.File
        write_attr = lambda handle, key, value: setattr(handle.attrs, key, value)

    with opener(path, "a") as handle:
        for key, value in attributes.items():
            write_attr(handle, key, value)


def update_metadata_from_table(
    config: MeopConfig,
    *,
    deployment: str = "",
    smru_name: str = "",
    modes: Iterable[str] = ("lr0", "hr0", "fr0"),
) -> list[Path]:
    """Patch file-level attributes from ``table_meta.csv``.

    Canonical metadata tables now live under ``data/tables`` and are synchronized with the
    the current pure-Python workflow can continue to run unchanged.
    """

    table_path = resolve_table_path(config, "table_meta.csv", required=False)
    rows = read_csv_rows(table_path)
    if not rows:
        return []

    selected_rows: list[dict[str, str]] = []
    if smru_name:
        selected_rows = [row for row in rows if row.get("smru_platform_code", "") == smru_name]
    elif deployment:
        selected_rows = [
            row
            for row in rows
            if row.get("smru_platform_code", "").split("-")[0] == deployment
        ]
    else:
        selected_rows = rows

    updated_files: list[Path] = []
    for row in selected_rows:
        tag = row.get("smru_platform_code", "")
        if not tag:
            continue
        attrs = non_empty_attributes(row, exclude=("smru_platform_code",))
        for qf in modes:
            path = fname_prof(tag, qf=qf, config=config)
            if not path.exists():
                continue
            _write_attrs_in_place(path, attrs)
            updated_files.append(path)
    return updated_files
