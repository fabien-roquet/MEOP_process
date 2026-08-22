"""Prepare CTD and Argo CORA source files as MEOP-compatible 10-degree tiles."""

from __future__ import annotations

import argparse
from collections import Counter, OrderedDict
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable, Sequence

import numpy as np

from .cora_download import CORA_DATASET_ID, EASYCORA_DATASET_ID, PRODUCT_ID


_PROFILE_NAME = re.compile(
    r"_(?P<date>\d{8})_PR_(?P<file_class>[A-Z0-9]{2})\.nc$",
    re.IGNORECASE,
)
_GOOD_QC = (1, 2)
_TIME_UNITS = "days since 1950-01-01 00:00:00 UTC"


class SourceFormatError(ValueError):
    """Raised when a source file cannot be classified or safely interpreted."""


@dataclass(frozen=True)
class SourcePolicy:
    reference_kind: str
    dataset_id: str
    file_class: str
    required_probe_type: int
    allow_missing_probe_type: bool


@dataclass
class ProfileBatch:
    latitude: list[float]
    longitude: list[float]
    juld: list[float]
    temp: list[np.ndarray]
    psal: list[np.ndarray]
    probe_type: list[int]
    source_profile_index: list[int]
    platform_number: list[str]
    cycle_number: list[str]
    dc_reference: list[str]
    data_mode: list[str]
    wmo_inst_type: list[str]
    source_file: list[str]
    source_file_class: list[str]
    reference_kind: list[str]
    temp_source: list[str]
    psal_source: list[str]
    pres_source: list[str]
    source_dataset: list[str]

    @classmethod
    def empty(cls) -> "ProfileBatch":
        return cls(*([] for _ in range(19)))

    def __len__(self) -> int:
        return len(self.latitude)


def _source_policy(path: Path) -> SourcePolicy:
    match = _PROFILE_NAME.search(path.name)
    if match is None:
        raise SourceFormatError(f"not a recognised CORA daily profile filename: {path.name}")
    file_class = match.group("file_class").upper()
    components = {item.casefold() for item in path.parts}
    if "easycora_argo" in components:
        if file_class != "PF":
            raise SourceFormatError(f"non-PF file found in easycora_argo: {path}")
        return SourcePolicy("ARGO", EASYCORA_DATASET_ID, file_class, 5, True)
    if "cora_argo" in components:
        if file_class != "PF":
            raise SourceFormatError(f"non-PF file found in cora_argo: {path}")
        return SourcePolicy("ARGO", CORA_DATASET_ID, file_class, 5, False)
    if "cora_ctd_candidates" in components:
        if file_class not in {"CT", "OC", "TE"}:
            raise SourceFormatError(f"unexpected CTD-candidate class {file_class}: {path}")
        return SourcePolicy("CTD", CORA_DATASET_ID, file_class, 2, False)

    # Permit direct use of a manually organised source tree, but retain the same
    # conservative classifications as the downloader.
    if file_class == "PF":
        return SourcePolicy("ARGO", CORA_DATASET_ID, file_class, 5, False)
    if file_class in {"CT", "OC", "TE"}:
        return SourcePolicy("CTD", CORA_DATASET_ID, file_class, 2, False)
    raise SourceFormatError(f"file class {file_class} is outside the CTD/Argo selection: {path}")


def _normalise_qc(values: Any) -> np.ndarray:
    array = np.asarray(values)
    result = np.full(array.shape, -127, dtype=np.int16)
    if array.dtype.kind in "iuf":
        finite = np.isfinite(array)
        result[finite] = array[finite].astype(np.int16)
        return result
    flattened = array.astype("U").reshape(-1)
    parsed = result.reshape(-1)
    for index, raw in enumerate(flattened):
        text = str(raw).strip().replace("b'", "").replace("'", "")
        try:
            parsed[index] = int(text)
        except ValueError:
            continue
    return result


def _profile_vector(data_array: Any, n_prof: int) -> np.ndarray:
    if "N_PROF" not in data_array.dims:
        raise SourceFormatError(f"{data_array.name} has no N_PROF dimension")
    axis = data_array.dims.index("N_PROF")
    values = np.asarray(data_array.values)
    values = np.moveaxis(values, axis, 0)
    if values.shape[0] != n_prof:
        raise SourceFormatError(f"{data_array.name} has an inconsistent N_PROF length")
    return values


def _numeric_profile_vector(data_array: Any, n_prof: int) -> np.ndarray:
    values = _profile_vector(data_array, n_prof)
    values = values.reshape(n_prof, -1)
    if values.shape[1] != 1:
        raise SourceFormatError(f"{data_array.name} is not one value per profile")
    return np.asarray(values[:, 0], dtype=np.float64)


def _decode_string(value: Any) -> str:
    array = np.asarray(value)
    if array.ndim == 0:
        item = array.item()
        if isinstance(item, (bytes, np.bytes_)):
            return bytes(item).decode("utf-8", errors="replace").strip(" \x00")
        return str(item).strip(" \x00")
    flat = array.reshape(-1)
    if flat.dtype.kind == "S":
        return b"".join(bytes(item) for item in flat).decode("utf-8", errors="replace").strip(" \x00")
    return "".join(str(item) for item in flat).strip(" \x00")


def _profile_strings(dataset: Any, names: Iterable[str], n_prof: int) -> list[str]:
    for name in names:
        if name not in dataset:
            continue
        values = _profile_vector(dataset[name], n_prof)
        return [_decode_string(values[index]) for index in range(n_prof)]
    return [""] * n_prof


def _profile_qc(dataset: Any, names: Iterable[str], n_prof: int) -> tuple[np.ndarray, str | None]:
    for name in names:
        if name not in dataset:
            continue
        values = _profile_vector(dataset[name], n_prof)
        values = values.reshape(n_prof, -1)
        if values.shape[1] != 1:
            raise SourceFormatError(f"{name} is not one QC value per profile")
        return _normalise_qc(values[:, 0]), name
    return np.full(n_prof, 1, dtype=np.int16), None


def _time_days_since_1950(data_array: Any, n_prof: int) -> np.ndarray:
    try:
        from netCDF4 import date2num, num2date  # type: ignore
    except ImportError as exc:  # pragma: no cover - required project dependency
        raise RuntimeError("netCDF4 is required to decode CORA time variables") from exc
    values = _numeric_profile_vector(data_array, n_prof)
    units = str(data_array.attrs.get("units", "")).strip()
    calendar = str(data_array.attrs.get("calendar", "standard")).strip() or "standard"
    if " since " not in units:
        raise SourceFormatError(f"{data_array.name} has no decodable time units: {units!r}")
    result = np.full(n_prof, np.nan, dtype=np.float64)
    finite = np.isfinite(values)
    if np.any(finite):
        dates = num2date(
            values[finite],
            units=units,
            calendar=calendar,
            only_use_cftime_datetimes=True,
        )
        result[finite] = np.asarray(
            date2num(dates, units=_TIME_UNITS, calendar="standard"), dtype=np.float64
        )
    return result


def _profile_level_values(dataset: Any, name: str, index: int) -> np.ndarray:
    if name not in dataset:
        return np.empty(0, dtype=np.float64)
    variable = dataset[name]
    if "N_PROF" not in variable.dims:
        raise SourceFormatError(f"{name} has no N_PROF dimension")
    values = np.asarray(variable.isel(N_PROF=index).values, dtype=np.float64)
    return np.squeeze(values).reshape(-1)


def _profile_level_qc(dataset: Any, name: str, index: int, length: int) -> np.ndarray:
    if name not in dataset:
        return np.ones(length, dtype=np.int16)
    variable = dataset[name]
    if "N_PROF" not in variable.dims:
        raise SourceFormatError(f"{name} has no N_PROF dimension")
    values = np.asarray(variable.isel(N_PROF=index).values)
    values = np.squeeze(values).reshape(-1)
    qc = _normalise_qc(values)
    if qc.size != length:
        raise SourceFormatError(f"{name} does not align with its measurement variable")
    return qc


def _best_parameter(dataset: Any, base: str, index: int) -> tuple[np.ndarray, np.ndarray, str]:
    adjusted = f"{base}_ADJUSTED"
    if adjusted in dataset:
        values = _profile_level_values(dataset, adjusted, index)
        if np.count_nonzero(np.isfinite(values)):
            return (
                values,
                _profile_level_qc(dataset, f"{adjusted}_QC", index, values.size),
                "adjusted",
            )
    values = _profile_level_values(dataset, base, index)
    return values, _profile_level_qc(dataset, f"{base}_QC", index, values.size), "raw"


def _pressure_parameter(dataset: Any, index: int, latitude: float) -> tuple[np.ndarray, np.ndarray, str]:
    for name in ("PRES",):
        values, qc, source = _best_parameter(dataset, name, index)
        if values.size and np.count_nonzero(np.isfinite(values)):
            return values, qc, source

    # Some CORA profiles contain depth but no pressure. Convert only when necessary
    # and mark that decision explicitly in PRES_SOURCE.
    for name in ("DEPH", "DEPTH"):
        if name not in dataset and f"{name}_ADJUSTED" not in dataset:
            continue
        depth, qc, source = _best_parameter(dataset, name, index)
        try:
            import gsw  # type: ignore
        except ImportError as exc:  # pragma: no cover - required project dependency
            raise RuntimeError("gsw is required to convert CORA depth to pressure") from exc
        pressure = np.asarray(gsw.p_from_z(-depth, latitude), dtype=np.float64)
        return pressure, qc, f"{name.lower()}_{source}_converted"
    return np.empty(0), np.empty(0, dtype=np.int16), "missing"


def _interpolate_profile(
    pressure: np.ndarray,
    pressure_qc: np.ndarray,
    values: np.ndarray,
    values_qc: np.ndarray,
    pressure_grid: np.ndarray,
) -> np.ndarray:
    output = np.full(pressure_grid.size, np.nan, dtype=np.float32)
    if pressure.size != values.size or pressure_qc.size != pressure.size or values_qc.size != values.size:
        return output
    valid = (
        np.isfinite(pressure)
        & np.isfinite(values)
        & (pressure >= 0.0)
        & np.isin(pressure_qc, _GOOD_QC)
        & np.isin(values_qc, _GOOD_QC)
    )
    if np.count_nonzero(valid) < 2:
        return output
    x = pressure[valid]
    y = values[valid]
    order = np.argsort(x, kind="stable")
    x = x[order]
    y = y[order]
    x, unique_indices = np.unique(x, return_index=True)
    y = y[unique_indices]
    if x.size < 2:
        return output
    output[:] = np.interp(pressure_grid, x, y, left=np.nan, right=np.nan).astype(np.float32)
    return output


def _tile_cell(longitude: float, latitude: float) -> tuple[int, int]:
    longitude = ((longitude + 180.0) % 360.0) - 180.0
    lon_cell = int(math.floor(longitude / 10.0)) * 10
    lat_cell = int(math.floor(latitude / 10.0)) * 10
    if latitude == 90.0:
        lat_cell = 80
    return lon_cell, lat_cell


def _tile_filename(cell: tuple[int, int]) -> str:
    longitude, latitude = cell
    lon_suffix = "E" if longitude >= 0 else "W"
    lat_suffix = "N" if latitude >= 0 else "S"
    return (
        f"CORA_lon{abs(longitude):02d}{lon_suffix}_"
        f"lat{abs(latitude):02d}{lat_suffix}.nc"
    )


def _source_manifest_fingerprint(source_root: Path, files: Sequence[Path]) -> str:
    digest = hashlib.sha256()
    for path in files:
        stat = path.stat()
        record = f"{path.relative_to(source_root)}\t{stat.st_size}\t{stat.st_mtime_ns}\n"
        digest.update(record.encode("utf-8"))
    return digest.hexdigest()


class _TileWriters:
    def __init__(self, root: Path, pressure_grid: np.ndarray, max_open: int = 32):
        self.root = root
        self.pressure_grid = pressure_grid
        self.max_open = max_open
        self.handles: OrderedDict[Path, Any] = OrderedDict()

    def close(self) -> None:
        while self.handles:
            _, handle = self.handles.popitem(last=False)
            handle.close()

    def _open(self, path: Path) -> Any:
        try:
            from netCDF4 import Dataset  # type: ignore
        except ImportError as exc:  # pragma: no cover - required project dependency
            raise RuntimeError("netCDF4 is required to write CORA tiles") from exc
        if path in self.handles:
            handle = self.handles.pop(path)
            self.handles[path] = handle
            return handle
        if len(self.handles) >= self.max_open:
            _, old_handle = self.handles.popitem(last=False)
            old_handle.close()
        exists = path.exists()
        handle = Dataset(path, "a" if exists else "w", format="NETCDF4")
        if not exists:
            self._initialise(handle)
        self.handles[path] = handle
        return handle

    def _initialise(self, handle: Any) -> None:
        handle.createDimension("N_PROF", None)
        handle.createDimension("PRES_GRID", self.pressure_grid.size)
        pressure = handle.createVariable("PRES_GRID", "f4", ("PRES_GRID",))
        pressure[:] = self.pressure_grid.astype(np.float32)
        pressure.units = "dbar"
        pressure.positive = "down"

        for name in ("LATITUDE", "LONGITUDE", "JULD"):
            variable = handle.createVariable(name, "f8", ("N_PROF",), zlib=True, complevel=4)
            if name == "LATITUDE":
                variable.units = "degree_north"
            elif name == "LONGITUDE":
                variable.units = "degree_east"
            else:
                variable.units = _TIME_UNITS
                variable.calendar = "standard"
        for name in ("TEMP", "PSAL"):
            variable = handle.createVariable(
                name,
                "f4",
                ("N_PROF", "PRES_GRID"),
                fill_value=np.float32(np.nan),
                zlib=True,
                complevel=4,
                chunksizes=(128, self.pressure_grid.size),
            )
            variable.units = "degree_C" if name == "TEMP" else "1e-3"
        handle.createVariable("PROBE_TYPE", "i1", ("N_PROF",), zlib=True, complevel=4)
        handle.createVariable("SOURCE_PROFILE_INDEX", "i4", ("N_PROF",), zlib=True, complevel=4)
        for name in (
            "PLATFORM_NUMBER",
            "CYCLE_NUMBER",
            "DC_REFERENCE",
            "DATA_MODE",
            "WMO_INST_TYPE",
            "SOURCE_FILE",
            "SOURCE_FILE_CLASS",
            "REFERENCE_KIND",
            "TEMP_SOURCE",
            "PSAL_SOURCE",
            "PRES_SOURCE",
            "SOURCE_DATASET",
        ):
            handle.createVariable(name, str, ("N_PROF",))
        handle.title = "CTD and Argo reference profiles prepared from CORA"
        handle.product_id = PRODUCT_ID
        handle.source_datasets = f"{CORA_DATASET_ID}, {EASYCORA_DATASET_ID}"
        handle.created_at_utc = datetime.now(timezone.utc).isoformat()
        handle.profile_selection = "CTD: PROBE_TYPE=2; Argo: EasyCORA PF or raw CORA PF with PROBE_TYPE=5"
        handle.value_policy = "adjusted parameter for whole profile when present, otherwise raw"
        handle.qc_policy = "level QC in {1,2}; position and time QC in {1,2} when present"

    def append(self, cell: tuple[int, int], batch: ProfileBatch) -> None:
        if not len(batch):
            return
        path = self.root / _tile_filename(cell)
        handle = self._open(path)
        start = handle.dimensions["N_PROF"].size
        stop = start + len(batch)
        section = slice(start, stop)
        handle.variables["LATITUDE"][section] = np.asarray(batch.latitude, dtype=np.float64)
        handle.variables["LONGITUDE"][section] = np.asarray(batch.longitude, dtype=np.float64)
        handle.variables["JULD"][section] = np.asarray(batch.juld, dtype=np.float64)
        handle.variables["TEMP"][section, :] = np.stack(batch.temp)
        handle.variables["PSAL"][section, :] = np.stack(batch.psal)
        handle.variables["PROBE_TYPE"][section] = np.asarray(batch.probe_type, dtype=np.int8)
        handle.variables["SOURCE_PROFILE_INDEX"][section] = np.asarray(
            batch.source_profile_index, dtype=np.int32
        )
        for name, values in (
            ("PLATFORM_NUMBER", batch.platform_number),
            ("CYCLE_NUMBER", batch.cycle_number),
            ("DC_REFERENCE", batch.dc_reference),
            ("DATA_MODE", batch.data_mode),
            ("WMO_INST_TYPE", batch.wmo_inst_type),
            ("SOURCE_FILE", batch.source_file),
            ("SOURCE_FILE_CLASS", batch.source_file_class),
            ("REFERENCE_KIND", batch.reference_kind),
            ("TEMP_SOURCE", batch.temp_source),
            ("PSAL_SOURCE", batch.psal_source),
            ("PRES_SOURCE", batch.pres_source),
            ("SOURCE_DATASET", batch.source_dataset),
        ):
            handle.variables[name][section] = np.asarray(values, dtype=object)
        handle.sync()


def _process_source_file(
    path: Path,
    source_root: Path,
    pressure_grid: np.ndarray,
    *,
    require_both: bool,
    missing_qc: set[str],
) -> tuple[dict[tuple[int, int], ProfileBatch], Counter[str]]:
    try:
        import xarray as xr  # type: ignore
    except ImportError as exc:  # pragma: no cover - required project dependency
        raise RuntimeError("xarray is required to read CORA source files") from exc
    policy = _source_policy(path)
    counts: Counter[str] = Counter()
    grouped: dict[tuple[int, int], ProfileBatch] = {}
    with xr.open_dataset(path, decode_times=False, mask_and_scale=True) as dataset:
        n_prof = int(dataset.sizes.get("N_PROF", 0))
        counts["profiles_seen"] += n_prof
        if not n_prof:
            return grouped, counts
        for required in ("LATITUDE", "LONGITUDE"):
            if required not in dataset:
                raise SourceFormatError(f"{path.name} is missing {required}")
        time_name = "TIME" if "TIME" in dataset else "JULD" if "JULD" in dataset else None
        if time_name is None:
            raise SourceFormatError(f"{path.name} has neither TIME nor JULD")

        latitude = _numeric_profile_vector(dataset["LATITUDE"], n_prof)
        longitude = _numeric_profile_vector(dataset["LONGITUDE"], n_prof)
        juld = _time_days_since_1950(dataset[time_name], n_prof)
        position_qc, position_qc_name = _profile_qc(dataset, ("POSITION_QC",), n_prof)
        time_qc, time_qc_name = _profile_qc(dataset, (f"{time_name}_QC", "TIME_QC", "JULD_QC"), n_prof)
        if position_qc_name is None:
            missing_qc.add(f"{path.name}:POSITION_QC")
        if time_qc_name is None:
            missing_qc.add(f"{path.name}:{time_name}_QC")

        if "PROBE_TYPE" in dataset:
            probe_values = _profile_vector(dataset["PROBE_TYPE"], n_prof).reshape(n_prof, -1)
            if probe_values.shape[1] != 1:
                raise SourceFormatError(f"{path.name} PROBE_TYPE is not one value per profile")
            probe_type = _normalise_qc(probe_values[:, 0])
            kind_mask = probe_type == policy.required_probe_type
        elif policy.allow_missing_probe_type:
            probe_type = np.full(n_prof, policy.required_probe_type, dtype=np.int16)
            kind_mask = np.ones(n_prof, dtype=bool)
            counts["profiles_selected_without_probe_type"] += n_prof
        else:
            raise SourceFormatError(
                f"{path.name} has no PROBE_TYPE; file class alone is not sufficient for "
                f"safe {policy.reference_kind} selection"
            )

        platform = _profile_strings(dataset, ("PLATFORM_CODE", "PLATFORM_NUMBER"), n_prof)
        cycle_number = _profile_strings(dataset, ("CYCLE_NUMBER",), n_prof)
        dc_reference = _profile_strings(dataset, ("DC_REFERENCE",), n_prof)
        data_mode = _profile_strings(dataset, ("DATA_MODE",), n_prof)
        wmo_inst_type = _profile_strings(dataset, ("WMO_INST_TYPE",), n_prof)
        candidate_indices = np.flatnonzero(kind_mask)
        counts["profiles_matching_type"] += candidate_indices.size

        for index in candidate_indices:
            if (
                not np.isfinite(latitude[index])
                or not -90.0 <= latitude[index] <= 90.0
                or not np.isfinite(longitude[index])
                or not np.isfinite(juld[index])
                or position_qc[index] not in _GOOD_QC
                or time_qc[index] not in _GOOD_QC
            ):
                counts["profiles_bad_time_or_position"] += 1
                continue
            pressure, pressure_qc, pressure_source = _pressure_parameter(
                dataset, int(index), float(latitude[index])
            )
            temp, temp_qc, temp_source = _best_parameter(dataset, "TEMP", int(index))
            psal, psal_qc, psal_source = _best_parameter(dataset, "PSAL", int(index))
            temp_grid = _interpolate_profile(
                pressure, pressure_qc, temp, temp_qc, pressure_grid
            )
            psal_grid = _interpolate_profile(
                pressure, pressure_qc, psal, psal_qc, pressure_grid
            )
            has_temp = bool(np.any(np.isfinite(temp_grid)))
            has_psal = bool(np.any(np.isfinite(psal_grid)))
            if (require_both and not (has_temp and has_psal)) or (
                not require_both and not (has_temp or has_psal)
            ):
                counts["profiles_insufficient_ts"] += 1
                continue

            normalised_lon = ((float(longitude[index]) + 180.0) % 360.0) - 180.0
            cell = _tile_cell(normalised_lon, float(latitude[index]))
            batch = grouped.setdefault(cell, ProfileBatch.empty())
            batch.latitude.append(float(latitude[index]))
            batch.longitude.append(normalised_lon)
            batch.juld.append(float(juld[index]))
            batch.temp.append(temp_grid)
            batch.psal.append(psal_grid)
            batch.probe_type.append(int(probe_type[index]))
            batch.source_profile_index.append(int(index))
            batch.platform_number.append(platform[index])
            batch.cycle_number.append(cycle_number[index])
            batch.dc_reference.append(dc_reference[index])
            batch.data_mode.append(data_mode[index])
            batch.wmo_inst_type.append(wmo_inst_type[index])
            batch.source_file.append(str(path.relative_to(source_root)))
            batch.source_file_class.append(policy.file_class)
            batch.reference_kind.append(policy.reference_kind)
            batch.temp_source.append(temp_source)
            batch.psal_source.append(psal_source)
            batch.pres_source.append(pressure_source)
            batch.source_dataset.append(policy.dataset_id)
            counts[f"profiles_written_{policy.reference_kind.lower()}"] += 1
    return grouped, counts


def prepare_cora_tiles(
    source_root: Path | str,
    output_dir: Path | str,
    *,
    pressure_start: float = 10.0,
    pressure_stop: float = 1000.0,
    pressure_step: float = 10.0,
    require_both: bool = True,
    overwrite: bool = False,
    skip_bad_files: bool = False,
) -> dict[str, Any]:
    """Build provenance-preserving 10-degree reference tiles from selected sources."""
    source = Path(source_root).expanduser().resolve()
    output = Path(output_dir).expanduser().resolve()
    if not source.is_dir():
        raise ValueError(f"Source directory does not exist: {source}")
    if pressure_step <= 0 or pressure_stop < pressure_start:
        raise ValueError("Invalid pressure grid")
    pressure_grid = np.arange(
        pressure_start, pressure_stop + pressure_step * 0.5, pressure_step, dtype=np.float64
    )
    files = sorted(path for path in (source / "source").rglob("*.nc"))
    if not files:
        files = sorted(source.rglob("*.nc"))
    files = [path for path in files if _PROFILE_NAME.search(path.name)]
    if not files:
        raise ValueError(f"No recognised CORA daily profile files found below {source}")

    output.mkdir(parents=True, exist_ok=True)
    existing = sorted(output.glob("CORA_lon*_lat*.nc"))
    if existing and not overwrite:
        raise ValueError(
            f"Output directory already contains {len(existing)} tiles: {output}. "
            "Use a new directory or pass --overwrite."
        )
    if overwrite:
        for path in existing:
            path.unlink()
        manifest_path = output / "PREPARATION_MANIFEST.json"
        if manifest_path.exists():
            manifest_path.unlink()

    counters: Counter[str] = Counter()
    missing_qc: set[str] = set()
    errors: list[dict[str, str]] = []
    writer = _TileWriters(output, pressure_grid)
    try:
        for number, path in enumerate(files, start=1):
            try:
                grouped, file_counts = _process_source_file(
                    path,
                    source,
                    pressure_grid,
                    require_both=require_both,
                    missing_qc=missing_qc,
                )
            except (OSError, SourceFormatError, ValueError) as exc:
                if not skip_bad_files:
                    raise SourceFormatError(f"Failed while processing {path}: {exc}") from exc
                errors.append({"file": str(path.relative_to(source)), "error": str(exc)})
                counters["source_files_skipped"] += 1
                continue
            counters.update(file_counts)
            counters["source_files_processed"] += 1
            for cell, batch in grouped.items():
                writer.append(cell, batch)
            if number % 250 == 0:
                print(f"processed {number}/{len(files)} source files")
    finally:
        writer.close()

    tiles = sorted(output.glob("CORA_lon*_lat*.nc"))
    written = counters["profiles_written_argo"] + counters["profiles_written_ctd"]
    if not written:
        raise RuntimeError("No CTD or Argo profiles passed the selection and QC policy")
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "product_id": PRODUCT_ID,
        "source_root": str(source),
        "source_download_manifest": (
            str(source / "download_manifest.json")
            if (source / "download_manifest.json").exists()
            else None
        ),
        "source_file_count": len(files),
        "source_inventory_fingerprint_sha256": _source_manifest_fingerprint(source, files),
        "output_dir": str(output),
        "tile_count": len(tiles),
        "pressure_grid_dbar": {
            "start": pressure_start,
            "stop": pressure_stop,
            "step": pressure_step,
            "count": pressure_grid.size,
        },
        "selection": {
            "ctd": "raw CORA candidate classes CT/OC/TE and PROBE_TYPE=2",
            "argo": "EasyCORA PF (default download) or raw CORA PF and PROBE_TYPE=5",
            "require_both_temperature_and_salinity": require_both,
        },
        "value_policy": "adjusted whole-profile parameter when present, otherwise raw",
        "accepted_level_qc": list(_GOOD_QC),
        "counters": dict(sorted(counters.items())),
        "missing_qc_variables": sorted(missing_qc),
        "skipped_file_errors": errors,
    }
    (output / "PREPARATION_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare downloaded CTD/Argo CORA files as MEOP 10-degree reference tiles."
    )
    parser.add_argument("--source-root", type=Path, required=True, help="Root created by meop-cora-download.")
    parser.add_argument("--output-dir", type=Path, required=True, help="New directory for prepared tiles.")
    parser.add_argument("--pressure-start", type=float, default=10.0)
    parser.add_argument("--pressure-stop", type=float, default=1000.0)
    parser.add_argument("--pressure-step", type=float, default=10.0)
    parser.add_argument(
        "--allow-single-variable",
        action="store_true",
        help="Keep profiles with usable temperature or salinity instead of requiring both.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing CORA tile files in the exact output directory.",
    )
    parser.add_argument(
        "--skip-bad-files",
        action="store_true",
        help="Record and skip malformed source files instead of stopping the build.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        manifest = prepare_cora_tiles(
            args.source_root,
            args.output_dir,
            pressure_start=args.pressure_start,
            pressure_stop=args.pressure_stop,
            pressure_step=args.pressure_step,
            require_both=not args.allow_single_variable,
            overwrite=args.overwrite,
            skip_bad_files=args.skip_bad_files,
        )
    except (RuntimeError, ValueError) as exc:
        print(f"error: {exc}")
        return 2
    counters = manifest["counters"]
    print(
        f"Prepared {manifest['tile_count']} tiles: "
        f"{counters.get('profiles_written_argo', 0)} Argo and "
        f"{counters.get('profiles_written_ctd', 0)} CTD profiles"
    )
    print(f"Manifest: {Path(manifest['output_dir']) / 'PREPARATION_MANIFEST.json'}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
