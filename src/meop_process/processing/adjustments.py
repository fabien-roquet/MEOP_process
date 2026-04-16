from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import numpy as np
import xarray as xr

from ..catalog.deployments import load_info_deployment
from ..catalog.filenames import fname_prof, list_fname_prof
from ..catalog.tables import read_csv_rows
from ..data.layout import resolve_table_path
from ..models import MeopConfig, Selection
from .hr import _open_dataset, create_fr1_python, create_hr1_python
from .netcdf import save_dataset_with_compression
from .qc import _load_coeff_row, _to_numeric_qc, ensure_processing_parameters


PROCESS_QFS = ("lr0", "lr1", "hr0", "hr1", "fr0", "fr1")


@dataclass(frozen=True)
class AdjustmentResult:
    written_files: tuple[Path, ...]
    processed_tags: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "written_files": [str(path) for path in self.written_files],
            "processed_tags": list(self.processed_tags),
        }


def _selected_tags(config: MeopConfig, selection: Selection) -> tuple[str, ...]:
    selection = selection.normalized()
    if selection.smru_name:
        return (selection.smru_name,)

    discovered: set[str] = set()
    for qf in PROCESS_QFS:
        discovered.update(path.name.split("_")[0] for path in list_fname_prof(deployment=selection.deployment, qf=qf, config=config))
    if discovered:
        return tuple(sorted(discovered))

    info = load_info_deployment(config, deployment=selection.deployment, smru_name=selection.smru_name)
    return tuple(sorted(info.list_smru_name))


def _ensure_default_coefficients(config: MeopConfig, smru_names: Iterable[str]) -> Path | None:
    path = resolve_table_path(config, "table_coeff.csv", required=False)
    rows = read_csv_rows(path)
    if not rows:
        return None

    existing = {row.get("smru_platform_code", "") or row.get("row_name", "") for row in rows}
    changed = False
    for smru_name in smru_names:
        if smru_name in existing:
            continue
        rows.append(
            {
                "smru_platform_code": smru_name,
                "T1": "0",
                "T2": "0",
                "S1": "0",
                "S2": "0",
                "remove": "0",
                "Sremove": "0",
                "comment": "no comment",
            }
        )
        changed = True
    if not changed:
        return None

    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key and key not in fieldnames:
                fieldnames.append(key)

    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return path


def _atomic_write_dataset(dataset: xr.Dataset, path: Path) -> None:
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    path.parent.mkdir(parents=True, exist_ok=True)
    save_dataset_with_compression(dataset, tmp_path)
    tmp_path.replace(path)


def _as_float(value: object, default: float = 0.0) -> float:
    if value in (None, ""):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _decode_text(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="ignore").strip()
    return str(value).strip()


def _parameter_names(dataset: xr.Dataset) -> list[str]:
    if "PARAMETER" in dataset:
        values = np.asarray(dataset["PARAMETER"].values, dtype=object)
        if values.ndim >= 3 and values.shape[0] > 0 and values.shape[1] > 0:
            return [_decode_text(item) for item in values[0, 0, :]]
    if "STATION_PARAMETERS" in dataset:
        values = np.asarray(dataset["STATION_PARAMETERS"].values, dtype=object)
        if values.ndim >= 2 and values.shape[0] > 0:
            return [_decode_text(item) for item in values[0, :]]
    names = ["PRES", "TEMP", "PSAL"]
    for optional in ("CHLA", "DOXY", "LIGHT"):
        if optional in dataset:
            names.append(optional)
    return names


def _copy_variable_metadata(target: xr.DataArray, source: xr.DataArray) -> None:
    target.attrs.update(source.attrs)
    encoding = dict(source.encoding)
    dtype_kind = np.asarray(target.values).dtype.kind
    if dtype_kind in {"U", "S", "O"}:
        for key in ("dtype", "char_dim_name"):
            encoding.pop(key, None)
    target.encoding.update(encoding)


def _format_coeff_string(kind: str, value1: float, value2: float) -> str:
    if kind == "PRES":
        return f"p1= {value1:8.6e} dbar/km, p2= {value2:8.6e} dbar"
    if kind == "TEMP":
        return f"t1= {value1:8.6e} degC/km, t2= {value2:8.6e} degC"
    if kind == "PSAL":
        return f"s1= {value1:8.6e}  psu/km, s2= {value2:8.6e}  psu"
    if kind == "CHLA":
        return "f1= 0.6, f2= 0.0"
    return " "


def _load_salinity_offset_series(config: MeopConfig, smru_name: str, n_prof: int) -> np.ndarray | None:
    path = resolve_table_path(config, "table_salinity_offsets.csv", required=False)
    rows = read_csv_rows(path)
    row = next((record for record in rows if record.get("smru_platform_code", "") == smru_name), None)
    if row is None or n_prof <= 0:
        return None

    indices: list[int] = []
    offsets: list[float] = []
    for idx in range(1, 5):
        index_value = int(_as_float(row.get(f"index_{idx}"), default=0.0))
        offset_value = _as_float(row.get(f"offset_{idx}"), default=0.0)
        indices.append(index_value)
        offsets.append(offset_value)

    zero_positions = [position for position, value in enumerate(indices) if value == 0]
    if zero_positions:
        indices[zero_positions[0]] = n_prof
        for position in reversed(zero_positions[1:]):
            del indices[position]
            del offsets[position]

    valid_pairs = [(index, offset) for index, offset in zip(indices, offsets, strict=False) if index > 0]
    if not valid_pairs:
        return None

    valid_pairs.sort(key=lambda item: item[0])
    xp = np.asarray([index for index, _ in valid_pairs], dtype=np.float64)
    fp = np.asarray([offset for _, offset in valid_pairs], dtype=np.float64)
    x = np.arange(1, n_prof + 1, dtype=np.float64)
    if xp.size == 1:
        return np.full(n_prof, fp[0], dtype=np.float64)
    return np.interp(x, xp, fp).astype(np.float64)


def _profile_good_counts_from_lr0(config: MeopConfig, smru_name: str, deployment: str) -> tuple[np.ndarray | None, np.ndarray | None]:
    path = fname_prof(smru_name, deployment=deployment, qf="lr0", config=config)
    if not path.exists():
        return None, None
    dataset = _open_dataset(path)
    try:
        temp_qc = _to_numeric_qc(dataset["TEMP_QC"].values) if "TEMP_QC" in dataset else None
        psal_qc = _to_numeric_qc(dataset["PSAL_QC"].values) if "PSAL_QC" in dataset else None
        temp_counts = np.sum(temp_qc <= 1, axis=1) if temp_qc is not None else None
        psal_counts = np.sum(psal_qc <= 1, axis=1) if psal_qc is not None else None
        return temp_counts, psal_counts
    finally:
        dataset.close()


def _apply_offset_and_calibration(dataset: xr.Dataset, *, config: MeopConfig, smru_name: str) -> xr.Dataset:
    result = dataset.copy(deep=True)
    coeff_row = _load_coeff_row(config, smru_name)

    p1 = 0.0
    p2 = 0.0
    t1 = _as_float(coeff_row.get("T1"))
    t2 = _as_float(coeff_row.get("T2"))
    s1 = _as_float(coeff_row.get("S1"))
    s2 = _as_float(coeff_row.get("S2"))

    pres = np.asarray(result["PRES"].values, dtype=np.float64)
    temp = np.asarray(result["TEMP"].values, dtype=np.float64)
    psal = np.asarray(result["PSAL"].values, dtype=np.float64)

    pres_adjusted = pres - p1 * 1e-3 * pres - p2
    temp_adjusted = temp - t1 * 1e-3 * pres_adjusted - t2
    psal_adjusted = psal - s1 * 1e-3 * pres_adjusted - s2

    result["PRES_ADJUSTED"] = (result["PRES_ADJUSTED"].dims, pres_adjusted.astype(np.float32))
    result["TEMP_ADJUSTED"] = (result["TEMP_ADJUSTED"].dims, temp_adjusted.astype(np.float32))
    result["PSAL_ADJUSTED"] = (result["PSAL_ADJUSTED"].dims, psal_adjusted.astype(np.float32))
    _copy_variable_metadata(result["PRES_ADJUSTED"], dataset["PRES_ADJUSTED"])
    _copy_variable_metadata(result["TEMP_ADJUSTED"], dataset["TEMP_ADJUSTED"])
    _copy_variable_metadata(result["PSAL_ADJUSTED"], dataset["PSAL_ADJUSTED"])

    if "CHLA" in result and "CHLA_ADJUSTED" in result:
        chla = np.asarray(result["CHLA"].values, dtype=np.float64)
        chla_adjusted = 0.6 * chla
        result["CHLA_ADJUSTED"] = (result["CHLA_ADJUSTED"].dims, chla_adjusted.astype(np.float32))
        _copy_variable_metadata(result["CHLA_ADJUSTED"], dataset["CHLA_ADJUSTED"])

    if "DOXY" in result and "DOXY_ADJUSTED" in result:
        result["DOXY_ADJUSTED"] = (result["DOXY_ADJUSTED"].dims, np.asarray(result["DOXY"].values).copy())
        _copy_variable_metadata(result["DOXY_ADJUSTED"], dataset["DOXY_ADJUSTED"])

    if "LIGHT" in result and "LIGHT_ADJUSTED" in result:
        result["LIGHT_ADJUSTED"] = (result["LIGHT_ADJUSTED"].dims, np.asarray(result["LIGHT"].values).copy())
        _copy_variable_metadata(result["LIGHT_ADJUSTED"], dataset["LIGHT_ADJUSTED"])

    salinity_offsets = _load_salinity_offset_series(config, smru_name, int(result.sizes.get("N_PROF", 0)))
    if salinity_offsets is not None and salinity_offsets.size == psal_adjusted.shape[0]:
        psal_adjusted = psal_adjusted + salinity_offsets[:, None]
        result["PSAL_ADJUSTED"] = (result["PSAL_ADJUSTED"].dims, psal_adjusted.astype(np.float32))
        _copy_variable_metadata(result["PSAL_ADJUSTED"], dataset["PSAL_ADJUSTED"])

    if "SCIENTIFIC_CALIB_COEFFICIENT" in result:
        coeff_values = np.asarray(result["SCIENTIFIC_CALIB_COEFFICIENT"].values, dtype=object).copy()
        parameter_names = _parameter_names(result)
        name_to_index = {name: index for index, name in enumerate(parameter_names)}
        n_prof = coeff_values.shape[0] if coeff_values.ndim >= 1 else int(result.sizes.get("N_PROF", 0))

        updates: dict[str, np.ndarray] = {
            "PRES": np.full(n_prof, _format_coeff_string("PRES", p1, p2), dtype=object),
            "TEMP": np.full(n_prof, _format_coeff_string("TEMP", t1, t2), dtype=object),
            "PSAL": np.full(n_prof, _format_coeff_string("PSAL", s1, s2), dtype=object),
        }
        if salinity_offsets is not None and salinity_offsets.size == n_prof:
            updates["PSAL"] = np.asarray(
                [_format_coeff_string("PSAL", s1, s2 - float(offset)) for offset in salinity_offsets],
                dtype=object,
            )
        if "CHLA" in name_to_index:
            updates["CHLA"] = np.full(n_prof, _format_coeff_string("CHLA", 0.6, 0.0), dtype=object)

        for parameter_name, strings in updates.items():
            index = name_to_index.get(parameter_name)
            if index is None:
                continue
            if coeff_values.ndim == 3:
                coeff_values[:, 0, index] = strings
            elif coeff_values.ndim == 2:
                coeff_values[:, index] = strings

        result["SCIENTIFIC_CALIB_COEFFICIENT"] = (
            result["SCIENTIFIC_CALIB_COEFFICIENT"].dims,
            coeff_values,
        )
        _copy_variable_metadata(result["SCIENTIFIC_CALIB_COEFFICIENT"], dataset["SCIENTIFIC_CALIB_COEFFICIENT"])

    return result


def _assign_error_estimates(dataset: xr.Dataset, *, config: MeopConfig, deployment: str, smru_name: str) -> xr.Dataset:
    result = dataset.copy(deep=True)
    params = ensure_processing_parameters(config, deployment)
    temp_error = _as_float(params.get("temp_error"), default=np.nan)
    psal_error = _as_float(params.get("psal_error"), default=np.nan)

    temp_qc = _to_numeric_qc(result["TEMP_QC"].values) if "TEMP_QC" in result else None
    psal_qc = _to_numeric_qc(result["PSAL_QC"].values) if "PSAL_QC" in result else None

    if temp_qc is not None and "TEMP_ADJUSTED_ERROR" in result:
        temp_errors = np.full(result["TEMP_ADJUSTED_ERROR"].shape, np.nan, dtype=np.float32)
        temp_errors[temp_qc < 9] = np.float32(temp_error)
        lr0_temp_counts, _ = _profile_good_counts_from_lr0(config, smru_name, deployment)
        if lr0_temp_counts is not None and lr0_temp_counts.shape[0] == temp_errors.shape[0]:
            low_counts = np.where(lr0_temp_counts < 10)[0]
            if low_counts.size:
                temp_errors[low_counts, :] = np.float32(temp_error * 2.0)
        temp_errors[temp_qc == 9] = np.nan
        result["TEMP_ADJUSTED_ERROR"] = (result["TEMP_ADJUSTED_ERROR"].dims, temp_errors)
        _copy_variable_metadata(result["TEMP_ADJUSTED_ERROR"], dataset["TEMP_ADJUSTED_ERROR"])

    if psal_qc is not None and "PSAL_ADJUSTED_ERROR" in result:
        psal_errors = np.full(result["PSAL_ADJUSTED_ERROR"].shape, np.nan, dtype=np.float32)
        psal_errors[psal_qc < 9] = np.float32(psal_error)
        _, lr0_psal_counts = _profile_good_counts_from_lr0(config, smru_name, deployment)
        if lr0_psal_counts is not None and lr0_psal_counts.shape[0] == psal_errors.shape[0]:
            low_counts = np.where(lr0_psal_counts < 10)[0]
            if low_counts.size:
                psal_errors[low_counts, :] = np.float32(psal_error * 2.0)
        psal_errors[psal_qc == 9] = np.nan
        result["PSAL_ADJUSTED_ERROR"] = (result["PSAL_ADJUSTED_ERROR"].dims, psal_errors)
        _copy_variable_metadata(result["PSAL_ADJUSTED_ERROR"], dataset["PSAL_ADJUSTED_ERROR"])

    return result


def _apply_adjustments_to_file(config: MeopConfig, path: Path, *, smru_name: str, deployment: str) -> bool:
    if not path.exists():
        return False
    dataset = _open_dataset(path)
    try:
        updated = _apply_offset_and_calibration(dataset, config=config, smru_name=smru_name)
        updated = _assign_error_estimates(updated, config=config, deployment=deployment, smru_name=smru_name)
        _atomic_write_dataset(updated, path)
        return True
    finally:
        dataset.close()


def apply_adjustments(config: MeopConfig, selection: Selection) -> AdjustmentResult:
    selection = selection.normalized()
    tags = _selected_tags(config, selection)
    _ensure_default_coefficients(config, tags)

    written: list[Path] = []
    processed_tags: list[str] = []
    for smru_name in tags:
        processed_any = False
        for qf in PROCESS_QFS:
            path = fname_prof(smru_name, deployment=selection.deployment, qf=qf, config=config)
            if _apply_adjustments_to_file(config, path, smru_name=smru_name, deployment=selection.deployment):
                written.append(path)
                processed_any = True
        if processed_any:
            processed_tags.append(smru_name)

    return AdjustmentResult(written_files=tuple(written), processed_tags=tuple(processed_tags))


def apply_tlc(config: MeopConfig, selection: Selection):
    return create_hr1_python(config, selection, thermal_lag=True)


def apply_tlc_fr(config: MeopConfig, selection: Selection):
    return create_fr1_python(config, selection, thermal_lag=True)


def apply_notlc(config: MeopConfig, selection: Selection):
    return create_hr1_python(config, selection, thermal_lag=False)


def apply_notlc_fr(config: MeopConfig, selection: Selection):
    return create_fr1_python(config, selection, thermal_lag=False)
