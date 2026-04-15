from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import xarray as xr

from ..catalog.tables import read_csv_rows, read_indexed_csv_rows, write_indexed_csv_rows
from ..data.layout import resolve_table_path
from ..models import DeploymentInfo, MeopConfig


DEFAULT_PARAM_ROW: dict[str, str] = {
    "temp_error": "0.1",
    "psal_error": "0.2",
    "minT": "-3",
    "maxT": "32",
    "minS": "4",
    "maxS": "40",
    "min_Nprof": "30",
    "pmax": "1000",
    "pmax_fluo": "200",
    "is_lon_centre_180": "0",
}

TABLE_PARAM_FIELD_ORDER = [
    "temp_error",
    "psal_error",
    "minT",
    "maxT",
    "minS",
    "maxS",
    "min_Nprof",
    "pmax",
    "pmax_fluo",
    "is_lon_centre_180",
]

TABLE_COEFF_FIELD_ORDER = ["smru_platform_code", "T1", "T2", "S1", "S2", "remove", "Sremove", "comment"]


@dataclass(frozen=True)
class QcFilterResult:
    full_profile_removals: int
    salinity_profile_removals: int
    dataset: xr.Dataset


class _WorkingState:
    def __init__(self, dataset: xr.Dataset) -> None:
        self.dataset = dataset
        self.pres = np.asarray(dataset["PRES"].values, dtype=np.float64)
        self.temp = np.asarray(dataset["TEMP"].values, dtype=np.float64)
        self.psal = np.asarray(dataset["PSAL"].values, dtype=np.float64)
        self.lat = np.asarray(dataset["LATITUDE"].values, dtype=np.float64)
        self.lon = np.asarray(dataset["LONGITUDE"].values, dtype=np.float64)
        self.juld = np.asarray(dataset["JULD"].values, dtype=np.float64)

        self.pres_qc = _to_numeric_qc(dataset["PRES_QC"].values)
        self.temp_qc = _to_numeric_qc(dataset["TEMP_QC"].values)
        self.psal_qc = _to_numeric_qc(dataset["PSAL_QC"].values)

        self.chla_qc = _optional_numeric_qc(dataset, "CHLA_QC")
        self.doxy_qc = _optional_numeric_qc(dataset, "DOXY_QC")
        self.light_qc = _optional_numeric_qc(dataset, "LIGHT_QC")

        self._sync_masks()

    def _sync_masks(self) -> None:
        self.pres_work = np.where(self.pres_qc <= 1, self.pres, np.nan)
        self.temp_work = np.where(self.temp_qc <= 1, self.temp, np.nan)
        self.psal_work = np.where(self.psal_qc <= 1, self.psal, np.nan)

    def remove_full_profiles(self, indices: np.ndarray) -> int:
        profile_index = _clean_profile_indices(indices, self.pres_qc.shape[0])
        if profile_index.size == 0:
            return 0
        _set_profile_qc(self.pres_qc, profile_index)
        _set_profile_qc(self.temp_qc, profile_index)
        _set_profile_qc(self.psal_qc, profile_index)
        if self.chla_qc is not None:
            _set_profile_qc(self.chla_qc, profile_index)
        if self.doxy_qc is not None:
            _set_profile_qc(self.doxy_qc, profile_index)
        if self.light_qc is not None:
            _set_profile_qc(self.light_qc, profile_index)
        self._sync_masks()
        return int(profile_index.size)

    def remove_salinity_profiles(self, indices: np.ndarray) -> int:
        profile_index = _clean_profile_indices(indices, self.psal_qc.shape[0])
        if profile_index.size == 0:
            return 0
        _set_profile_qc(self.psal_qc, profile_index)
        pres_extra = (self.psal_qc[profile_index] == 4) & (self.temp_qc[profile_index] == 9) & (self.pres_qc[profile_index] != 9)
        self.pres_qc[profile_index] = np.where(pres_extra, 4, self.pres_qc[profile_index])
        self._sync_masks()
        return int(profile_index.size)

    def remove_tag(self) -> int:
        valid_before = np.any(self.pres_qc <= 1)
        count = self.remove_full_profiles(np.arange(self.pres_qc.shape[0]))
        return count if valid_before else 0

    def to_dataset(self) -> xr.Dataset:
        pres_qc = _to_string_qc(self.pres_qc)
        temp_qc = _to_string_qc(self.temp_qc)
        psal_qc = _to_string_qc(self.psal_qc)
        self.dataset["PRES_QC"] = (("N_PROF", "N_LEVELS"), pres_qc)
        self.dataset["TEMP_QC"] = (("N_PROF", "N_LEVELS"), temp_qc)
        self.dataset["PSAL_QC"] = (("N_PROF", "N_LEVELS"), psal_qc)
        self.dataset["PRES_ADJUSTED_QC"] = (("N_PROF", "N_LEVELS"), pres_qc.copy())
        self.dataset["TEMP_ADJUSTED_QC"] = (("N_PROF", "N_LEVELS"), temp_qc.copy())
        self.dataset["PSAL_ADJUSTED_QC"] = (("N_PROF", "N_LEVELS"), psal_qc.copy())

        if self.chla_qc is not None and "CHLA_QC" in self.dataset:
            chla_qc = _to_string_qc(self.chla_qc)
            self.dataset["CHLA_QC"] = (("N_PROF", "N_LEVELS"), chla_qc)
            if "CHLA_ADJUSTED_QC" in self.dataset:
                self.dataset["CHLA_ADJUSTED_QC"] = (("N_PROF", "N_LEVELS"), chla_qc.copy())
        if self.doxy_qc is not None and "DOXY_QC" in self.dataset:
            doxy_qc = _to_string_qc(self.doxy_qc)
            self.dataset["DOXY_QC"] = (("N_PROF", "N_LEVELS"), doxy_qc)
            if "DOXY_ADJUSTED_QC" in self.dataset:
                self.dataset["DOXY_ADJUSTED_QC"] = (("N_PROF", "N_LEVELS"), doxy_qc.copy())
        if self.light_qc is not None and "LIGHT_QC" in self.dataset:
            light_qc = _to_string_qc(self.light_qc)
            self.dataset["LIGHT_QC"] = (("N_PROF", "N_LEVELS"), light_qc)
            if "LIGHT_ADJUSTED_QC" in self.dataset:
                self.dataset["LIGHT_ADJUSTED_QC"] = (("N_PROF", "N_LEVELS"), light_qc.copy())

        _update_count_and_geo_attributes(self.dataset, self)
        return self.dataset


try:  # pragma: no cover - exercised only when gsw is installed
    import gsw as _gsw  # type: ignore
except Exception:  # pragma: no cover - used in the current execution environment
    _gsw = None


def apply_lr0_qc_filters(config: MeopConfig, info: DeploymentInfo, smru_name: str, dataset: xr.Dataset) -> QcFilterResult:
    """Apply the historical seal-QC logic to one ``*_lr0_prof.nc`` dataset.

    The Python port intentionally matches the historical behavior closely:
    - raw data values remain untouched;
    - QC arrays are updated in place;
    - count/geospatial attributes are recomputed from QC arrays.

    Density-dependent criteria prefer the Python ``gsw`` package when available. When the
    runtime environment does not provide it, a lightweight EOS-80 fallback is used so the
    workflow remains importable and testable.
    """

    state = _WorkingState(dataset)
    params = ensure_processing_parameters(config, info.EXP)
    coeff_row = _load_coeff_row(config, smru_name)
    filter_rows = _load_filter_rows(config, info.EXP, smru_name)

    full_removed = 0
    sal_removed = 0

    # lat/lon/date
    bad_llj = np.where(~np.isfinite(state.lat * state.lon * state.juld))[0]
    full_removed += state.remove_full_profiles(bad_llj)

    # full-profile outliers
    n_temp = np.sum(state.temp_qc <= 1, axis=1)
    full_removed += state.remove_full_profiles(np.where((n_temp < 3) & (n_temp > 0))[0])
    full_removed += state.remove_full_profiles(_criterion_indices(state, "Tmin", [_as_float(params.get("minT"))], sal_only=False))
    full_removed += state.remove_full_profiles(_criterion_indices(state, "Tmax", [_as_float(params.get("maxT"))], sal_only=False))
    temp_diff = np.diff(state.temp_work, axis=1)
    constant_temp = np.where(np.all(temp_diff == 0, axis=1))[0]
    full_removed += state.remove_full_profiles(constant_temp)
    full_removed += state.remove_full_profiles(np.where(state.lat == 0)[0])
    lat_mean = np.nanmean(state.lat)
    lat_std = np.nanstd(state.lat)
    if np.isfinite(lat_std) and lat_std > 0:
        full_removed += state.remove_full_profiles(np.where(np.abs(state.lat - lat_mean) > 5 * lat_std)[0])

    # salinity-profile outliers
    n_psal = np.sum(state.psal_qc < 2, axis=1)
    sal_removed += state.remove_salinity_profiles(np.where((n_psal <= 5) & (n_psal > 0))[0])
    sal_removed += state.remove_salinity_profiles(_criterion_indices(state, "Smin", [_as_float(params.get("minS"))], sal_only=True))
    sal_removed += state.remove_salinity_profiles(_criterion_indices(state, "Smax", [_as_float(params.get("maxS"))], sal_only=True))

    # density-related salinity filters
    n_psal = np.sum(state.psal_qc < 2, axis=1)
    max_pressure = _safe_nanmax(state.pres_work, axis=1)
    pressure_step = np.diff(state.pres_work, axis=1)
    vertical_resolution = np.where(np.any(pressure_step > (max_pressure[:, None] / 3.0), axis=1))[0]
    sal_removed += state.remove_salinity_profiles(vertical_resolution)

    inversion_profiles = _top_to_bottom_inversion_indices(state, n_psal)
    sal_removed += state.remove_salinity_profiles(inversion_profiles)

    # manual editing from table_coeff
    if _as_bool(coeff_row.get("remove")):
        full_removed += state.remove_tag()
    if _as_bool(coeff_row.get("Sremove")):
        sal_removed += state.remove_salinity_profiles(np.arange(state.psal_qc.shape[0]))

    # deployment/tag-specific filters from table_filter
    for row in filter_rows:
        criterion = row.get("filter", "")
        value = _filter_values(row)
        if not criterion:
            continue
        if _as_bool(row.get("Sonly")):
            sal_removed += state.remove_salinity_profiles(_criterion_indices(state, criterion, value, sal_only=True))
        else:
            full_removed += state.remove_full_profiles(_criterion_indices(state, criterion, value, sal_only=False))

    # minimum valid profile count
    n_temp_valid = np.sum(state.temp_qc <= 1, axis=1)
    valid_temp_profiles = int(np.sum(n_temp_valid > 5))
    minimum_profiles = int(_as_float(params.get("min_Nprof"), 30))
    if 0 < valid_temp_profiles < minimum_profiles:
        full_removed += state.remove_tag()
        _set_coeff_flag(config, smru_name, "remove", "1")

    return QcFilterResult(
        full_profile_removals=full_removed,
        salinity_profile_removals=sal_removed,
        dataset=state.to_dataset(),
    )


def ensure_processing_parameters(config: MeopConfig, deployment: str) -> dict[str, str]:
    rows = _load_indexed_runtime_rows(config, "table_param.csv")
    for row in rows:
        if row.get("row_name") == deployment:
            return row

    new_row = {"row_name": deployment, **DEFAULT_PARAM_ROW}
    rows.append(new_row)
    _write_runtime_indexed_rows(config, "table_param.csv", rows, field_order=TABLE_PARAM_FIELD_ORDER)
    return new_row


def _load_filter_rows(config: MeopConfig, deployment: str, smru_name: str) -> list[dict[str, str]]:
    path = resolve_table_path(config, "table_filter.csv", required=False)
    rows = read_csv_rows(path)
    deployment_rows = [row for row in rows if row.get("smru_platform_name") == deployment]
    tag_rows = [row for row in rows if row.get("smru_platform_name") == smru_name]
    return deployment_rows + tag_rows


def _load_coeff_row(config: MeopConfig, smru_name: str) -> dict[str, str]:
    rows = _load_indexed_runtime_rows(config, "table_coeff.csv")
    for row in rows:
        if row.get("row_name") == smru_name or row.get("smru_platform_code") == smru_name:
            return row
    return {"row_name": smru_name, "smru_platform_code": smru_name, "remove": "0", "Sremove": "0"}


def _set_coeff_flag(config: MeopConfig, smru_name: str, field: str, value: str) -> None:
    rows = _load_indexed_runtime_rows(config, "table_coeff.csv")
    target: dict[str, str] | None = None
    for row in rows:
        if row.get("row_name") == smru_name or row.get("smru_platform_code") == smru_name:
            target = row
            break
    if target is None:
        target = {
            "row_name": smru_name,
            "smru_platform_code": smru_name,
            "T1": "0",
            "T2": "0",
            "S1": "0",
            "S2": "0",
            "remove": "0",
            "Sremove": "0",
            "comment": "no comment",
        }
        rows.append(target)
    target[field] = value
    _write_runtime_indexed_rows(config, "table_coeff.csv", rows, field_order=TABLE_COEFF_FIELD_ORDER)


def _load_indexed_runtime_rows(config: MeopConfig, name: str) -> list[dict[str, str]]:
    path = resolve_table_path(config, name, required=False)
    if not path.exists():
        return []
    return read_indexed_csv_rows(path)


def _write_runtime_indexed_rows(config: MeopConfig, name: str, rows: list[dict[str, str]], *, field_order: Iterable[str]) -> Path:
    path = resolve_table_path(config, name, required=False)
    destination = write_indexed_csv_rows(path, rows, field_order=field_order)
    legacy = config.processdir / name
    if legacy != destination:
        legacy.write_text(destination.read_text(encoding="utf-8"), encoding="utf-8")
    return destination


def _optional_numeric_qc(dataset: xr.Dataset, name: str) -> np.ndarray | None:
    if name not in dataset:
        return None
    return _to_numeric_qc(dataset[name].values)


def _to_numeric_qc(values: np.ndarray) -> np.ndarray:
    text = np.asarray(values).astype("U4")
    text = np.char.strip(text)
    numeric = np.full(text.shape, 9, dtype=np.int16)
    numeric[text == "0"] = 0
    numeric[text == "1"] = 1
    numeric[text == "4"] = 4
    numeric[text == "9"] = 9
    numeric[numeric == 0] = 1
    return numeric


def _to_string_qc(values: np.ndarray) -> np.ndarray:
    return values.astype("U1").astype(object)


def _clean_profile_indices(indices: Iterable[int] | np.ndarray, size: int) -> np.ndarray:
    array = np.asarray(list(indices) if not isinstance(indices, np.ndarray) else indices, dtype=np.int64).ravel()
    if array.size == 0:
        return np.array([], dtype=np.int64)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return np.array([], dtype=np.int64)
    array = np.unique(array.astype(np.int64))
    array = array[(array >= 0) & (array < size)]
    return array


def _set_profile_qc(qc: np.ndarray, indices: np.ndarray) -> None:
    if indices.size == 0:
        return
    current = qc[indices]
    qc[indices] = np.where(current != 9, 4, current)


def _filter_values(row: dict[str, str]) -> list[float]:
    values = [_as_float(row.get("x1")), _as_float(row.get("x2"))]
    if len(values) == 2 and np.isnan(values[1]):
        return [values[0]]
    return values


def _criterion_indices(state: _WorkingState, criterion: str, value: list[float], *, sal_only: bool) -> np.ndarray:
    if not value:
        return np.array([], dtype=np.int64)
    if len(value) == 2 and np.isnan(value[1]):
        value = [value[0]]

    pres = state.pres
    temp = state.temp
    psal = state.psal
    sigma0 = None

    def dens() -> np.ndarray:
        nonlocal sigma0
        if sigma0 is None:
            sigma0 = _sigma0_matrix(state, masked=False)
        return sigma0

    def any_profile(mask: np.ndarray) -> np.ndarray:
        return np.where(np.any(mask, axis=1))[0]

    match criterion:
        case "all":
            return np.arange(pres.shape[0], dtype=np.int64)
        case "index":
            # table_filter uses 1-based profile indices
            raw = np.asarray(value, dtype=float)
            raw = raw[np.isfinite(raw)]
            return np.unique(raw.astype(np.int64) - 1)
        case "Pmin":
            return any_profile(pres < value[0])
        case "Pmax":
            return any_profile(pres > value[0])
        case "Tmin":
            return any_profile(temp < value[0])
        case "Tmax":
            return any_profile(temp > value[0])
        case "Smin":
            return any_profile(psal < value[0])
        case "Smax":
            return any_profile(psal > value[0])
        case "Dmin":
            return any_profile(dens() < value[0])
        case "Dmax":
            return any_profile(dens() > value[0])
        case "P+D-":
            return any_profile((pres > value[0]) & (dens() < value[1]))
        case "P+D+":
            return any_profile((pres > value[0]) & (dens() > value[1]))
        case "P+S-":
            return any_profile((pres > value[0]) & (psal < value[1]))
        case "P+S+":
            return any_profile((pres > value[0]) & (psal > value[1]))
        case "P+T-":
            return any_profile((pres > value[0]) & (temp < value[1]))
        case "P+T+":
            return any_profile((pres > value[0]) & (temp > value[1]))
        case "P-S-":
            return any_profile((pres < value[0]) & (psal < value[1]))
        case "P-S+":
            return any_profile((pres < value[0]) & (psal > value[1]))
        case "P-T-":
            return any_profile((pres < value[0]) & (temp < value[1]))
        case "P-T+":
            return any_profile((pres < value[0]) & (temp > value[1]))
        case "P-D-":
            return any_profile((pres < value[0]) & (dens() < value[1]))
        case "P-D+":
            return any_profile((pres < value[0]) & (dens() > value[1]))
        case "T+S-":
            threshold = value[1] if sal_only and len(value) > 1 else value[0]
            return any_profile((temp > threshold) & (psal < value[1]))
        case "T+S+":
            threshold = value[1] if sal_only and len(value) > 1 else value[0]
            return any_profile((temp > threshold) & (psal > value[1]))
        case "T-S-":
            threshold = value[1] if sal_only and len(value) > 1 else value[0]
            return any_profile((temp < threshold) & (psal < value[1]))
        case "T-S+":
            threshold = value[1] if sal_only and len(value) > 1 else value[0]
            return any_profile((temp < threshold) & (psal > value[1]))
        case "S-D-":
            return any_profile((psal < value[0]) & (dens() < value[1]))
        case "S+D+":
            return any_profile((psal > value[0]) & (dens() > value[1]))
        case "S-D+":
            return any_profile((psal < value[0]) & (dens() > value[1]))
        case "S+D-":
            return any_profile((psal > value[0]) & (dens() < value[1]))
        case "date_min":
            relative = state.juld - np.nanmin(state.juld)
            return np.where(relative < value[0])[0]
        case "date_max":
            relative = state.juld - np.nanmin(state.juld)
            return np.where(relative > value[0])[0]
        case "lat_max":
            return np.where(state.lat < value[0])[0]
        case "lat_min":
            return np.where(state.lat > value[0])[0]
        case _:
            return np.array([], dtype=np.int64)


def _top_to_bottom_inversion_indices(state: _WorkingState, n_psal: np.ndarray) -> np.ndarray:
    sigma0 = _sigma0_matrix(state, masked=True)
    max_pressure_index = _safe_nanargmax(state.pres_work, axis=1)
    bottom_index = np.maximum(max_pressure_index - 1, 0)
    profile_index = np.arange(state.pres_work.shape[0])
    second_temp = state.temp_work[:, 1] if state.temp_work.shape[1] > 1 else np.full(state.temp_work.shape[0], np.nan)
    bottom_sigma = sigma0[profile_index, bottom_index]
    second_sigma = sigma0[:, 1] if sigma0.shape[1] > 1 else np.full(state.pres_work.shape[0], np.nan)
    mask = (second_temp < 7) & (n_psal > 5) & np.isfinite(bottom_sigma) & np.isfinite(second_sigma) & ((bottom_sigma - second_sigma) < -0.03)
    return np.where(mask)[0]


def _safe_nanmax(values: np.ndarray, axis: int) -> np.ndarray:
    work = np.where(np.isfinite(values), values, -np.inf)
    out = np.max(work, axis=axis)
    out = out.astype(np.float64)
    out[out == -np.inf] = np.nan
    return out


def _safe_nanargmax(values: np.ndarray, axis: int) -> np.ndarray:
    work = np.where(np.isfinite(values), values, -np.inf)
    index = np.argmax(work, axis=axis)
    return index.astype(np.int64)


def _sigma0_matrix(state: _WorkingState, *, masked: bool) -> np.ndarray:
    psal = state.psal_work if masked else state.psal
    temp_field = state.temp_work if masked else state.temp
    pres_field = state.pres_work if masked else state.pres
    result = np.full(psal.shape, np.nan, dtype=np.float64)
    for idx in range(psal.shape[0]):
        mask = np.isfinite(psal[idx]) & np.isfinite(temp_field[idx]) & np.isfinite(pres_field[idx])
        if not np.any(mask):
            continue
        sp = psal[idx, mask]
        temp = temp_field[idx, mask]
        pres = pres_field[idx, mask]
        lon = float(state.lon[idx]) if np.isfinite(state.lon[idx]) else 0.0
        lat = float(state.lat[idx]) if np.isfinite(state.lat[idx]) else 0.0
        result[idx, mask] = _sigma0_profile(sp, temp, pres, lon=lon, lat=lat)
    return result


def _sigma0_profile(sp: np.ndarray, temp: np.ndarray, pres: np.ndarray, *, lon: float, lat: float) -> np.ndarray:
    if _gsw is not None:  # pragma: no branch - deterministic branch when dependency is installed
        sa = _gsw.SA_from_SP(sp, pres, lon, lat)
        ct = _gsw.CT_from_t(sa, temp, pres)
        return np.asarray(_gsw.sigma0(sa, ct), dtype=np.float64)
    return _sigma0_eos80_fallback(sp, temp)


def _sigma0_eos80_fallback(sp: np.ndarray, temp: np.ndarray) -> np.ndarray:
    # EOS-80 density at atmospheric pressure, used only when the gsw package is unavailable.
    rho_w = ((((6.536332e-9 * temp - 1.120083e-6) * temp + 1.001685e-4) * temp - 9.095290e-3) * temp + 6.793952e-2) * temp + 999.842594
    a = (((5.3875e-9 * temp - 8.2467e-7) * temp + 7.6438e-5) * temp - 4.0899e-3) * temp + 0.824493
    b = (-1.6546e-6 * temp + 1.0227e-4) * temp - 5.72466e-3
    c = 4.8314e-4
    rho = rho_w + a * sp + b * np.power(sp, 1.5) + c * sp * sp
    return rho - 1000.0


def _update_count_and_geo_attributes(dataset: xr.Dataset, state: _WorkingState) -> None:
    temp_good = state.temp_qc == 1
    psal_good = state.psal_qc == 1
    dataset.attrs["number_of_ts_profiles"] = float(np.sum(np.any(temp_good & psal_good, axis=1)))
    dataset.attrs["number_of_t_profiles"] = float(np.sum(np.any(temp_good, axis=1)))
    if state.chla_qc is not None:
        dataset.attrs["number_chla_profiles"] = float(np.sum(np.any(state.chla_qc == 1, axis=1)))
    if state.doxy_qc is not None:
        dataset.attrs["number_doxy_profiles"] = float(np.sum(np.any(state.doxy_qc == 1, axis=1)))
    if state.light_qc is not None:
        dataset.attrs["number_light_profiles"] = float(np.sum(np.any(state.light_qc == 1, axis=1)))

    valid_profiles = np.where(np.any(temp_good | psal_good, axis=1))[0]
    if valid_profiles.size == 0:
        dataset.attrs["geospatial_lat_min"] = " "
        dataset.attrs["geospatial_lat_max"] = " "
        dataset.attrs["geospatial_lon_min"] = " "
        dataset.attrs["geospatial_lon_max"] = " "
        return

    lat = state.lat[valid_profiles]
    lon = state.lon[valid_profiles].copy()
    lon = np.where(lon < 0, lon + 360.0, lon)
    dataset.attrs["geospatial_lat_min"] = float(np.nanmin(lat))
    dataset.attrs["geospatial_lat_max"] = float(np.nanmax(lat))
    dataset.attrs["geospatial_lon_min"] = float(np.nanmin(lon))
    dataset.attrs["geospatial_lon_max"] = float(np.nanmax(lon))


def _as_float(value: str | None, default: float = np.nan) -> float:
    if value is None:
        return float(default)
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return float(default)
    try:
        return float(text)
    except ValueError:
        return float(default)


def _as_bool(value: str | None) -> bool:
    text = (value or "").strip().lower()
    return text in {"1", "true", "y", "yes"}
