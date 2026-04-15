from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

from ..catalog.filenames import fname_prof, list_fname_prof
from ..models import MeopConfig, Selection
from .conductivity import smooth_profile, thermal_cell_correction
from .qc import _WorkingState, _to_numeric_qc
from .stabilise_sa_const_ct import SA_SCALE, stabilise_SP_const_CT


try:  # pragma: no cover - exercised only when gsw is installed
    import gsw as _gsw  # type: ignore
except Exception:  # pragma: no cover - current container does not provide gsw
    _gsw = None


STANDARD_HR_LEVELS = np.arange(1, 1001, dtype=np.float32)


@dataclass(frozen=True)
class HrResult:
    written_files: tuple[Path, ...]
    processed_tags: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "written_files": [str(path) for path in self.written_files],
            "processed_tags": list(self.processed_tags),
        }


class Hr1PortPending(RuntimeError):
    """Raised when the requested hr1 pathway still depends on a not-yet-ported logic branch."""



def _open_dataset(path: Path) -> xr.Dataset:
    try:
        return xr.load_dataset(path, engine="h5netcdf", decode_times=False)
    except Exception:
        return xr.load_dataset(path, decode_times=False)



def _selected_tags(config: MeopConfig, selection: Selection, *, source_qf: str) -> tuple[str, ...]:
    if selection.smru_name:
        return (selection.smru_name,)
    files = list_fname_prof(deployment=selection.deployment, qf=source_qf, config=config)
    return tuple(path.name.split("_")[0] for path in files)



def _as_utc(now: datetime | None) -> datetime:
    return (now or datetime.now(timezone.utc)).astimezone(timezone.utc)



def _update_update_timestamps(dataset: xr.Dataset, now: datetime) -> None:
    if "DATE_UPDATE" in dataset:
        dataset["DATE_UPDATE"] = xr.DataArray(np.asarray(now.strftime("%Y%m%d%H%M%S"), dtype=object))
    dataset.attrs["date_update"] = now.strftime("%Y-%m-%dT%H:%M:%SZ")



def _prepare_profile_input(values: np.ndarray, pressure: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mask = np.isfinite(values) & np.isfinite(pressure)
    if np.count_nonzero(mask) <= 5:
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)
    z = pressure[mask].astype(np.float64, copy=False)
    x = values[mask].astype(np.float64, copy=False)
    order = np.argsort(z)
    z = z[order]
    x = x[order]
    keep = np.r_[np.diff(z) != 0, True]
    z = z[keep]
    x = x[keep]
    if z.size <= 5:
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)
    return z, x



def _interp_linear(depth: np.ndarray, values: np.ndarray, levels: np.ndarray) -> np.ndarray:
    if depth.size <= 5:
        return np.full(levels.shape, np.nan, dtype=np.float64)
    return np.interp(levels, depth, values, left=np.nan, right=np.nan).astype(np.float64)



def _interp_temp_psal(
    pressure: np.ndarray,
    temp: np.ndarray,
    psal: np.ndarray,
    *,
    lon: float,
    lat: float,
    levels: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    z_t, t = _prepare_profile_input(temp, pressure)
    z_s, s = _prepare_profile_input(psal, pressure)

    t_std = _interp_linear(z_t, t, levels)
    s_std = _interp_linear(z_s, s, levels)
    # Keep the default HR interpolation aligned with the historical workflow: no extrapolation
    # above the shallowest observed level and no TEOS-10 interpolation at this stage. The GSW
    # package is still used downstream for CT/sigma0-based stabilisation.
    return t_std, s_std



def _interp_sensor(pressure: np.ndarray, values: np.ndarray, levels: np.ndarray) -> np.ndarray:
    z, x = _prepare_profile_input(values, pressure)
    return _interp_linear(z, x, levels)



def _profile_qc(level_qc: np.ndarray) -> np.ndarray:
    return np.where(np.any(level_qc == b"1", axis=1), b"A", b"F").astype(object)



def _full_like_fill(source: xr.DataArray, shape: tuple[int, ...]) -> np.ndarray:
    fill_value = source.attrs.get("_FillValue")
    if fill_value is None:
        if source.dtype.kind in {"U", "S", "O"}:
            fill_value = " "
        else:
            fill_value = np.nan
    dtype = object if source.dtype.kind in {"U", "S", "O"} else source.dtype
    return np.full(shape, fill_value, dtype=dtype)



def _copy_non_level_variables(dataset: xr.Dataset) -> xr.Dataset:
    result = xr.Dataset(attrs=dict(dataset.attrs))
    for name, coord in dataset.coords.items():
        if "N_LEVELS" not in coord.dims:
            result = result.assign_coords({name: coord.copy(deep=True)})
    for name, variable in dataset.data_vars.items():
        if "N_LEVELS" not in variable.dims:
            result[name] = variable.copy(deep=True)
    return result



def _build_hr0_dataset(source: xr.Dataset, *, now: datetime) -> xr.Dataset:
    levels = STANDARD_HR_LEVELS.astype(np.float64)
    result = _copy_non_level_variables(source)

    n_prof = int(source.sizes.get("N_PROF", 0))
    n_levels = int(levels.size)
    shape = (n_prof, n_levels)

    pressure = np.broadcast_to(levels[None, :], shape).astype(np.float32)
    pres_qc = np.full(shape, b"1", dtype=object)

    pres_src = np.asarray(source["PRES"].values, dtype=np.float64)
    temp_src = np.asarray(source["TEMP"].values, dtype=np.float64)
    psal_src = np.asarray(source["PSAL"].values, dtype=np.float64)

    pres_mask = _to_numeric_qc(source["PRES_QC"].values) == 1
    temp_mask = _to_numeric_qc(source["TEMP_QC"].values) == 1
    psal_mask = _to_numeric_qc(source["PSAL_QC"].values) == 1

    pres_work = np.where(pres_mask, pres_src, np.nan)
    temp_work = np.where(temp_mask, temp_src, np.nan)
    psal_work = np.where(psal_mask, psal_src, np.nan)

    temp = np.full(shape, np.nan, dtype=np.float32)
    psal = np.full(shape, np.nan, dtype=np.float32)
    temp_qc = np.full(shape, b"9", dtype=object)
    psal_qc = np.full(shape, b"9", dtype=object)

    optional_values: dict[str, np.ndarray] = {}
    optional_qc: dict[str, np.ndarray] = {}
    for optional in ("CHLA", "DOXY", "LIGHT"):
        if optional not in source:
            continue
        optional_values[optional] = np.full(shape, np.nan, dtype=np.float32)
        optional_qc[optional] = np.full(shape, b"9", dtype=object)

    lat = np.asarray(source["LATITUDE"].values, dtype=np.float64)
    lon = np.asarray(source["LONGITUDE"].values, dtype=np.float64)

    for idx in range(n_prof):
        p = pres_work[idx]
        t = temp_work[idx]
        s = psal_work[idx]
        t_i, s_i = _interp_temp_psal(p, t, s, lon=float(lon[idx]), lat=float(lat[idx]), levels=levels)
        temp[idx, :] = t_i.astype(np.float32)
        psal[idx, :] = s_i.astype(np.float32)
        temp_qc[idx, np.isfinite(t_i)] = b"1"
        psal_qc[idx, np.isfinite(s_i)] = b"1"

        for optional in optional_values:
            source_values = np.asarray(source[optional].values, dtype=np.float64)
            source_qc = _to_numeric_qc(source[f"{optional}_QC"].values) == 1
            source_work = np.where(source_qc, source_values, np.nan)
            value_i = _interp_sensor(p, source_work[idx], levels)
            optional_values[optional][idx, :] = value_i.astype(np.float32)
            optional_qc[optional][idx, np.isfinite(value_i)] = b"1"

    arrays: dict[str, np.ndarray] = {
        "PRES": pressure,
        "PRES_ADJUSTED": pressure.copy(),
        "PRES_QC": pres_qc,
        "PRES_ADJUSTED_QC": pres_qc.copy(),
        "TEMP": temp,
        "TEMP_ADJUSTED": temp.copy(),
        "TEMP_QC": temp_qc,
        "TEMP_ADJUSTED_QC": temp_qc.copy(),
        "PSAL": psal,
        "PSAL_ADJUSTED": psal.copy(),
        "PSAL_QC": psal_qc,
        "PSAL_ADJUSTED_QC": psal_qc.copy(),
        "PROFILE_PRES_QC": _profile_qc(pres_qc),
        "PROFILE_TEMP_QC": _profile_qc(temp_qc),
        "PROFILE_PSAL_QC": _profile_qc(psal_qc),
    }
    for optional in optional_values:
        arrays[optional] = optional_values[optional]
        arrays[f"{optional}_ADJUSTED"] = optional_values[optional].copy()
        arrays[f"{optional}_QC"] = optional_qc[optional]
        arrays[f"{optional}_ADJUSTED_QC"] = optional_qc[optional].copy()

    for name, variable in source.data_vars.items():
        if "N_LEVELS" not in variable.dims:
            continue
        values = arrays.get(name)
        if values is None:
            values = _full_like_fill(variable, shape)
        result[name] = (variable.dims, values)
        result[name].attrs.update(variable.attrs)

    state = _WorkingState(result)
    result = state.to_dataset()
    result.attrs["thermal_lag_adjustment"] = "no"
    _update_update_timestamps(result, now)
    return result



def _choose_adjusted_source(dataset: xr.Dataset, name: str) -> np.ndarray:
    adjusted = f"{name}_ADJUSTED"
    if adjusted in dataset:
        values = np.asarray(dataset[adjusted].values, dtype=np.float64)
        if np.isfinite(values).any():
            return values
    return np.asarray(dataset[name].values, dtype=np.float64)



def _choose_adjusted_qc(dataset: xr.Dataset, name: str) -> np.ndarray:
    adjusted_qc = f"{name}_ADJUSTED_QC"
    if adjusted_qc in dataset:
        return _to_numeric_qc(dataset[adjusted_qc].values)
    return _to_numeric_qc(dataset[f"{name}_QC"].values)



def _compute_ct_from_sp(sp: np.ndarray, temp: np.ndarray, pres: np.ndarray) -> tuple[np.ndarray, str]:
    if _gsw is not None:
        try:  # pragma: no cover - only exercised when gsw is installed
            sa = np.asarray(sp, dtype=np.float64) * SA_SCALE
            ct = np.asarray(_gsw.CT_from_t(sa, temp, pres), dtype=np.float64)
            return ct, "gsw"
        except Exception:
            pass
    return np.asarray(temp, dtype=np.float64), "temperature-fallback"



def _smooth_valid_profiles(values: np.ndarray) -> np.ndarray:
    out = np.asarray(values, dtype=np.float64).copy()
    for idx in range(out.shape[0]):
        out[idx, :] = smooth_profile(out[idx, :])
    return out



def _build_hr1_dataset(source: xr.Dataset, *, thermal_lag: bool, now: datetime) -> xr.Dataset:
    result = source.copy(deep=True)
    pres_source = _choose_adjusted_source(source, "PRES")
    temp_source = _choose_adjusted_source(source, "TEMP")
    psal_source = _choose_adjusted_source(source, "PSAL")

    pres_qc = _choose_adjusted_qc(source, "PRES")
    temp_qc = _choose_adjusted_qc(source, "TEMP")
    psal_qc = _choose_adjusted_qc(source, "PSAL")

    pres_adjusted = np.where(pres_qc == 1, pres_source, np.nan).astype(np.float64)
    temp_adjusted = np.where(temp_qc == 1, temp_source, np.nan).astype(np.float64)
    psal_adjusted = np.where(psal_qc == 1, psal_source, np.nan).astype(np.float64)

    if thermal_lag and _gsw is not None:
        corrected = thermal_cell_correction(psal_adjusted, temp_adjusted, pres_adjusted)
        temp_corrected = corrected.temperature
        psal_corrected = corrected.salinity
    else:
        temp_corrected = temp_adjusted.copy()
        psal_corrected = psal_adjusted.copy()
        corrected = None

    ct, ct_method = _compute_ct_from_sp(psal_corrected, temp_corrected, pres_adjusted)
    psal_stable_levels, metadata = stabilise_SP_const_CT(
        psal_corrected.T,
        ct.T,
        pres_adjusted.T,
        return_metadata=True,
    )
    psal_stable = np.asarray(psal_stable_levels, dtype=np.float64).T

    temp_smoothed = _smooth_valid_profiles(temp_corrected)
    psal_smoothed = _smooth_valid_profiles(np.asarray(psal_stable, dtype=np.float64))

    result["TEMP_ADJUSTED"] = (source["TEMP_ADJUSTED"].dims, temp_smoothed.astype(np.float32))
    result["PSAL_ADJUSTED"] = (source["PSAL_ADJUSTED"].dims, psal_smoothed.astype(np.float32))
    result["TEMP_ADJUSTED"].attrs.update(source["TEMP_ADJUSTED"].attrs)
    result["PSAL_ADJUSTED"].attrs.update(source["PSAL_ADJUSTED"].attrs)
    result.attrs["thermal_lag_adjustment"] = "yes" if thermal_lag else "no"
    if thermal_lag:
        result.attrs["density_inversion_adjustment"] = "python-gsw" if _gsw is not None else "python-fallback"
        result.attrs["thermal_lag_method"] = getattr(corrected, "method", "disabled-without-gsw")
    else:
        result.attrs["density_inversion_adjustment"] = "python-gsw" if _gsw is not None else "python-fallback"
    result.attrs["ct_method"] = ct_method
    successful = sum(1 for item in metadata if item is not None and item.success)
    if successful:
        result.attrs["stabilised_profiles"] = int(successful)
    _update_update_timestamps(result, now)
    return result



def _nearest_profile_projection(source_pressure: np.ndarray, source_values: np.ndarray, target_pressure: np.ndarray) -> np.ndarray:
    out = np.full(target_pressure.shape, np.nan, dtype=np.float64)
    mask = np.isfinite(source_pressure) & np.isfinite(source_values)
    if np.count_nonzero(mask) <= 1:
        return out
    p_src = source_pressure[mask].astype(np.float64, copy=False)
    x_src = source_values[mask].astype(np.float64, copy=False)
    order = np.argsort(p_src)
    p_src = p_src[order]
    x_src = x_src[order]
    finite_target = np.isfinite(target_pressure)
    if not np.any(finite_target):
        return out
    p_tgt = target_pressure[finite_target].astype(np.float64, copy=False)
    inside = (p_tgt >= p_src[0]) & (p_tgt <= p_src[-1])
    if not np.any(inside):
        return out
    p_inside = p_tgt[inside]
    insert = np.searchsorted(p_src, p_inside, side="left")
    insert = np.clip(insert, 1, p_src.size - 1)
    left = insert - 1
    right = insert
    choose_right = np.abs(p_src[right] - p_inside) <= np.abs(p_inside - p_src[left])
    nearest = np.where(choose_right, right, left)
    values = np.full(p_tgt.shape, np.nan, dtype=np.float64)
    values[inside] = x_src[nearest]
    out[finite_target] = values
    return out



def _build_lr1_dataset(lr0: xr.Dataset, hr1: xr.Dataset, *, thermal_lag: bool, now: datetime) -> xr.Dataset:
    result = lr0.copy(deep=True)
    n_prof = int(result.sizes.get("N_PROF", 0))

    lr_pres = np.where(_choose_adjusted_qc(lr0, "PRES") == 1, _choose_adjusted_source(lr0, "PRES"), np.nan)
    lr_temp = _choose_adjusted_source(lr0, "TEMP")
    lr_psal = _choose_adjusted_source(lr0, "PSAL")
    hr_pres = np.where(_choose_adjusted_qc(hr1, "PRES") == 1, _choose_adjusted_source(hr1, "PRES"), np.nan)
    hr_temp = np.where(_choose_adjusted_qc(hr1, "TEMP") == 1, _choose_adjusted_source(hr1, "TEMP"), np.nan)
    hr_psal = np.where(_choose_adjusted_qc(hr1, "PSAL") == 1, _choose_adjusted_source(hr1, "PSAL"), np.nan)

    temp_adjusted = np.asarray(lr_temp, dtype=np.float64).copy()
    psal_adjusted = np.asarray(lr_psal, dtype=np.float64).copy()

    for idx in range(n_prof):
        projected_temp = _nearest_profile_projection(hr_pres[idx], hr_temp[idx], lr_pres[idx])
        projected_psal = _nearest_profile_projection(hr_pres[idx], hr_psal[idx], lr_pres[idx])

        replace_t = np.isfinite(projected_temp)
        replace_s = np.isfinite(projected_psal)
        temp_adjusted[idx, replace_t] = projected_temp[replace_t]
        psal_adjusted[idx, replace_s] = projected_psal[replace_s]

    result["TEMP_ADJUSTED"] = (lr0["TEMP_ADJUSTED"].dims, temp_adjusted.astype(np.float32))
    result["PSAL_ADJUSTED"] = (lr0["PSAL_ADJUSTED"].dims, psal_adjusted.astype(np.float32))
    result["TEMP_ADJUSTED"].attrs.update(lr0["TEMP_ADJUSTED"].attrs)
    result["PSAL_ADJUSTED"].attrs.update(lr0["PSAL_ADJUSTED"].attrs)
    result.attrs["thermal_lag_adjustment"] = "yes" if thermal_lag else "no"
    _update_update_timestamps(result, now)
    return result



def create_hr0_python(config: MeopConfig, selection: Selection, *, now: datetime | None = None) -> HrResult:
    tags = _selected_tags(config, selection, source_qf="lr0")
    written: list[Path] = []
    timestamp = _as_utc(now)
    for smru_name in tags:
        source_path = fname_prof(smru_name, deployment=selection.deployment, qf="lr0", config=config)
        if not source_path.exists():
            continue
        dataset = _open_dataset(source_path)
        try:
            target_dataset = _build_hr0_dataset(dataset, now=timestamp)
            target_path = fname_prof(smru_name, deployment=selection.deployment, qf="hr0", config=config)
            target_path.parent.mkdir(parents=True, exist_ok=True)
            target_dataset.to_netcdf(target_path, engine="h5netcdf")
            written.append(target_path)
        finally:
            dataset.close()
    return HrResult(
        written_files=tuple(written),
        processed_tags=tuple(
            tag for tag in tags if fname_prof(tag, deployment=selection.deployment, qf="hr0", config=config).exists()
        ),
    )



def _create_adjusted_profile_dataset(
    config: MeopConfig,
    selection: Selection,
    *,
    source_qf: str,
    target_qf: str,
    thermal_lag: bool,
    now: datetime | None = None,
    projected_source_qf: str | None = None,
    projected_target_qf: str | None = None,
) -> HrResult:
    tags = _selected_tags(config, selection, source_qf=source_qf)
    written: list[Path] = []
    processed: list[str] = []
    timestamp = _as_utc(now)

    for smru_name in tags:
        source_path = fname_prof(smru_name, deployment=selection.deployment, qf=source_qf, config=config)
        if not source_path.exists():
            continue

        source_dataset = _open_dataset(source_path)
        projected_dataset = None
        if projected_source_qf and projected_target_qf:
            projected_path = fname_prof(smru_name, deployment=selection.deployment, qf=projected_source_qf, config=config)
            if not projected_path.exists():
                source_dataset.close()
                continue
            projected_dataset = _open_dataset(projected_path)

        try:
            target_dataset = _build_hr1_dataset(source_dataset, thermal_lag=thermal_lag, now=timestamp)
            target_path = fname_prof(smru_name, deployment=selection.deployment, qf=target_qf, config=config)
            target_path.parent.mkdir(parents=True, exist_ok=True)
            target_dataset.to_netcdf(target_path, engine="h5netcdf")
            written.append(target_path)

            if projected_dataset is not None and projected_target_qf is not None:
                projected_target_dataset = _build_lr1_dataset(projected_dataset, target_dataset, thermal_lag=thermal_lag, now=timestamp)
                projected_target_path = fname_prof(smru_name, deployment=selection.deployment, qf=projected_target_qf, config=config)
                projected_target_path.parent.mkdir(parents=True, exist_ok=True)
                projected_target_dataset.to_netcdf(projected_target_path, engine="h5netcdf")
                written.append(projected_target_path)

            processed.append(smru_name)
        finally:
            source_dataset.close()
            if projected_dataset is not None:
                projected_dataset.close()

    return HrResult(
        written_files=tuple(written),
        processed_tags=tuple(processed),
    )



def create_hr1_python(
    config: MeopConfig,
    selection: Selection,
    *,
    thermal_lag: bool = False,
    now: datetime | None = None,
) -> HrResult:
    return _create_adjusted_profile_dataset(
        config,
        selection,
        source_qf="hr0",
        target_qf="hr1",
        projected_source_qf="lr0",
        projected_target_qf="lr1",
        thermal_lag=thermal_lag,
        now=now,
    )



def create_fr1_python(
    config: MeopConfig,
    selection: Selection,
    *,
    thermal_lag: bool = False,
    now: datetime | None = None,
) -> HrResult:
    return _create_adjusted_profile_dataset(
        config,
        selection,
        source_qf="fr0",
        target_qf="fr1",
        thermal_lag=thermal_lag,
        now=now,
    )
