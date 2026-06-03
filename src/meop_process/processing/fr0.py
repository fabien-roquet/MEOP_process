from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import xarray as xr

from ..catalog.deployments import load_info_deployment
from ..catalog.filenames import fname_prof, fname_traj
from ..io.hr_ctd import resolve_hr_ctd_path
from ..io.raw_odv import OdvProfile
from ..models import MeopConfig, Selection
from .hr import STANDARD_HR_LEVELS, _interp_sensor, _prepare_profile_input
from .netcdf import save_dataset_with_compression
from .ncargo import NcargoResult, _build_ncargo_dataset
from .qc import apply_lr0_qc_filters


try:  # pragma: no cover - exercised only when gsw is available
    import gsw as _gsw  # type: ignore
except Exception:  # pragma: no cover - fallback path remains testable without gsw
    _gsw = None


DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
JULD_REF = datetime(1950, 1, 1, tzinfo=timezone.utc)
PSEUIL_DEPTH = 15.0
MIN_DIVE_DURATION_SECONDS = 300.0
SEGMENT_GAP_MIN_SECONDS = 60.0
SURFACE_END_THRESHOLD = 6.0
MIN_SEGMENT_SAMPLES = 30


@dataclass(frozen=True)
class HrRawData:
    timestamp: pd.DatetimeIndex
    juld: np.ndarray
    pressure: np.ndarray
    temperature: np.ndarray
    salinity: np.ndarray
    conductivity: np.ndarray
    fluorescence: np.ndarray
    oxygen: np.ndarray
    light: np.ndarray
    continuous: bool
    isfluo: bool
    isoxy: bool
    islight: bool

    @property
    def n_samples(self) -> int:
        return int(self.pressure.size)


@dataclass(frozen=True)
class DetectedProfiles:
    ascent_start_index: np.ndarray
    profile_end_index: np.ndarray

    @property
    def n_profiles(self) -> int:
        return int(self.ascent_start_index.size)


def _normalize_header(name: str) -> str:
    text = str(name).strip().lower()
    for token in ("[", "]", "(", ")", "/"):
        text = text.replace(token, " ")
    return " ".join(text.split())


def _pick_column(columns: list[str], kind: str) -> str | None:
    normalized = {column: _normalize_header(column) for column in columns}

    if kind == "date":
        for column, name in normalized.items():
            if name.startswith("date"):
                return column
        return None

    if kind == "pressure":
        exact = [column for column, name in normalized.items() if name in {"pressure dbar", "corrected depth", "pressure_re_surface dbar"}]
        if exact:
            return exact[0]
        general = [
            column
            for column, name in normalized.items()
            if name.startswith("pressure") and "surface" not in name and "re surface" not in name
        ]
        if general:
            return general[0]
        return None

    if kind == "temperature":
        for column, name in normalized.items():
            if name.startswith("temp") or name.startswith("temperature"):
                return column
        return None

    if kind == "salinity":
        preferred = [column for column, name in normalized.items() if name.startswith("sal lagged")]
        if preferred:
            return preferred[0]
        for column, name in normalized.items():
            if name.startswith("salinity") or name.startswith("sal"):
                return column
        return None

    if kind == "conductivity":
        for column, name in normalized.items():
            if name.startswith("cond") or name.startswith("conductivity"):
                return column
        return None

    if kind == "fluorescence":
        for column, name in normalized.items():
            if "fluo" in name or "fluorescence" in name:
                return column
        return None

    if kind == "oxygen":
        for column, name in normalized.items():
            if "oxy" in name or "oxygen" in name:
                return column
        return None

    if kind == "light":
        for column, name in normalized.items():
            if "ppfd" in name or "light" in name:
                return column
        return None

    raise ValueError(f"Unknown HR column kind: {kind}")


def _to_numeric(series: pd.Series | None, size: int) -> np.ndarray:
    if series is None:
        return np.full(size, np.nan, dtype=np.float64)
    values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=np.float64, copy=True)
    values[np.isclose(values, 999.0, equal_nan=False)] = np.nan
    return values


def _parse_hr_datetime(series: pd.Series) -> pd.DatetimeIndex:
    values = pd.to_datetime(series, errors="coerce", utc=True)
    if values.isna().any():
        raise ValueError("Failed to parse one or more HR timestamps")
    return pd.DatetimeIndex(values)


def load_hr_data(path: str | Path, *, continuous: bool) -> HrRawData:
    source = Path(path)
    try:
        dataframe = pd.read_csv(source, sep="\t")
    except pd.errors.ParserError as exc:
        warnings.warn(
            f"Malformed HR rows skipped while reading {source.name}: {exc}",
            stacklevel=2,
        )
        dataframe = pd.read_csv(source, sep="\t", on_bad_lines="skip")
    columns = list(dataframe.columns)

    date_column = _pick_column(columns, "date")
    pressure_column = _pick_column(columns, "pressure")
    temp_column = _pick_column(columns, "temperature")
    salinity_column = _pick_column(columns, "salinity")
    conductivity_column = _pick_column(columns, "conductivity")
    fluorescence_column = _pick_column(columns, "fluorescence")
    oxygen_column = _pick_column(columns, "oxygen")
    light_column = _pick_column(columns, "light")

    if date_column is None or pressure_column is None or temp_column is None:
        raise ValueError(f"Missing required HR columns in {source.name}")

    timestamp = _parse_hr_datetime(dataframe[date_column])
    pressure = _to_numeric(dataframe[pressure_column], len(dataframe))
    temperature = _to_numeric(dataframe[temp_column], len(dataframe))
    salinity = _to_numeric(dataframe[salinity_column] if salinity_column else None, len(dataframe))
    conductivity = _to_numeric(dataframe[conductivity_column] if conductivity_column else None, len(dataframe))
    fluorescence = _to_numeric(dataframe[fluorescence_column] if fluorescence_column else None, len(dataframe))
    oxygen = _to_numeric(dataframe[oxygen_column] if oxygen_column else None, len(dataframe))
    light = _to_numeric(dataframe[light_column] if light_column else None, len(dataframe))

    unix_seconds = timestamp.astype("int64") / 1e9
    keep_monotonic = np.r_[True, np.diff(unix_seconds) > 0]

    timestamp = timestamp[keep_monotonic]
    pressure = pressure[keep_monotonic]
    temperature = temperature[keep_monotonic]
    salinity = salinity[keep_monotonic]
    conductivity = conductivity[keep_monotonic]
    fluorescence = fluorescence[keep_monotonic]
    oxygen = oxygen[keep_monotonic]
    light = light[keep_monotonic]

    salinity[salinity < 3.0] = np.nan
    ref_ns = pd.Timestamp(JULD_REF).value
    timestamp_ns = timestamp.astype("int64")  # already in nanoseconds
    juld = (timestamp_ns - ref_ns).astype(np.float64) / 1e9 / 86400.0

    return HrRawData(
        timestamp=timestamp,
        juld=juld,
        pressure=pressure,
        temperature=temperature,
        salinity=salinity,
        conductivity=conductivity,
        fluorescence=fluorescence,
        oxygen=oxygen,
        light=light,
        continuous=continuous,
        isfluo=bool(np.isfinite(fluorescence).any()),
        isoxy=bool(np.isfinite(oxygen).any()),
        islight=bool(np.isfinite(light).any()),
    )


def _resolution_seconds(juld: np.ndarray) -> float:
    if juld.size < 2:
        return 1.0
    diffs = np.diff(juld) * 86400.0
    diffs = diffs[np.isfinite(diffs) & (diffs > 0)]
    if diffs.size == 0:
        return 1.0
    return float(np.nanmedian(diffs))


def _moving_mean_legacy(values: np.ndarray, window_half_width: int = 6) -> np.ndarray:
    smoothed = np.asarray(values, dtype=np.float64).copy()
    start = 1 + window_half_width
    stop = smoothed.size - window_half_width
    for idx in range(start, stop):
        smoothed[idx] = np.nanmean(values[idx - window_half_width : idx + window_half_width + 1])
    return smoothed


def _id_dives(juld: np.ndarray, pressure: np.ndarray, *, pseuil: float = PSEUIL_DEPTH, dureed: float = MIN_DIVE_DURATION_SECONDS) -> tuple[np.ndarray, np.ndarray, float]:
    resolution = _resolution_seconds(juld)
    states = np.where(np.asarray(pressure, dtype=np.float64) >= pseuil, -1, 1)

    starts: list[int] = []
    ends: list[int] = []
    armed = False
    for idx in range(states.size - 1):
        delta = int(states[idx] - states[idx + 1])
        if delta == 2 and not armed:
            starts.append(idx)
            armed = True
        elif delta == -2 and armed:
            ends.append(idx)
            armed = False

    if len(starts) > len(ends):
        starts = starts[: len(ends)]

    filtered_starts: list[int] = []
    filtered_ends: list[int] = []
    for start, end in zip(starts, ends):
        if (end - start) * resolution >= dureed:
            filtered_starts.append(start)
            filtered_ends.append(end)

    return np.asarray(filtered_starts, dtype=np.int64), np.asarray(filtered_ends, dtype=np.int64), resolution


def _poly_4_delim_fabien(pressure: np.ndarray, starts: np.ndarray, ends: np.ndarray, *, resolution: float) -> np.ndarray:
    vert_speed = np.zeros_like(pressure, dtype=np.float64)
    if pressure.size > 1:
        vert_speed[1:] = np.diff(pressure) / max(resolution, 1e-6)
    vert_speed = _moving_mean_legacy(vert_speed)

    crit = 0.75
    delim = np.zeros((2, starts.size), dtype=np.int64)
    for profile_index, (start, end) in enumerate(zip(starts, ends)):
        local_speed = vert_speed[start : end + 1]
        local_pressure = pressure[start : end + 1]
        max_local_index = int(np.nanargmax(local_pressure))
        max_depth = float(local_pressure[max_local_index])
        depthmax_index = start + max_local_index

        x = np.arange(1, local_speed.size + 1, dtype=np.float64)
        ok = np.isfinite(local_speed)
        if np.count_nonzero(ok) >= 2:
            degree = min(4, np.count_nonzero(ok) - 1)
            coeffs = np.polyfit(x[ok], local_speed[ok], deg=degree)
            vf = np.polyval(coeffs, x)
        else:
            vf = np.zeros_like(x)

        phase = np.full(local_speed.shape, -1.0, dtype=np.float64)
        phase[vf >= crit] = -0.5
        phase[vf < -crit] = 0.5

        bottom = np.flatnonzero(phase == -1)
        ascent = np.flatnonzero(phase == 0.5)

        delim_start = max(depthmax_index - 1, start)
        if bottom.size:
            candidates = start + bottom
            deep_enough = candidates[pressure[candidates] / max(max_depth, 1e-6) >= 0.4]
            if deep_enough.size:
                delim_start = int(deep_enough[0])

        delim_end = min(depthmax_index + 1, end)
        if ascent.size:
            candidate = int(start + ascent[0])
            if pressure[candidate] / max(max_depth, 1e-6) >= 0.4:
                delim_end = candidate
            else:
                backtracked = candidate
                while backtracked > delim_start and pressure[backtracked] / max(max_depth, 1e-6) < 0.4:
                    backtracked -= 1
                delim_end = backtracked if backtracked > delim_start else min(depthmax_index + 1, end)

        delim[0, profile_index] = int(np.clip(delim_start, start, end))
        delim[1, profile_index] = int(np.clip(delim_end, start, end))

    return delim


def _gap_threshold_seconds(raw: HrRawData) -> float:
    return max(SEGMENT_GAP_MIN_SECONDS, 30.0 * _resolution_seconds(raw.juld))


def _split_by_time_gaps(raw: HrRawData) -> list[tuple[int, int]]:
    if raw.n_samples == 0:
        return []
    threshold = _gap_threshold_seconds(raw)
    gap_seconds = np.diff(raw.juld) * 86400.0
    cut = np.flatnonzero(gap_seconds > threshold)
    starts = np.r_[0, cut + 1]
    ends = np.r_[cut, raw.n_samples - 1]
    return [(int(start), int(end)) for start, end in zip(starts, ends)]


def _segment_looks_complete(raw: HrRawData, start: int, end: int) -> bool:
    pressure = raw.pressure[start : end + 1]
    if pressure.size < MIN_SEGMENT_SAMPLES:
        return False
    finite = np.isfinite(pressure)
    if np.count_nonzero(finite) < MIN_SEGMENT_SAMPLES:
        return False
    pressure = pressure[finite]
    if np.nanmax(pressure) < PSEUIL_DEPTH:
        return False
    if pressure[-1] > SURFACE_END_THRESHOLD:
        return False
    return True


def detect_profiles(raw: HrRawData) -> DetectedProfiles:
    if raw.n_samples == 0:
        return DetectedProfiles(np.array([], dtype=np.int64), np.array([], dtype=np.int64))

    segmented = _split_by_time_gaps(raw)
    if len(segmented) > 1:
        ascent_start = []
        profile_end = []
        for start, end in segmented:
            if _segment_looks_complete(raw, start, end):
                ascent_start.append(start)
                profile_end.append(end)
        return DetectedProfiles(np.asarray(ascent_start, dtype=np.int64), np.asarray(profile_end, dtype=np.int64))

    if raw.continuous:
        starts, ends, resolution = _id_dives(raw.juld, raw.pressure)
        if starts.size == 0:
            return DetectedProfiles(np.array([], dtype=np.int64), np.array([], dtype=np.int64))

        delim = _poly_4_delim_fabien(raw.pressure, starts, ends, resolution=resolution)
        ascent_start = delim[1, :].astype(np.int64)

        profile_end: list[int] = []
        threshold = 2.0
        for idx in range(starts.size - 1):
            segment = raw.pressure[ends[idx] : starts[idx + 1] + 1]
            shallow = np.flatnonzero(segment < threshold)
            profile_end.append(int(ends[idx] + shallow[0]) if shallow.size else int(ends[idx]))
        profile_end.append(int(ends[-1]))
        profile_end_array = np.asarray(profile_end, dtype=np.int64)
    else:
        breaks = np.flatnonzero(np.abs(np.diff(raw.pressure)) > 10.0)
        ascent_start = np.r_[0, breaks + 1].astype(np.int64)
        profile_end_array = np.r_[breaks, raw.pressure.size - 1].astype(np.int64)

    if profile_end_array.size > 1:
        bad = np.flatnonzero(np.diff(raw.juld[profile_end_array]) < 0)
        if bad.size:
            keep = int(bad[0] + 1)
            ascent_start = ascent_start[:keep]
            profile_end_array = profile_end_array[:keep]

    valid = profile_end_array > ascent_start
    return DetectedProfiles(ascent_start[valid], profile_end_array[valid])


def _interp_linear_surface_fill(depth: np.ndarray, values: np.ndarray, levels: np.ndarray) -> np.ndarray:
    if depth.size <= 5:
        return np.full(levels.shape, np.nan, dtype=np.float64)
    out = np.interp(levels, depth, values, left=values[0], right=np.nan).astype(np.float64)
    out[levels >= depth[-1]] = np.nan
    return out


def _interp_temp_psal_fr0(pressure: np.ndarray, temp: np.ndarray, psal: np.ndarray, levels: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    z_t, t = _prepare_profile_input(temp, pressure)
    z_s, s = _prepare_profile_input(psal, pressure)

    t_std = _interp_linear_surface_fill(z_t, t, levels)
    s_std = _interp_linear_surface_fill(z_s, s, levels)

    if _gsw is not None and z_t.size > 5 and z_s.size > 5:
        union = np.unique(np.concatenate([z_t, z_s]).astype(np.float64))
        t_union = np.interp(union, z_t, t, left=t[0], right=t[-1]).astype(np.float64)
        s_union = np.interp(union, z_s, s, left=s[0], right=s[-1]).astype(np.float64)
        valid = np.isfinite(union) & np.isfinite(t_union) & np.isfinite(s_union)
        if np.count_nonzero(valid) > 3 and np.all(np.diff(union[valid]) > 0):
            try:
                s_i, t_i = _gsw.sa_ct_interp(s_union[valid], t_union[valid], union[valid], levels)
                t_std = np.asarray(t_i, dtype=np.float64)
                s_std = np.asarray(s_i, dtype=np.float64)
                t_std[levels >= z_t[-1]] = np.nan
                s_std[levels >= z_s[-1]] = np.nan
                t_std[levels < z_t[0]] = t[0]
                s_std[levels < z_s[0]] = s[0]
            except Exception:
                pass

    return t_std, s_std


def _linear_interp_extrap(x: np.ndarray, xp: np.ndarray, fp: np.ndarray) -> np.ndarray:
    valid = np.isfinite(xp) & np.isfinite(fp)
    if np.count_nonzero(valid) == 0:
        return np.full(x.shape, np.nan, dtype=np.float64)
    xp = xp[valid]
    fp = fp[valid]
    order = np.argsort(xp)
    xp = xp[order]
    fp = fp[order]
    keep = np.r_[np.diff(xp) > 0, True]
    xp = xp[keep]
    fp = fp[keep]
    if xp.size == 1:
        return np.full(x.shape, fp[0], dtype=np.float64)
    out = np.interp(x, xp, fp).astype(np.float64)
    left = x < xp[0]
    if np.any(left):
        slope = (fp[1] - fp[0]) / max(xp[1] - xp[0], 1e-12)
        out[left] = fp[0] + slope * (x[left] - xp[0])
    right = x > xp[-1]
    if np.any(right):
        slope = (fp[-1] - fp[-2]) / max(xp[-1] - xp[-2], 1e-12)
        out[right] = fp[-1] + slope * (x[right] - xp[-1])
    return out


def _interpolate_locations(lr0_path: Path, juld_hr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    dataset = xr.load_dataset(lr0_path, decode_times=False)
    try:
        juld = np.asarray(dataset["JULD"].values, dtype=np.float64)
        lat = np.asarray(dataset["LATITUDE"].values, dtype=np.float64)
        lon = np.asarray(dataset["LONGITUDE"].values, dtype=np.float64)
    finally:
        dataset.close()

    valid = np.isfinite(juld) & np.isfinite(lat) & np.isfinite(lon)
    if np.count_nonzero(valid) < 2:
        return np.full(juld_hr.shape, np.nan, dtype=np.float64), np.full(juld_hr.shape, np.nan, dtype=np.float64)

    return (
        _linear_interp_extrap(juld_hr, juld[valid], lat[valid]),
        _linear_interp_extrap(juld_hr, juld[valid], lon[valid]),
    )


def _format_profile_timestamp(timestamp: pd.Timestamp) -> str:
    return timestamp.tz_convert("UTC").strftime(DATE_FORMAT)


def _build_full_resolution_profiles(
    smru_name: str,
    raw: HrRawData,
    bounds: DetectedProfiles,
    *,
    lat_hr: np.ndarray,
    lon_hr: np.ndarray,
) -> list[OdvProfile]:
    levels = STANDARD_HR_LEVELS.astype(np.float64)
    profiles: list[OdvProfile] = []

    for station_index, (start, end) in enumerate(zip(bounds.ascent_start_index, bounds.profile_end_index), start=1):
        segment_pressure = raw.pressure[start : end + 1]
        if segment_pressure.size == 0:
            continue

        max_index = int(np.nanargmax(segment_pressure))
        cut = min(max_index + 6, segment_pressure.size)
        if cut >= segment_pressure.size:
            continue

        sample_slice = slice(start + cut, end + 1)
        pressure = raw.pressure[sample_slice]
        temperature = raw.temperature[sample_slice]
        salinity = raw.salinity[sample_slice]
        fluorescence = raw.fluorescence[sample_slice]
        oxygen = raw.oxygen[sample_slice]
        light = raw.light[sample_slice]

        lat_value = float(lat_hr[station_index - 1])
        lon_value = float(lon_hr[station_index - 1])
        if not np.isfinite(lat_value) or not np.isfinite(lon_value) or lat_value == 0.0:
            continue

        temp_interp, sal_interp = _interp_temp_psal_fr0(pressure, temperature, salinity, levels)
        fluo_interp = _interp_sensor(pressure, fluorescence, levels)
        oxy_interp = _interp_sensor(pressure, oxygen, levels)
        light_interp = _interp_sensor(pressure, light, levels)

        if not (np.isfinite(temp_interp).any() or np.isfinite(sal_interp).any()):
            continue

        profiles.append(
            OdvProfile(
                smru_name=smru_name,
                station=station_index,
                timestamp=_format_profile_timestamp(raw.timestamp[end]),
                longitude=lon_value,
                latitude=lat_value,
                pressure=tuple(levels.astype(float)),
                temperature=tuple(temp_interp.astype(float)),
                salinity=tuple(sal_interp.astype(float)),
                fluorescence=tuple(fluo_interp.astype(float)),
                oxygen=tuple(oxy_interp.astype(float)),
                light=tuple(light_interp.astype(float)),
                conductivity=tuple(np.full(levels.shape, np.nan, dtype=np.float64)),
            )
        )

    return profiles


def _write_traj_file(
    config: MeopConfig,
    info,
    smru_name: str,
    raw: HrRawData,
    *,
    juld_hr: np.ndarray,
    lat_hr: np.ndarray,
    lon_hr: np.ndarray,
    lr0_path: Path,
    now: datetime,
) -> Path:
    lr0 = xr.load_dataset(lr0_path, decode_times=False)
    try:
        attrs = {key: value for key, value in lr0.attrs.items() if "number" not in key.lower()}
    finally:
        lr0.close()

    lat_interp = _linear_interp_extrap(raw.juld, juld_hr, lat_hr)
    lon_interp = _linear_interp_extrap(raw.juld, juld_hr, lon_hr)

    data_vars: dict[str, tuple[tuple[str, ...], np.ndarray, dict[str, object]]] = {
        "TIME": (("TIME",), raw.juld.astype(np.float32), {"units": "days since 1950-01-01 00:00:00 UTC"}),
        "LATITUDE": (("TIME",), lat_interp.astype(np.float32), {"units": "degree_north"}),
        "LONGITUDE": (("TIME",), lon_interp.astype(np.float32), {"units": "degree_east"}),
        "PRES": (("TIME",), raw.pressure.astype(np.float32), {"units": "SEA PRESSURE"}),
        "TEMP": (("TIME",), raw.temperature.astype(np.float32), {"units": "degree_Celsius"}),
        "PSAL": (("TIME",), raw.salinity.astype(np.float32), {"units": "PRACTICAL SALINITY"}),
    }
    if raw.isfluo:
        data_vars["CHLA"] = (("TIME",), raw.fluorescence.astype(np.float32), {"units": "mg/m3"})
    if raw.isoxy:
        data_vars["DOXY"] = (("TIME",), raw.oxygen.astype(np.float32), {"units": "micromole/kg"})

    dataset = xr.Dataset(attrs=attrs)
    dataset = dataset.assign_coords(TIME=("TIME", np.arange(raw.n_samples, dtype=np.int32)))
    for name, (dims, values, var_attrs) in data_vars.items():
        dataset[name] = (dims, values)
        dataset[name].attrs.update(var_attrs)
        dataset[name].encoding["_FillValue"] = np.float32(99999)

    dataset.attrs["DATE_CREATION"] = now.strftime("%Y%m%d%H%M%S")
    dataset.attrs["date_update"] = now.strftime("%Y-%m-%dT%H:%M:%SZ")
    dataset.attrs["profile_source"] = "full-resolution-traj"

    target = fname_traj(smru_name, deployment=info.EXP, config=config)
    target.parent.mkdir(parents=True, exist_ok=True)
    save_dataset_with_compression(dataset, target)
    return target


def create_fr0_python(
    config: MeopConfig,
    selection: Selection,
    *,
    now: datetime | None = None,
) -> NcargoResult:
    info = load_info_deployment(config, deployment=selection.deployment, smru_name=selection.smru_name)
    timestamp = now or datetime.now(timezone.utc)

    tags = (selection.smru_name,) if selection.smru_name else tuple(info.hr_platform_codes)
    written: list[Path] = []
    processed: list[str] = []

    for smru_name in tags:
        resolved = resolve_hr_ctd_path(config, smru_name)
        if resolved is None or not resolved.exists:
            continue

        lr0_path = fname_prof(smru_name, deployment=info.EXP, qf="lr0", config=config)
        if not lr0_path.exists():
            continue

        raw = load_hr_data(resolved.expected_path, continuous=resolved.continuous)
        bounds = detect_profiles(raw)
        if bounds.n_profiles == 0:
            continue

        juld_hr = raw.juld[bounds.profile_end_index]
        lat_hr, lon_hr = _interpolate_locations(lr0_path, juld_hr)
        profiles = _build_full_resolution_profiles(smru_name, raw, bounds, lat_hr=lat_hr, lon_hr=lon_hr)
        if not profiles:
            continue

        dataset = _build_ncargo_dataset(config, info, smru_name, profiles, now=timestamp)
        dataset = apply_lr0_qc_filters(config, info, smru_name, dataset).dataset
        dataset.attrs["profile_source"] = "full-resolution"
        target = fname_prof(smru_name, deployment=info.EXP, qf="fr0", config=config)
        target.parent.mkdir(parents=True, exist_ok=True)
        save_dataset_with_compression(dataset, target)
        written.append(target)

        traj_path = _write_traj_file(
            config,
            info,
            smru_name,
            raw,
            juld_hr=juld_hr,
            lat_hr=lat_hr,
            lon_hr=lon_hr,
            lr0_path=lr0_path,
            now=timestamp,
        )
        written.append(traj_path)
        processed.append(smru_name)

    return NcargoResult(written_files=tuple(written), processed_tags=tuple(processed))


def create_fr0(config: MeopConfig, selection: Selection, *, now: datetime | None = None) -> NcargoResult:
    return create_fr0_python(config, selection, now=now)


def create_fr0_without_lr0(*args, **kwargs):
    raise NotImplementedError("create_fr0_without_lr0 still depends on external location sources and is deferred.")
