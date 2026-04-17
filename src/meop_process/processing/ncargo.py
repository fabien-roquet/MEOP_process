from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import xarray as xr

from ..catalog.deployments import load_info_deployment, load_platform_catalog
from ..catalog.filenames import fname_prof
from ..catalog.tables import read_csv_rows
from ..data.layout import resolve_table_path
from ..io.raw_odv import OdvProfile, discover_raw_odv_files, load_raw_odv_profiles
from ..models import DeploymentInfo, MeopConfig, Selection
from .netcdf import DEFAULT_FORMAT, save_dataset_with_compression
from .qc import apply_lr0_qc_filters


def prepare_ncargo_inputs(config: MeopConfig, info: DeploymentInfo) -> bool:
    """Compatibility stub kept for callers from the migration phase.

    The pure-Python package reads low-resolution ODV inputs directly and no longer requires any
    additional preprocessing step before creating LR0.
    """

    _ = config
    _ = info
    return True


@dataclass(frozen=True)
class NcargoResult:
    written_files: tuple[Path, ...]
    processed_tags: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "written_files": [str(path) for path in self.written_files],
            "processed_tags": list(self.processed_tags),
        }


def _load_table_coeff_rows(config: MeopConfig) -> list[dict[str, str]]:
    path = resolve_table_path(config, "table_coeff.csv", required=False)
    return read_csv_rows(path)


def _write_table_coeff_rows(config: MeopConfig, rows: list[dict[str, str]]) -> Path:
    path = resolve_table_path(config, "table_coeff.csv", required=False)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    import csv

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def ensure_default_coefficients(config: MeopConfig, smru_names: Iterable[str]) -> Path | None:
    rows = _load_table_coeff_rows(config)
    if not rows:
        return None
    existing = {row.get("smru_platform_code", "") for row in rows}
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
    return _write_table_coeff_rows(config, rows)


def _parse_timestamp(text: str) -> datetime:
    value = str(text).strip()
    for fmt in (
        "%Y-%m-%d %H:%M",
        "%Y-%m-%d %H:%M:%S",
        "%Y/%m/%d %H:%M:%S",
        "%Y/%m/%d %H:%M:%S.%f",
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%dT%H:%M:%SZ",
    ):
        try:
            return datetime.strptime(value, fmt).replace(tzinfo=timezone.utc)
        except ValueError:
            continue
    raise ValueError(f"Unsupported timestamp format: {text!r}")


def _juld_from_timestamp(text: str) -> float:
    dt = _parse_timestamp(text)
    ref = datetime(1950, 1, 1, tzinfo=timezone.utc)
    return (dt - ref).total_seconds() / 86400.0


def _format_argo_datetime(dt: datetime) -> str:
    return dt.astimezone(timezone.utc).strftime("%Y%m%d%H%M%S")


def _format_global_datetime(dt: datetime) -> str:
    return dt.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _char_array(shape: tuple[int, ...], value: str) -> np.ndarray:
    width = shape[-1]
    fill = value.encode("ascii", "replace")[:width].ljust(width, b" ")
    return np.full(shape, fill, dtype="S1")


def _char_array_from_strings(values: list[str], width: int) -> np.ndarray:
    out = np.full((len(values), width), b" ", dtype="S1")
    for index, value in enumerate(values):
        encoded = value.encode("ascii", "replace")[:width].ljust(width, b" ")
        out[index, :] = np.frombuffer(encoded, dtype="S1")
    return out


def _char_array_from_string_matrix(values: list[list[str]], width: int) -> np.ndarray:
    n_rows = len(values)
    n_cols = len(values[0]) if values else 0
    out = np.full((n_rows, n_cols, width), b" ", dtype="S1")
    for row, row_values in enumerate(values):
        for col, value in enumerate(row_values):
            encoded = value.encode("ascii", "replace")[:width].ljust(width, b" ")
            out[row, col, :] = np.frombuffer(encoded, dtype="S1")
    return out


def _numeric_matrix(profiles: list[OdvProfile], field: str, *, fill: float = np.nan) -> np.ndarray:
    n_prof = len(profiles)
    n_levels = max((profile.n_levels for profile in profiles), default=0)
    data = np.full((n_prof, n_levels), fill, dtype=np.float32)
    for row, profile in enumerate(profiles):
        values = np.asarray(getattr(profile, field), dtype=np.float32)
        data[row, : len(values)] = values
    return data


def _qc_matrix(data: np.ndarray) -> np.ndarray:
    valid = np.isfinite(data)
    return np.where(valid, "1", "9").astype('U1')


def _profile_qc(level_qc: np.ndarray) -> np.ndarray:
    return np.where(np.any(level_qc == "1", axis=1), "A", "F").astype('U1')


def _valid_profile_count(*matrices: np.ndarray) -> int:
    if not matrices:
        return 0
    valid = np.ones(matrices[0].shape[0], dtype=bool)
    for matrix in matrices:
        valid &= np.any(matrix == "1", axis=1)
    return int(valid.sum())


def _valid_geospatial_mask(temp_qc: np.ndarray, psal_qc: np.ndarray) -> np.ndarray:
    return np.any((temp_qc == "1") | (psal_qc == "1"), axis=1)


def _format_wmo_platform_code(raw: str | None) -> str:
    if raw is None:
        return " "
    value = str(raw).strip()
    if not value:
        return " "
    if value.startswith("Q99"):
        return value
    digits = "".join(character for character in value if character.isdigit())
    if not digits:
        return value
    return f"Q99{int(digits):05d}"


def _loc_algorithm_label(raw: str | None) -> str:
    value = (raw or "").strip().upper()
    if value == "L":
        return "CLS LEAST SQUARES"
    if value == "K":
        return "CLS KALMAN"
    if value == "S":
        return "SMRU"
    return value or "ARGOS"


def _positioning_system_value(raw: str | None) -> str:
    value = (raw or "").strip().upper()
    if value == "L":
        return "LS"
    if value == "K":
        return "K"
    return "ARGOS"


def _platform_metadata_for_tag(config: MeopConfig, smru_name: str) -> dict[str, object]:
    import json

    records: list[dict[str, object]] = []
    for root in config.config_files_search_dirs:
        for name in ("platform3.json", "platform2_patch.json", "platform2.json"):
            path = root / name
            if not path.exists():
                continue
            payload = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(payload, dict):
                payload = [payload]
            if isinstance(payload, list):
                records.extend(record for record in payload if isinstance(record, dict))
    prefix = smru_name.split("-N")[0]
    for record in records:
        if str(record.get("smru_platform_code", "")).strip() == prefix:
            return record
    return {}


def _build_station_parameter_array(n_prof: int, n_param: int, parameter_names: list[str]) -> np.ndarray:
    rows = [parameter_names for _ in range(n_prof)]
    return _char_array_from_string_matrix(rows, 16)


def _calibration_equations(parameter_names: list[str]) -> list[str]:
    equations = {
        "PRES": "Pc = P - ( p1 [dbar/km] * P  * 1e-3 + p2 [dbar] )",
        "TEMP": "Tc = T - ( t1 [degC/km] * Pc * 1e-3 + t2 [degC] )",
        "PSAL": "Sc = S - ( s1 [ psu/km] * Pc * 1e-3 + s2 [ psu] )",
        "CHLA": "Fc = f1 * F + f2",
        "DOXY": " ",
        "LIGHT": " ",
    }
    return [equations.get(name, " ") for name in parameter_names]


def _build_ncargo_dataset(
    config: MeopConfig,
    info: DeploymentInfo,
    smru_name: str,
    profiles: list[OdvProfile],
    *,
    now: datetime,
) -> xr.Dataset:
    profiles = sorted(profiles, key=lambda profile: (_parse_timestamp(profile.timestamp), profile.station))
    n_prof = len(profiles)
    n_levels = max((profile.n_levels for profile in profiles), default=0)

    pressure = _numeric_matrix(profiles, "pressure")
    temperature = _numeric_matrix(profiles, "temperature")
    salinity = _numeric_matrix(profiles, "salinity")
    fluorescence = _numeric_matrix(profiles, "fluorescence")
    light = _numeric_matrix(profiles, "light")
    oxygen = _numeric_matrix(profiles, "oxygen")

    has_fluorescence = bool(np.isfinite(fluorescence).any())
    has_light = bool(np.isfinite(light).any())
    has_oxygen = bool(np.isfinite(oxygen).any())

    parameter_names = ["PRES", "TEMP", "PSAL"]
    if has_fluorescence:
        parameter_names.append("CHLA")
    if has_oxygen:
        parameter_names.append("DOXY")
    if has_light:
        parameter_names.append("LIGHT")
    n_param = len(parameter_names)

    pres_qc = _qc_matrix(pressure)
    temp_qc = _qc_matrix(temperature)
    psal_qc = _qc_matrix(salinity)
    fl_qc = _qc_matrix(fluorescence)
    oxy_qc = _qc_matrix(oxygen)
    light_qc = _qc_matrix(light)

    platform = _platform_metadata_for_tag(config, smru_name)
    platform_code = str(platform.get("platform_code", "")).strip()
    if platform_code:
        platform_number_value = f"{int(platform_code):08d}"
    else:
        platform_number_value = "00000000"

    station_parameters = _build_station_parameter_array(n_prof, n_param, parameter_names)
    parameter_array = station_parameters[:, np.newaxis, :, :]
    calib_equations = np.full((n_prof, 1, n_param, 256), b" ", dtype="S1")
    calib_coefficients = np.full((n_prof, 1, n_param, 256), b" ", dtype="S1")
    for index, equation in enumerate(_calibration_equations(parameter_names)):
        encoded = equation.encode("ascii", "replace")[:256].ljust(256, b" ")
        calib_equations[:, 0, index, :] = np.frombuffer(encoded, dtype="S1")

    juld = np.array([_juld_from_timestamp(profile.timestamp) for profile in profiles], dtype=np.float64)
    latitude = np.array([profile.latitude for profile in profiles], dtype=np.float64)
    longitude = np.array([profile.longitude for profile in profiles], dtype=np.float64)
    timestamps = [_parse_timestamp(profile.timestamp) for profile in profiles]

    data_vars: dict[str, tuple[tuple[str, ...], np.ndarray]] = {
        "DATA_TYPE": (("STRING16",), _char_array((16,), "Argo profile")),
        "FORMAT_VERSION": (("STRING4",), _char_array((4,), "3.0 ")),
        "HANDBOOK_VERSION": (("STRING4",), _char_array((4,), "3.0 ")),
        "REFERENCE_DATE_TIME": (("DATE_TIME",), _char_array((14,), "19500101000000")),
        "DATE_CREATION": (("DATE_TIME",), _char_array((14,), _format_argo_datetime(now))),
        "DATE_UPDATE": (("DATE_TIME",), _char_array((14,), _format_argo_datetime(now))),
        "PLATFORM_NUMBER": (("N_PROF", "STRING8"), _char_array((n_prof, 8), platform_number_value)),
        "PROJECT_NAME": (("N_PROF", "STRING64"), _char_array((n_prof, 64), "MEOP")),
        "PI_NAME": (("N_PROF", "STRING64"), _char_array((n_prof, 64), info.PI or str(platform.get("pi_code", "")).strip())),
        "STATION_PARAMETERS": (("N_PROF", "N_PARAM", "STRING16"), station_parameters),
        "CYCLE_NUMBER": (("N_PROF",), np.arange(1, n_prof + 1, dtype=np.int32)),
        "DIRECTION": (("N_PROF",), np.full((n_prof,), b"A", dtype="S1")),
        "DATA_CENTRE": (("N_PROF", "STRING2"), _char_array((n_prof, 2), "IF")),
        "DC_REFERENCE": (("N_PROF", "STRING32"), _char_array((n_prof, 32), f"{info.EXP[:24]:>24}{index + 1:08d}")),
        "DATA_STATE_INDICATOR": (("N_PROF", "STRING4"), _char_array((n_prof, 4), " ")),
        "DATA_MODE": (("N_PROF",), np.full((n_prof,), b"D", dtype="S1")),
        "INST_REFERENCE": (("N_PROF", "STRING64"), _char_array((n_prof, 64), " ")),
        "WMO_INST_TYPE": (("N_PROF", "STRING4"), _char_array((n_prof, 4), "995 ")),
        "JULD": (("N_PROF",), juld),
        "JULD_QC": (("N_PROF",), np.full((n_prof,), b"1", dtype="S1")),
        "JULD_LOCATION": (("N_PROF",), juld),
        "LATITUDE": (("N_PROF",), latitude),
        "LONGITUDE": (("N_PROF",), longitude),
        "POSITION_QC": (("N_PROF",), np.full((n_prof,), b"1", dtype="S1")),
        "POSITIONING_SYSTEM": (("N_PROF", "STRING8"), _char_array((n_prof, 8), _positioning_system_value(platform.get("loc_algorithm")))),
        "PROFILE_PRES_QC": (("N_PROF",), _profile_qc(pres_qc)),
        "PROFILE_PSAL_QC": (("N_PROF",), _profile_qc(psal_qc)),
        "PROFILE_TEMP_QC": (("N_PROF",), _profile_qc(temp_qc)),
        "PRES": (("N_PROF", "N_LEVELS"), pressure),
        "PRES_QC": (("N_PROF", "N_LEVELS"), pres_qc),
        "PRES_ADJUSTED": (("N_PROF", "N_LEVELS"), pressure.copy()),
        "PRES_ADJUSTED_QC": (("N_PROF", "N_LEVELS"), pres_qc.copy()),
        "PRES_ADJUSTED_ERROR": (("N_PROF", "N_LEVELS"), np.full((n_prof, n_levels), np.nan, dtype=np.float32)),
        "TEMP": (("N_PROF", "N_LEVELS"), temperature),
        "TEMP_QC": (("N_PROF", "N_LEVELS"), temp_qc),
        "TEMP_ADJUSTED": (("N_PROF", "N_LEVELS"), temperature.copy()),
        "TEMP_ADJUSTED_QC": (("N_PROF", "N_LEVELS"), temp_qc.copy()),
        "TEMP_ADJUSTED_ERROR": (("N_PROF", "N_LEVELS"), np.full((n_prof, n_levels), np.nan, dtype=np.float32)),
        "PSAL": (("N_PROF", "N_LEVELS"), salinity),
        "PSAL_QC": (("N_PROF", "N_LEVELS"), psal_qc),
        "PSAL_ADJUSTED": (("N_PROF", "N_LEVELS"), salinity.copy()),
        "PSAL_ADJUSTED_QC": (("N_PROF", "N_LEVELS"), psal_qc.copy()),
        "PSAL_ADJUSTED_ERROR": (("N_PROF", "N_LEVELS"), np.full((n_prof, n_levels), np.nan, dtype=np.float32)),
        "PARAMETER": (("N_PROF", "N_CALIB", "N_PARAM", "STRING16"), parameter_array),
        "SCIENTIFIC_CALIB_EQUATION": (("N_PROF", "N_CALIB", "N_PARAM", "STRING256"), calib_equations),
        "SCIENTIFIC_CALIB_COEFFICIENT": (("N_PROF", "N_CALIB", "N_PARAM", "STRING256"), calib_coefficients),
    }

    if has_fluorescence:
        data_vars["PROFILE_CHLA_QC"] = (("N_PROF",), _profile_qc(fl_qc))
        data_vars["CHLA"] = (("N_PROF", "N_LEVELS"), fluorescence)
        data_vars["CHLA_QC"] = (("N_PROF", "N_LEVELS"), fl_qc)
        data_vars["CHLA_ADJUSTED"] = (("N_PROF", "N_LEVELS"), fluorescence.copy())
        data_vars["CHLA_ADJUSTED_QC"] = (("N_PROF", "N_LEVELS"), fl_qc.copy())
        data_vars["CHLA_ADJUSTED_ERROR"] = (("N_PROF", "N_LEVELS"), np.full((n_prof, n_levels), np.nan, dtype=np.float32))
    if has_oxygen:
        data_vars["PROFILE_DOXY_QC"] = (("N_PROF",), _profile_qc(oxy_qc))
        data_vars["DOXY"] = (("N_PROF", "N_LEVELS"), oxygen)
        data_vars["DOXY_QC"] = (("N_PROF", "N_LEVELS"), oxy_qc)
        data_vars["DOXY_ADJUSTED"] = (("N_PROF", "N_LEVELS"), oxygen.copy())
        data_vars["DOXY_ADJUSTED_QC"] = (("N_PROF", "N_LEVELS"), oxy_qc.copy())
        data_vars["DOXY_ADJUSTED_ERROR"] = (("N_PROF", "N_LEVELS"), np.full((n_prof, n_levels), np.nan, dtype=np.float32))
    if has_light:
        data_vars["PROFILE_LIGHT_QC"] = (("N_PROF",), _profile_qc(light_qc))
        data_vars["LIGHT"] = (("N_PROF", "N_LEVELS"), light)
        data_vars["LIGHT_QC"] = (("N_PROF", "N_LEVELS"), light_qc)
        data_vars["LIGHT_ADJUSTED"] = (("N_PROF", "N_LEVELS"), light.copy())
        data_vars["LIGHT_ADJUSTED_QC"] = (("N_PROF", "N_LEVELS"), light_qc.copy())
        data_vars["LIGHT_ADJUSTED_ERROR"] = (("N_PROF", "N_LEVELS"), np.full((n_prof, n_levels), np.nan, dtype=np.float32))

    dataset = xr.Dataset(
        data_vars=data_vars,
        coords={
            "N_PROF": np.arange(n_prof),
            "N_LEVELS": np.arange(n_levels),
            "N_PARAM": np.arange(n_param),
            "N_CALIB": np.arange(1),
        }
    )

    fill_float = 99999.0
    for name in [
        "PRES",
        "PRES_ADJUSTED",
        "PRES_ADJUSTED_ERROR",
        "TEMP",
        "TEMP_ADJUSTED",
        "TEMP_ADJUSTED_ERROR",
        "PSAL",
        "PSAL_ADJUSTED",
        "PSAL_ADJUSTED_ERROR",
        "CHLA",
        "CHLA_ADJUSTED",
        "CHLA_ADJUSTED_ERROR",
        "DOXY",
        "DOXY_ADJUSTED",
        "DOXY_ADJUSTED_ERROR",
        "LIGHT",
        "LIGHT_ADJUSTED",
        "LIGHT_ADJUSTED_ERROR",
    ]:
        if name in dataset:
            dataset[name].encoding["_FillValue"] = fill_float

    dataset["CYCLE_NUMBER"].encoding["_FillValue"] = np.int32(99999)
    for name in ("JULD", "JULD_LOCATION", "LATITUDE", "LONGITUDE"):
        dataset[name].encoding["_FillValue"] = float(fill_float)

    valid_mask = _valid_geospatial_mask(temp_qc, psal_qc)
    valid_lat = latitude[valid_mask]
    valid_lon = longitude[valid_mask]
    attrs: dict[str, object] = {
        "comment": " ",
        "pi_name": info.PI,
        "data_type": "Marine mammals time-series data",
        "format_version": "1.1",
        "date_update": _format_global_datetime(now),
        "version_database": config.version,
        "PI": info.PI,
        "reference_file_name": f"{smru_name}_prof.nc",
        "is_the_data_public": (
            "yes: data can be used freely providing that data owner is properly cited (see meop.net for citing information)"
            if str(info.public).strip() in {"1", "true", "True", "Y", "y"}
            else "no: data cannot be used without the prior consent of the data owner"
        ),
        "nation": info.NATION,
        "deployment_code": info.EXP,
        "source": "Marine mammal observation",
        "data_mode": "D",
        "references": "http://www.meop.net",
        "reference_doi": " ",
        "Conventions": "CF-1.6 Sea-mammals-1.1",
        "Netcdf_version": "4",
        "naming_authority": "MEOP consortium (Marine Mammals Exploring the Oceans Pole to Pole)",
        "cdm_data_type": "Station",
        "geospatial_vertical_min": 0.0,
        "geospatial_vertical_max": 2000.0,
        "data_assembly_center": "MEOP/Fabien Roquet (MISU)",
        "distribution_statement": "Follow MEOP data policy standards, cf. http://www.meop.net/the-dataset/data-access.html. Data available free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data",
        "citation": "The marine mammal data were collected and made freely available by the International MEOP Consortium and the national programs that contribute to it (http://www.meop.net).",
        "thermal_lag_adjustment": "no",
        "platform_code": platform_code or " ",
        "wmo_platform_code": _format_wmo_platform_code(platform.get("wmo_platform_code")),
        "smru_platform_code": smru_name,
        "species": str(platform.get("species", "")).strip(),
        "time_coverage_start": str(platform.get("time_coverage_start", timestamps[0].strftime("%Y-%m-%dT00:00:00Z"))).strip() if timestamps else "",
        "time_coverage_end": str(platform.get("time_coverage_end", timestamps[-1].strftime("%Y-%m-%dT00:00:00Z"))).strip() if timestamps else "",
        "location": str(platform.get("location", info.EXP)).replace(",", ".").strip() or info.EXP,
        "loc_algorithm": _loc_algorithm_label(platform.get("loc_algorithm")),
        "firmware_version": str(platform.get("firmware_version", "")).strip(),
        "firmware_parameters": str(platform.get("firmware_parameters", "")).strip(),
        "instr_id": str(platform.get("instr_id", "")).strip(),
        "ptt": str(platform.get("ptt", "")).strip(),
        "number_of_ts_profiles": float(_valid_profile_count(temp_qc, psal_qc)),
        "number_of_t_profiles": float(_valid_profile_count(temp_qc)),
        "number_chla_profiles": float(_valid_profile_count(fl_qc)) if has_fluorescence else 0.0,
        "number_doxy_profiles": float(_valid_profile_count(oxy_qc)) if has_oxygen else 0.0,
        "number_light_profiles": float(_valid_profile_count(light_qc)) if has_light else 0.0,
        "geospatial_lat_min": float(np.nanmin(valid_lat)) if valid_lat.size else np.nan,
        "geospatial_lat_max": float(np.nanmax(valid_lat)) if valid_lat.size else np.nan,
        "geospatial_lon_min": float(np.nanmin(np.where(valid_lon < 0, valid_lon + 360, valid_lon))) if valid_lon.size else np.nan,
        "geospatial_lon_max": float(np.nanmax(np.where(valid_lon < 0, valid_lon + 360, valid_lon))) if valid_lon.size else np.nan,
    }
    dataset.attrs.update(attrs)

    dataset["DATA_TYPE"].attrs.update({"comment": "Data type"})
    dataset["FORMAT_VERSION"].attrs.update({"comment": "File format version"})
    dataset["HANDBOOK_VERSION"].attrs.update({"comment": "Data handbook version"})
    dataset["REFERENCE_DATE_TIME"].attrs.update({"comment": "Date of reference for Julian days", "conventions": "YYYYMMDDHHMISS"})
    dataset["DATE_CREATION"].attrs.update({"comment": "Date of file creation", "conventions": "YYYYMMDDHHMISS"})
    dataset["DATE_UPDATE"].attrs.update({"long_name": "Date of update of this file", "conventions": "YYYYMMDDHHMISS"})
    dataset["PLATFORM_NUMBER"].attrs.update({"long_name": "Float unique identifier", "conventions": "WMO float identifier : A9IIIII"})
    dataset["PROJECT_NAME"].attrs.update({"comment": "Name of the project"})
    dataset["PI_NAME"].attrs.update({"comment": "Name of the principal investigator"})
    dataset["STATION_PARAMETERS"].attrs.update({"long_name": "List of available parameters for the station", "conventions": "Argo reference table 3"})
    dataset["CYCLE_NUMBER"].attrs.update({"long_name": "Float cycle number", "units": "1", "conventions": "0..N, 0 : launch cycle (if exists), 1 : first complete cycle"})
    dataset["DIRECTION"].attrs.update({"long_name": "Direction of the station profiles", "conventions": "A: ascending profiles, D: descending profiles"})
    dataset["DATA_CENTRE"].attrs.update({"long_name": "Data centre in charge of float data processing", "conventions": "Argo reference table 4"})
    dataset["DC_REFERENCE"].attrs.update({"long_name": "Station unique identifier in data centre", "conventions": "Data centre convention"})
    dataset["DATA_STATE_INDICATOR"].attrs.update({"long_name": "Degree of processing the data have passed through", "conventions": "Argo reference table 6"})
    dataset["DATA_MODE"].attrs.update({"long_name": "Delayed mode or real time data", "conventions": "R : real time; D : delayed mode; A : real time with adjustment"})
    dataset["INST_REFERENCE"].attrs.update({"long_name": "Instrument type", "conventions": "Brand, type, serial number"})
    dataset["WMO_INST_TYPE"].attrs.update({"long_name": "Coded instrument type", "conventions": "Argo reference table 8"})
    dataset["JULD"].attrs.update({"long_name": "Julian day (UTC) of the station relative to REFERENCE_DATE_TIME", "units": "days since 1950-01-01 00:00:00 UTC", "calendar": "julian", "conventions": "Relative julian days with decimal part (as parts of day)"})
    dataset["JULD_QC"].attrs.update({"long_name": "Quality on Date and Time", "conventions": "Argo reference table 2"})
    dataset["JULD_LOCATION"].attrs.update({"long_name": "Julian day (UTC) of the location relative to REFERENCE_DATE_TIME", "units": "days since 1950-01-01 00:00:00 UTC", "calendar": "julian", "conventions": "Relative julian days with decimal part (as parts of day)"})
    dataset["LATITUDE"].attrs.update({"long_name": "Latitude of the station, best estimate", "units": "degree_north", "valid_min": -90.0, "valid_max": 90.0})
    dataset["LONGITUDE"].attrs.update({"long_name": "Longitude of the station, best estimate", "units": "degree_east", "valid_min": -180.0, "valid_max": 180.0})
    dataset["POSITION_QC"].attrs.update({"long_name": "Quality on position (latitude and longitude)", "conventions": "Argo reference table 2"})
    dataset["POSITIONING_SYSTEM"].attrs.update({"long_name": "Positioning system"})
    dataset["PROFILE_PRES_QC"].attrs.update({"long_name": "Global quality flag of PRES profile", "conventions": "Argo reference table 2a"})
    dataset["PROFILE_PSAL_QC"].attrs.update({"long_name": "Global quality flag of PSAL profile", "conventions": "Argo reference table 2a"})
    dataset["PROFILE_TEMP_QC"].attrs.update({"long_name": "Global quality flag of TEMP profile", "conventions": "Argo reference table 2a"})
    dataset["PRES"].attrs.update({"long_name": "SEA PRESSURE", "units": "decibar", "valid_min": 0.0, "valid_max": 12000.0, "comment": "In situ measurement, sea surface = 0"})
    dataset["PRES_QC"].attrs.update({"long_name": "quality flag", "conventions": "Argo reference table 2"})
    dataset["PRES_ADJUSTED"].attrs.update({"long_name": "SEA PRESSURE", "units": "decibar", "valid_min": 0.0, "valid_max": 12000.0, "comment": "In situ measurement, sea surface = 0"})
    dataset["PRES_ADJUSTED_QC"].attrs.update({"long_name": "quality flag", "conventions": "Argo reference table 2"})
    dataset["PRES_ADJUSTED_ERROR"].attrs.update({"long_name": "SEA PRESSURE", "units": "decibar", "comment": "Contains the error on the adjusted values as determined by the delayed mode QC process."})
    dataset["TEMP"].attrs.update({"long_name": "SEA TEMPERATURE IN SITU ITS-90 SCALE", "units": "degree_Celsius", "valid_min": -2.0, "valid_max": 40.0, "comment": "In situ measurement"})
    dataset["TEMP_QC"].attrs.update({"long_name": "quality flag", "conventions": "Argo reference table 2"})
    dataset["TEMP_ADJUSTED"].attrs.update({"long_name": "SEA TEMPERATURE IN SITU ITS-90 SCALE", "units": "degree_Celsius", "valid_min": -2.0, "valid_max": 40.0, "comment": "In situ measurement"})
    dataset["TEMP_ADJUSTED_QC"].attrs.update({"long_name": "quality flag", "conventions": "Argo reference table 2"})
    dataset["TEMP_ADJUSTED_ERROR"].attrs.update({"long_name": "SEA TEMPERATURE ERROR IN SITU ITS-90 SCALE", "units": "degree_Celsius", "comment": "Contains the error on the adjusted values as determined by the delayed mode QC process."})
    dataset["PSAL"].attrs.update({"long_name": "PRACTICAL SALINITY", "units": "1e-3", "valid_min": 0.0, "valid_max": 42.0, "comment": "In situ measurement"})
    dataset["PSAL_QC"].attrs.update({"long_name": "quality flag", "conventions": "Argo reference table 2"})
    dataset["PSAL_ADJUSTED"].attrs.update({"long_name": "ADJUSTED PRACTICAL SALINITY", "units": "1e-3", "valid_min": 0.0, "valid_max": 42.0, "comment": "In situ measurement"})
    dataset["PSAL_ADJUSTED_QC"].attrs.update({"long_name": "quality flag", "conventions": "Argo reference table 2"})
    dataset["PSAL_ADJUSTED_ERROR"].attrs.update({"long_name": "PRACTICAL SALINITY ERROR", "units": "1e-3", "comment": "Contains the error on the adjusted values as determined by the delayed mode QC process."})
    dataset["PARAMETER"].attrs.update({"long_name": "List of parameters with calibration information", "conventions": "Argo reference table 3"})
    dataset["SCIENTIFIC_CALIB_EQUATION"].attrs.update({"long_name": "Calibration equation for this parameter"})
    dataset["SCIENTIFIC_CALIB_COEFFICIENT"].attrs.update({"long_name": "SRDL identifier in the SMRU database"})

    if has_fluorescence:
        dataset["PROFILE_CHLA_QC"].attrs.update({"long_name": "Global quality flag of CHLA profile", "conventions": "Argo reference table 2a"})
        dataset["CHLA"].attrs.update({"long_name": "CHLOROPHYLL-A", "units": "mg/m3", "valid_min": 0.0, "valid_max": 10.0, "comment": "In situ measurement"})
        dataset["CHLA_QC"].attrs.update({"long_name": "quality flag", "conventions": "Argo reference table 2"})
        dataset["CHLA_ADJUSTED"].attrs.update({"long_name": "CHLOROPHYLL-A", "units": "mg/m3", "valid_min": 0.0, "valid_max": 10.0, "comment": "In situ measurement"})
        dataset["CHLA_ADJUSTED_QC"].attrs.update({"long_name": "quality flag", "conventions": "Argo reference table 2"})
        dataset["CHLA_ADJUSTED_ERROR"].attrs.update({"long_name": "CHLOROPHYLL-A", "units": "mg/m3", "comment": "Contains the error on the adjusted values as determined by the delayed mode QC process."})
    if has_oxygen:
        dataset["PROFILE_DOXY_QC"].attrs.update({"long_name": "Global quality flag of DOXY profile", "conventions": "Argo reference table 2a"})
        dataset["DOXY"].attrs.update({"long_name": "DISSOLVED OXYGEN", "units": "micromole/kg", "valid_min": 0.0, "valid_max": 600.0, "comment": "In situ measurement"})
        dataset["DOXY_QC"].attrs.update({"long_name": "quality flag", "conventions": "Argo reference table 2"})
        dataset["DOXY_ADJUSTED"].attrs.update({"long_name": "DISSOLVED OXYGEN", "units": "micromole/kg", "valid_min": 0.0, "valid_max": 600.0, "comment": "In situ measurement"})
        dataset["DOXY_ADJUSTED_QC"].attrs.update({"long_name": "quality flag", "conventions": "Argo reference table 2"})
        dataset["DOXY_ADJUSTED_ERROR"].attrs.update({"long_name": "DISSOLVED OXYGEN", "units": "micromole/kg", "comment": "Contains the error on the adjusted values as determined by the delayed mode QC process."})
    if has_light:
        dataset["PROFILE_LIGHT_QC"].attrs.update({"long_name": "Global quality flag of LIGHT profile", "conventions": "Argo reference table 2a"})
        dataset["LIGHT"].attrs.update({"long_name": "ln(PPFD)", "units": "ln(μmol/m²/s)", "valid_min": 0.0, "valid_max": 600.0, "comment": "In situ measurement"})
        dataset["LIGHT_QC"].attrs.update({"long_name": "quality flag", "conventions": "Argo reference table 2"})
        dataset["LIGHT_ADJUSTED"].attrs.update({"long_name": "ln(PPFD)", "units": "ln(μmol/m²/s)", "valid_min": 0.0, "valid_max": 600.0, "comment": "In situ measurement"})
        dataset["LIGHT_ADJUSTED_QC"].attrs.update({"long_name": "quality flag", "conventions": "Argo reference table 2"})
        dataset["LIGHT_ADJUSTED_ERROR"].attrs.update({"long_name": "ln(PPFD)", "units": "ln(μmol/m²/s)", "comment": "Contains the error on the adjusted values as determined by the delayed mode QC process."})

    return dataset


def create_ncargo_python(
    config: MeopConfig,
    selection: Selection,
    *,
    now: datetime | None = None,
    format: str = DEFAULT_FORMAT,
) -> NcargoResult:
    """Create core ``*_lr0_prof.nc`` files directly from ODV raw text.

    This implementation covers raw ODV discovery, parsing, tag grouping, lr0 netCDF creation,
    and default coefficient-table seeding.
    """

    info = load_info_deployment(config, deployment=selection.deployment, smru_name=selection.smru_name)
    raw_files = discover_raw_odv_files(config, info.EXP)
    raw_profiles = load_raw_odv_profiles(raw_files)
    if not raw_profiles:
        return NcargoResult(written_files=(), processed_tags=())

    grouped: dict[str, list[OdvProfile]] = defaultdict(list)
    for profile in raw_profiles:
        grouped[profile.smru_name].append(profile)

    selected_tags = tuple(sorted(grouped))
    if selection.smru_name:
        selected_tags = tuple(tag for tag in selected_tags if tag == selection.smru_name)

    output_dir = config.final_dataset_dir / info.EXP
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    timestamp = now or datetime.now(timezone.utc)

    ensure_default_coefficients(config, selected_tags)

    for smru_name in selected_tags:
        dataset = _build_ncargo_dataset(config, info, smru_name, grouped[smru_name], now=timestamp)
        dataset = apply_lr0_qc_filters(config, info, smru_name, dataset).dataset
        target = fname_prof(smru_name, deployment=info.EXP, qf="lr0", config=config)
        save_dataset_with_compression(dataset, target, format=format)
        written.append(target)

    return NcargoResult(written_files=tuple(written), processed_tags=selected_tags)
