from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import xarray as xr

from ..catalog.deployments import load_deployment_catalog, load_hr_catalog
from ..catalog.filenames import deployment_from_smru_name
from ..io.netcdf import decode_text, juld_to_datetime
from ..models import MeopConfig
from ..plotting.regions import label_region


PREFERRED_QF_ORDER: tuple[str, ...] = ("hr2", "fr1", "hr1", "lr1", "fr0", "hr0", "lr0")
TAG_FIELD_ORDER: tuple[str, ...] = (
    "SMRU_PLATFORM_CODE",
    "DEPLOYMENT_CODE",
    "JULD",
    "LATITUDE",
    "LONGITUDE",
    "MASK",
    "JULD_END",
    "N_PROF_TEMP",
    "N_PROF_PSAL",
    "N_PROF_CHLA",
    "PUBLIC",
    "T1",
    "T2",
    "S1",
    "S2",
    "remove",
    "Sremove",
    "comments",
    "variable_offset",
    "instr_id",
    "year",
    "period",
    "continuous",
    "prefix",
)
DEPLOYMENT_FIELD_ORDER: tuple[str, ...] = (
    "DEPLOYMENT_CODE",
    "JULD",
    "LATITUDE",
    "LONGITUDE",
    "N_PROF_TEMP",
    "N_PROF_PSAL",
    "N_PROF_CHLA",
    "N_TAGS",
    "PI_CODE",
    "PROCESS",
    "PUBLIC",
    "COUNTRY",
    "FIRST_VERSION",
    "LAST_VERSION",
)

@dataclass(frozen=True)
class SummaryUpdateResult:
    output_dir: Path
    list_tags_path: Path
    list_deployments_path: Path
    impacted_deployments: tuple[str, ...]
    written: bool

    def as_dict(self) -> dict[str, object]:
        return {
            "output_dir": str(self.output_dir),
            "list_tags_path": str(self.list_tags_path),
            "list_deployments_path": str(self.list_deployments_path),
            "impacted_deployments": list(self.impacted_deployments),
            "written": self.written,
        }


def _open_dataset_lazy(path: Path) -> xr.Dataset:
    errors: list[Exception] = []
    for engine in ("h5netcdf", "scipy", None):
        try:
            if engine is None:
                return xr.open_dataset(path, decode_times=False)
            return xr.open_dataset(path, engine=engine, decode_times=False)
        except Exception as exc:  # pragma: no cover - backend availability varies by environment
            errors.append(exc)
    raise OSError(f"Unable to open MEOP netCDF file: {path}") from errors[-1]


def _safe_numeric(array: object) -> np.ndarray:
    return np.asarray(array, dtype=np.float64)


def _valid_numeric(values: np.ndarray, *, upper_fill: float = 9e4) -> np.ndarray:
    return np.isfinite(values) & (np.abs(values) < upper_fill)


def _decode_juld(values: np.ndarray) -> np.ndarray:
    numeric = _safe_numeric(values)
    valid = _valid_numeric(numeric)
    out = np.full(numeric.shape, None, dtype=object)
    if np.any(valid):
        out[valid] = juld_to_datetime(numeric[valid])
    return out


def _format_datetime(value: datetime | None) -> str:
    if value is None:
        return ""
    return value.strftime("%Y-%m-%d %H:%M:%S")


def _decode_scalar_attr(dataset: xr.Dataset, key: str, default: str = "") -> str:
    return decode_text(dataset.attrs.get(key, default)) or default


def _preferred_field_name(dataset: xr.Dataset, name: str) -> str | None:
    if f"{name}_ADJUSTED" in dataset:
        return f"{name}_ADJUSTED"
    if name in dataset:
        return name
    return None


def _preferred_qc_name(dataset: xr.Dataset, name: str) -> str | None:
    if f"{name}_ADJUSTED_QC" in dataset:
        return f"{name}_ADJUSTED_QC"
    if f"{name}_QC" in dataset:
        return f"{name}_QC"
    return None


def _valid_profile_count(dataset: xr.Dataset, name: str) -> int:
    qc_name = _preferred_qc_name(dataset, name)
    if qc_name is not None:
        text = np.asarray(dataset[qc_name].values).astype("U4")
        text = np.char.strip(text)
        valid = np.isin(text, ["1", "8"])
        return int(np.sum(np.any(valid, axis=1)))

    field_name = _preferred_field_name(dataset, name)
    if field_name is None:
        return 0
    values = _safe_numeric(dataset[field_name].values)
    return int(np.sum(np.any(_valid_numeric(values), axis=1)))


def _profile_start_index(times: np.ndarray) -> int | None:
    for index, item in enumerate(times.tolist()):
        if item is not None:
            return int(index)
    return None


def _profile_end_index(times: np.ndarray) -> int | None:
    for index in range(times.size - 1, -1, -1):
        if times[index] is not None:
            return int(index)
    return None


def _first_valid_location(dataset: xr.Dataset, index: int | None) -> tuple[float | None, float | None]:
    if index is None or "LATITUDE" not in dataset or "LONGITUDE" not in dataset:
        return None, None
    lat = _safe_numeric(dataset["LATITUDE"].values)
    lon = _safe_numeric(dataset["LONGITUDE"].values)
    lat_value = float(lat[index]) if index < lat.size and _valid_numeric(lat[index:index + 1]).item() else None
    lon_value = float(lon[index]) if index < lon.size and _valid_numeric(lon[index:index + 1]).item() else None
    return lat_value, lon_value


def _resolve_region_label(lat: float | None, lon: float | None) -> str:
    if lat is None or lon is None or not np.isfinite(lat) or not np.isfinite(lon):
        return ""
    region = label_region(float(lon), float(lat))
    return "" if region == "Unknown" else region


def _truthy_flag(value: object) -> int:
    text = str(value).strip().lower()
    return 1 if text in {"1", "y", "yes", "true"} else 0


def _float_or_none(value: float | None) -> float | None:
    if value is None or not np.isfinite(value):
        return None
    return float(value)


def _coefficient_row(config: MeopConfig) -> dict[str, dict[str, object]]:
    path = config.tablesdir / "table_coeff.csv"
    if not path.exists():
        return {}
    frame = pd.read_csv(path)
    if "comment" in frame.columns and "comments" not in frame.columns:
        frame = frame.rename(columns={"comment": "comments"})
    frame = frame.fillna("")
    records: dict[str, dict[str, object]] = {}
    for row in frame.to_dict(orient="records"):
        smru = str(row.get("smru_platform_code", "")).strip()
        if smru:
            records[smru] = row
    return records


def _variable_offset_map(config: MeopConfig) -> dict[str, int]:
    path = config.tablesdir / "table_salinity_offsets.csv"
    if not path.exists():
        return {}
    frame = pd.read_csv(path)
    mapping: dict[str, int] = {}
    for smru in frame.get("smru_platform_code", pd.Series(dtype=str)).astype(str):
        smru_name = smru.strip()
        if smru_name:
            mapping[smru_name] = 1
    return mapping


def _discover_processed_products(config: MeopConfig) -> dict[str, dict[str, Path]]:
    inventory: dict[str, dict[str, Path]] = {}
    root = config.final_dataset_dir
    if not root.exists():
        return inventory
    for deployment_dir in sorted(path for path in root.iterdir() if path.is_dir()):
        deployment = deployment_dir.name
        by_smru: dict[str, Path] = {}
        for qf in PREFERRED_QF_ORDER:
            for path in sorted(deployment_dir.glob(f"*_{qf}_prof.nc")):
                smru = path.name.split("_")[0]
                by_smru.setdefault(smru, path)
        if by_smru:
            inventory[deployment] = by_smru
    return inventory


def _resolve_output_dir(config: MeopConfig, output_dir: str | Path | None = None) -> Path:
    if output_dir is not None:
        resolved = Path(output_dir)
        resolved.mkdir(parents=True, exist_ok=True)
        return resolved
    if not config.has_publish_version:
        dev_dir = config.datadir / "batch" / "latest" / "metadata_summaries"
        dev_dir.mkdir(parents=True, exist_ok=True)
        return dev_dir
    candidates = [config.publicdir_ctd, config.publicdir]
    for candidate in candidates:
        if (candidate / "list_tags.csv").exists() or (candidate / "list_deployments.csv").exists():
            candidate.mkdir(parents=True, exist_ok=True)
            return candidate
    config.publicdir_ctd.mkdir(parents=True, exist_ok=True)
    return config.publicdir_ctd


def _existing_membership(path: Path) -> dict[str, set[str]]:
    if not path.exists():
        return {}
    frame = pd.read_csv(path)
    if frame.empty or "DEPLOYMENT_CODE" not in frame.columns or "SMRU_PLATFORM_CODE" not in frame.columns:
        return {}
    membership: dict[str, set[str]] = {}
    for deployment, smru in zip(frame["DEPLOYMENT_CODE"], frame["SMRU_PLATFORM_CODE"], strict=False):
        dep = str(deployment).strip()
        tag = str(smru).strip()
        if dep and tag:
            membership.setdefault(dep, set()).add(tag)
    return membership


def _tag_inventory_membership(inventory: dict[str, dict[str, Path]]) -> dict[str, set[str]]:
    return {deployment: set(paths.keys()) for deployment, paths in inventory.items()}


def _summarize_tag_file(
    config: MeopConfig,
    path: Path,
    *,
    deployment_catalog: dict[str, object],
    hr_catalog: dict[str, dict[str, str]],
    coeff_rows: dict[str, dict[str, object]],
    variable_offsets: dict[str, int],
) -> dict[str, object]:
    with _open_dataset_lazy(path) as dataset:
        smru_name = _decode_scalar_attr(dataset, "smru_platform_code", path.name.split("_")[0]) or path.name.split("_")[0]
        deployment_code = _decode_scalar_attr(dataset, "deployment_code", deployment_from_smru_name(smru_name)) or deployment_from_smru_name(smru_name)

        times = _decode_juld(dataset["JULD"].values) if "JULD" in dataset else np.full((0,), None, dtype=object)
        start_idx = _profile_start_index(times)
        end_idx = _profile_end_index(times)
        start_time = times[start_idx] if start_idx is not None else None
        end_time = times[end_idx] if end_idx is not None else None
        lat, lon = _first_valid_location(dataset, start_idx)

        record = deployment_catalog.get(deployment_code)
        public_value = _truthy_flag(getattr(record, "public", "")) if record is not None else 0
        coeff = coeff_rows.get(smru_name, {})
        hr_row = hr_catalog.get(smru_name, {})

        year_value = str(hr_row.get("year", "")).strip()
        if not year_value and isinstance(start_time, datetime):
            year_value = str(start_time.year)

        row = {
            "SMRU_PLATFORM_CODE": smru_name,
            "DEPLOYMENT_CODE": deployment_code,
            "JULD": _format_datetime(start_time if isinstance(start_time, datetime) else None),
            "LATITUDE": _float_or_none(lat),
            "LONGITUDE": _float_or_none(lon),
            "MASK": _resolve_region_label(lat, lon),
            "JULD_END": _format_datetime(end_time if isinstance(end_time, datetime) else None),
            "N_PROF_TEMP": _valid_profile_count(dataset, "TEMP"),
            "N_PROF_PSAL": _valid_profile_count(dataset, "PSAL"),
            "N_PROF_CHLA": _valid_profile_count(dataset, "CHLA"),
            "PUBLIC": public_value,
            "T1": coeff.get("T1", ""),
            "T2": coeff.get("T2", ""),
            "S1": coeff.get("S1", ""),
            "S2": coeff.get("S2", ""),
            "remove": coeff.get("remove", ""),
            "Sremove": coeff.get("Sremove", ""),
            "comments": coeff.get("comments", coeff.get("comment", "OK")) or "OK",
            "variable_offset": variable_offsets.get(smru_name, ""),
            "instr_id": _decode_scalar_attr(dataset, "instr_id", str(hr_row.get("instr_id", ""))).strip(),
            "year": year_value,
            "period": str(hr_row.get("period", "")).strip(),
            "continuous": str(hr_row.get("continuous", "")).strip(),
            "prefix": str(hr_row.get("prefix", "")).strip(),
        }
        return row


def _build_deployments_table(tags: pd.DataFrame, config: MeopConfig) -> pd.DataFrame:
    if tags.empty:
        return pd.DataFrame(columns=list(DEPLOYMENT_FIELD_ORDER))

    deployment_catalog = load_deployment_catalog(config, persist=False)
    frame = tags.copy()
    frame["JULD_TS"] = pd.to_datetime(frame["JULD"], errors="coerce")
    frame["LATITUDE_NUM"] = pd.to_numeric(frame["LATITUDE"], errors="coerce")
    frame["LONGITUDE_NUM"] = pd.to_numeric(frame["LONGITUDE"], errors="coerce")
    frame["N_PROF_TEMP_NUM"] = pd.to_numeric(frame["N_PROF_TEMP"], errors="coerce").fillna(0)
    frame["N_PROF_PSAL_NUM"] = pd.to_numeric(frame["N_PROF_PSAL"], errors="coerce").fillna(0)
    frame["N_PROF_CHLA_NUM"] = pd.to_numeric(frame["N_PROF_CHLA"], errors="coerce").fillna(0)

    rows: list[dict[str, object]] = []
    for deployment, subset in frame.groupby("DEPLOYMENT_CODE", sort=True):
        record = deployment_catalog.get(str(deployment))
        start = subset["JULD_TS"].dropna()
        row = {
            "DEPLOYMENT_CODE": deployment,
            "JULD": _format_datetime(start.min().to_pydatetime() if not start.empty else None),
            "LATITUDE": float(subset["LATITUDE_NUM"].mean()) if subset["LATITUDE_NUM"].notna().any() else None,
            "LONGITUDE": float(subset["LONGITUDE_NUM"].mean()) if subset["LONGITUDE_NUM"].notna().any() else None,
            "N_PROF_TEMP": int(subset["N_PROF_TEMP_NUM"].sum()),
            "N_PROF_PSAL": int(subset["N_PROF_PSAL_NUM"].sum()),
            "N_PROF_CHLA": int(subset["N_PROF_CHLA_NUM"].sum()),
            "N_TAGS": int(len(subset)),
            "PI_CODE": getattr(record, "pi_code", "") if record is not None else "",
            "PROCESS": _truthy_flag(getattr(record, "process", "")) if record is not None else 0,
            "PUBLIC": _truthy_flag(getattr(record, "public", "")) if record is not None else 0,
            "COUNTRY": getattr(record, "country", "") if record is not None else "",
            "FIRST_VERSION": getattr(record, "extra", {}).get("first_version", "") if record is not None else "",
            "LAST_VERSION": getattr(record, "extra", {}).get("last_version", "") if record is not None else "",
        }
        rows.append(row)

    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(columns=list(DEPLOYMENT_FIELD_ORDER))
    return result.loc[:, list(DEPLOYMENT_FIELD_ORDER)].sort_values(["DEPLOYMENT_CODE"]).reset_index(drop=True)


def _ordered_tags_frame(rows: Iterable[dict[str, object]]) -> pd.DataFrame:
    frame = pd.DataFrame(list(rows))
    if frame.empty:
        return pd.DataFrame(columns=list(TAG_FIELD_ORDER))
    for column in TAG_FIELD_ORDER:
        if column not in frame.columns:
            frame[column] = ""
    frame = frame.loc[:, list(TAG_FIELD_ORDER)]
    return frame.sort_values(["DEPLOYMENT_CODE", "SMRU_PLATFORM_CODE"]).reset_index(drop=True)


def _frames_equal(left: pd.DataFrame, right: pd.DataFrame) -> bool:
    try:
        return left.fillna("").equals(right.fillna(""))
    except Exception:
        return False


def update_metadata_summaries(
    config: MeopConfig,
    *,
    processed_deployments: Iterable[str] | None = None,
    force: bool = False,
    output_dir: str | Path | None = None,
) -> SummaryUpdateResult:
    destination = _resolve_output_dir(config, output_dir=output_dir)
    list_tags_path = destination / "list_tags.csv"
    list_deployments_path = destination / "list_deployments.csv"

    inventory = _discover_processed_products(config)
    current_membership = _tag_inventory_membership(inventory)
    existing_membership = _existing_membership(list_tags_path)

    impacted = {str(item).strip() for item in (processed_deployments or []) if str(item).strip()}
    impacted.update({deployment for deployment in set(current_membership) | set(existing_membership) if current_membership.get(deployment, set()) != existing_membership.get(deployment, set())})
    if force:
        impacted.update(current_membership.keys())
        impacted.update(existing_membership.keys())

    existing_tags = pd.read_csv(list_tags_path) if list_tags_path.exists() else pd.DataFrame(columns=list(TAG_FIELD_ORDER))
    if existing_tags.empty and inventory and not impacted:
        impacted.update(current_membership.keys())

    if not impacted and list_tags_path.exists() and list_deployments_path.exists():
        return SummaryUpdateResult(
            output_dir=destination,
            list_tags_path=list_tags_path,
            list_deployments_path=list_deployments_path,
            impacted_deployments=tuple(),
            written=False,
        )

    deployment_catalog = load_deployment_catalog(config, persist=False)
    hr_catalog = load_hr_catalog(config, persist=False)
    coeff_rows = _coefficient_row(config)
    variable_offsets = _variable_offset_map(config)

    preserved = existing_tags[~existing_tags["DEPLOYMENT_CODE"].astype(str).isin(impacted)].copy() if not existing_tags.empty else pd.DataFrame(columns=list(TAG_FIELD_ORDER))
    refreshed_rows: list[dict[str, object]] = []
    for deployment in sorted(impacted):
        for path in inventory.get(deployment, {}).values():
            refreshed_rows.append(
                _summarize_tag_file(
                    config,
                    path,
                    deployment_catalog=deployment_catalog,
                    hr_catalog=hr_catalog,
                    coeff_rows=coeff_rows,
                    variable_offsets=variable_offsets,
                )
            )

    pieces: list[pd.DataFrame] = []
    if not preserved.empty:
        pieces.append(preserved)
    refreshed_frame = _ordered_tags_frame(refreshed_rows)
    if not refreshed_frame.empty:
        pieces.append(refreshed_frame)
    combined = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=list(TAG_FIELD_ORDER))
    combined = _ordered_tags_frame(combined.to_dict(orient="records"))
    deployments = _build_deployments_table(combined, config)

    previous_tags = existing_tags.loc[:, [column for column in TAG_FIELD_ORDER if column in existing_tags.columns]] if not existing_tags.empty else pd.DataFrame(columns=list(TAG_FIELD_ORDER))
    previous_tags = _ordered_tags_frame(previous_tags.to_dict(orient="records")) if not previous_tags.empty else pd.DataFrame(columns=list(TAG_FIELD_ORDER))
    previous_deployments = pd.read_csv(list_deployments_path) if list_deployments_path.exists() else pd.DataFrame(columns=list(DEPLOYMENT_FIELD_ORDER))
    if not previous_deployments.empty:
        for column in DEPLOYMENT_FIELD_ORDER:
            if column not in previous_deployments.columns:
                previous_deployments[column] = ""
        previous_deployments = previous_deployments.loc[:, list(DEPLOYMENT_FIELD_ORDER)].sort_values(["DEPLOYMENT_CODE"]).reset_index(drop=True)
    else:
        previous_deployments = pd.DataFrame(columns=list(DEPLOYMENT_FIELD_ORDER))

    written = False
    if force or not _frames_equal(combined, previous_tags):
        destination.mkdir(parents=True, exist_ok=True)
        combined.to_csv(list_tags_path, index=False)
        written = True
    if force or not _frames_equal(deployments, previous_deployments):
        destination.mkdir(parents=True, exist_ok=True)
        deployments.to_csv(list_deployments_path, index=False)
        written = True

    return SummaryUpdateResult(
        output_dir=destination,
        list_tags_path=list_tags_path,
        list_deployments_path=list_deployments_path,
        impacted_deployments=tuple(sorted(impacted)),
        written=written,
    )
