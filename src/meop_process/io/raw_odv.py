from __future__ import annotations

import csv
import filecmp
import logging
import math
import shutil
import warnings
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np

from ..catalog.tables import read_csv_rows
from ..data.layout import resolve_table_path
from ..config.paths import ensure_runtime_directories
from ..models import MeopConfig

logger = logging.getLogger(__name__)


RAW_TEXT_SUFFIXES = (
    "_ODV.txt",
    "_CTD_ODV.txt",
    "_FL_ODV.txt",
)


@dataclass(frozen=True)
class RawOdvFiles:
    """Resolved raw ODV assets for one deployment.

    Search and staging are intentionally package-owned so callers no longer need to know
    whether the operator provided a zip archive, plain text ODV files, or separate CTD/FL
    files. All supported inputs are mirrored into ``config.raw_odv_dir``.
    """

    deployment: str
    raw_root: Path
    source_root: Path
    archive: Path
    combined_text: Path | None
    ctd_text: Path | None
    fl_text: Path | None

    @property
    def preferred_ctd_text(self) -> Path | None:
        if self.combined_text is not None and self.combined_text.exists():
            return self.combined_text
        if self.ctd_text is not None and self.ctd_text.exists():
            return self.ctd_text
        return None

    @property
    def has_combined_text(self) -> bool:
        return self.combined_text is not None and self.combined_text.exists()

    @property
    def has_ctd_text(self) -> bool:
        return self.ctd_text is not None and self.ctd_text.exists()

    @property
    def has_fl_text(self) -> bool:
        return self.fl_text is not None and self.fl_text.exists()

    def as_dict(self) -> dict[str, str | None]:
        return {
            "deployment": self.deployment,
            "raw_root": str(self.raw_root),
            "source_root": str(self.source_root),
            "archive": str(self.archive),
            "combined_text": str(self.combined_text) if self.combined_text else None,
            "ctd_text": str(self.ctd_text) if self.ctd_text else None,
            "fl_text": str(self.fl_text) if self.fl_text else None,
            "preferred_ctd_text": str(self.preferred_ctd_text) if self.preferred_ctd_text else None,
        }


@dataclass(frozen=True)
class OdvProfile:
    smru_name: str
    station: int
    timestamp: str
    longitude: float
    latitude: float
    pressure: tuple[float, ...]
    temperature: tuple[float, ...]
    salinity: tuple[float, ...]
    fluorescence: tuple[float, ...] = ()
    light: tuple[float, ...] = ()
    oxygen: tuple[float, ...] = ()
    conductivity: tuple[float, ...] = ()

    @property
    def key(self) -> tuple[str, int]:
        return (self.smru_name, self.station)

    @property
    def n_levels(self) -> int:
        return len(self.pressure)


@dataclass(frozen=True)
class OdvProfileIndex:
    deployment: str
    smru_names: tuple[str, ...]
    profile_count_by_smru: dict[str, int]
    total_profiles: int
    max_levels: int


def _same_file(source: Path, destination: Path) -> bool:
    return destination.exists() and filecmp.cmp(source, destination, shallow=False)


def _copy_if_needed(source: Path, destination: Path) -> bool:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if _same_file(source, destination):
        return False
    shutil.copy2(source, destination)
    return True


def discover_raw_odv_files(config: MeopConfig, deployment: str, *, root: Path | None = None) -> RawOdvFiles:
    ensure_runtime_directories(config)
    raw_root = root or config.raw_odv_dir
    archive = raw_root / f"{deployment}_ODV.zip"
    combined = raw_root / f"{deployment}_ODV.txt"
    ctd = raw_root / f"{deployment}_CTD_ODV.txt"
    fl = raw_root / f"{deployment}_FL_ODV.txt"
    return RawOdvFiles(
        deployment=deployment,
        raw_root=raw_root,
        source_root=raw_root,
        archive=archive,
        combined_text=combined if combined.exists() else None,
        ctd_text=ctd if ctd.exists() else None,
        fl_text=fl if fl.exists() else None,
    )


def import_raw_data_zip(config: MeopConfig, deployment: str) -> bool:
    """Prepare low-resolution raw inputs directly from ``data/data_raw/raw_smru_data_odv``.

    Raw archives and extracted text files are now staged directly under ``config.raw_odv_dir``.
    If a ``<deployment>_ODV.zip`` archive is present there, it is extracted in place. Plain text
    ODV files already staged in that directory are used as-is. No ``input/`` mirror is required.
    """

    if not deployment:
        return False

    ensure_runtime_directories(config)
    target_dir = config.raw_odv_dir
    target_dir.mkdir(parents=True, exist_ok=True)

    archive_path = target_dir / f"{deployment}_ODV.zip"
    changed = False
    if archive_path.exists():
        with zipfile.ZipFile(archive_path) as archive:
            archive.extractall(path=target_dir)
        changed = True

    combined = target_dir / f"{deployment}_ODV.txt"
    ctd = target_dir / f"{deployment}_CTD_ODV.txt"
    fl = target_dir / f"{deployment}_FL_ODV.txt"

    if not combined.exists() and ctd.exists():
        shutil.copy2(ctd, combined)
        changed = True

    return archive_path.exists() or combined.exists() or ctd.exists() or fl.exists() or changed


def _normalize_header(name: str) -> str:
    normalized = name.strip().lower()
    normalized = normalized.replace("[", " ").replace("]", " ")
    normalized = normalized.replace("(", " ").replace(")", " ")
    normalized = normalized.replace("/", " ")
    normalized = " ".join(normalized.split())
    return normalized


def _to_float(value: str) -> float:
    text = value.strip()
    if not text:
        return math.nan
    number = float(text)
    if number == 999 or number == 999.0:
        return math.nan
    return number


def _parse_timestamp(text: str) -> datetime:
    value = str(text).strip()
    # Detect placeholder/template dates that were never filled in
    if value in ("yyyy-mm-dd hh:mm", "Loc. date=yyyy/MM/dd HH:mm:ss", "yyyy/mm/dd hh:mm:ss"):
        warnings.warn(
            f"Detected placeholder/template date format (not a real timestamp): {value!r}. "
            f"Using current UTC time as fallback.",
            stacklevel=3,
        )
        logger.warning(f"Template date format detected: {value!r}")
        return datetime.now(timezone.utc)
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
    warnings.warn(
        f"Unsupported timestamp format: {text!r}. Using current UTC time as fallback.",
        stacklevel=3,
    )
    logger.warning(f"Could not parse timestamp: {text!r}")
    return datetime.now(timezone.utc)


def _sensor_key(column_name: str) -> str | None:
    name = _normalize_header(column_name)
    if name.startswith("pressure"):
        return "pressure"
    if name.startswith("temperature"):
        return "temperature"
    if name.startswith("salinity"):
        return "salinity"
    if name.startswith("fluorescence"):
        return "fluorescence"
    if name.startswith("light"):
        return "light"
    if name.startswith("oxygen"):
        return "oxygen"
    if name.startswith("conductivity"):
        return "conductivity"
    return None


def _iter_odv_rows(path: Path) -> tuple[list[str], Iterable[list[str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter=";")
        try:
            first = next(reader)
        except StopIteration:
            return [], []
        if first and first[0].startswith("//"):
            try:
                header = next(reader)
            except StopIteration:
                return [], []
        else:
            header = first
        rows = []
        for line_num, row in enumerate(reader, start=3 if first and first[0].startswith("//") else 2):
            # Skip empty rows
            if not any(cell.strip() for cell in row):
                continue
            # Validate row has reasonable number of fields
            if row and len(header) > 0:
                # Allow some flexibility in field count (±2 fields)
                if len(row) < len(header) - 2 or len(row) > len(header) + 2:
                    warnings.warn(
                        f"Malformed ODV row at line {line_num}: expected ~{len(header)} fields, got {len(row)}. "
                        f"Row will be skipped: {row[:5]}...",
                        stacklevel=3,
                    )
                    logger.debug(f"Full malformed row: {row}")
                    continue
            rows.append(row)
    return header, rows


def read_odv_profiles(path: str | Path) -> list[OdvProfile]:
    """Parse a generic ODV text file into profile records.

    Each new non-empty ``Cruise`` field starts a new profile; continuation rows inherit the
    current profile identity and only contribute additional vertical levels.
    """

    source = Path(path)
    if not source.exists():
        return []

    header, rows = _iter_odv_rows(source)
    if not header:
        return []

    sensor_columns = {
        _sensor_key(name): index
        for index, name in enumerate(header)
        if _sensor_key(name) is not None
    }

    profiles: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    for row in rows:
        values = list(row) + [""] * max(0, len(header) - len(row))
        cruise = values[0].strip() if len(values) > 0 else ""
        station = values[1].strip() if len(values) > 1 else ""
        timestamp = values[3].strip() if len(values) > 3 else ""
        longitude = values[4].strip() if len(values) > 4 else ""
        latitude = values[5].strip() if len(values) > 5 else ""

        if cruise:
            if current is not None:
                profiles.append(current)
            current = {
                "smru_name": cruise,
                "station": int(float(station)) if station else len(profiles) + 1,
                "timestamp": timestamp,
                "longitude": _to_float(longitude),
                "latitude": _to_float(latitude),
                "pressure": [],
                "temperature": [],
                "salinity": [],
                "fluorescence": [],
                "light": [],
                "oxygen": [],
                "conductivity": [],
            }

        if current is None:
            continue

        for sensor in ("pressure", "temperature", "salinity", "fluorescence", "light", "oxygen", "conductivity"):
            index = sensor_columns.get(sensor)
            if index is None or index >= len(values):
                value = math.nan
            else:
                value = _to_float(values[index])
            current[sensor].append(value)

    if current is not None:
        profiles.append(current)

    return [
        OdvProfile(
            smru_name=str(profile["smru_name"]),
            station=int(profile["station"]),
            timestamp=str(profile["timestamp"]),
            longitude=float(profile["longitude"]),
            latitude=float(profile["latitude"]),
            pressure=tuple(float(value) for value in profile["pressure"]),
            temperature=tuple(float(value) for value in profile["temperature"]),
            salinity=tuple(float(value) for value in profile["salinity"]),
            fluorescence=tuple(float(value) for value in profile["fluorescence"]),
            light=tuple(float(value) for value in profile["light"]),
            oxygen=tuple(float(value) for value in profile["oxygen"]),
            conductivity=tuple(float(value) for value in profile["conductivity"]),
        )
        for profile in profiles
    ]


def build_odv_profile_index(files: RawOdvFiles, *, config: MeopConfig | None = None) -> OdvProfileIndex:
    profiles = load_raw_odv_profiles(files, config=config)
    if not profiles:
        return OdvProfileIndex(
            deployment=files.deployment,
            smru_names=(),
            profile_count_by_smru={},
            total_profiles=0,
            max_levels=0,
        )

    counts: dict[str, int] = {}
    max_levels = 0
    for profile in profiles:
        counts[profile.smru_name] = counts.get(profile.smru_name, 0) + 1
        if profile.n_levels > max_levels:
            max_levels = profile.n_levels

    return OdvProfileIndex(
        deployment=files.deployment,
        smru_names=tuple(sorted(counts)),
        profile_count_by_smru=dict(sorted(counts.items())),
        total_profiles=sum(counts.values()),
        max_levels=max_levels,
    )


def _resample_to_pressure(pressure: tuple[float, ...], values: tuple[float, ...], target_pressure: tuple[float, ...]) -> tuple[float, ...]:
    """Resample sensor values to a target pressure grid using linear interpolation.
    
    If source pressure grid differs from target, interpolate values. Otherwise return as-is.
    """
    if not values or not pressure or not target_pressure:
        return values
    if pressure == target_pressure:
        return values
    
    # Use numpy to interpolate
    p_src = np.asarray(pressure, dtype=np.float64)
    v_src = np.asarray(values, dtype=np.float64)
    p_tgt = np.asarray(target_pressure, dtype=np.float64)
    
    # Only interpolate if we have valid data
    valid_mask = np.isfinite(p_src) & np.isfinite(v_src)
    if not np.any(valid_mask):
        # No valid data, return NaNs for target grid
        return tuple(np.full_like(p_tgt, np.nan, dtype=float))
    
    p_valid = p_src[valid_mask]
    v_valid = v_src[valid_mask]
    
    # Sort by pressure
    sort_idx = np.argsort(p_valid)
    p_valid = p_valid[sort_idx]
    v_valid = v_valid[sort_idx]
    
    # Interpolate to target pressure grid
    v_interp = np.interp(p_tgt, p_valid, v_valid, left=np.nan, right=np.nan)
    return tuple(float(v) for v in v_interp)


def merge_sensor_profiles(ctd_profiles: Iterable[OdvProfile], fl_profiles: Iterable[OdvProfile]) -> list[OdvProfile]:
    index = {profile.key: profile for profile in fl_profiles}
    merged: list[OdvProfile] = []
    for ctd_profile in ctd_profiles:
        extra = index.get(ctd_profile.key)
        if extra is None:
            merged.append(ctd_profile)
            continue
        
        # Resample FL sensor values to match CTD pressure grid
        fluorescence = _resample_to_pressure(extra.pressure, extra.fluorescence, ctd_profile.pressure)
        light = _resample_to_pressure(extra.pressure, extra.light, ctd_profile.pressure)
        oxygen = _resample_to_pressure(extra.pressure, extra.oxygen, ctd_profile.pressure)
        
        # Warn if resampling was necessary
        if extra.pressure != ctd_profile.pressure and any((extra.fluorescence, extra.light, extra.oxygen)):
            logger.debug(
                f"{ctd_profile.smru_name} station {ctd_profile.station}: "
                f"Resampled FL sensors from {len(extra.pressure)} to {len(ctd_profile.pressure)} levels"
            )
        
        merged.append(
            OdvProfile(
                smru_name=ctd_profile.smru_name,
                station=ctd_profile.station,
                timestamp=ctd_profile.timestamp,
                longitude=ctd_profile.longitude,
                latitude=ctd_profile.latitude,
                pressure=ctd_profile.pressure,
                temperature=ctd_profile.temperature,
                salinity=ctd_profile.salinity,
                fluorescence=fluorescence,
                light=light,
                oxygen=oxygen,
                conductivity=ctd_profile.conductivity,
            )
        )
    return merged


def _load_split_tag_rules(config: MeopConfig | None) -> dict[str, int]:
    if config is None:
        return {}
    path = resolve_table_path(config, "table_split_tags.csv", required=False)
    rows = read_csv_rows(path)
    rules: dict[str, int] = {}
    for row in rows:
        smru_name = str(row.get("smru_platform_name", "")).strip()
        if not smru_name:
            continue
        try:
            nsplit = int(float(row.get("nsplit", "0") or 0))
        except ValueError:
            continue
        if nsplit > 1:
            rules[smru_name] = nsplit
    return rules


def _split_profiles_by_gap(profiles: list[OdvProfile], nsplit: int) -> list[OdvProfile]:
    if nsplit <= 1 or len(profiles) <= 1:
        return profiles

    timestamps = np.asarray([_parse_timestamp(profile.timestamp).timestamp() for profile in profiles], dtype=np.float64)
    if timestamps.size < 2:
        return profiles

    gaps = np.diff(timestamps)
    split_count = min(nsplit - 1, gaps.size)
    if split_count <= 0:
        return profiles
    breakpoints = np.argsort(gaps)[::-1][:split_count] + 1
    boundaries = [0] + sorted(int(point) for point in breakpoints.tolist()) + [len(profiles)]

    split_profiles: list[OdvProfile] = []
    for part_index, (start, end) in enumerate(zip(boundaries[:-1], boundaries[1:], strict=False), start=1):
        for profile in profiles[start:end]:
            split_profiles.append(
                OdvProfile(
                    smru_name=f"{profile.smru_name}-N{part_index}",
                    station=profile.station,
                    timestamp=profile.timestamp,
                    longitude=profile.longitude,
                    latitude=profile.latitude,
                    pressure=profile.pressure,
                    temperature=profile.temperature,
                    salinity=profile.salinity,
                    fluorescence=profile.fluorescence,
                    light=profile.light,
                    oxygen=profile.oxygen,
                    conductivity=profile.conductivity,
                )
            )
    return split_profiles


def _apply_split_tag_rules(profiles: list[OdvProfile], config: MeopConfig | None) -> list[OdvProfile]:
    rules = _load_split_tag_rules(config)
    if not rules:
        return profiles

    grouped: dict[str, list[OdvProfile]] = {}
    order: list[str] = []
    for profile in profiles:
        if profile.smru_name not in grouped:
            grouped[profile.smru_name] = []
            order.append(profile.smru_name)
        grouped[profile.smru_name].append(profile)

    result: list[OdvProfile] = []
    for smru_name in order:
        group = grouped[smru_name]
        nsplit = rules.get(smru_name, 1)
        if nsplit > 1:
            group = [profile for profile in group if np.isfinite(profile.longitude) or np.isfinite(profile.latitude)]
        if nsplit > 1:
            result.extend(_split_profiles_by_gap(group, nsplit))
        else:
            result.extend(group)
    return result


def load_raw_odv_profiles(files: RawOdvFiles, *, config: MeopConfig | None = None) -> list[OdvProfile]:
    ctd_source = files.preferred_ctd_text
    if ctd_source is None:
        return []
    ctd_profiles = read_odv_profiles(ctd_source)
    # If a separate FL file exists, merge it regardless of whether a combined CTD file is present.
    # Some deployments have {dep}_ODV.txt (CTD-only) alongside {dep}_FL_ODV.txt with sensor data.
    if not files.has_fl_text:
        return _apply_split_tag_rules(ctd_profiles, config)
    fl_profiles = read_odv_profiles(files.fl_text)
    return _apply_split_tag_rules(merge_sensor_profiles(ctd_profiles, fl_profiles), config)
