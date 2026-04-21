from __future__ import annotations

import csv
import filecmp
import math
import shutil
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
        rows = [row for row in reader if any(cell.strip() for cell in row)]
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


def merge_sensor_profiles(ctd_profiles: Iterable[OdvProfile], fl_profiles: Iterable[OdvProfile]) -> list[OdvProfile]:
    index = {profile.key: profile for profile in fl_profiles}
    merged: list[OdvProfile] = []
    for ctd_profile in ctd_profiles:
        extra = index.get(ctd_profile.key)
        if extra is None:
            merged.append(ctd_profile)
            continue
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
                fluorescence=extra.fluorescence,
                light=extra.light,
                oxygen=extra.oxygen,
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
    if files.has_combined_text:
        return _apply_split_tag_rules(ctd_profiles, config)
    if not files.has_fl_text:
        return _apply_split_tag_rules(ctd_profiles, config)
    fl_profiles = read_odv_profiles(files.fl_text)
    return _apply_split_tag_rules(merge_sensor_profiles(ctd_profiles, fl_profiles), config)
