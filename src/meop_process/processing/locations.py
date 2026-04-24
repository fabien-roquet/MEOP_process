from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from ..catalog.deployments import load_info_deployment
from ..catalog.filenames import fname_prof
from ..io.netcdf import decode_text, open_meop_netcdf
from ..models import MeopConfig, Selection
from .netcdf import save_dataset_with_compression


DEFAULT_REFERENCE_DATE = datetime(1950, 1, 1, tzinfo=timezone.utc)


@dataclass(frozen=True)
class LocationSource:
    type: str
    jul: np.ndarray
    lat: np.ndarray
    lon: np.ndarray


@dataclass(frozen=True)
class LocationAdjustmentResult:
    written_files: tuple[Path, ...] = ()
    placeholder: bool = True

    def as_dict(self) -> dict[str, object]:
        return {
            "written_files": [str(path) for path in self.written_files],
            "placeholder": self.placeholder,
        }


def _smru_prefix(smru_name: str) -> str:
    return smru_name.split("-N")[0]


def _to_juld(datetimes: pd.Series | pd.DatetimeIndex) -> np.ndarray:
    parsed = pd.to_datetime(datetimes, utc=True, errors="coerce")
    values = parsed.to_numpy() if isinstance(parsed, pd.Series) else np.asarray(parsed)
    out = np.full(values.shape, np.nan, dtype=np.float64)
    for index, value in np.ndenumerate(values):
        if pd.isna(value):
            continue
        moment = pd.Timestamp(value).to_pydatetime()
        out[index] = (moment - DEFAULT_REFERENCE_DATE).total_seconds() / 86400.0
    return out


def _normalize_longitudes(values: np.ndarray) -> np.ndarray:
    out = np.asarray(values, dtype=np.float64).copy()
    mask = out > 180.0
    out[mask] = out[mask] - 360.0
    return out


def _haversine_km(lat1: np.ndarray, lon1: np.ndarray, lat2: np.ndarray, lon2: np.ndarray) -> np.ndarray:
    earth_radius_km = 6371.0
    lat1_rad = np.deg2rad(lat1)
    lon1_rad = np.deg2rad(lon1)
    lat2_rad = np.deg2rad(lat2)
    lon2_rad = np.deg2rad(lon2)
    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2.0) ** 2
    return 2.0 * earth_radius_km * np.arcsin(np.sqrt(a))


def _cumulative_track_distance_km(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    out = np.zeros(lat.shape, dtype=np.float64)
    if lat.size <= 1:
        return out
    step = _haversine_km(lat[:-1], lon[:-1], lat[1:], lon[1:])
    out[1:] = np.cumsum(step)
    return out


def _positioning_system_code(loc_algorithm: str) -> str | None:
    text = loc_algorithm.strip().upper()
    if "LEAST SQUARES" in text:
        return "LS"
    if "KALMAN" in text:
        return "K"
    if text == "SMRU":
        return "ARGOS"
    return None


def _positioning_system_array(length: int, value: str) -> np.ndarray:
    width = 8
    encoded = value.encode("ascii", "replace")[:width].ljust(width, b" ")
    out = np.full((length, width), b" ", dtype="S1")
    out[:, :] = np.frombuffer(encoded, dtype="S1")
    return out


def _load_crawl_data(config: MeopConfig, smru_name: str, ptt: str, jul: np.ndarray) -> LocationSource | None:
    candidates = sorted(config.crawl_locdir.glob(f"*_{_smru_prefix(smru_name)}_{ptt}_*_crawl.csv"))
    if not candidates:
        return None

    datemin = float(np.nanmin(jul))
    datemax = float(np.nanmax(jul))
    for path in candidates:
        frame = pd.read_csv(path)
        columns = {str(name).strip().strip('"'): name for name in frame.columns}
        date_column = columns.get("GMT") or columns.get("date")
        lon_column = columns.get("mu.x") or columns.get("mu_x") or columns.get("lon")
        lat_column = columns.get("mu.y") or columns.get("mu_y") or columns.get("lat")
        if date_column is None or lon_column is None or lat_column is None:
            continue

        loc_jul = _to_juld(frame[date_column])
        loc_lon = _normalize_longitudes(np.asarray(frame[lon_column], dtype=np.float64))
        loc_lat = np.asarray(frame[lat_column], dtype=np.float64)
        mask = (
            np.isfinite(loc_jul)
            & np.isfinite(loc_lon)
            & np.isfinite(loc_lat)
            & (loc_jul > datemin - 10.0)
            & (loc_jul < datemax + 10.0)
        )
        if not np.any(mask):
            continue
        return LocationSource(type="crawl", jul=loc_jul[mask], lat=loc_lat[mask], lon=loc_lon[mask])
    return None


def _load_cls_data(config: MeopConfig, ptt: str, jul: np.ndarray) -> LocationSource | None:
    datemin = float(np.nanmin(jul))
    datemax = float(np.nanmax(jul))
    chosen: Path | None = None
    for path in sorted(config.cls_locdir.glob(f"{ptt}_*.csv")):
        parts = path.name.split("_")
        if len(parts) < 3:
            continue
        start = pd.to_datetime(parts[1], utc=True, errors="coerce")
        end = pd.to_datetime(parts[2].split(".")[0], utc=True, errors="coerce")
        if pd.isna(start) or pd.isna(end):
            continue
        start_jul = (start.to_pydatetime() - DEFAULT_REFERENCE_DATE).total_seconds() / 86400.0
        end_jul = (end.to_pydatetime() - DEFAULT_REFERENCE_DATE).total_seconds() / 86400.0
        if datemin > start_jul - 300.0 and datemax < end_jul + 300.0:
            chosen = path
    if chosen is None:
        return None

    frame = pd.read_csv(chosen, sep=";")
    columns = {str(name).strip(): name for name in frame.columns}
    lat = np.asarray(frame[columns["Latitude"]], dtype=np.float64)
    lon = _normalize_longitudes(np.asarray(frame[columns["Longitude"]], dtype=np.float64))
    lat2 = np.asarray(frame[columns["Latitude 2"]], dtype=np.float64)
    lon2 = _normalize_longitudes(np.asarray(frame[columns["Longitude 2"]], dtype=np.float64))
    quality = frame[columns["Loc. quality"]].astype(str).str.strip().to_numpy()
    loc_jul = _to_juld(frame[columns["Loc. date=yyyy/MM/dd HH:mm:ss"]])

    order = np.argsort(loc_jul)
    loc_jul = loc_jul[order]
    lat = lat[order]
    lon = lon[order]
    lat2 = lat2[order]
    lon2 = lon2[order]
    quality = quality[order]

    _, unique_indices = np.unique(loc_jul, return_index=True)
    unique_indices.sort()
    loc_jul = loc_jul[unique_indices]
    lat = lat[unique_indices]
    lon = lon[unique_indices]
    lat2 = lat2[unique_indices]
    lon2 = lon2[unique_indices]
    quality = quality[unique_indices]

    valid_quality = quality != "Z"
    loc_jul = loc_jul[valid_quality]
    lat = lat[valid_quality]
    lon = lon[valid_quality]
    lat2 = lat2[valid_quality]
    lon2 = lon2[valid_quality]
    if loc_jul.size == 0:
        return None

    diff1 = np.abs(lon2 - lon) + np.abs(lat2 - lat)
    replace = np.flatnonzero(diff1 != 0.0)
    replace = replace[replace > 0]
    if replace.size:
        dist1 = np.abs(lon[replace] - lon[replace - 1]) + np.abs(lat[replace] - lat[replace - 1])
        dist2 = np.abs(lon2[replace] - lon[replace - 1]) + np.abs(lat2[replace] - lat[replace - 1])
        better = replace[dist2 < dist1]
        lon[better] = lon2[better]
        lat[better] = lat2[better]

    finite = np.isfinite(loc_jul) & np.isfinite(lat) & np.isfinite(lon)
    if not np.any(finite):
        return None
    return LocationSource(type="cls", jul=loc_jul[finite], lat=lat[finite], lon=lon[finite])


def _load_smru_locations(jul: np.ndarray, lat: np.ndarray, lon: np.ndarray) -> LocationSource | None:
    mask = np.isfinite(jul) & np.isfinite(lat) & np.isfinite(lon)
    if not np.any(mask):
        return None
    base_jul = jul[mask]
    base_lat = lat[mask]
    base_lon = lon[mask]
    order = np.argsort(base_jul)
    base_jul = base_jul[order]
    base_lat = base_lat[order]
    base_lon = base_lon[order]
    _, unique_indices = np.unique(base_jul, return_index=True)
    unique_indices.sort()
    return LocationSource(type="smru", jul=base_jul[unique_indices], lat=base_lat[unique_indices], lon=base_lon[unique_indices])


def _load_location_source(config: MeopConfig, smru_name: str, ptt: str, jul: np.ndarray, lat: np.ndarray, lon: np.ndarray) -> LocationSource | None:
    for loader in (
        lambda: _load_crawl_data(config, smru_name, ptt, jul),
        lambda: _load_cls_data(config, ptt, jul),
        lambda: _load_smru_locations(jul, lat, lon),
    ):
        source = loader()
        if source is not None and source.jul.size:
            return source
    return None


def _interpolate_locations(source: LocationSource, jul: np.ndarray, lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    newlat = np.interp(jul, source.jul, source.lat)
    newlon = np.interp(jul, source.jul, source.lon)
    outside = (jul < np.nanmin(source.jul) - 2.0) | (jul > np.nanmax(source.jul) + 2.0)
    newlat[outside] = lat[outside]
    newlon[outside] = lon[outside]

    cumulative = _cumulative_track_distance_km(newlat, newlon)
    delta_jul = np.diff(jul)
    with np.errstate(divide="ignore", invalid="ignore"):
        velocity = np.diff(cumulative) / delta_jul / 86.4
    fast = np.flatnonzero(velocity > 3.0)
    bad = np.intersect1d(fast, fast + 1)
    bad = bad[(bad > 0) & (bad < newlat.size - 1)]
    duplicate = np.flatnonzero(delta_jul == 0.0)
    bad = np.union1d(bad, duplicate)
    if bad.size:
        good_mask = np.ones(newlat.shape, dtype=bool)
        good_mask[bad] = False
        if np.count_nonzero(good_mask) >= 3:
            newlat[bad] = np.interp(jul[bad], jul[good_mask], newlat[good_mask])
            newlon[bad] = np.interp(jul[bad], jul[good_mask], newlon[good_mask])
    return newlat, newlon


def _atomic_save_dataset(dataset, path: Path) -> None:
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    save_dataset_with_compression(dataset, tmp_path)
    tmp_path.replace(path)


def _apply_location_adjustment_to_file(config: MeopConfig, smru_name: str, path: Path) -> bool:
    if not path.exists():
        return False
    dataset = open_meop_netcdf(path, decode_times=False)
    try:
        if "JULD" not in dataset or "LATITUDE" not in dataset or "LONGITUDE" not in dataset:
            return False
        jul = np.asarray(dataset["JULD"].values, dtype=np.float64)
        lat = np.asarray(dataset["LATITUDE"].values, dtype=np.float64)
        lon = np.asarray(dataset["LONGITUDE"].values, dtype=np.float64)
        if lat.size < 3:
            return False
        ptt = decode_text(dataset.attrs.get("ptt", ""))
        if not ptt:
            return False

        source = _load_location_source(config, smru_name, ptt, jul, lat, lon)
        if source is None or source.type not in {"cls", "crawl"}:
            return False

        updated = dataset.copy(deep=True)
        newlat, newlon = _interpolate_locations(source, jul, lat, lon)
        updated["LATITUDE"] = (dataset["LATITUDE"].dims, newlat.astype(dataset["LATITUDE"].dtype, copy=False))
        updated["LONGITUDE"] = (dataset["LONGITUDE"].dims, newlon.astype(dataset["LONGITUDE"].dtype, copy=False))
        updated["LATITUDE"].attrs.update(dataset["LATITUDE"].attrs)
        updated["LONGITUDE"].attrs.update(dataset["LONGITUDE"].attrs)

        original_algorithm = decode_text(dataset.attrs.get("loc_algorithm", ""))
        if "POSITIONING_SYSTEM" in updated:
            system = _positioning_system_code(original_algorithm)
            if system:
                variable = updated["POSITIONING_SYSTEM"]
                dims = variable.dims
                if len(dims) == 1:
                    dims = (dims[0], "STRING8")
                updated["POSITIONING_SYSTEM"] = (
                    dims,
                    _positioning_system_array(int(updated.sizes.get("N_PROF", lat.size)), system),
                )
                updated["POSITIONING_SYSTEM"].attrs.update(variable.attrs)

        updated.attrs["loc_algorithm"] = "CLS SMOOTH KALMAN" if source.type == "cls" else "CRAWL"
        print(f"location update smru_name={smru_name} ptt={ptt} source={source.type}")
        _atomic_save_dataset(updated, path)
        return True
    finally:
        dataset.close()


def apply_location_adjustment_placeholder(config: MeopConfig, selection: Selection) -> LocationAdjustmentResult:
    """Adjust LR0 positions from crawl/CLS post-processed locations when available."""

    info = load_info_deployment(config, deployment=selection.deployment, smru_name=selection.smru_name)
    if selection.smru_name:
        tags = (selection.smru_name,)
    elif info.list_tag_lr0:
        tags = tuple(path.name.split("_")[0] for path in info.list_tag_lr0)
    else:
        tags = info.list_smru_name

    written: list[Path] = []
    for smru_name in tags:
        path = fname_prof(smru_name, deployment=info.EXP, qf="lr0", config=config)
        if _apply_location_adjustment_to_file(config, smru_name, path):
            written.append(path)
    return LocationAdjustmentResult(written_files=tuple(written), placeholder=False)
