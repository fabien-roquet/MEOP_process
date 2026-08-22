from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import fnmatch
import sys
from pathlib import Path
from typing import Iterable

import numpy as np

from .workflows.compare import ComparisonReport, compare_netcdf_outputs
from .config.loader import load_config
from .reference.cora import cora_cells_for_track, load_cora_track
from .plotting.calibration import plot_ts_calibration
from .catalog.filenames import deployment_from_smru_name, list_fname_prof

CALIBRATION_QF_PREFERENCE = ("hr1", "lr1", "hr2", "hr0", "lr0", "fr1", "fr0")
MIN_CALIBRATION_PROFILES = 10
MIN_CALIBRATION_LEVELS = 3
CALIBRATION_COVERAGE_MIN_DBAR = 200.0
CALIBRATION_COVERAGE_MAX_DBAR = 600.0
NULL_ISLAND_TOLERANCE_DEGREES = 1e-6
DEPLOYMENT_SEGMENT_LINK_KM = 2500.0
MINOR_SEGMENT_FRACTION = 0.10
EARTH_RADIUS_KM = 6371.0088


class CalibrationStatusError(RuntimeError):
    """Expected data-availability outcome from a calibration-plot run."""

    status = "failed"

    def __init__(self, message: str, *, written_files: Iterable[Path] = ()) -> None:
        super().__init__(message)
        self.written_files = tuple(Path(path) for path in written_files)


class NoReferenceDataError(CalibrationStatusError):
    status = "no_reference"


class InsufficientReferenceDataError(CalibrationStatusError):
    status = "insufficient_reference"


class InsufficientTargetDataError(CalibrationStatusError):
    status = "insufficient_target"


class InvalidTargetDataError(CalibrationStatusError):
    status = "invalid_target"


@dataclass(frozen=True)
class TargetCalibrationSelection:
    """Accepted target profiles and their target-only calibration support."""

    indices: np.ndarray
    latitudes: np.ndarray
    longitudes: np.ndarray
    juld: np.ndarray
    component_sizes: tuple[int, ...]
    dropped_invalid: int
    dropped_minor_segments: int
    support: dict[str, dict[str, float | int | None]]


def _geographic_components(
    latitudes: np.ndarray,
    longitudes: np.ndarray,
    *,
    max_link_km: float = DEPLOYMENT_SEGMENT_LINK_KM,
) -> tuple[np.ndarray, ...]:
    """Return spatially connected profile components on a sphere."""
    from scipy.spatial import cKDTree

    lat = np.asarray(latitudes, dtype=np.float64).reshape(-1)
    lon = np.asarray(longitudes, dtype=np.float64).reshape(-1)
    if lat.shape != lon.shape:
        raise ValueError("latitude and longitude arrays must have the same shape")
    if lat.size == 0:
        return ()
    if max_link_km <= 0.0:
        raise ValueError("component link distance must be positive")

    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)
    cos_lat = np.cos(lat_rad)
    xyz = np.column_stack(
        (cos_lat * np.cos(lon_rad), cos_lat * np.sin(lon_rad), np.sin(lat_rad))
    )
    angular_radius = min(np.pi, max_link_km / EARTH_RADIUS_KM)
    chord_radius = 2.0 * np.sin(0.5 * angular_radius)
    pairs = cKDTree(xyz).query_pairs(chord_radius, output_type="ndarray")

    parent = np.arange(lat.size, dtype=np.int64)

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = int(parent[index])
        return index

    for left, right in pairs:
        left_root = find(int(left))
        right_root = find(int(right))
        if left_root != right_root:
            parent[right_root] = left_root

    grouped: dict[int, list[int]] = {}
    for index in range(lat.size):
        grouped.setdefault(find(index), []).append(index)
    components = tuple(
        np.asarray(indices, dtype=np.int64)
        for indices in sorted(grouped.values(), key=lambda item: (-len(item), item[0]))
    )
    return components


def _select_deployment_component_indices(
    latitudes: np.ndarray,
    longitudes: np.ndarray,
    juld: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, tuple[int, ...], int, int]:
    """QC positions and discard small components disconnected from the deployment."""
    lat = np.asarray(latitudes, dtype=np.float64).reshape(-1)
    lon = np.asarray(longitudes, dtype=np.float64).reshape(-1)
    time = np.asarray(juld, dtype=np.float64).reshape(-1)
    if lat.shape != lon.shape or lat.shape != time.shape:
        raise ValueError("target LATITUDE, LONGITUDE and JULD must have matching shapes")

    normalised_lon = ((lon + 180.0) % 360.0) - 180.0
    null_island = (
        np.abs(lat) <= NULL_ISLAND_TOLERANCE_DEGREES
    ) & (np.abs(normalised_lon) <= NULL_ISLAND_TOLERANCE_DEGREES)
    valid = (
        np.isfinite(lat)
        & np.isfinite(lon)
        & np.isfinite(time)
        & (lat >= -90.0)
        & (lat <= 90.0)
        & (np.abs(lon) <= 360.0)
        & (np.abs(time) < 99999.0)
        & ~null_island
    )
    valid_indices = np.flatnonzero(valid)
    dropped_invalid = int(lat.size - valid_indices.size)
    if valid_indices.size == 0:
        return valid_indices, normalised_lon, (), dropped_invalid, 0

    components = _geographic_components(
        lat[valid_indices], normalised_lon[valid_indices]
    )
    component_sizes = tuple(int(component.size) for component in components)
    largest = component_sizes[0]
    minimum_component = max(
        MIN_CALIBRATION_PROFILES,
        int(np.ceil(largest * MINOR_SEGMENT_FRACTION)),
    )
    retained_components = [
        component for component in components if component.size >= minimum_component
    ]
    if not retained_components:
        retained_components = [components[0]]
    retained_relative = np.sort(np.concatenate(retained_components))
    retained = valid_indices[retained_relative]
    dropped_minor = int(valid_indices.size - retained.size)
    return retained, normalised_lon, component_sizes, dropped_invalid, dropped_minor


def _selected_target_values(dataset, names: tuple[str, ...], indices: np.ndarray) -> np.ndarray:
    name = next((candidate for candidate in names if candidate in dataset), None)
    if name is None:
        raise ValueError(f"target file is missing {' or '.join(names)}")
    data_array = dataset[name]
    if "N_PROF" in data_array.dims:
        data_array = data_array.isel(N_PROF=indices)
    values = np.asarray(data_array.values, dtype=np.float64)
    values[np.abs(values) >= 9999.0] = np.nan
    return values


def _target_variable_support(
    pressure: np.ndarray,
    values: np.ndarray,
) -> dict[str, float | int | None]:
    data = np.asarray(values, dtype=np.float64)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.ndim != 2:
        raise ValueError("target T/S variables must be one- or two-dimensional")
    pres = np.asarray(pressure, dtype=np.float64)
    if pres.ndim == 1:
        pres = np.broadcast_to(pres.reshape(1, -1), data.shape)
    if pres.shape != data.shape:
        raise ValueError("target pressure and T/S arrays have incompatible shapes")

    finite = np.isfinite(pres) & np.isfinite(data)
    linear = finite & (pres >= 200.0) & (pres <= 1000.0)
    band = finite & (pres >= 400.0) & (pres <= 600.0)
    valid_pressure = pres[finite]
    return {
        "n_profiles_linear": int(np.count_nonzero(np.any(linear, axis=1))),
        "n_levels_linear": int(np.count_nonzero(np.any(linear, axis=0))),
        "n_profiles_band": int(np.count_nonzero(np.any(band, axis=1))),
        "n_levels_band": int(np.count_nonzero(np.any(band, axis=0))),
        "pressure_min": (
            float(np.nanmin(valid_pressure)) if valid_pressure.size else None
        ),
        "pressure_max": (
            float(np.nanmax(valid_pressure)) if valid_pressure.size else None
        ),
    }


def _target_support_is_sufficient(support: dict[str, float | int | None]) -> bool:
    pressure_min = support["pressure_min"]
    pressure_max = support["pressure_max"]
    return bool(
        int(support["n_profiles_linear"] or 0) >= MIN_CALIBRATION_PROFILES
        and int(support["n_levels_linear"] or 0) >= MIN_CALIBRATION_LEVELS
        and int(support["n_profiles_band"] or 0) >= MIN_CALIBRATION_PROFILES
        and int(support["n_levels_band"] or 0) >= MIN_CALIBRATION_LEVELS
        and pressure_min is not None
        and float(pressure_min) <= CALIBRATION_COVERAGE_MIN_DBAR
        and pressure_max is not None
        and float(pressure_max) >= CALIBRATION_COVERAGE_MAX_DBAR
    )


def _preflight_calibration_target(path: Path) -> TargetCalibrationSelection:
    """Select the deployment component and establish target-only T/S support."""
    import xarray as xr

    with xr.open_dataset(path, decode_times=False) as dataset:
        for required in ("LATITUDE", "LONGITUDE", "JULD"):
            if required not in dataset:
                raise ValueError(f"target file is missing {required}")
        latitudes = np.asarray(dataset["LATITUDE"].values, dtype=np.float64)
        longitudes = np.asarray(dataset["LONGITUDE"].values, dtype=np.float64)
        juld = np.asarray(dataset["JULD"].values, dtype=np.float64)
        (
            indices,
            normalised_longitudes,
            component_sizes,
            dropped_invalid,
            dropped_minor,
        ) = _select_deployment_component_indices(latitudes, longitudes, juld)
        if indices.size == 0:
            raise ValueError("no valid deployment positions remain after coordinate QC")
        pressure = _selected_target_values(
            dataset, ("PRES_ADJUSTED", "PRES"), indices
        )
        temperature = _selected_target_values(
            dataset, ("TEMP_ADJUSTED", "TEMP"), indices
        )
        salinity = _selected_target_values(
            dataset, ("PSAL_ADJUSTED", "PSAL"), indices
        )

    support = {
        "TEMP": _target_variable_support(pressure, temperature),
        "PSAL": _target_variable_support(pressure, salinity),
    }
    return TargetCalibrationSelection(
        indices=indices,
        latitudes=np.asarray(latitudes, dtype=np.float64)[indices],
        longitudes=normalised_longitudes[indices],
        juld=np.asarray(juld, dtype=np.float64)[indices],
        component_sizes=component_sizes,
        dropped_invalid=dropped_invalid,
        dropped_minor_segments=dropped_minor,
        support=support,
    )


def _load_cora_for_track(
    cora_dir: Path,
    *,
    longitudes: np.ndarray,
    latitudes: np.ndarray,
    margin: float = 5.0,
) -> tuple[dict[str, np.ndarray], tuple[tuple[int, int], ...]]:
    """Load only the buffered cells and profile rows along an accepted track."""
    cells = cora_cells_for_track(latitudes, longitudes, margin=margin)
    data = load_cora_track(
        cora_dir,
        latitudes=latitudes,
        longitudes=longitudes,
        margin=margin,
    )
    return data, cells


def _open_dataset(path: Path):
    import xarray as xr

    last_error: Exception | None = None
    for engine in (None, "h5netcdf", "scipy"):
        try:
            kwargs = {"decode_times": False}
            if engine is not None:
                kwargs["engine"] = engine
            return xr.open_dataset(path, **kwargs)
        except Exception as error:  # pragma: no cover
            last_error = error
            continue
    raise OSError(f"Unable to open comparison dataset {path}: {last_error}")


def _iter_netcdf_files(root: Path, patterns: tuple[str, ...]) -> dict[str, Path]:
    files: dict[str, Path] = {}
    for path in sorted(root.rglob("*.nc")):
        relative = path.relative_to(root).as_posix()
        if patterns and not any(fnmatch.fnmatch(relative, pattern) for pattern in patterns):
            continue
        files[relative] = path
    return files


def _format_report(label: str, report: ComparisonReport) -> list[str]:
    lines = [f"[{label}] {'OK' if report.is_equal else 'DIFF'}"]
    for item in report.notes:
        lines.append(f"note: {item}")
    for item in report.dimension_differences:
        lines.append(f"dimension: {item}")
    for item in report.variable_differences:
        lines.append(f"variable: {item}")
    for item in report.attribute_differences:
        lines.append(f"attribute: {item}")
    return lines


def _profile_tuple(dataset, index: int) -> tuple[float, float, float]:
    return (
        float(np.asarray(dataset["JULD"].values, dtype=np.float64)[index]),
        float(np.asarray(dataset["LATITUDE"].values, dtype=np.float64)[index]),
        float(np.asarray(dataset["LONGITUDE"].values, dtype=np.float64)[index]),
    )


def _profiles_close(left: tuple[float, float, float], right: tuple[float, float, float], *, atol: float) -> bool:
    return bool(np.allclose(np.asarray(left), np.asarray(right), equal_nan=True, atol=atol))


def _format_profile(profile: tuple[float, float, float]) -> str:
    return f"JULD={profile[0]:.6f} LAT={profile[1]:.6f} LON={profile[2]:.6f}"


def _time_location_notes(reference: Path, candidate: Path, *, atol: float) -> list[str]:
    with _open_dataset(reference) as ref, _open_dataset(candidate) as cand:
        required = ("JULD", "LATITUDE", "LONGITUDE")
        missing = [name for name in required if name not in ref or name not in cand]
        if missing:
            return [f"time/location summary unavailable: missing variables {', '.join(missing)}"]

        ref_count = int(ref.sizes.get("N_PROF", 0))
        cand_count = int(cand.sizes.get("N_PROF", 0))
        notes = [f"time/location profiles reference={ref_count} candidate={cand_count}"]
        if ref_count == 0 or cand_count == 0:
            return notes

        prefix = 0
        while prefix < min(ref_count, cand_count):
            if not _profiles_close(_profile_tuple(ref, prefix), _profile_tuple(cand, prefix), atol=atol):
                break
            prefix += 1

        suffix = 0
        while suffix < min(ref_count - prefix, cand_count - prefix):
            if not _profiles_close(_profile_tuple(ref, ref_count - 1 - suffix), _profile_tuple(cand, cand_count - 1 - suffix), atol=atol):
                break
            suffix += 1

        notes.append(f"time/location aligned prefix={prefix} suffix={suffix}")
        ref_unmatched = ref_count - prefix - suffix
        cand_unmatched = cand_count - prefix - suffix
        notes.append(f"time/location unmatched reference={ref_unmatched} candidate={cand_unmatched}")

        if ref_unmatched == 0 and cand_unmatched > 0 and suffix == ref_count:
            start = _profile_tuple(cand, 0)
            end = _profile_tuple(cand, cand_unmatched - 1)
            notes.append(f"candidate has {cand_unmatched} leading extra profiles: {_format_profile(start)} -> {_format_profile(end)}")
        elif ref_unmatched == 0 and cand_unmatched > 0 and prefix == ref_count:
            start = _profile_tuple(cand, prefix)
            end = _profile_tuple(cand, cand_count - 1)
            notes.append(f"candidate has {cand_unmatched} trailing extra profiles: {_format_profile(start)} -> {_format_profile(end)}")
        elif cand_unmatched == 0 and ref_unmatched > 0 and suffix == cand_count:
            start = _profile_tuple(ref, 0)
            end = _profile_tuple(ref, ref_unmatched - 1)
            notes.append(f"reference has {ref_unmatched} leading extra profiles: {_format_profile(start)} -> {_format_profile(end)}")
        elif cand_unmatched == 0 and ref_unmatched > 0 and prefix == cand_count:
            start = _profile_tuple(ref, prefix)
            end = _profile_tuple(ref, ref_count - 1)
            notes.append(f"reference has {ref_unmatched} trailing extra profiles: {_format_profile(start)} -> {_format_profile(end)}")
        elif ref_unmatched > 0 or cand_unmatched > 0:
            if ref_unmatched > 0:
                notes.append(
                    "reference unmatched window: "
                    f"{_format_profile(_profile_tuple(ref, prefix))} -> {_format_profile(_profile_tuple(ref, ref_count - suffix - 1))}"
                )
            if cand_unmatched > 0:
                notes.append(
                    "candidate unmatched window: "
                    f"{_format_profile(_profile_tuple(cand, prefix))} -> {_format_profile(_profile_tuple(cand, cand_count - suffix - 1))}"
                )

        notes.append(f"reference first={_format_profile(_profile_tuple(ref, 0))} last={_format_profile(_profile_tuple(ref, ref_count - 1))}")
        notes.append(f"candidate first={_format_profile(_profile_tuple(cand, 0))} last={_format_profile(_profile_tuple(cand, cand_count - 1))}")
        return notes


def _compare_pair(
    reference: Path,
    candidate: Path,
    *,
    variables: tuple[str, ...] | None,
    attributes: tuple[str, ...] | None,
    atol: float,
) -> ComparisonReport:
    return compare_netcdf_outputs(
        reference,
        candidate,
        variables=variables,
        attributes=attributes,
        atol=atol,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare MEOP NetCDF outputs or directories of outputs.")
    parser.add_argument("reference", nargs="?", default=None, help="Reference NetCDF file or directory.")
    parser.add_argument("candidate", nargs="?", default=None, help="Candidate NetCDF file or directory.")
    parser.add_argument(
        "--include",
        action="append",
        default=None,
        help="Glob pattern for directory mode, relative to the compared roots. May be supplied multiple times.",
    )
    parser.add_argument("--variable", action="append", default=None, help="Restrict comparison to one or more variables.")
    parser.add_argument("--attribute", action="append", default=None, help="Restrict comparison to one or more global attributes.")
    parser.add_argument("--time-location", action="store_true", help="Add a JULD/LATITUDE/LONGITUDE profile alignment summary.")
    parser.add_argument("--atol", type=float, default=1e-6, help="Absolute tolerance for numeric comparisons.")
    parser.add_argument(
        "--plot1",
        metavar="SMRU_NAME",
        default=None,
        help=(
            "Generate CORA-based T/S calibration plots for the given tag.  "
            "Requires that references.cora_dir is configured in configs.json.  "
            "Plots are saved to config.plotdir / deployment / {smru_name}_calibration*.png."
        ),
    )
    parser.add_argument(
        "--config",
        metavar="FILE",
        default=None,
        help="Path to root configs.json (overrides default discovery).",
    )
    return parser


def _select_calibration_target(smru_name: str, *, config) -> tuple[Path | None, str]:
    for qf in CALIBRATION_QF_PREFERENCE:
        paths = list_fname_prof(smru_name, qf=qf, config=config)
        if paths:
            return paths[-1], qf
    paths = list_fname_prof(smru_name, config=config)
    return (paths[-1], "") if paths else (None, "")


def _diagnostic_support(written: Iterable[Path]) -> dict[str, bool]:
    """Return whether TEMP and PSAL meet the minimum calibration support gate."""
    paths = tuple(Path(path) for path in written)
    csv_path = next((path for path in paths if path.name.endswith("_calibration_offsets.csv")), None)
    support = {"TEMP": False, "PSAL": False}
    if csv_path is None or not csv_path.exists():
        return support
    with csv_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    def parse_int(row: dict[str, str] | None, field: str) -> int | None:
        if row is None:
            return None
        try:
            return int(row.get(field, "") or "")
        except (TypeError, ValueError):
            return None

    def parse_float(row: dict[str, str] | None, field: str) -> float | None:
        if row is None:
            return None
        try:
            value = float(row.get(field, "") or "")
        except (TypeError, ValueError):
            return None
        return value if np.isfinite(value) else None

    for variable in support:
        selected = {
            row.get("diagnostic", ""): row
            for row in rows
            if row.get("variable") == variable
        }
        linear = selected.get("linear_200_1000")
        band = selected.get("band_400_600")
        all_valid = selected.get("all_valid_profile_median")
        coefficient_prefix = "T" if variable == "TEMP" else "S"
        observed_min = parse_float(all_valid, "pressure_min")
        observed_max = parse_float(all_valid, "pressure_max")
        support[variable] = all(
            (
                (parse_int(linear, "n_profiles") or 0) >= MIN_CALIBRATION_PROFILES,
                (parse_int(linear, "n_levels") or 0) >= MIN_CALIBRATION_LEVELS,
                parse_float(linear, f"suggested_{coefficient_prefix}1") is not None,
                parse_float(linear, f"suggested_{coefficient_prefix}2") is not None,
                (parse_int(band, "n_profiles") or 0) >= MIN_CALIBRATION_PROFILES,
                (parse_int(band, "n_levels") or 0) >= MIN_CALIBRATION_LEVELS,
                parse_float(band, f"suggested_{coefficient_prefix}2") is not None,
                observed_min is not None
                and observed_min <= CALIBRATION_COVERAGE_MIN_DBAR,
                observed_max is not None
                and observed_max >= CALIBRATION_COVERAGE_MAX_DBAR,
            )
        )
    return support


def _ensure_sufficient_diagnostics(smru_name: str, written: Iterable[Path]) -> None:
    paths = tuple(Path(path) for path in written)
    support = _diagnostic_support(paths)
    missing = [variable for variable, usable in support.items() if not usable]
    if missing:
        raise InsufficientReferenceDataError(
            f"reference profiles were found for {smru_name}, but {', '.join(missing)} "
            f"did not meet the calibration support gate (at least "
            f"{MIN_CALIBRATION_PROFILES} profiles, {MIN_CALIBRATION_LEVELS} levels in "
            f"400-600 dbar, and observed coverage from "
            f"{CALIBRATION_COVERAGE_MIN_DBAR:.0f} to "
            f"{CALIBRATION_COVERAGE_MAX_DBAR:.0f} dbar)",
            written_files=paths,
        )


def generate_calibration_plots(smru_name: str, *, config) -> list[Path]:
    """Generate CORA calibration plots for one tag and return written paths.

    Expected target/reference availability outcomes raise ``CalibrationStatusError``
    subclasses so callers can report them separately from software failures.
    """
    if config.cora_dir is None:
        raise ValueError(
            'references.cora_dir is not set in configs.json. Add "references": '
            '{"cora_dir": "/path/to/CORA_ncfiles"} to enable calibration plots.'
        )

    target_path, target_qf = _select_calibration_target(smru_name, config=config)
    if target_path is None:
        raise InvalidTargetDataError(f"no profile files found for {smru_name!r}")

    try:
        target = _preflight_calibration_target(target_path)
    except Exception as exc:
        raise InvalidTargetDataError(f"error reading target path: {exc}") from exc

    unsupported = [
        variable
        for variable, support in target.support.items()
        if not _target_support_is_sufficient(support)
    ]
    if unsupported:
        raise InsufficientTargetDataError(
            f"{smru_name} has insufficient target support for {', '.join(unsupported)} "
            f"before CORA loading: retained {target.indices.size} profiles from spatial "
            f"components {target.component_sizes}, dropped {target.dropped_invalid} invalid "
            f"positions and {target.dropped_minor_segments} minor-component profiles; "
            f"requires at least {MIN_CALIBRATION_PROFILES} profiles, "
            f"{MIN_CALIBRATION_LEVELS} levels in 400-600 dbar, and coverage from "
            f"{CALIBRATION_COVERAGE_MIN_DBAR:.0f} to "
            f"{CALIBRATION_COVERAGE_MAX_DBAR:.0f} dbar"
        )

    cora_data, candidate_cells = _load_cora_for_track(
        config.cora_dir,
        longitudes=target.longitudes,
        latitudes=target.latitudes,
        margin=5.0,
    )
    if np.asarray(cora_data.get("lat", ())).size == 0:
        raise NoReferenceDataError(
            f"no CORA profiles found for {smru_name} within the 5-degree corridor "
            f"around {target.indices.size} accepted profiles in "
            f"{len(candidate_cells)} candidate cells"
        )

    deployment = deployment_from_smru_name(smru_name)
    context_qf = target_qf or "*"
    other_paths = [
        p for p in list_fname_prof(deployment=deployment, qf=context_qf, config=config)
        if p != target_path
    ]

    output_dir = config.plotdir / deployment
    written = plot_ts_calibration(
        smru_name,
        cora_data=cora_data,
        target_path=target_path,
        target_profile_indices=target.indices,
        other_paths=other_paths,
        output_dir=output_dir,
        write_diagnostics=True,
    )
    _ensure_sufficient_diagnostics(smru_name, written)
    return written


def _run_calibration_plots(smru_name: str, config_file: str | None) -> int:
    """Generate CORA calibration plots for *smru_name* and return an exit code."""
    try:
        config = load_config(config_file=config_file, require_config=True)
    except (FileNotFoundError, ValueError) as exc:
        print(f"WARNING: {exc}", file=sys.stderr)
        return 2
    try:
        written = generate_calibration_plots(smru_name, config=config)
    except CalibrationStatusError as exc:
        for path in exc.written_files:
            print(f"wrote: {path}")
        print(f"{exc.status}: {exc}", file=sys.stderr)
        return 1
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    for path in written:
        print(f"wrote: {path}")
    return 0


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.plot1:
        return _run_calibration_plots(args.plot1, getattr(args, "config", None))

    if args.reference is None or args.candidate is None:
        parser.error("reference and candidate are required unless --plot1 is used")

    reference = Path(args.reference)
    candidate = Path(args.candidate)
    patterns = tuple(args.include or ())
    variables = tuple(args.variable) if args.variable else None
    attributes = tuple(args.attribute) if args.attribute else (() if variables is not None else None)

    if reference.is_file() and candidate.is_file():
        report = _compare_pair(reference, candidate, variables=variables, attributes=attributes, atol=args.atol)
        if args.time_location:
            report.notes.extend(_time_location_notes(reference, candidate, atol=args.atol))
        print("\n".join(_format_report(reference.name, report)))
        return 0 if report.is_equal else 1

    if reference.is_dir() and candidate.is_dir():
        ref_files = _iter_netcdf_files(reference, patterns)
        cand_files = _iter_netcdf_files(candidate, patterns)
        labels = sorted(set(ref_files).union(cand_files))
        had_difference = False
        for label in labels:
            ref_path = ref_files.get(label)
            cand_path = cand_files.get(label)
            if ref_path is None:
                print(f"[{label}] DIFF")
                print(f"note: missing from reference: {label}")
                had_difference = True
                continue
            if cand_path is None:
                print(f"[{label}] DIFF")
                print(f"note: missing from candidate: {label}")
                had_difference = True
                continue
            report = _compare_pair(ref_path, cand_path, variables=variables, attributes=attributes, atol=args.atol)
            if args.time_location:
                report.notes.extend(_time_location_notes(ref_path, cand_path, atol=args.atol))
            print("\n".join(_format_report(label, report)))
            had_difference = had_difference or not report.is_equal
        return 1 if had_difference else 0

    parser.error("reference and candidate must both be files or both be directories")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main(sys.argv[1:]))
