from __future__ import annotations

import argparse
import fnmatch
import sys
from pathlib import Path
from typing import Iterable

import numpy as np

from .workflows.compare import ComparisonReport, compare_netcdf_outputs
from .config.loader import load_config
from .reference.cora import load_cora_tiles
from .plotting.calibration import plot_ts_calibration
from .catalog.filenames import deployment_from_smru_name, list_fname_prof


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


def _run_calibration_plots(smru_name: str, config_file: str | None) -> int:
    """Generate CORA calibration plots for *smru_name* and return an exit code."""
    try:
        config = load_config(config_file=config_file, require_config=True)
    except (FileNotFoundError, ValueError) as exc:
        print(f"WARNING: {exc}", file=sys.stderr)
        return 2
    if config.cora_dir is None:
        print(
            f"WARNING: references.cora_dir is not set in configs.json.  "
            f"Add \"references\": {{\"cora_dir\": \"/path/to/CORA_ncfiles\"}} to enable --plot1.",
            file=sys.stderr,
        )
        return 2

    target_paths = list_fname_prof(smru_name, config=config)
    if not target_paths:
        print(f"error: no profile files found for {smru_name!r}", file=sys.stderr)
        return 2
    target_path = target_paths[-1]  # prefer lr1 if multiple qf levels exist

    # Determine bounding box from the target tag
    try:
        import xarray as xr

        with xr.open_dataset(target_path, decode_times=False) as ds:
            lats = np.asarray(ds["LATITUDE"].values, dtype=np.float64)
            lons = np.asarray(ds["LONGITUDE"].values, dtype=np.float64)
        valid_lat = lats[~np.isnan(lats)]
        valid_lon = lons[~np.isnan(lons)]
        if valid_lat.size == 0 or valid_lon.size == 0:
            print(f"error: no valid lat/lon in {target_path}", file=sys.stderr)
            return 2
        margin = 5.0
        lon_min, lon_max = float(valid_lon.min()) - margin, float(valid_lon.max()) + margin
        lat_min, lat_max = float(valid_lat.min()) - margin, float(valid_lat.max()) + margin
    except Exception as exc:
        print(f"error reading target path: {exc}", file=sys.stderr)
        return 2

    cora_data = load_cora_tiles(
        config.cora_dir,
        lon_min=lon_min,
        lon_max=lon_max,
        lat_min=lat_min,
        lat_max=lat_max,
    )

    deployment = deployment_from_smru_name(smru_name)
    other_paths = [
        p for p in list_fname_prof(deployment=deployment, config=config)
        if p != target_path
    ]

    output_dir = config.plotdir / deployment
    written = plot_ts_calibration(
        smru_name,
        cora_data=cora_data,
        target_path=target_path,
        other_paths=other_paths,
        output_dir=output_dir,
    )
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
