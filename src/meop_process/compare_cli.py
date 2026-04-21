from __future__ import annotations

import argparse
import fnmatch
import sys
from pathlib import Path
from typing import Iterable

from .workflows.compare import ComparisonReport, compare_netcdf_outputs


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
    parser.add_argument("reference", help="Reference NetCDF file or directory.")
    parser.add_argument("candidate", help="Candidate NetCDF file or directory.")
    parser.add_argument(
        "--include",
        action="append",
        default=None,
        help="Glob pattern for directory mode, relative to the compared roots. May be supplied multiple times.",
    )
    parser.add_argument("--variable", action="append", default=None, help="Restrict comparison to one or more variables.")
    parser.add_argument("--attribute", action="append", default=None, help="Restrict comparison to one or more global attributes.")
    parser.add_argument("--atol", type=float, default=1e-6, help="Absolute tolerance for numeric comparisons.")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)

    reference = Path(args.reference)
    candidate = Path(args.candidate)
    patterns = tuple(args.include or ())
    variables = tuple(args.variable) if args.variable else None
    attributes = tuple(args.attribute) if args.attribute else (() if variables is not None else None)

    if reference.is_file() and candidate.is_file():
        report = _compare_pair(reference, candidate, variables=variables, attributes=attributes, atol=args.atol)
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
            print("\n".join(_format_report(label, report)))
            had_difference = had_difference or not report.is_equal
        return 1 if had_difference else 0

    parser.error("reference and candidate must both be files or both be directories")
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main(sys.argv[1:]))
