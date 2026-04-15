from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np
import xarray as xr


@dataclass
class ComparisonReport:
    reference_files: list[Path] = field(default_factory=list)
    candidate_files: list[Path] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    dimension_differences: list[str] = field(default_factory=list)
    variable_differences: list[str] = field(default_factory=list)
    attribute_differences: list[str] = field(default_factory=list)

    @property
    def is_equal(self) -> bool:
        return not (self.dimension_differences or self.variable_differences or self.attribute_differences)


_DEF_OPEN_ENGINES = (None, "h5netcdf", "scipy")


def _open_dataset(path: Path) -> xr.Dataset:
    last_error: Exception | None = None
    for engine in _DEF_OPEN_ENGINES:
        try:
            kwargs = {"decode_times": False}
            if engine is not None:
                kwargs["engine"] = engine
            return xr.open_dataset(path, **kwargs)
        except Exception as error:  # pragma: no cover - exercised indirectly by multiple file formats
            last_error = error
            continue
    raise OSError(f"Unable to open comparison dataset {path}: {last_error}")


def _normalise_value(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="ignore").strip()
    if isinstance(value, str):
        return value.strip()
    return value


def _values_equal(left, right, *, atol: float = 1e-6) -> bool:
    left_array = np.asarray(left)
    right_array = np.asarray(right)
    if left_array.shape != right_array.shape:
        return False
    if np.issubdtype(left_array.dtype, np.number) and np.issubdtype(right_array.dtype, np.number):
        return bool(np.allclose(left_array, right_array, equal_nan=True, atol=atol))
    return np.array_equal(np.vectorize(_normalise_value)(left_array), np.vectorize(_normalise_value)(right_array))


def compare_netcdf_outputs(
    reference: str | Path,
    candidate: str | Path,
    *,
    variables: Iterable[str] | None = None,
    attributes: Iterable[str] | None = None,
    atol: float = 1e-6,
) -> ComparisonReport:
    report = ComparisonReport(reference_files=[Path(reference)], candidate_files=[Path(candidate)])
    with _open_dataset(Path(reference)) as ref, _open_dataset(Path(candidate)) as cand:
        ref_dims = {name: int(size) for name, size in ref.sizes.items()}
        cand_dims = {name: int(size) for name, size in cand.sizes.items()}
        if variables is None:
            selected_variables = sorted(set(ref.data_vars).intersection(cand.data_vars))
        else:
            selected_variables = list(variables)
        if attributes is None:
            selected_attributes = sorted(set(ref.attrs).intersection(cand.attrs))
        else:
            selected_attributes = list(attributes)

        for name in sorted(set(ref_dims).union(cand_dims)):
            if ref_dims.get(name) != cand_dims.get(name):
                report.dimension_differences.append(
                    f"dimension {name}: reference={ref_dims.get(name)} candidate={cand_dims.get(name)}"
                )

        for name in selected_variables:
            if name not in ref:
                report.variable_differences.append(f"variable {name}: missing from reference")
                continue
            if name not in cand:
                report.variable_differences.append(f"variable {name}: missing from candidate")
                continue
            if ref[name].shape != cand[name].shape:
                report.variable_differences.append(
                    f"variable {name}: shape reference={ref[name].shape} candidate={cand[name].shape}"
                )
                continue
            if not _values_equal(ref[name].values, cand[name].values, atol=atol):
                report.variable_differences.append(f"variable {name}: values differ")

        for name in selected_attributes:
            if name not in ref.attrs:
                report.attribute_differences.append(f"attribute {name}: missing from reference")
                continue
            if name not in cand.attrs:
                report.attribute_differences.append(f"attribute {name}: missing from candidate")
                continue
            left = _normalise_value(ref.attrs[name])
            right = _normalise_value(cand.attrs[name])
            if isinstance(left, (float, int)) and isinstance(right, (float, int)):
                if not np.isclose(left, right, equal_nan=True, atol=atol):
                    report.attribute_differences.append(f"attribute {name}: reference={left!r} candidate={right!r}")
            elif left != right:
                report.attribute_differences.append(f"attribute {name}: reference={left!r} candidate={right!r}")

    return report


def compare_reference_outputs(*args, **kwargs) -> ComparisonReport:
    return compare_netcdf_outputs(*args, **kwargs)
