"""CORA-based TS calibration plots.

For a given target tag, produce a two-panel figure per chunk of ≤200 profiles:

  Left panel  — T/S diagram
    * CORA background profiles (light grey)
    * Other MEOP profiles for the same deployment (light blue)
    * Target profiles coloured by time (viridis)

  Right panel — Salinity anomaly vs pressure
    * Anomaly of each target PSAL profile against the per-pressure-level CORA
      median, coloured by time (viridis)

Output files are saved alongside the existing diagnostics plots under
``config.plotdir / deployment / {smru_name}_calibration*.png``.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Sequence
import warnings

import numpy as np


DIAGNOSTIC_FIELDNAMES = (
    "smru_platform_code",
    "variable",
    "diagnostic",
    "pressure_min",
    "pressure_max",
    "n_points",
    "n_profiles",
    "n_levels",
    "median_anomaly",
    "mean_anomaly",
    "mad_anomaly",
    "q25_anomaly",
    "q75_anomaly",
    "slope_per_dbar",
    "intercept",
    "suggested_T1",
    "suggested_T2",
    "suggested_S1",
    "suggested_S2",
    "note",
)

DIAGNOSTIC_BANDS = (
    ("band_400_600", 400.0, 600.0, "primary deep comparison band"),
    ("band_200_400", 200.0, 400.0, "upper comparison band"),
    ("band_600_1000", 600.0, 1000.0, "deep comparison band"),
    ("band_0_200", 0.0, 200.0, "shallow comparison band"),
    ("band_1000_2000", 1000.0, 2000.0, "very deep comparison band"),
)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _open_any(path: Path):
    """Open a NetCDF file using xarray, trying multiple engines."""
    import xarray as xr

    last_error: Exception | None = None
    for engine in (None, "h5netcdf", "scipy"):
        try:
            kwargs = {"decode_times": False}
            if engine is not None:
                kwargs["engine"] = engine
            return xr.open_dataset(path, **kwargs)
        except Exception as exc:
            last_error = exc
    raise OSError(f"Cannot open {path}: {last_error}")


def _extract_profiles(
    path: Path,
    *,
    profile_indices: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (pres_grid, temp, psal, juld) from a MEOP lr0/lr1 prof file.

    All output arrays are float64.  ``pres_grid`` shape is (N_LEVELS,),
    ``temp`` / ``psal`` shape is (N_PROF, N_LEVELS), ``juld`` shape is (N_PROF,).
    """
    with _open_any(path) as ds:
        def selected_values(name: str) -> np.ndarray:
            data_array = ds[name]
            if profile_indices is not None and "N_PROF" in data_array.dims:
                data_array = data_array.isel(N_PROF=profile_indices)
            return np.asarray(data_array.values, dtype=np.float64)

        pres = np.asarray(
            selected_values("PRES_ADJUSTED")
            if "PRES_ADJUSTED" in ds
            else selected_values("PRES"),
            dtype=np.float64,
        )
        temp = np.asarray(
            selected_values("TEMP_ADJUSTED")
            if "TEMP_ADJUSTED" in ds
            else selected_values("TEMP"),
            dtype=np.float64,
        )
        psal = np.asarray(
            selected_values("PSAL_ADJUSTED")
            if "PSAL_ADJUSTED" in ds
            else selected_values("PSAL"),
            dtype=np.float64,
        )
        juld = selected_values("JULD")
        # Replace fill values (≥9999) with NaN
        for arr in (pres, temp, psal):
            arr[arr >= 9999.0] = np.nan

        # pres may be (N_PROF, N_LEVELS) → take median as shared grid
        if pres.ndim == 2:
            pres_grid = _nanmedian_no_warning(pres, axis=0)
        else:
            pres_grid = pres
        return pres_grid, temp, psal, juld


def _same_pressure_grid(left: np.ndarray, right: np.ndarray) -> bool:
    left = np.asarray(left, dtype=np.float64).reshape(-1)
    right = np.asarray(right, dtype=np.float64).reshape(-1)
    return left.shape == right.shape and bool(np.allclose(left, right, equal_nan=True))


def _project_matrix_to_grid(source_pres: np.ndarray, values: np.ndarray, target_pres: np.ndarray) -> np.ndarray | None:
    source_pres = np.asarray(source_pres, dtype=np.float64).reshape(-1)
    target_pres = np.asarray(target_pres, dtype=np.float64).reshape(-1)
    values = np.asarray(values, dtype=np.float64)
    if values.ndim == 1:
        values = values.reshape(1, -1)
    if values.ndim != 2 or values.shape[1] != source_pres.size:
        return None
    if _same_pressure_grid(source_pres, target_pres):
        return values

    projected = np.full((values.shape[0], target_pres.size), np.nan, dtype=np.float64)
    finite_pres = np.isfinite(source_pres)
    for row in range(values.shape[0]):
        valid = finite_pres & np.isfinite(values[row])
        if np.count_nonzero(valid) < 2:
            continue
        order = np.argsort(source_pres[valid])
        src_pres = source_pres[valid][order]
        src_values = values[row][valid][order]
        projected[row] = np.interp(target_pres, src_pres, src_values, left=np.nan, right=np.nan)
    return projected


def _nanmedian_no_warning(values: np.ndarray, axis: int = 0) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(values, axis=axis)


def _nanpercentile_no_warning(values: np.ndarray, q: float, axis: int = 0) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanpercentile(values, q, axis=axis)


def _finite_sample(values: np.ndarray, *, max_points: int = 50_000) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64).reshape(-1)
    finite = array[np.isfinite(array)]
    if finite.size <= max_points:
        return finite
    step = int(np.ceil(finite.size / max_points))
    return finite[::step]


def _robust_range(
    arrays: Sequence[np.ndarray],
    *,
    min_span: float,
    lower_bound: float,
    upper_bound: float,
) -> tuple[float, float] | None:
    samples = [_finite_sample(array) for array in arrays]
    samples = [sample for sample in samples if sample.size > 0]
    if not samples:
        return None
    values = np.concatenate(samples)
    values = values[(values >= lower_bound) & (values <= upper_bound)]
    if values.size < 2:
        return None
    lo, hi = np.nanpercentile(values, [1.0, 99.0])
    if not np.isfinite(lo) or not np.isfinite(hi):
        return None
    if hi - lo < min_span:
        center = 0.5 * (lo + hi)
        lo = center - 0.5 * min_span
        hi = center + 0.5 * min_span
    margin = max(0.05 * (hi - lo), 0.05 * min_span)
    return float(lo - margin), float(hi + margin)


def _add_potential_density_contours(
    ax,
    *,
    psal_arrays: Sequence[np.ndarray],
    temp_arrays: Sequence[np.ndarray],
) -> None:
    """Draw sigma0 contours behind the T/S curves when gsw is available."""
    sal_range = _robust_range(psal_arrays, min_span=0.2, lower_bound=20.0, upper_bound=45.0)
    temp_range = _robust_range(temp_arrays, min_span=0.5, lower_bound=-5.0, upper_bound=35.0)
    if sal_range is None or temp_range is None:
        return
    try:
        import gsw  # type: ignore
    except ImportError:
        return

    salinity = np.linspace(sal_range[0], sal_range[1], 80)
    temperature = np.linspace(temp_range[0], temp_range[1], 80)
    salinity_grid, temperature_grid = np.meshgrid(salinity, temperature)
    try:
        absolute_salinity = gsw.SA_from_SP(salinity_grid, 0.0, 0.0, 0.0)
        conservative_temperature = gsw.CT_from_t(absolute_salinity, temperature_grid, 0.0)
        sigma0 = gsw.sigma0(absolute_salinity, conservative_temperature)
    except Exception:
        return
    if not np.isfinite(sigma0).any():
        return

    sigma_min = float(np.nanmin(sigma0))
    sigma_max = float(np.nanmax(sigma0))
    level_start = np.floor(sigma_min * 2.0) / 2.0
    level_end = np.ceil(sigma_max * 2.0) / 2.0
    levels = np.arange(level_start, level_end + 0.25, 0.5)
    if levels.size < 2:
        return
    contours = ax.contour(
        salinity_grid,
        temperature_grid,
        sigma0,
        levels=levels,
        colors="0.55",
        linewidths=0.45,
        alpha=0.55,
        zorder=0,
    )
    ax.clabel(contours, inline=True, fontsize=6, fmt=lambda value: f"{value:g}")


def _fmt_float(value: float | np.floating | None, *, digits: int = 6) -> str:
    if value is None:
        return ""
    value = float(value)
    if not np.isfinite(value):
        return ""
    return f"{value:.{digits}g}"


def _valid_counts(anomaly: np.ndarray, level_mask: np.ndarray) -> tuple[int, int, int]:
    if not np.any(level_mask):
        return 0, 0, 0
    values = anomaly[:, level_mask]
    finite = np.isfinite(values)
    return (
        int(np.count_nonzero(finite)),
        int(np.count_nonzero(np.any(finite, axis=1))),
        int(np.count_nonzero(np.any(finite, axis=0))),
    )


def _profile_median_values(anomaly: np.ndarray, level_mask: np.ndarray) -> np.ndarray:
    if not np.any(level_mask):
        return np.empty(0, dtype=np.float64)
    values = _nanmedian_no_warning(anomaly[:, level_mask], axis=1)
    return values[np.isfinite(values)]


def _summary_stats(values: np.ndarray) -> dict[str, float | None]:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {
            "median_anomaly": None,
            "mean_anomaly": None,
            "mad_anomaly": None,
            "q25_anomaly": None,
            "q75_anomaly": None,
        }
    median = float(np.nanmedian(values))
    return {
        "median_anomaly": median,
        "mean_anomaly": float(np.nanmean(values)),
        "mad_anomaly": float(np.nanmedian(np.abs(values - median))),
        "q25_anomaly": float(np.nanpercentile(values, 25.0)),
        "q75_anomaly": float(np.nanpercentile(values, 75.0)),
    }


def _blank_diagnostic_row(
    smru_name: str,
    *,
    variable: str,
    diagnostic: str,
    pressure_min: float | None,
    pressure_max: float | None,
    n_points: int,
    n_profiles: int,
    n_levels: int,
    note: str,
) -> dict[str, str]:
    return {
        "smru_platform_code": smru_name,
        "variable": variable,
        "diagnostic": diagnostic,
        "pressure_min": _fmt_float(pressure_min),
        "pressure_max": _fmt_float(pressure_max),
        "n_points": str(n_points),
        "n_profiles": str(n_profiles),
        "n_levels": str(n_levels),
        "median_anomaly": "",
        "mean_anomaly": "",
        "mad_anomaly": "",
        "q25_anomaly": "",
        "q75_anomaly": "",
        "slope_per_dbar": "",
        "intercept": "",
        "suggested_T1": "",
        "suggested_T2": "",
        "suggested_S1": "",
        "suggested_S2": "",
        "note": note,
    }


def _coefficient_fields(variable: str, *, slope: float | None = None, intercept: float | None = None) -> dict[str, str]:
    fields = {"suggested_T1": "", "suggested_T2": "", "suggested_S1": "", "suggested_S2": ""}
    if variable == "PSAL":
        fields["suggested_S1"] = _fmt_float(None if slope is None else slope * 1000.0)
        fields["suggested_S2"] = _fmt_float(intercept)
    elif variable == "TEMP":
        fields["suggested_T1"] = _fmt_float(None if slope is None else slope * 1000.0)
        fields["suggested_T2"] = _fmt_float(intercept)
    return fields


def _band_diagnostic_row(
    smru_name: str,
    *,
    variable: str,
    anomaly: np.ndarray,
    pres_grid: np.ndarray,
    diagnostic: str,
    pressure_min: float | None,
    pressure_max: float | None,
    note: str,
) -> dict[str, str]:
    finite_pres = np.isfinite(pres_grid)
    level_mask = finite_pres.copy()
    if pressure_min is not None:
        level_mask &= pres_grid >= pressure_min
    if pressure_max is not None:
        level_mask &= pres_grid <= pressure_max
    n_points, n_profiles, n_levels = _valid_counts(anomaly, level_mask)
    profile_values = _profile_median_values(anomaly, level_mask)
    row = _blank_diagnostic_row(
        smru_name,
        variable=variable,
        diagnostic=diagnostic,
        pressure_min=pressure_min,
        pressure_max=pressure_max,
        n_points=n_points,
        n_profiles=n_profiles,
        n_levels=n_levels,
        note=f"{note}; profile-median tag-minus-CORA anomaly; suggested coefficient is the value to subtract",
    )
    stats = _summary_stats(profile_values)
    for key, value in stats.items():
        row[key] = _fmt_float(value)
    median = stats["median_anomaly"]
    row.update(_coefficient_fields(variable, intercept=median))
    return row


def _linear_diagnostic_row(
    smru_name: str,
    *,
    variable: str,
    anomaly: np.ndarray,
    pres_grid: np.ndarray,
    pressure_min: float = 200.0,
    pressure_max: float = 1000.0,
) -> dict[str, str]:
    finite_pres = np.isfinite(pres_grid)
    level_mask = finite_pres & (pres_grid >= pressure_min) & (pres_grid <= pressure_max)
    n_points, n_profiles, n_levels = _valid_counts(anomaly, level_mask)
    level_median = _nanmedian_no_warning(anomaly, axis=0)
    fit_mask = level_mask & np.isfinite(level_median)
    row = _blank_diagnostic_row(
        smru_name,
        variable=variable,
        diagnostic=f"linear_{int(pressure_min)}_{int(pressure_max)}",
        pressure_min=pressure_min,
        pressure_max=pressure_max,
        n_points=n_points,
        n_profiles=n_profiles,
        n_levels=n_levels,
        note="level-median tag-minus-CORA anomaly fit as intercept + slope * pressure",
    )
    y = level_median[fit_mask]
    x = pres_grid[fit_mask]
    stats = _summary_stats(y)
    for key, value in stats.items():
        row[key] = _fmt_float(value)
    if x.size < 3:
        return row

    slope, intercept = np.polyfit(x, y, 1)
    residual = y - (slope * x + intercept)
    residual_median = float(np.nanmedian(residual))
    residual_mad = float(np.nanmedian(np.abs(residual - residual_median)))
    if np.isfinite(residual_mad) and residual_mad > 0:
        keep = np.abs(residual - residual_median) <= 3.0 * 1.4826 * residual_mad
        if np.count_nonzero(keep) >= 3 and np.count_nonzero(keep) < keep.size:
            slope, intercept = np.polyfit(x[keep], y[keep], 1)

    row["slope_per_dbar"] = _fmt_float(float(slope))
    row["intercept"] = _fmt_float(float(intercept))
    row.update(_coefficient_fields(variable, slope=float(slope), intercept=float(intercept)))
    return row


def _all_valid_pressure_span(anomaly: np.ndarray, pres_grid: np.ndarray) -> tuple[float | None, float | None]:
    finite_level = np.any(np.isfinite(anomaly), axis=0) & np.isfinite(pres_grid)
    if not np.any(finite_level):
        return None, None
    values = pres_grid[finite_level]
    return float(np.nanmin(values)), float(np.nanmax(values))


def _calibration_diagnostic_rows(
    smru_name: str,
    *,
    pres_grid: np.ndarray,
    temp_anomaly: np.ndarray,
    psal_anomaly: np.ndarray,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for variable, anomaly in (("PSAL", psal_anomaly), ("TEMP", temp_anomaly)):
        for diagnostic, pmin, pmax, note in DIAGNOSTIC_BANDS:
            rows.append(
                _band_diagnostic_row(
                    smru_name,
                    variable=variable,
                    anomaly=anomaly,
                    pres_grid=pres_grid,
                    diagnostic=diagnostic,
                    pressure_min=pmin,
                    pressure_max=pmax,
                    note=note,
                )
            )
        all_min, all_max = _all_valid_pressure_span(anomaly, pres_grid)
        rows.append(
            _band_diagnostic_row(
                smru_name,
                variable=variable,
                anomaly=anomaly,
                pres_grid=pres_grid,
                diagnostic="all_valid_profile_median",
                pressure_min=all_min,
                pressure_max=all_max,
                note="all valid target/CORA comparison levels",
            )
        )
        rows.append(
            _linear_diagnostic_row(
                smru_name,
                variable=variable,
                anomaly=anomaly,
                pres_grid=pres_grid,
            )
        )
    return rows


def _write_diagnostics_csv(path: Path, rows: Sequence[dict[str, str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=DIAGNOSTIC_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _row_float(row: dict[str, str] | None, key: str) -> float | None:
    if row is None:
        return None
    value = row.get(key, "")
    if value == "":
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if np.isfinite(parsed) else None


def _find_diagnostic_row(rows: Sequence[dict[str, str]], *, variable: str, diagnostic: str) -> dict[str, str] | None:
    for row in rows:
        if row.get("variable") == variable and row.get("diagnostic") == diagnostic:
            return row
    return None


def _annotation_text(rows: Sequence[dict[str, str]], *, variable: str) -> str:
    band = _find_diagnostic_row(rows, variable=variable, diagnostic="band_400_600")
    linear = _find_diagnostic_row(rows, variable=variable, diagnostic="linear_200_1000")
    if variable == "PSAL":
        band_value = _row_float(band, "suggested_S2")
        linear_1 = _row_float(linear, "suggested_S1")
        linear_2 = _row_float(linear, "suggested_S2")
        prefix = "PSAL"
        linear_label = "S1/S2"
    else:
        band_value = _row_float(band, "suggested_T2")
        linear_1 = _row_float(linear, "suggested_T1")
        linear_2 = _row_float(linear, "suggested_T2")
        prefix = "TEMP"
        linear_label = "T1/T2"
    band_text = "n/a" if band_value is None else f"{band_value:+.4g}"
    if linear_1 is None or linear_2 is None:
        linear_text = "n/a"
    else:
        linear_text = f"{linear_1:+.4g} / {linear_2:+.4g}"
    return f"{prefix} 400-600 offset: {band_text}\n{linear_label} linear 200-1000: {linear_text}"


def _write_offset_diagnostic_plot(
    smru_name: str,
    *,
    pres_grid: np.ndarray,
    temp_anomaly: np.ndarray,
    psal_anomaly: np.ndarray,
    rows: Sequence[dict[str, str]],
    output_dir: Path,
    plt,
) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(11, 7), sharey=True)
    panels = (
        (axes[0], psal_anomaly, "PSAL", "PSAL anomaly (tag - CORA, PSU)", "#1565C0"),
        (axes[1], temp_anomaly, "TEMP", "TEMP anomaly (tag - CORA, degC)", "#C62828"),
    )
    finite_pressure = np.isfinite(pres_grid)
    for ax, anomaly, variable, xlabel, color in panels:
        median = _nanmedian_no_warning(anomaly, axis=0)
        q25 = _nanpercentile_no_warning(anomaly, 25.0, axis=0)
        q75 = _nanpercentile_no_warning(anomaly, 75.0, axis=0)
        valid = finite_pressure & np.isfinite(median)
        if np.any(valid):
            ax.fill_betweenx(pres_grid[valid], q25[valid], q75[valid], color=color, alpha=0.18, linewidth=0)
            ax.plot(median[valid], pres_grid[valid], color=color, lw=1.4, label="median")
        ax.axvline(0.0, color="k", lw=0.8, ls="--")
        ax.axhspan(400.0, 600.0, color="#F9A825", alpha=0.15, linewidth=0)
        ax.grid(True, alpha=0.22, lw=0.6)
        ax.set_xlabel(xlabel)
        ax.set_title(f"{variable} offset diagnostics")
        ax.text(
            0.03,
            0.03,
            _annotation_text(rows, variable=variable),
            transform=ax.transAxes,
            va="bottom",
            ha="left",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "0.75", "alpha": 0.86},
        )
    axes[0].set_ylabel("Pressure (dbar)")
    if np.any(finite_pressure):
        finite_level = finite_pressure & (
            np.any(np.isfinite(psal_anomaly), axis=0) | np.any(np.isfinite(temp_anomaly), axis=0)
        )
        if np.any(finite_level):
            y_min = float(np.nanmin(pres_grid[finite_level]))
            y_max = float(np.nanmax(pres_grid[finite_level]))
            axes[0].set_ylim(y_max, y_min)
        else:
            axes[0].invert_yaxis()
    else:
        axes[0].invert_yaxis()
    fig.suptitle(f"{smru_name} — CORA anomaly coefficient diagnostics")
    fig.tight_layout()
    out_path = output_dir / f"{smru_name}_calibration_offsets.png"
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def plot_ts_calibration(
    smru_name: str,
    *,
    cora_data: dict[str, np.ndarray],
    target_path: Path,
    target_profile_indices: np.ndarray | None = None,
    other_paths: Sequence[Path] = (),
    output_dir: Path,
    chunk_size: int = 200,
    write_diagnostics: bool = False,
) -> list[Path]:
    """Generate CORA-based T/S calibration plots for *smru_name*.

    Parameters
    ----------
    smru_name:
        SMRU platform code of the target tag.
    cora_data:
        Output of ``load_cora_tiles(...)`` — dict with keys
        ``lat``, ``lon``, ``temp``, ``psal``, ``pres``.
    target_path:
        Path to the target tag's lr0 or lr1 prof NetCDF file.
    target_profile_indices:
        Optional accepted ``N_PROF`` rows after position/component QC. Profiles
        outside this selection do not contribute to figures or diagnostics.
    other_paths:
        Paths to other tags in the same deployment (used for context cloud).
    output_dir:
        Directory where PNG files are written.
    chunk_size:
        Maximum number of target profiles per figure.
    write_diagnostics:
        When true, also write a CSV and PNG summarising tag-minus-CORA
        PSAL/TEMP offsets in coefficient units compatible with table_coeff.csv.

    Returns
    -------
    list[Path]
        Paths to the written PNG files.
    """
    try:
        import matplotlib  # type: ignore

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
        import matplotlib.cm as cm  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise ImportError("matplotlib is required for calibration plots") from exc

    pres_grid, tgt_temp, tgt_psal, tgt_juld = _extract_profiles(
        target_path, profile_indices=target_profile_indices
    )
    n_tgt = tgt_temp.shape[0]

    # Build CORA median salinity per pressure level for anomaly panel
    cora_pres = cora_data.get("pres", np.array([]))
    cora_temp = cora_data.get("temp", np.empty((0, len(pres_grid))))
    cora_psal = cora_data.get("psal", np.empty((0, len(pres_grid))))

    # Interpolate CORA to the target pressure grid if needed
    if len(cora_pres) > 0:
        projected_temp = _project_matrix_to_grid(cora_pres, cora_temp, pres_grid)
        projected_psal = _project_matrix_to_grid(cora_pres, cora_psal, pres_grid)
        if projected_temp is not None and projected_psal is not None:
            cora_temp = projected_temp
            cora_psal = projected_psal
            cora_pres = pres_grid

    cora_temp_median = _nanmedian_no_warning(cora_temp, axis=0) if cora_temp.shape[0] > 0 else np.full_like(pres_grid, np.nan)
    cora_psal_median = _nanmedian_no_warning(cora_psal, axis=0) if cora_psal.shape[0] > 0 else np.full_like(pres_grid, np.nan)
    temp_anomaly = tgt_temp - cora_temp_median.reshape(1, -1)
    psal_anomaly = tgt_psal - cora_psal_median.reshape(1, -1)

    # Build "other" profiles cloud
    other_temp_list: list[np.ndarray] = []
    other_psal_list: list[np.ndarray] = []
    for op in other_paths:
        if not op.exists():
            continue
        try:
            other_pres, ot, op_psal, _ = _extract_profiles(op)
            projected_temp = _project_matrix_to_grid(other_pres, ot, pres_grid)
            projected_psal = _project_matrix_to_grid(other_pres, op_psal, pres_grid)
            if projected_temp is None or projected_psal is None:
                continue
            other_temp_list.append(projected_temp)
            other_psal_list.append(projected_psal)
        except Exception:
            pass

    other_temp = np.concatenate(other_temp_list, axis=0) if other_temp_list else np.empty((0, len(pres_grid)))
    other_psal = np.concatenate(other_psal_list, axis=0) if other_psal_list else np.empty((0, len(pres_grid)))

    # Colour target profiles by time
    if n_tgt > 1:
        valid_juld = np.where(np.isnan(tgt_juld), 0.0, tgt_juld)
        juld_norm = (valid_juld - valid_juld.min()) / max(valid_juld.max() - valid_juld.min(), 1.0)
    else:
        juld_norm = np.zeros(n_tgt)

    cmap = cm.viridis

    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    diagnostic_written: list[Path] = []
    if write_diagnostics:
        rows = _calibration_diagnostic_rows(
            smru_name,
            pres_grid=pres_grid,
            temp_anomaly=temp_anomaly,
            psal_anomaly=psal_anomaly,
        )
        diagnostic_written.append(
            _write_diagnostics_csv(output_dir / f"{smru_name}_calibration_offsets.csv", rows)
        )
        diagnostic_written.append(
            _write_offset_diagnostic_plot(
                smru_name,
                pres_grid=pres_grid,
                temp_anomaly=temp_anomaly,
                psal_anomaly=psal_anomaly,
                rows=rows,
                output_dir=output_dir,
                plt=plt,
            )
        )
    n_chunks = max(1, int(np.ceil(n_tgt / chunk_size)))

    for chunk_idx in range(n_chunks):
        start = chunk_idx * chunk_size
        end = min(start + chunk_size, n_tgt)
        c_temp = tgt_temp[start:end]
        c_psal = tgt_psal[start:end]
        c_norm = juld_norm[start:end]

        fig, axes = plt.subplots(1, 2, figsize=(12, 7))
        ax_ts, ax_anom = axes

        # ── Left panel: T/S diagram ──────────────────────────────────────────
        _add_potential_density_contours(
            ax_ts,
            psal_arrays=(cora_psal, other_psal, c_psal),
            temp_arrays=(cora_temp, other_temp, c_temp),
        )

        # CORA background
        for i in range(min(cora_temp.shape[0], 2000)):
            ax_ts.plot(cora_psal[i], cora_temp[i], color="0.80", lw=0.3, alpha=0.5, rasterized=True)

        # Other MEOP context
        for i in range(min(other_temp.shape[0], 500)):
            ax_ts.plot(other_psal[i], other_temp[i], color="#2196F3", lw=0.4, alpha=0.4, rasterized=True)

        # Target profiles coloured by time
        for i in range(c_temp.shape[0]):
            ax_ts.plot(c_psal[i], c_temp[i], color=cmap(c_norm[i]), lw=0.8, alpha=0.8)

        # Colour bar
        sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax_ts, fraction=0.03, pad=0.02)
        cbar.set_label("Time (normalised)")
        ax_ts.set_xlabel("Salinity (PSU)")
        ax_ts.set_ylabel("Temperature (°C)")
        ax_ts.set_title(f"{smru_name} — T/S diagram (profiles {start + 1}–{end})")

        # ── Right panel: salinity anomaly vs pressure ────────────────────────
        c_psal_anomaly = psal_anomaly[start:end]
        for i in range(c_psal.shape[0]):
            ax_anom.plot(c_psal_anomaly[i], pres_grid, color=cmap(c_norm[i]), lw=0.8, alpha=0.8)

        ax_anom.axvline(0, color="k", lw=0.8, ls="--")
        ax_anom.invert_yaxis()
        ax_anom.set_xlabel("PSAL anomaly (PSU)")
        ax_anom.set_ylabel("Pressure (dbar)")
        ax_anom.set_title(f"{smru_name} — PSAL anomaly vs CORA median")

        fig.tight_layout()

        if n_chunks == 1:
            out_path = output_dir / f"{smru_name}_calibration.png"
        else:
            out_path = output_dir / f"{smru_name}_calibration_part{chunk_idx}.png"

        fig.savefig(out_path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        written.append(out_path)

    written.extend(diagnostic_written)
    return written
