from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import re
from typing import Iterable

import numpy as np
import xarray as xr

from ..catalog.filenames import fname_plots, list_fname_prof
from ..io.netcdf import decode_text, juld_to_datetime, open_meop_netcdf
from ..models import MeopConfig, Selection
from ..processing.qc import _sigma0_profile, _to_numeric_qc


try:  # pragma: no cover - optional dependency
    import cartopy.crs as ccrs  # type: ignore
    import cartopy.feature as cfeature  # type: ignore
except Exception:  # pragma: no cover - exercised when cartopy is unavailable
    ccrs = None
    cfeature = None


@dataclass(frozen=True)
class DiagnosticResult:
    written_files: tuple[Path, ...]
    processed_tags: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "written_files": [str(path) for path in self.written_files],
            "processed_tags": list(self.processed_tags),
        }


@dataclass(frozen=True)
class TagDiagnosticData:
    smru_name: str
    deployment: str
    adjusted: bool
    pressure: np.ndarray
    temp: np.ndarray
    psal: np.ndarray
    sigma0: np.ndarray
    lon: np.ndarray
    lat: np.ndarray
    times: np.ndarray
    summary_lines: tuple[str, ...]


@dataclass(frozen=True)
class DeploymentDiagnosticSummary:
    deployment: str
    adjusted: bool
    smru_names: tuple[str, ...]
    n_profiles: int
    n_temp_profiles: int
    n_psal_profiles: int
    start_time: datetime | None
    end_time: datetime | None
    lon: np.ndarray
    lat: np.ndarray


def _import_matplotlib():
    import matplotlib

    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    from matplotlib.gridspec import GridSpec

    return plt, Normalize, ScalarMappable, GridSpec


def _qc_mask(dataset: xr.Dataset, name: str, *, adjusted: bool) -> np.ndarray:
    suffix = "_ADJUSTED_QC" if adjusted and f"{name}_ADJUSTED_QC" in dataset else "_QC"
    numeric = _to_numeric_qc(dataset[f"{name}{suffix}"].values)
    return numeric <= 1


def _field(dataset: xr.Dataset, name: str, *, adjusted: bool) -> np.ndarray:
    variable = f"{name}_ADJUSTED" if adjusted and f"{name}_ADJUSTED" in dataset else name
    return np.asarray(dataset[variable].values, dtype=np.float64)


def _pressure(dataset: xr.Dataset, *, adjusted: bool) -> np.ndarray:
    return _field(dataset, "PRES", adjusted=adjusted)


def _masked_field(dataset: xr.Dataset, name: str, *, adjusted: bool) -> np.ndarray:
    values = _field(dataset, name, adjusted=adjusted)
    mask = _qc_mask(dataset, name, adjusted=adjusted)
    return np.where(mask, values, np.nan)


def _smru_name(dataset: xr.Dataset, fallback: str) -> str:
    return decode_text(dataset.attrs.get("smru_platform_code", "")) or fallback


def _profile_times(dataset: xr.Dataset) -> np.ndarray:
    if "JULD" not in dataset:
        return np.asarray([None] * int(dataset.sizes.get("N_PROF", 0)), dtype=object)
    return juld_to_datetime(dataset["JULD"].values)


def _elapsed_days(times: np.ndarray) -> tuple[np.ndarray, str]:
    valid = [item for item in times.ravel().tolist() if item is not None]
    if not valid:
        n_prof = int(times.size)
        return np.arange(n_prof, dtype=float), "profile"
    start = min(valid)
    days = np.full(times.shape, np.nan, dtype=float)
    for idx, item in enumerate(times.ravel().tolist()):
        if item is not None:
            days.ravel()[idx] = (item - start).total_seconds() / 86400.0
    finite = np.isfinite(days)
    if np.count_nonzero(finite) >= 2 and np.all(np.diff(days[finite]) > 0):
        return days.astype(float), "days"
    fallback = np.arange(times.size, dtype=float)
    return fallback.reshape(times.shape), "profile"


def _profile_colors(times: np.ndarray) -> tuple[np.ndarray, object, object]:
    plt, Normalize, _, _ = _import_matplotlib()

    days, axis_label = _elapsed_days(times)
    finite = np.isfinite(days)
    if not np.any(finite):
        days = np.arange(times.size, dtype=float)
        axis_label = "profile"
        finite = np.isfinite(days)
    vmin = float(np.nanmin(days[finite]))
    vmax = float(np.nanmax(days[finite]))
    if vmax <= vmin:
        vmax = vmin + 1.0
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap("turbo")
    colors = cmap(norm(days.astype(float)))
    return colors, norm, cmap, axis_label, days.astype(float)


def _decode_datetime(item) -> str:
    if item is None:
        return "unknown"
    return item.strftime("%Y-%m-%d")


def _compute_sigma0(dataset: xr.Dataset, *, adjusted: bool) -> np.ndarray:
    pres = _pressure(dataset, adjusted=adjusted)
    temp = _masked_field(dataset, "TEMP", adjusted=adjusted)
    psal = _masked_field(dataset, "PSAL", adjusted=adjusted)
    lon = np.asarray(dataset["LONGITUDE"].values, dtype=np.float64)
    lat = np.asarray(dataset["LATITUDE"].values, dtype=np.float64)

    sigma = np.full(temp.shape, np.nan, dtype=np.float64)
    for idx in range(temp.shape[0]):
        mask = np.isfinite(temp[idx]) & np.isfinite(psal[idx]) & np.isfinite(pres[idx])
        if not np.any(mask):
            continue
        lon_i = float(lon[idx]) if np.isfinite(lon[idx]) else 0.0
        lat_i = float(lat[idx]) if np.isfinite(lat[idx]) else 0.0
        sigma[idx, mask] = _sigma0_profile(psal[idx, mask], temp[idx, mask], pres[idx, mask], lon=lon_i, lat=lat_i)
    return sigma


def _value_range(values: np.ndarray, *, pad_fraction: float = 0.03) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 0.0, 1.0
    lo, hi = np.nanpercentile(finite, [1.0, 99.0])
    if hi <= lo:
        hi = lo + 1.0
    pad = (hi - lo) * pad_fraction
    return float(lo - pad), float(hi + pad)


def _section_grid(pressure: np.ndarray, values: np.ndarray, *, pmax: int = 1000, step: int = 5) -> tuple[np.ndarray, np.ndarray]:
    depth = np.arange(0, pmax + step, step, dtype=float)
    section = np.full((values.shape[0], depth.size), np.nan, dtype=np.float64)
    for idx in range(values.shape[0]):
        mask = np.isfinite(pressure[idx]) & np.isfinite(values[idx])
        if np.count_nonzero(mask) < 2:
            continue
        p = pressure[idx, mask].astype(np.float64, copy=False)
        v = values[idx, mask].astype(np.float64, copy=False)
        order = np.argsort(p)
        p = p[order]
        v = v[order]
        keep = np.r_[True, np.diff(p) > 0]
        p = p[keep]
        v = v[keep]
        if p.size < 2:
            continue
        inside = (depth >= p[0]) & (depth <= p[-1])
        if not np.any(inside):
            continue
        section[idx, inside] = np.interp(depth[inside], p, v)
    return depth, section


def _centers_to_edges(values: Iterable[float]) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.size == 1:
        return np.asarray([arr[0] - 0.5, arr[0] + 0.5], dtype=float)
    mids = 0.5 * (arr[1:] + arr[:-1])
    left = arr[0] - (mids[0] - arr[0])
    right = arr[-1] + (arr[-1] - mids[-1])
    return np.concatenate([[left], mids, [right]])


def _extract_adjustment_lines(dataset: xr.Dataset) -> list[str]:
    lines: list[str] = []
    coeff = dataset.get("SCIENTIFIC_CALIB_COEFFICIENT")
    if coeff is not None:
        values = np.asarray(coeff.values, dtype=object)
        if values.ndim >= 3 and values.shape[0] > 0:
            for item in values[0, 0, :]:
                text = decode_text(item)
                if text and text != " ":
                    lines.append(text)
    if not lines:
        lines.append("No calibration strings found")
    return lines


def _find_first_float(text: str) -> float | None:
    match = re.search(r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", text)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def _summary_lines(dataset: xr.Dataset, *, adjusted: bool, sigma0: np.ndarray) -> list[str]:
    smru = decode_text(dataset.attrs.get("smru_platform_code", ""))
    deployment = decode_text(dataset.attrs.get("deployment_code", ""))
    location = decode_text(dataset.attrs.get("location", ""))
    species = decode_text(dataset.attrs.get("species", ""))
    tlc = decode_text(dataset.attrs.get("thermal_lag_adjustment", "")) or ("yes" if adjusted else "no")
    times = _profile_times(dataset)
    valid_times = [item for item in times.ravel().tolist() if item is not None]
    start = _decode_datetime(min(valid_times)) if valid_times else "unknown"
    end = _decode_datetime(max(valid_times)) if valid_times else "unknown"
    span = "?"
    if len(valid_times) >= 2:
        span = str(int(round((max(valid_times) - min(valid_times)).total_seconds() / 86400.0)))

    temp = _masked_field(dataset, "TEMP", adjusted=adjusted)
    psal = _masked_field(dataset, "PSAL", adjusted=adjusted)
    n_temp = int(np.sum(np.any(np.isfinite(temp), axis=1)))
    n_psal = int(np.sum(np.any(np.isfinite(psal), axis=1)))

    temp_error = np.asarray(dataset.get("TEMP_ADJUSTED_ERROR", np.nan), dtype=float)
    psal_error = np.asarray(dataset.get("PSAL_ADJUSTED_ERROR", np.nan), dtype=float)
    sigma_error = np.nan
    if np.isfinite(temp_error).any() and np.isfinite(psal_error).any():
        sigma_error = float(np.nanstd(sigma0)) if np.isfinite(sigma0).any() else np.nan

    lines = [
        smru or deployment or "MEOP tag",
        f"deployment: {deployment}" if deployment else "",
        f"species: {species}" if species else "",
        f"location: {location}" if location else "",
        f"period: {start} to {end}",
        f"span: {span} days",
        f"thermal lag: {tlc}",
        f"TEMP profiles: {n_temp}",
        f"PSAL profiles: {n_psal}",
    ]

    coeff_lines = _extract_adjustment_lines(dataset)
    pretty_coeffs: list[str] = []
    for item in coeff_lines:
        pretty_coeffs.append(item)
    lines.extend(pretty_coeffs)

    if np.isfinite(temp_error).any():
        lines.append(f"TEMP err median: {float(np.nanmedian(temp_error)):.2f}")
    if np.isfinite(psal_error).any():
        lines.append(f"PSAL err median: {float(np.nanmedian(psal_error)):.2f}")
    if np.isfinite(sigma_error):
        lines.append(f"sigma0 spread: {sigma_error:.2f}")
    return [line for line in lines if line]


def _tag_diagnostic_data(dataset: xr.Dataset, source_name: str, *, adjusted: bool) -> TagDiagnosticData:
    pressure = _pressure(dataset, adjusted=adjusted)
    temp = _masked_field(dataset, "TEMP", adjusted=adjusted)
    psal = _masked_field(dataset, "PSAL", adjusted=adjusted)
    sigma0 = _compute_sigma0(dataset, adjusted=adjusted)
    lon = np.asarray(dataset["LONGITUDE"].values, dtype=float)
    lat = np.asarray(dataset["LATITUDE"].values, dtype=float)
    times = _profile_times(dataset)
    smru_name = _smru_name(dataset, source_name)
    deployment = decode_text(dataset.attrs.get("deployment_code", "")) or source_name.split("-")[0]
    return TagDiagnosticData(
        smru_name=smru_name,
        deployment=deployment,
        adjusted=adjusted,
        pressure=pressure,
        temp=temp,
        psal=psal,
        sigma0=sigma0,
        lon=lon,
        lat=lat,
        times=times,
        summary_lines=tuple(_summary_lines(dataset, adjusted=adjusted, sigma0=sigma0)),
    )


def _profile_count(values: np.ndarray) -> int:
    return int(np.sum(np.any(np.isfinite(values), axis=1)))


def _deployment_summary_lines(tags: tuple[TagDiagnosticData, ...]) -> list[str]:
    if not tags:
        return ["No diagnostics inputs found"]

    deployment = tags[0].deployment or "deployment"
    valid_times = [item for tag in tags for item in tag.times.ravel().tolist() if item is not None]
    start = _decode_datetime(min(valid_times)) if valid_times else "unknown"
    end = _decode_datetime(max(valid_times)) if valid_times else "unknown"
    span = "?"
    if len(valid_times) >= 2:
        span = str(int(round((max(valid_times) - min(valid_times)).total_seconds() / 86400.0)))

    total_prof = sum(int(tag.pressure.shape[0]) for tag in tags)
    total_temp = sum(_profile_count(tag.temp) for tag in tags)
    total_psal = sum(_profile_count(tag.psal) for tag in tags)
    labels = ", ".join(tag.smru_name for tag in tags[:8])
    if len(tags) > 8:
        labels = f"{labels}, ..."

    return [
        f"deployment: {deployment}",
        f"tags: {len(tags)}",
        f"profiles: {total_prof}",
        f"TEMP profiles: {total_temp}",
        f"PSAL profiles: {total_psal}",
        f"period: {start} to {end}",
        f"span: {span} days",
        f"variables: {'adjusted' if tags[0].adjusted else 'raw'}",
        f"tags included: {labels}",
    ]


def _deployment_summary(tags: tuple[TagDiagnosticData, ...]) -> DeploymentDiagnosticSummary:
    valid_times = [item for tag in tags for item in tag.times.ravel().tolist() if item is not None]
    lon_parts = [tag.lon[np.isfinite(tag.lon)] for tag in tags if np.isfinite(tag.lon).any()]
    lat_parts = [tag.lat[np.isfinite(tag.lat)] for tag in tags if np.isfinite(tag.lat).any()]
    lon = np.concatenate(lon_parts) if lon_parts else np.asarray([], dtype=float)
    lat = np.concatenate(lat_parts) if lat_parts else np.asarray([], dtype=float)
    return DeploymentDiagnosticSummary(
        deployment=tags[0].deployment,
        adjusted=tags[0].adjusted,
        smru_names=tuple(tag.smru_name for tag in tags),
        n_profiles=sum(int(tag.pressure.shape[0]) for tag in tags),
        n_temp_profiles=sum(_profile_count(tag.temp) for tag in tags),
        n_psal_profiles=sum(_profile_count(tag.psal) for tag in tags),
        start_time=min(valid_times) if valid_times else None,
        end_time=max(valid_times) if valid_times else None,
        lon=lon,
        lat=lat,
    )


def _overview_summary_lines(summaries: tuple[DeploymentDiagnosticSummary, ...], *, qf: str) -> list[str]:
    total_tags = sum(len(summary.smru_names) for summary in summaries)
    total_profiles = sum(summary.n_profiles for summary in summaries)
    total_temp = sum(summary.n_temp_profiles for summary in summaries)
    total_psal = sum(summary.n_psal_profiles for summary in summaries)
    valid_starts = [summary.start_time for summary in summaries if summary.start_time is not None]
    valid_ends = [summary.end_time for summary in summaries if summary.end_time is not None]
    start = _decode_datetime(min(valid_starts)) if valid_starts else "unknown"
    end = _decode_datetime(max(valid_ends)) if valid_ends else "unknown"
    return [
        "MEOP diagnostics overview",
        f"deployments: {len(summaries)}",
        f"tags: {total_tags}",
        f"profiles: {total_profiles}",
        f"TEMP profiles: {total_temp}",
        f"PSAL profiles: {total_psal}",
        f"period: {start} to {end}",
        f"product: {qf}",
        f"variables: {'adjusted' if summaries[0].adjusted else 'raw'}",
    ]


def _plot_overview_map(ax, summaries: tuple[DeploymentDiagnosticSummary, ...], *, colors: np.ndarray, title: str) -> None:
    if not summaries:
        ax.axis("off")
        return
    mask_any = any(summary.lon.size and summary.lat.size for summary in summaries)
    if not mask_any:
        ax.set_title(title)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.grid(True, alpha=0.25)
        ax.text(0.5, 0.5, "No valid track coordinates", ha="center", va="center", transform=ax.transAxes)
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        return

    lon_all = np.concatenate([summary.lon for summary in summaries if summary.lon.size])
    central_longitude = 180 if lon_all.size and float(np.nanmax(lon_all) - np.nanmin(lon_all)) > 180 else 0

    if ccrs is not None:
        projection = ccrs.PlateCarree(central_longitude=central_longitude)
        ax.remove()
        fig = ax.figure
        ax = fig.add_subplot(ax.get_subplotspec(), projection=projection)
        for idx, summary in enumerate(summaries):
            if not summary.lon.size or not summary.lat.size:
                continue
            ax.plot(summary.lon, summary.lat, color=colors[idx], linewidth=0.9, alpha=0.8, transform=ccrs.PlateCarree(), label=summary.deployment)
        ax.coastlines(resolution="110m", linewidth=0.7)
        ax.add_feature(cfeature.LAND, facecolor="0.88")
        gl = ax.gridlines(draw_labels=True, linewidth=0.4, alpha=0.4)
        gl.top_labels = False
        gl.right_labels = False
        ax.set_title(title)
    else:
        for idx, summary in enumerate(summaries):
            if not summary.lon.size or not summary.lat.size:
                continue
            ax.plot(summary.lon, summary.lat, color=colors[idx], linewidth=0.9, alpha=0.8, label=summary.deployment)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(title)
        ax.grid(True, alpha=0.25)
        lon_pad = max(0.5, 0.05 * max(1.0, float(np.nanmax(lon_all) - np.nanmin(lon_all))))
        lat_all = np.concatenate([summary.lat for summary in summaries if summary.lat.size])
        lat_pad = max(0.5, 0.05 * max(1.0, float(np.nanmax(lat_all) - np.nanmin(lat_all))))
        ax.set_xlim(float(np.nanmin(lon_all) - lon_pad), float(np.nanmax(lon_all) + lon_pad))
        ax.set_ylim(float(np.nanmin(lat_all) - lat_pad), float(np.nanmax(lat_all) + lat_pad))

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="upper left", fontsize=8, frameon=False)


def _plot_overview_counts(ax, summaries: tuple[DeploymentDiagnosticSummary, ...], *, colors: np.ndarray) -> None:
    labels = [summary.deployment for summary in summaries]
    profiles = np.asarray([summary.n_profiles for summary in summaries], dtype=float)
    tags = np.asarray([len(summary.smru_names) for summary in summaries], dtype=float)
    x = np.arange(len(labels), dtype=float)
    ax.bar(x - 0.18, profiles, width=0.36, color=colors[: len(labels)], alpha=0.8, label="profiles")
    ax.bar(x + 0.18, tags, width=0.36, color="0.35", alpha=0.6, label="tags")
    ax.set_xticks(x, labels, rotation=30, ha="right")
    ax.set_ylabel("count")
    ax.set_title("Deployment counts")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)


def _plot_overview_timeline(ax, summaries: tuple[DeploymentDiagnosticSummary, ...], *, colors: np.ndarray) -> None:
    labels = [summary.deployment for summary in summaries]
    y = np.arange(len(labels), dtype=float)
    for idx, summary in enumerate(summaries):
        if summary.start_time is None or summary.end_time is None:
            continue
        ax.plot([summary.start_time, summary.end_time], [y[idx], y[idx]], color=colors[idx], linewidth=3.0, solid_capstyle="round")
        ax.scatter([summary.start_time, summary.end_time], [y[idx], y[idx]], color=colors[idx], s=18)
    ax.set_yticks(y, labels)
    ax.set_title("Deployment time span")
    ax.grid(True, axis="x", alpha=0.25)


def _save_global_overview_figure(summaries: tuple[DeploymentDiagnosticSummary, ...], target: Path, *, qf: str) -> None:
    plt, _, _, GridSpec = _import_matplotlib()
    cmap = plt.get_cmap("tab20", max(len(summaries), 1))
    colors = cmap(np.arange(max(len(summaries), 1)))

    fig = plt.figure(figsize=(16, 10), constrained_layout=False)
    gs = GridSpec(2, 3, figure=fig, width_ratios=[1.0, 1.15, 1.0], height_ratios=[1.0, 1.0])
    ax_info = fig.add_subplot(gs[:, 0])
    ax_map = fig.add_subplot(gs[:, 1])
    ax_counts = fig.add_subplot(gs[0, 2])
    ax_timeline = fig.add_subplot(gs[1, 1:])

    _plot_info(ax_info, _overview_summary_lines(summaries, qf=qf))
    _plot_overview_map(ax_map, summaries, colors=colors, title="Cross-deployment tracks")
    _plot_overview_counts(ax_counts, summaries, colors=colors)
    _plot_overview_timeline(ax_timeline, summaries, colors=colors)

    adjusted = summaries[0].adjusted if summaries else True
    fig.suptitle(f"MEOP overview diagnostics ({qf}, {'adjusted' if adjusted else 'raw'})", fontsize=16, y=0.985)
    fig.subplots_adjust(left=0.05, right=0.98, bottom=0.08, top=0.92, wspace=0.25, hspace=0.28)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_labelled_profiles(ax, tags: tuple[TagDiagnosticData, ...], *, field_name: str, colors: np.ndarray, pmax: int, xlabel: str) -> None:
    total = 0
    for idx, tag in enumerate(tags):
        values = getattr(tag, field_name)
        first = True
        for prof in range(values.shape[0]):
            mask = np.isfinite(values[prof]) & np.isfinite(tag.pressure[prof])
            if np.count_nonzero(mask) < 2:
                continue
            ax.plot(
                values[prof, mask],
                tag.pressure[prof, mask],
                color=colors[idx],
                linewidth=0.8,
                alpha=0.65,
                label=tag.smru_name if first else None,
            )
            total += 1
            first = False
    ax.set_ylim(pmax, 0)
    ax.set_xlabel(xlabel)
    ax.grid(True, alpha=0.25)
    ax.set_title(f"{field_name.upper()}: {total} profiles")


def _plot_tag_tracks(ax, tags: tuple[TagDiagnosticData, ...], *, colors: np.ndarray, title: str) -> None:
    mask_any = any(np.isfinite(tag.lon).any() and np.isfinite(tag.lat).any() for tag in tags)
    if not mask_any:
        ax.set_title(title)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.grid(True, alpha=0.25)
        ax.text(0.5, 0.5, "No valid track coordinates", ha="center", va="center", transform=ax.transAxes)
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        return

    all_lon = np.concatenate([tag.lon[np.isfinite(tag.lon)] for tag in tags if np.isfinite(tag.lon).any()])
    central_longitude = 180 if all_lon.size and float(np.nanmax(all_lon) - np.nanmin(all_lon)) > 180 else 0

    if ccrs is not None:
        projection = ccrs.PlateCarree(central_longitude=central_longitude)
        ax.remove()
        fig = ax.figure
        ax = fig.add_subplot(ax.get_subplotspec(), projection=projection)
        for idx, tag in enumerate(tags):
            mask = np.isfinite(tag.lon) & np.isfinite(tag.lat)
            if not np.any(mask):
                continue
            ax.plot(tag.lon[mask], tag.lat[mask], color=colors[idx], linewidth=0.8, alpha=0.8, transform=ccrs.PlateCarree(), label=tag.smru_name)
            ax.scatter(tag.lon[mask], tag.lat[mask], color=colors[idx], s=8, alpha=0.8, transform=ccrs.PlateCarree(), edgecolors="none")
        ax.coastlines(resolution="110m", linewidth=0.7)
        ax.add_feature(cfeature.LAND, facecolor="0.88")
        gl = ax.gridlines(draw_labels=True, linewidth=0.4, alpha=0.4)
        gl.top_labels = False
        gl.right_labels = False
        ax.set_title(title)
    else:
        lon_values: list[np.ndarray] = []
        lat_values: list[np.ndarray] = []
        for idx, tag in enumerate(tags):
            mask = np.isfinite(tag.lon) & np.isfinite(tag.lat)
            if not np.any(mask):
                continue
            lon = tag.lon[mask]
            lat = tag.lat[mask]
            lon_values.append(lon)
            lat_values.append(lat)
            ax.plot(lon, lat, color=colors[idx], linewidth=0.8, alpha=0.8, label=tag.smru_name)
            ax.scatter(lon, lat, color=colors[idx], s=12, alpha=0.8, edgecolors="none")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(title)
        ax.grid(True, alpha=0.25)
        if lon_values and lat_values:
            lon_all = np.concatenate(lon_values)
            lat_all = np.concatenate(lat_values)
            lon_pad = max(0.5, 0.05 * max(1.0, float(np.nanmax(lon_all) - np.nanmin(lon_all))))
            lat_pad = max(0.5, 0.05 * max(1.0, float(np.nanmax(lat_all) - np.nanmin(lat_all))))
            ax.set_xlim(float(np.nanmin(lon_all) - lon_pad), float(np.nanmax(lon_all) + lon_pad))
            ax.set_ylim(float(np.nanmin(lat_all) - lat_pad), float(np.nanmax(lat_all) + lat_pad))

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="upper left", fontsize=8, frameon=False)


def _plot_deployment_ts(ax, tags: tuple[TagDiagnosticData, ...], *, colors: np.ndarray) -> None:
    for idx, tag in enumerate(tags):
        first = True
        for prof in range(tag.temp.shape[0]):
            mask = np.isfinite(tag.temp[prof]) & np.isfinite(tag.psal[prof])
            if np.count_nonzero(mask) < 2:
                continue
            ax.plot(
                tag.psal[prof, mask],
                tag.temp[prof, mask],
                color=colors[idx],
                linewidth=0.8,
                alpha=0.65,
                label=tag.smru_name if first else None,
            )
            first = False
    all_temp = np.concatenate([tag.temp.ravel() for tag in tags]) if tags else np.asarray([], dtype=float)
    all_psal = np.concatenate([tag.psal.ravel() for tag in tags]) if tags else np.asarray([], dtype=float)
    ax.set_xlabel("Salinity")
    ax.set_ylabel("In-situ temperature [°C]")
    ax.grid(True, alpha=0.25)
    if np.isfinite(all_psal).any() and np.isfinite(all_temp).any():
        smin, smax = _value_range(all_psal)
        tmin, tmax = _value_range(all_temp)
        sx = np.linspace(smin, smax, 80)
        ty = np.linspace(tmin, tmax, 80)
        S, T = np.meshgrid(sx, ty)
        try:
            contour_sigma = _sigma0_profile(S.ravel(), T.ravel(), np.zeros(S.size), lon=0.0, lat=0.0).reshape(S.shape)
            ax.contour(S, T, contour_sigma, colors="0.35", linewidths=0.5, levels=8)
        except Exception:
            pass
    ax.set_title("Deployment TS envelope")
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="best", fontsize=8, frameon=False)


def _save_deployment_overview_figure(tags: tuple[TagDiagnosticData, ...], target: Path, *, qf: str, pmax: int) -> None:
    plt, _, _, GridSpec = _import_matplotlib()
    cmap = plt.get_cmap("tab20", max(len(tags), 1))
    colors = cmap(np.arange(max(len(tags), 1)))

    fig = plt.figure(figsize=(16, 10), constrained_layout=False)
    gs = GridSpec(2, 3, figure=fig, width_ratios=[1.1, 1.0, 1.0], height_ratios=[1.0, 1.0])
    ax_info = fig.add_subplot(gs[:, 0])
    ax_map = fig.add_subplot(gs[0, 1])
    ax_ts = fig.add_subplot(gs[0, 2])
    ax_t = fig.add_subplot(gs[1, 1])
    ax_s = fig.add_subplot(gs[1, 2])

    _plot_info(ax_info, _deployment_summary_lines(tags))
    _plot_tag_tracks(ax_map, tags, colors=colors, title="Deployment tracks")
    _plot_deployment_ts(ax_ts, tags, colors=colors)
    _plot_labelled_profiles(ax_t, tags, field_name="temp", colors=colors, pmax=pmax, xlabel="Temperature [°C]")
    _plot_labelled_profiles(ax_s, tags, field_name="psal", colors=colors, pmax=pmax, xlabel="Salinity")

    deployment = tags[0].deployment if tags else target.stem
    adjusted = tags[0].adjusted if tags else True
    fig.suptitle(f"{deployment} deployment diagnostics ({qf}, {'adjusted' if adjusted else 'raw'})", fontsize=16, y=0.985)
    fig.subplots_adjust(left=0.05, right=0.98, bottom=0.06, top=0.92, wspace=0.24, hspace=0.24)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_profiles(ax, values: np.ndarray, pressure: np.ndarray, colors: np.ndarray, *, label: str, pmax: int) -> None:
    valid_profiles = 0
    for idx in range(values.shape[0]):
        mask = np.isfinite(values[idx]) & np.isfinite(pressure[idx])
        if np.count_nonzero(mask) < 2:
            continue
        ax.plot(values[idx, mask], pressure[idx, mask], color=colors[idx], linewidth=0.9, alpha=0.9)
        valid_profiles += 1
    ax.set_ylim(pmax, 0)
    ax.set_title(f"{label}: {valid_profiles} profiles")
    ax.grid(True, alpha=0.25)


def _plot_ts(ax, temp: np.ndarray, psal: np.ndarray, colors: np.ndarray, sigma0: np.ndarray | None = None) -> None:
    plt, _, _, _ = _import_matplotlib()

    for idx in range(temp.shape[0]):
        mask = np.isfinite(temp[idx]) & np.isfinite(psal[idx])
        if np.count_nonzero(mask) < 2:
            continue
        ax.plot(psal[idx, mask], temp[idx, mask], color=colors[idx], linewidth=0.9, alpha=0.9)
    ax.set_xlabel("Salinity")
    ax.set_ylabel("In-situ temperature [°C]")
    ax.grid(True, alpha=0.25)
    if sigma0 is not None and np.isfinite(psal).any() and np.isfinite(temp).any():
        smin, smax = _value_range(psal)
        tmin, tmax = _value_range(temp)
        sx = np.linspace(smin, smax, 80)
        ty = np.linspace(tmin, tmax, 80)
        S, T = np.meshgrid(sx, ty)
        try:
            contour_sigma = _sigma0_profile(S.ravel(), T.ravel(), np.zeros(S.size), lon=0.0, lat=0.0).reshape(S.shape)
            ax.contour(S, T, contour_sigma, colors="0.3", linewidths=0.6, levels=8)
        except Exception:
            pass


def _plot_map(ax, lon: np.ndarray, lat: np.ndarray, colors: np.ndarray, *, title: str) -> None:
    mask = np.isfinite(lon) & np.isfinite(lat)
    if not np.any(mask):
        ax.set_title(title)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.grid(True, alpha=0.25)
        ax.text(0.5, 0.5, "No valid track coordinates", ha="center", va="center", transform=ax.transAxes)
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        return ax

    lon = lon[mask]
    lat = lat[mask]
    colors = colors[mask]

    if ccrs is not None:
        projection = ccrs.PlateCarree(central_longitude=180 if np.nanmax(lon) - np.nanmin(lon) > 180 else 0)
        ax.remove()
        fig = ax.figure
        new_ax = fig.add_subplot(ax.get_subplotspec(), projection=projection)
        new_ax.scatter(lon, lat, c=colors, s=10, transform=ccrs.PlateCarree(), edgecolors="none")
        new_ax.plot(lon, lat, color="0.5", linewidth=0.5, alpha=0.3, transform=ccrs.PlateCarree())
        new_ax.coastlines(resolution="110m", linewidth=0.7)
        new_ax.add_feature(cfeature.LAND, facecolor="0.88")
        gl = new_ax.gridlines(draw_labels=True, linewidth=0.4, alpha=0.4)
        gl.top_labels = False
        gl.right_labels = False
        new_ax.set_title(title)
        return new_ax

    ax.scatter(lon, lat, c=colors, s=14, edgecolors="none")
    ax.plot(lon, lat, color="0.5", linewidth=0.5, alpha=0.3)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    lon_pad = max(0.5, 0.05 * (np.nanmax(lon) - np.nanmin(lon) if np.isfinite(lon).any() else 1.0))
    lat_pad = max(0.5, 0.05 * (np.nanmax(lat) - np.nanmin(lat) if np.isfinite(lat).any() else 1.0))
    ax.set_xlim(np.nanmin(lon) - lon_pad, np.nanmax(lon) + lon_pad)
    ax.set_ylim(np.nanmin(lat) - lat_pad, np.nanmax(lat) + lat_pad)
    return ax


def _plot_info(ax, lines: list[str]) -> None:
    ax.axis("off")
    ax.text(0.02, 0.98, "\n".join(lines), va="top", ha="left", fontsize=11, linespacing=1.45)


def _section_panel(ax, x: np.ndarray, depth: np.ndarray, section: np.ndarray, *, cmap: str, title: str, value_label: str, pressure_obs: np.ndarray, colorbar: bool = False):
    plt, _, _, _ = _import_matplotlib()

    x_edges = _centers_to_edges(x)
    depth_edges = _centers_to_edges(depth)
    mesh_x, mesh_y = np.meshgrid(x_edges, depth_edges)
    data = section.T
    vmin, vmax = _value_range(data)
    mesh = ax.pcolormesh(mesh_x, mesh_y, data, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    if np.isfinite(data).any():
        try:
            levels = np.linspace(np.nanpercentile(data, 10), np.nanpercentile(data, 90), 6)
            if np.unique(np.round(levels, 6)).size >= 3:
                contour = ax.contour(x, depth, data, levels=levels, colors="0.15", linewidths=0.5)
                ax.clabel(contour, fmt="%.2f", fontsize=8)
        except Exception:
            pass

    obs_x = np.broadcast_to(x[:, None], pressure_obs.shape)
    mask = np.isfinite(pressure_obs)
    ax.scatter(obs_x[mask], pressure_obs[mask], s=3, color="k", alpha=0.7)
    ax.set_ylim(depth[-1], depth[0])
    ax.set_ylabel("depth [dbar]")
    ax.set_title(title)
    ax.grid(True, alpha=0.15)
    return mesh, value_label


def _save_overview_figure(dataset: xr.Dataset, target: Path, *, adjusted: bool, pmax: int) -> None:
    plt, _, ScalarMappable, GridSpec = _import_matplotlib()

    tag = _tag_diagnostic_data(dataset, target.stem, adjusted=adjusted)
    smru = tag.smru_name
    pressure = tag.pressure
    temp = tag.temp
    psal = tag.psal
    sigma0 = tag.sigma0
    lon = tag.lon
    lat = tag.lat
    times = tag.times
    colors, norm, cmap, axis_label, days = _profile_colors(times)

    fig = plt.figure(figsize=(15, 10), constrained_layout=False)
    gs = GridSpec(2, 4, figure=fig, width_ratios=[1.1, 1.0, 1.0, 1.15], height_ratios=[1.0, 1.05])

    ax_info = fig.add_subplot(gs[:, 0])
    ax_t = fig.add_subplot(gs[0, 1])
    ax_s = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[0, 3])
    ax_map = fig.add_subplot(gs[1, 1])
    ax_ts = fig.add_subplot(gs[1, 2:])

    _plot_info(ax_info, list(tag.summary_lines))
    _plot_profiles(ax_t, temp, pressure, colors, label="TEMP", pmax=pmax)
    _plot_profiles(ax_s, psal, pressure, colors, label="PSAL", pmax=pmax)
    _plot_profiles(ax_d, sigma0, pressure, colors, label="SIGMA0", pmax=pmax)
    ax_t.set_xlabel("Temperature [°C]")
    ax_s.set_xlabel("Salinity")
    ax_d.set_xlabel("sigma0 [kg m$^{-3}$]")

    new_map_ax = _plot_map(ax_map, lon, lat, colors, title="Track")
    if new_map_ax is not None:
        ax_map = new_map_ax
    _plot_ts(ax_ts, temp, psal, colors, sigma0=sigma0)
    ax_ts.set_title("TS diagram")

    color_ax = fig.add_axes([0.18, 0.93, 0.62, 0.02])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, cax=color_ax, orientation="horizontal")
    cbar.set_label(f"Profile colour key [{axis_label}]")

    fig.suptitle(f"{smru} diagnostics ({'adjusted' if adjusted else 'raw'})", fontsize=16, y=0.985)
    fig.subplots_adjust(left=0.05, right=0.98, bottom=0.06, top=0.90, wspace=0.28, hspace=0.28)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _save_section_figure(dataset: xr.Dataset, target: Path, *, adjusted: bool, pmax: int) -> None:
    plt, _, ScalarMappable, GridSpec = _import_matplotlib()

    pressure = _pressure(dataset, adjusted=adjusted)
    temp = _masked_field(dataset, "TEMP", adjusted=adjusted)
    psal = _masked_field(dataset, "PSAL", adjusted=adjusted)
    sigma0 = _compute_sigma0(dataset, adjusted=adjusted)
    times = _profile_times(dataset)
    _, norm, cmap, axis_label, x = _profile_colors(times)

    depth, temp_section = _section_grid(pressure, temp, pmax=pmax)
    _, psal_section = _section_grid(pressure, psal, pmax=pmax)
    _, sigma_section = _section_grid(pressure, sigma0, pmax=pmax)

    fig = plt.figure(figsize=(14, 10), constrained_layout=False)
    gs = GridSpec(3, 1, figure=fig, height_ratios=[1.0, 1.0, 1.0])
    axes = [fig.add_subplot(gs[idx, 0]) for idx in range(3)]

    m0, _ = _section_panel(axes[0], x, depth, temp_section, cmap="coolwarm", title="Temperature section", value_label="°C", pressure_obs=pressure)
    m1, _ = _section_panel(axes[1], x, depth, psal_section, cmap="viridis", title="Salinity section", value_label="psu", pressure_obs=pressure)
    m2, _ = _section_panel(axes[2], x, depth, sigma_section, cmap="cividis", title="Sigma0 section", value_label="kg m$^{-3}$", pressure_obs=pressure)
    axes[2].set_xlabel(axis_label)

    top_ax = fig.add_axes([0.12, 0.94, 0.76, 0.02])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, cax=top_ax, orientation="horizontal")
    cbar.set_label(f"Profile colour key [{axis_label}]")

    for ax, mesh, label in zip(axes, (m0, m1, m2), ("°C", "psu", "kg m$^{-3}$"), strict=False):
        cbar_local = fig.colorbar(mesh, ax=ax, orientation="vertical", pad=0.012, fraction=0.03)
        cbar_local.set_label(label)

    smru = _smru_name(dataset, target.stem)
    fig.suptitle(f"{smru} sections ({'adjusted' if adjusted else 'raw'})", fontsize=16, y=0.985)
    fig.subplots_adjust(left=0.07, right=0.93, bottom=0.06, top=0.90, hspace=0.24)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _source_file_paths(config: MeopConfig, selection: Selection, *, qf: str) -> tuple[Path, ...]:
    selection = selection.normalized()
    if selection.smru_name:
        path = list_fname_prof(smru_name=selection.smru_name, deployment=selection.deployment, qf=qf, config=config)
        return tuple(path)
    if not selection.deployment:
        paths: list[Path] = []
        for deployment_dir in sorted(path for path in config.final_dataset_dir.iterdir() if path.is_dir()):
            paths.extend(sorted(deployment_dir.glob(f"*_{qf}_prof.nc")))
        return tuple(paths)
    return tuple(list_fname_prof(deployment=selection.deployment, qf=qf, config=config))


def generate_diagnostics(
    config: MeopConfig,
    selection: Selection,
    *,
    qf: str = "lr1",
    adjusted: bool = True,
    pmax: int = 1000,
) -> DiagnosticResult:
    written: list[Path] = []
    processed: list[str] = []
    deployment_tags: dict[str, list[TagDiagnosticData]] = {}

    for source_path in _source_file_paths(config, selection, qf=qf):
        dataset = open_meop_netcdf(source_path)
        try:
            smru_name = source_path.name.split("_")[0]
            tag = _tag_diagnostic_data(dataset, smru_name, adjusted=adjusted)
            deployment_tags.setdefault(tag.deployment, []).append(tag)
            suffix = "adj" if adjusted else "raw"
            overview_path = fname_plots(smru_name, deployment=tag.deployment, qf=qf, suffix=f"diags_TS_{suffix}", config=config)
            section_path = fname_plots(smru_name, deployment=tag.deployment, qf=qf, suffix=f"transect_{suffix}", config=config)
            _save_overview_figure(dataset, overview_path, adjusted=adjusted, pmax=pmax)
            _save_section_figure(dataset, section_path, adjusted=adjusted, pmax=pmax)
            written.extend([overview_path, section_path])
            processed.append(smru_name)
        finally:
            dataset.close()

    overview_summaries: list[DeploymentDiagnosticSummary] = []
    if deployment_tags:
        suffix = "adj" if adjusted else "raw"
        for deployment, tags in sorted(deployment_tags.items()):
            deployment_path = config.plots_by_deployment_dir / deployment / f"{deployment}_{qf}_deployment_overview_{suffix}.png"
            _save_deployment_overview_figure(tuple(tags), deployment_path, qf=qf, pmax=pmax)
            written.append(deployment_path)
            overview_summaries.append(_deployment_summary(tuple(tags)))

    if len(overview_summaries) >= 2:
        suffix = "adj" if adjusted else "raw"
        overview_path = config.plots_overview_dir / f"all_deployments_{qf}_overview_{suffix}.png"
        _save_global_overview_figure(tuple(overview_summaries), overview_path, qf=qf)
        written.append(overview_path)

    return DiagnosticResult(written_files=tuple(written), processed_tags=tuple(processed))
