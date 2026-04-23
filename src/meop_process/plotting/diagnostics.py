from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import json
from pathlib import Path
import re
from typing import Iterable

import numpy as np
import xarray as xr

from ..catalog.filenames import fname_plots, list_fname_prof
from ..io.netcdf import decode_text, juld_to_datetime, open_meop_netcdf
from ..models import MeopConfig, Selection
from ..processing.qc import _sigma0_profile, _to_numeric_qc
from .regions import label_region


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
    region: str = "Unknown"
    has_chla: bool = False
    has_doxy: bool = False
    has_hr: bool = False


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
    regions: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, object]:
        return {
            "deployment": self.deployment,
            "adjusted": self.adjusted,
            "smru_names": list(self.smru_names),
            "n_profiles": self.n_profiles,
            "n_temp_profiles": self.n_temp_profiles,
            "n_psal_profiles": self.n_psal_profiles,
            "start_time": self.start_time.isoformat() if self.start_time is not None else None,
            "end_time": self.end_time.isoformat() if self.end_time is not None else None,
            "lon": self.lon.astype(float).tolist(),
            "lat": self.lat.astype(float).tolist(),
            "regions": list(self.regions),
        }

    @classmethod
    def from_dict(cls, payload: dict[str, object]) -> "DeploymentDiagnosticSummary":
        start_raw = payload.get("start_time")
        end_raw = payload.get("end_time")
        return cls(
            deployment=str(payload.get("deployment", "")),
            adjusted=bool(payload.get("adjusted", True)),
            smru_names=tuple(str(item) for item in payload.get("smru_names", [])),
            n_profiles=int(payload.get("n_profiles", 0)),
            n_temp_profiles=int(payload.get("n_temp_profiles", 0)),
            n_psal_profiles=int(payload.get("n_psal_profiles", 0)),
            start_time=datetime.fromisoformat(str(start_raw)) if start_raw else None,
            end_time=datetime.fromisoformat(str(end_raw)) if end_raw else None,
            lon=np.asarray(payload.get("lon", []), dtype=float),
            lat=np.asarray(payload.get("lat", []), dtype=float),
            regions=tuple(str(r) for r in payload.get("regions", [])),
        )


DIAGNOSTIC_PARTS = ("tag", "deployment", "overview")


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


def _tag_diagnostic_data(dataset: xr.Dataset, source_name: str, *, adjusted: bool, config: MeopConfig | None = None) -> TagDiagnosticData:
    pressure = _pressure(dataset, adjusted=adjusted)
    temp = _masked_field(dataset, "TEMP", adjusted=adjusted)
    psal = _masked_field(dataset, "PSAL", adjusted=adjusted)
    sigma0 = _compute_sigma0(dataset, adjusted=adjusted)
    lon = np.asarray(dataset["LONGITUDE"].values, dtype=float)
    lat = np.asarray(dataset["LATITUDE"].values, dtype=float)
    times = _profile_times(dataset)
    smru_name = _smru_name(dataset, source_name)
    deployment = decode_text(dataset.attrs.get("deployment_code", "")) or source_name.split("-")[0]

    finite_mask = np.isfinite(lon) & np.isfinite(lat)
    if finite_mask.any():
        region = label_region(float(np.nanmedian(lon[finite_mask])), float(np.nanmedian(lat[finite_mask])))
    else:
        region = "Unknown"

    has_chla = "CHLA" in dataset or "CHLA_ADJUSTED" in dataset
    has_doxy = "DOXY" in dataset or "DOXY_ADJUSTED" in dataset
    has_hr = False
    if config is not None:
        from ..catalog.filenames import fname_prof as _fname_prof
        has_hr = _fname_prof(smru_name, deployment=deployment, qf="hr2", config=config).exists()

    summary = list(_summary_lines(dataset, adjusted=adjusted, sigma0=sigma0))
    if region and region != "Unknown":
        summary.append(f"region: {region}")

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
        summary_lines=tuple(summary),
        region=region,
        has_chla=has_chla,
        has_doxy=has_doxy,
        has_hr=has_hr,
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

    unique_regions = sorted({tag.region for tag in tags if tag.region and tag.region != "Unknown"})
    lines = [
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
    if unique_regions:
        lines.append(f"regions: {', '.join(unique_regions)}")
    return lines


def _deployment_summary(tags: tuple[TagDiagnosticData, ...]) -> DeploymentDiagnosticSummary:
    valid_times = [item for tag in tags for item in tag.times.ravel().tolist() if item is not None]
    lon_parts = [tag.lon[np.isfinite(tag.lon)] for tag in tags if np.isfinite(tag.lon).any()]
    lat_parts = [tag.lat[np.isfinite(tag.lat)] for tag in tags if np.isfinite(tag.lat).any()]
    lon = np.concatenate(lon_parts) if lon_parts else np.asarray([], dtype=float)
    lat = np.concatenate(lat_parts) if lat_parts else np.asarray([], dtype=float)
    unique_regions = tuple(sorted({tag.region for tag in tags if tag.region and tag.region != "Unknown"}))
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
        regions=unique_regions,
    )


def _overview_summary_lines(summaries: tuple[DeploymentDiagnosticSummary, ...], *, qf: str) -> list[str]:
    if not summaries:
        return [
            "MEOP diagnostics overview",
            "deployments: 0",
            "tags: 0",
            "profiles: 0",
            "TEMP profiles: 0",
            "PSAL profiles: 0",
            "period: unknown to unknown",
            f"product: {qf}",
            "variables: adjusted",
        ]
    total_tags = sum(len(summary.smru_names) for summary in summaries)
    total_profiles = sum(summary.n_profiles for summary in summaries)
    total_temp = sum(summary.n_temp_profiles for summary in summaries)
    total_psal = sum(summary.n_psal_profiles for summary in summaries)
    valid_starts = [summary.start_time for summary in summaries if summary.start_time is not None]
    valid_ends = [summary.end_time for summary in summaries if summary.end_time is not None]
    start = _decode_datetime(min(valid_starts)) if valid_starts else "unknown"
    end = _decode_datetime(max(valid_ends)) if valid_ends else "unknown"
    unique_regions = sorted({r for summary in summaries for r in summary.regions})
    lines = [
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
    if unique_regions:
        lines.append(f"regions: {', '.join(unique_regions)}")
    return lines


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


def _make_track_fig(figsize=(10, 7), central_longitude=0.0):
    """Create a figure+axes for a track map, using cartopy if available."""
    plt, _, _, _ = _import_matplotlib()
    if ccrs is not None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree(central_longitude=central_longitude))
        ax.coastlines(resolution="110m", linewidth=0.7)
        ax.add_feature(cfeature.LAND, facecolor="0.88")
        gl = ax.gridlines(draw_labels=True, linewidth=0.4, alpha=0.4)
        gl.top_labels = False
        gl.right_labels = False
        transform = ccrs.PlateCarree()
    else:
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.grid(True, alpha=0.25)
        transform = None
    return fig, ax, transform


def _save_section_figure(dataset: xr.Dataset, target: Path, *, adjusted: bool, pmax: int) -> None:
    """Time-section colorfill with CHLA/DOXY panels if present."""
    plt, _, ScalarMappable, GridSpec = _import_matplotlib()

    pressure = _pressure(dataset, adjusted=adjusted)
    times = _profile_times(dataset)
    _, norm, cmap, axis_label, x = _profile_colors(times)

    panels: list[tuple[np.ndarray, str, str, str]] = []
    temp = _masked_field(dataset, "TEMP", adjusted=adjusted)
    psal = _masked_field(dataset, "PSAL", adjusted=adjusted)
    sigma0 = _compute_sigma0(dataset, adjusted=adjusted)
    depth, temp_section = _section_grid(pressure, temp, pmax=pmax)
    _, psal_section = _section_grid(pressure, psal, pmax=pmax)
    _, sigma_section = _section_grid(pressure, sigma0, pmax=pmax)
    panels.append((temp_section, "coolwarm", "Temperature section", "°C"))
    panels.append((psal_section, "viridis", "Salinity section", "psu"))
    panels.append((sigma_section, "cividis", "Sigma0 section", "kg m$^{-3}$"))

    for opt_name, opt_cmap, opt_title, opt_label in (
        ("CHLA", "YlGn", "Chlorophyll-a section", "mg m$^{-3}$"),
        ("DOXY", "Blues", "Dissolved oxygen section", "µmol kg$^{-1}$"),
    ):
        if opt_name in dataset or f"{opt_name}_ADJUSTED" in dataset:
            opt_field = _masked_field(dataset, opt_name, adjusted=adjusted)
            _, opt_section = _section_grid(pressure, opt_field, pmax=pmax)
            panels.append((opt_section, opt_cmap, opt_title, opt_label))

    n_panels = len(panels)
    fig = plt.figure(figsize=(14, 3.5 * n_panels + 1), constrained_layout=False)
    gs = GridSpec(n_panels, 1, figure=fig, height_ratios=[1.0] * n_panels)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(n_panels)]

    meshes = []
    for ax, (section, pcmap, title, label) in zip(axes, panels, strict=False):
        mesh, _ = _section_panel(ax, x, depth, section, cmap=pcmap, title=title, value_label=label, pressure_obs=pressure)
        meshes.append((mesh, label))
    axes[-1].set_xlabel(axis_label)

    for ax, (mesh, label) in zip(axes, meshes, strict=False):
        cbar_local = fig.colorbar(mesh, ax=ax, orientation="vertical", pad=0.012, fraction=0.03)
        cbar_local.set_label(label)

    top_ax = fig.add_axes([0.12, 0.965, 0.76, 0.018])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, cax=top_ax, orientation="horizontal")
    cbar.set_label(f"Profile colour key [{axis_label}]")

    smru = _smru_name(dataset, target.stem)
    fig.suptitle(f"{smru} sections ({'adjusted' if adjusted else 'raw'})", fontsize=14, y=0.99)
    fig.subplots_adjust(left=0.07, right=0.93, bottom=0.06, top=0.93, hspace=0.30)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _save_ts_figure(dataset: xr.Dataset, target: Path, *, adjusted: bool) -> None:
    """TS diagram, profiles colored by time, sigma0 isopycnals."""
    plt, _, ScalarMappable, _ = _import_matplotlib()

    temp = _masked_field(dataset, "TEMP", adjusted=adjusted)
    psal = _masked_field(dataset, "PSAL", adjusted=adjusted)
    sigma0 = _compute_sigma0(dataset, adjusted=adjusted)
    times = _profile_times(dataset)
    colors, norm, cmap, axis_label, _ = _profile_colors(times)

    fig, ax = plt.subplots(figsize=(9, 7))

    for idx in range(temp.shape[0]):
        mask = np.isfinite(temp[idx]) & np.isfinite(psal[idx])
        if np.count_nonzero(mask) < 2:
            continue
        ax.plot(psal[idx, mask], temp[idx, mask], color=colors[idx], linewidth=0.9, alpha=0.9)

    ax.set_xlabel("Salinity")
    ax.set_ylabel("In-situ temperature [°C]")
    ax.grid(True, alpha=0.25)

    if np.isfinite(psal).any() and np.isfinite(temp).any():
        smin, smax = _value_range(psal)
        tmin, tmax = _value_range(temp)
        sx = np.linspace(smin, smax, 80)
        ty = np.linspace(tmin, tmax, 80)
        S, T = np.meshgrid(sx, ty)
        try:
            contour_sigma = _sigma0_profile(S.ravel(), T.ravel(), np.zeros(S.size), lon=0.0, lat=0.0).reshape(S.shape)
            cs = ax.contour(S, T, contour_sigma, colors="0.3", linewidths=0.6, levels=8)
            ax.clabel(cs, fmt="%.1f", fontsize=7)
        except Exception:
            pass

    smru = _smru_name(dataset, target.stem)
    ax.set_title(f"{smru} TS diagram ({'adjusted' if adjusted else 'raw'})")

    color_ax = fig.add_axes([0.12, 0.93, 0.76, 0.02])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, cax=color_ax, orientation="horizontal")
    cbar.set_label(f"Profile colour [{axis_label}]")

    fig.subplots_adjust(left=0.10, right=0.95, bottom=0.09, top=0.88)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _save_map_figure(dataset: xr.Dataset, target: Path, *, adjusted: bool) -> None:
    """Single-tag track map, positions colored by time."""
    plt, _, ScalarMappable, _ = _import_matplotlib()

    lon = np.asarray(dataset["LONGITUDE"].values, dtype=float)
    lat = np.asarray(dataset["LATITUDE"].values, dtype=float)
    times = _profile_times(dataset)
    colors, norm, cmap, axis_label, _ = _profile_colors(times)

    mask = np.isfinite(lon) & np.isfinite(lat)
    lon_valid = lon[mask]
    lat_valid = lat[mask]
    central_longitude = 180.0 if (lon_valid.size > 0 and float(np.nanmax(lon_valid) - np.nanmin(lon_valid)) > 180) else 0.0

    fig, ax, transform = _make_track_fig(figsize=(10, 7), central_longitude=central_longitude)

    if not np.any(mask):
        ax.text(0.5, 0.5, "No valid track coordinates", ha="center", va="center", transform=ax.transAxes)
    elif transform is not None:
        ax.plot(lon_valid, lat_valid, color="0.5", linewidth=0.5, alpha=0.5, transform=transform)
        ax.scatter(lon_valid, lat_valid, c=colors[mask], s=14, transform=transform, edgecolors="none")
    else:
        ax.plot(lon_valid, lat_valid, color="0.5", linewidth=0.5, alpha=0.5)
        ax.scatter(lon_valid, lat_valid, c=colors[mask], s=14, edgecolors="none")
        lon_pad = max(0.5, 0.05 * max(1.0, float(np.nanmax(lon_valid) - np.nanmin(lon_valid))))
        lat_pad = max(0.5, 0.05 * max(1.0, float(np.nanmax(lat_valid) - np.nanmin(lat_valid))))
        ax.set_xlim(float(np.nanmin(lon_valid) - lon_pad), float(np.nanmax(lon_valid) + lon_pad))
        ax.set_ylim(float(np.nanmin(lat_valid) - lat_pad), float(np.nanmax(lat_valid) + lat_pad))

    smru = _smru_name(dataset, target.stem)
    ax.set_title(f"{smru} track ({'adjusted' if adjusted else 'raw'})")

    color_ax = fig.add_axes([0.12, 0.03, 0.76, 0.025])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, cax=color_ax, orientation="horizontal")
    cbar.set_label(f"Profile colour [{axis_label}]")

    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _save_profiles_figure(dataset: xr.Dataset, target: Path, *, adjusted: bool, pmax: int) -> None:
    """Vertical profiles: TEMP, PSAL, SIG0 side by side, colored by time."""
    plt, _, ScalarMappable, _ = _import_matplotlib()

    pressure = _pressure(dataset, adjusted=adjusted)
    temp = _masked_field(dataset, "TEMP", adjusted=adjusted)
    psal = _masked_field(dataset, "PSAL", adjusted=adjusted)
    sigma0 = _compute_sigma0(dataset, adjusted=adjusted)
    times = _profile_times(dataset)
    colors, norm, cmap, axis_label, _ = _profile_colors(times)

    fig, axes = plt.subplots(1, 3, figsize=(13, 7), sharey=True)

    for ax, values, xlabel, title in zip(
        axes,
        (temp, psal, sigma0),
        ("Temperature [°C]", "Salinity", "sigma0 [kg m$^{-3}$]"),
        ("TEMP", "PSAL", "SIGMA0"),
        strict=False,
    ):
        valid = 0
        for idx in range(values.shape[0]):
            mask = np.isfinite(values[idx]) & np.isfinite(pressure[idx])
            if np.count_nonzero(mask) < 2:
                continue
            ax.plot(values[idx, mask], pressure[idx, mask], color=colors[idx], linewidth=0.9, alpha=0.9)
            valid += 1
        ax.set_ylim(pmax, 0)
        ax.set_xlabel(xlabel)
        ax.set_title(f"{title}: {valid} profiles")
        ax.grid(True, alpha=0.25)

    axes[0].set_ylabel("Pressure [dbar]")

    color_ax = fig.add_axes([0.12, 0.93, 0.76, 0.02])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, cax=color_ax, orientation="horizontal")
    cbar.set_label(f"Profile colour [{axis_label}]")

    smru = _smru_name(dataset, target.stem)
    fig.suptitle(f"{smru} profiles ({'adjusted' if adjusted else 'raw'})", fontsize=14, y=0.99)
    fig.subplots_adjust(left=0.07, right=0.97, bottom=0.09, top=0.88, wspace=0.12)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _save_flags_figure(dataset: xr.Dataset, target: Path, *, adjusted: bool) -> None:
    """QC flag overview: raw in grey/red, adjusted good in steelblue."""
    plt, _, _, _ = _import_matplotlib()

    pressure_raw = _field(dataset, "PRES", adjusted=False)
    temp_raw = _field(dataset, "TEMP", adjusted=False)
    psal_raw = _field(dataset, "PSAL", adjusted=False)
    temp_qc = _to_numeric_qc(dataset["TEMP_QC"].values if "TEMP_QC" in dataset else np.ones(temp_raw.shape, dtype=float))
    psal_qc = _to_numeric_qc(dataset["PSAL_QC"].values if "PSAL_QC" in dataset else np.ones(psal_raw.shape, dtype=float))

    fig, axes = plt.subplots(1, 2, figsize=(13, 7), sharey=True)

    pmax = int(np.nanmax(pressure_raw[np.isfinite(pressure_raw)])) if np.isfinite(pressure_raw).any() else 1000
    pmax = min(pmax + 50, 2000)

    for ax, raw_vals, raw_qc, varname in zip(
        axes, (temp_raw, psal_raw), (temp_qc, psal_qc), ("TEMP", "PSAL"), strict=False
    ):
        n_levels = raw_vals.shape[1] if raw_vals.ndim > 1 else 1
        for idx in range(raw_vals.shape[0]):
            mask = np.isfinite(raw_vals[idx]) & np.isfinite(pressure_raw[idx])
            if np.count_nonzero(mask) < 2:
                continue
            bad_frac = np.mean(raw_qc[idx] >= 3) if raw_qc.ndim > 1 else float(raw_qc[idx] >= 3)
            color = "tomato" if bad_frac > 0.30 else "0.7"
            ax.plot(raw_vals[idx, mask], pressure_raw[idx, mask], color=color, linewidth=0.6, alpha=0.55)

        if adjusted:
            temp_adj = _masked_field(dataset, varname, adjusted=True)
            pressure_adj = _pressure(dataset, adjusted=True)
            for idx in range(temp_adj.shape[0]):
                mask = np.isfinite(temp_adj[idx]) & np.isfinite(pressure_adj[idx])
                if np.count_nonzero(mask) < 2:
                    continue
                ax.plot(temp_adj[idx, mask], pressure_adj[idx, mask], color="steelblue", linewidth=0.7, alpha=0.7)

        ax.set_ylim(pmax, 0)
        ax.set_xlabel(varname)
        ax.set_title(f"{varname} QC flags")
        ax.grid(True, alpha=0.25)

    axes[0].set_ylabel("Pressure [dbar]")

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color="0.7", linewidth=1.5, label="raw accepted"),
        Line2D([0], [0], color="tomato", linewidth=1.5, label="raw rejected (>30% bad)"),
    ]
    if adjusted:
        legend_elements.append(Line2D([0], [0], color="steelblue", linewidth=1.5, label="adjusted good"))
    axes[0].legend(handles=legend_elements, loc="lower left", fontsize=9, frameon=True)

    smru = _smru_name(dataset, target.stem)
    fig.suptitle(f"{smru} QC flags ({'adjusted' if adjusted else 'raw'})", fontsize=14)
    fig.subplots_adjust(left=0.07, right=0.97, bottom=0.09, top=0.92, wspace=0.12)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _save_info_text(dataset: xr.Dataset, target: Path, *, adjusted: bool, has_hr: bool = False) -> None:
    """Write plain-text metadata."""
    sigma0 = _compute_sigma0(dataset, adjusted=adjusted)
    lines = _summary_lines(dataset, adjusted=adjusted, sigma0=sigma0)

    smru = _smru_name(dataset, target.stem)
    times = _profile_times(dataset)
    valid_times = [t for t in times.ravel().tolist() if t is not None]

    finite_mask = np.isfinite(np.asarray(dataset["LONGITUDE"].values, dtype=float)) & \
                  np.isfinite(np.asarray(dataset["LATITUDE"].values, dtype=float))
    lon_arr = np.asarray(dataset["LONGITUDE"].values, dtype=float)
    lat_arr = np.asarray(dataset["LATITUDE"].values, dtype=float)
    region = "Unknown"
    if finite_mask.any():
        region = label_region(float(np.nanmedian(lon_arr[finite_mask])), float(np.nanmedian(lat_arr[finite_mask])))

    has_chla = "CHLA" in dataset or "CHLA_ADJUSTED" in dataset
    has_doxy = "DOXY" in dataset or "DOXY_ADJUSTED" in dataset
    sensors = ["TEMP", "PSAL"]
    if has_chla:
        sensors.append("CHLA")
    if has_doxy:
        sensors.append("DOXY")

    output_lines = list(lines)
    if region and region != "Unknown":
        output_lines.append(f"region: {region}")
    output_lines.append(f"sensors: {'/'.join(sensors)}")
    if has_hr:
        output_lines.append("HR data: yes")
    output_lines.append("")
    output_lines.append("--- calibration ---")
    output_lines.extend(_extract_adjustment_lines(dataset))

    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(output_lines) + "\n", encoding="utf-8")


def _info_text_path(smru_name: str, deployment: str, qf: str, *, suffix: str, config: MeopConfig) -> Path:
    base = fname_plots(smru_name, deployment=deployment, qf=qf, suffix=suffix, config=config)
    return base.with_suffix(".txt")


def _save_deployment_map_figure(tags: tuple[TagDiagnosticData, ...], target: Path, *, qf: str) -> None:
    """Track map of all tags, one color per tag (tab20 colormap)."""
    plt, _, _, _ = _import_matplotlib()

    cmap = plt.get_cmap("tab20", max(len(tags), 1))
    colors = cmap(np.arange(max(len(tags), 1)))

    lon_parts = [tag.lon[np.isfinite(tag.lon)] for tag in tags if np.isfinite(tag.lon).any()]
    all_lon = np.concatenate(lon_parts) if lon_parts else np.array([])
    central_longitude = 180.0 if (all_lon.size > 0 and float(np.nanmax(all_lon) - np.nanmin(all_lon)) > 180) else 0.0

    fig, ax, transform = _make_track_fig(figsize=(10, 7), central_longitude=central_longitude)

    for idx, tag in enumerate(tags):
        mask = np.isfinite(tag.lon) & np.isfinite(tag.lat)
        if not np.any(mask):
            continue
        lon_v = tag.lon[mask]
        lat_v = tag.lat[mask]
        if transform is not None:
            ax.plot(lon_v, lat_v, color=colors[idx], linewidth=0.9, alpha=0.8, transform=transform, label=tag.smru_name)
            ax.scatter(lon_v, lat_v, color=colors[idx], s=8, alpha=0.8, transform=transform, edgecolors="none")
        else:
            ax.plot(lon_v, lat_v, color=colors[idx], linewidth=0.9, alpha=0.8, label=tag.smru_name)
            ax.scatter(lon_v, lat_v, color=colors[idx], s=8, alpha=0.8, edgecolors="none")

    if not tags or not any(np.isfinite(tag.lon).any() and np.isfinite(tag.lat).any() for tag in tags):
        ax.text(0.5, 0.5, "No valid track coordinates", ha="center", va="center", transform=ax.transAxes)
    elif transform is None and all_lon.size:
        all_lat = np.concatenate([tag.lat[np.isfinite(tag.lat)] for tag in tags if np.isfinite(tag.lat).any()])
        lon_pad = max(0.5, 0.05 * max(1.0, float(np.nanmax(all_lon) - np.nanmin(all_lon))))
        lat_pad = max(0.5, 0.05 * max(1.0, float(np.nanmax(all_lat) - np.nanmin(all_lat))))
        ax.set_xlim(float(np.nanmin(all_lon) - lon_pad), float(np.nanmax(all_lon) + lon_pad))
        ax.set_ylim(float(np.nanmin(all_lat) - lat_pad), float(np.nanmax(all_lat) + lat_pad))

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="upper left", fontsize=8, frameon=False)

    deployment = tags[0].deployment if tags else ""
    adjusted = tags[0].adjusted if tags else True
    ax.set_title(f"{deployment} tracks ({qf}, {'adjusted' if adjusted else 'raw'})")

    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _save_deployment_ts_figure(tags: tuple[TagDiagnosticData, ...], target: Path, *, qf: str) -> None:
    """TS diagram for all tags, one color per tag, sigma0 isopycnals."""
    plt, _, _, _ = _import_matplotlib()

    cmap = plt.get_cmap("tab20", max(len(tags), 1))
    colors = cmap(np.arange(max(len(tags), 1)))

    fig, ax = plt.subplots(figsize=(9, 7))

    for idx, tag in enumerate(tags):
        first = True
        for prof in range(tag.temp.shape[0]):
            mask = np.isfinite(tag.temp[prof]) & np.isfinite(tag.psal[prof])
            if np.count_nonzero(mask) < 2:
                continue
            ax.plot(
                tag.psal[prof, mask], tag.temp[prof, mask],
                color=colors[idx], linewidth=0.8, alpha=0.65,
                label=tag.smru_name if first else None,
            )
            first = False

    all_temp = np.concatenate([tag.temp.ravel() for tag in tags]) if tags else np.asarray([], dtype=float)
    all_psal = np.concatenate([tag.psal.ravel() for tag in tags]) if tags else np.asarray([], dtype=float)

    if np.isfinite(all_psal).any() and np.isfinite(all_temp).any():
        smin, smax = _value_range(all_psal)
        tmin, tmax = _value_range(all_temp)
        sx = np.linspace(smin, smax, 80)
        ty = np.linspace(tmin, tmax, 80)
        S, T = np.meshgrid(sx, ty)
        try:
            contour_sigma = _sigma0_profile(S.ravel(), T.ravel(), np.zeros(S.size), lon=0.0, lat=0.0).reshape(S.shape)
            cs = ax.contour(S, T, contour_sigma, colors="0.35", linewidths=0.5, levels=8)
            ax.clabel(cs, fmt="%.1f", fontsize=7)
        except Exception:
            pass

    ax.set_xlabel("Salinity")
    ax.set_ylabel("In-situ temperature [°C]")
    ax.grid(True, alpha=0.25)

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="best", fontsize=8, frameon=False)

    deployment = tags[0].deployment if tags else ""
    adjusted = tags[0].adjusted if tags else True
    ax.set_title(f"{deployment} TS diagram ({qf}, {'adjusted' if adjusted else 'raw'})")

    fig.tight_layout()
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _save_deployment_histograms_figure(tags: tuple[TagDiagnosticData, ...], target: Path, *, qf: str) -> None:
    """3-panel: max pressure histogram, profiles-per-tag bar, profiles-per-month bar."""
    plt, _, _, _ = _import_matplotlib()

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    # Panel 1: max pressure histogram
    max_pres = []
    for tag in tags:
        for idx in range(tag.pressure.shape[0]):
            p = tag.pressure[idx]
            finite_p = p[np.isfinite(p)]
            if finite_p.size > 0:
                max_pres.append(float(np.max(finite_p)))
    if max_pres:
        axes[0].hist(max_pres, bins=30, color="steelblue", edgecolor="white", linewidth=0.5)
    axes[0].set_xlabel("Max pressure [dbar]")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Max pressure distribution")
    axes[0].grid(True, axis="y", alpha=0.25)

    # Panel 2: profiles per tag bar
    tag_names = [tag.smru_name for tag in tags]
    n_profs = [int(tag.pressure.shape[0]) for tag in tags]
    x = np.arange(len(tag_names))
    axes[1].bar(x, n_profs, color="steelblue", edgecolor="white", linewidth=0.5)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(tag_names, rotation=30, ha="right", fontsize=8)
    axes[1].set_ylabel("Profiles")
    axes[1].set_title("Profiles per tag")
    axes[1].grid(True, axis="y", alpha=0.25)

    # Panel 3: profiles per month
    from collections import Counter
    month_counts: Counter = Counter()
    for tag in tags:
        for t in tag.times.ravel().tolist():
            if t is not None:
                month_counts[t.strftime("%Y-%m")] += 1
    if month_counts:
        sorted_months = sorted(month_counts.keys())
        counts = [month_counts[m] for m in sorted_months]
        x_m = np.arange(len(sorted_months))
        axes[2].bar(x_m, counts, color="steelblue", edgecolor="white", linewidth=0.5)
        axes[2].set_xticks(x_m[::max(1, len(sorted_months) // 8)])
        axes[2].set_xticklabels(sorted_months[::max(1, len(sorted_months) // 8)], rotation=30, ha="right", fontsize=8)
    axes[2].set_ylabel("Profiles")
    axes[2].set_title("Profiles per month")
    axes[2].grid(True, axis="y", alpha=0.25)

    deployment = tags[0].deployment if tags else ""
    adjusted = tags[0].adjusted if tags else True
    fig.suptitle(f"{deployment} histograms ({qf}, {'adjusted' if adjusted else 'raw'})", fontsize=13)
    fig.tight_layout()
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _save_deployment_timing_figure(tags: tuple[TagDiagnosticData, ...], target: Path, *, qf: str) -> None:
    """Gantt chart: horizontal bar per tag with date range."""
    import matplotlib.dates as mdates
    plt, _, _, _ = _import_matplotlib()

    n = len(tags)
    fig, ax = plt.subplots(figsize=(11, max(3.5, 0.45 * n + 1.5)))

    for idx, tag in enumerate(tags):
        valid_times = [t for t in tag.times.ravel().tolist() if t is not None]
        if len(valid_times) < 2:
            continue
        start = min(valid_times)
        end = max(valid_times)
        ax.barh(idx, (end - start).total_seconds() / 86400.0, left=mdates.date2num(start), height=0.6, color="steelblue", alpha=0.8)
        ax.text(mdates.date2num(end) + 0.5, idx, tag.smru_name, va="center", fontsize=8)

    ax.set_yticks(range(n))
    ax.set_yticklabels([tag.smru_name for tag in tags], fontsize=8)
    ax.xaxis_date()
    fig.autofmt_xdate()
    ax.set_xlabel("Date")
    ax.set_title(f"{tags[0].deployment if tags else ''} timing ({qf})")
    ax.grid(True, axis="x", alpha=0.25)
    ax.invert_yaxis()

    fig.tight_layout()
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _save_overview_map_figure(summaries: tuple[DeploymentDiagnosticSummary, ...], target: Path, *, qf: str) -> None:
    """Global track map, one color per deployment."""
    plt, _, _, _ = _import_matplotlib()

    cmap = plt.get_cmap("tab20", max(len(summaries), 1))
    colors = cmap(np.arange(max(len(summaries), 1)))

    all_lon = np.concatenate([s.lon for s in summaries if s.lon.size]) if summaries else np.array([])
    central_longitude = 180.0 if (all_lon.size > 0 and float(np.nanmax(all_lon) - np.nanmin(all_lon)) > 180) else 0.0

    fig, ax, transform = _make_track_fig(figsize=(12, 7), central_longitude=central_longitude)

    for idx, summary in enumerate(summaries):
        if not summary.lon.size or not summary.lat.size:
            continue
        if transform is not None:
            ax.plot(summary.lon, summary.lat, color=colors[idx], linewidth=0.9, alpha=0.8, transform=transform, label=summary.deployment)
            ax.scatter(summary.lon, summary.lat, color=colors[idx], s=6, alpha=0.6, transform=transform, edgecolors="none")
        else:
            ax.plot(summary.lon, summary.lat, color=colors[idx], linewidth=0.9, alpha=0.8, label=summary.deployment)
            ax.scatter(summary.lon, summary.lat, color=colors[idx], s=6, alpha=0.6, edgecolors="none")

    if not summaries or not any(s.lon.size and s.lat.size for s in summaries):
        ax.text(0.5, 0.5, "No valid track coordinates", ha="center", va="center", transform=ax.transAxes)
    elif transform is None and all_lon.size:
        all_lat = np.concatenate([s.lat for s in summaries if s.lat.size])
        lon_pad = max(0.5, 0.05 * max(1.0, float(np.nanmax(all_lon) - np.nanmin(all_lon))))
        lat_pad = max(0.5, 0.05 * max(1.0, float(np.nanmax(all_lat) - np.nanmin(all_lat))))
        ax.set_xlim(float(np.nanmin(all_lon) - lon_pad), float(np.nanmax(all_lon) + lon_pad))
        ax.set_ylim(float(np.nanmin(all_lat) - lat_pad), float(np.nanmax(all_lat) + lat_pad))

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="upper left", fontsize=8, frameon=False)

    ax.set_title(f"MEOP deployment tracks ({qf})")

    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _save_overview_histograms_figure(summaries: tuple[DeploymentDiagnosticSummary, ...], target: Path, *, qf: str) -> None:
    """2-panel: profiles per deployment bar chart, profiles per year bar chart."""
    plt, _, _, _ = _import_matplotlib()

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel 1: profiles per deployment
    dep_names = [s.deployment for s in summaries]
    prof_counts = [s.n_profiles for s in summaries]
    x = np.arange(len(dep_names))
    axes[0].bar(x, prof_counts, color="steelblue", edgecolor="white", linewidth=0.5)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(dep_names, rotation=30, ha="right", fontsize=8)
    axes[0].set_ylabel("Profiles")
    axes[0].set_title("Profiles per deployment")
    axes[0].grid(True, axis="y", alpha=0.25)

    # Panel 2: profiles per year
    from collections import Counter
    year_counts: Counter = Counter()
    for s in summaries:
        if s.start_time is not None and s.end_time is not None:
            for yr in range(s.start_time.year, s.end_time.year + 1):
                year_counts[yr] += s.n_profiles // max(1, s.end_time.year - s.start_time.year + 1)
    if year_counts:
        sorted_years = sorted(year_counts.keys())
        counts = [year_counts[y] for y in sorted_years]
        x_y = np.arange(len(sorted_years))
        axes[1].bar(x_y, counts, color="steelblue", edgecolor="white", linewidth=0.5)
        axes[1].set_xticks(x_y)
        axes[1].set_xticklabels([str(y) for y in sorted_years], rotation=30, ha="right")
    axes[1].set_ylabel("Profiles (approx)")
    axes[1].set_title("Profiles per year")
    axes[1].grid(True, axis="y", alpha=0.25)

    fig.suptitle(f"MEOP overview histograms ({qf})", fontsize=13)
    fig.tight_layout()
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _save_overview_timing_figure(summaries: tuple[DeploymentDiagnosticSummary, ...], target: Path, *, qf: str) -> None:
    """Gantt chart: one bar per deployment."""
    import matplotlib.dates as mdates
    plt, _, _, _ = _import_matplotlib()

    n = len(summaries)
    fig, ax = plt.subplots(figsize=(12, max(4, 0.35 * n + 2)))

    for idx, summary in enumerate(summaries):
        if summary.start_time is None or summary.end_time is None:
            continue
        span = (summary.end_time - summary.start_time).total_seconds() / 86400.0
        ax.barh(idx, span, left=mdates.date2num(summary.start_time), height=0.6, color="steelblue", alpha=0.8)
        ax.text(mdates.date2num(summary.end_time) + 0.5, idx, summary.deployment, va="center", fontsize=8)

    ax.set_yticks(range(n))
    ax.set_yticklabels([s.deployment for s in summaries], fontsize=8)
    ax.xaxis_date()
    fig.autofmt_xdate()
    ax.set_xlabel("Date")
    ax.set_title(f"MEOP deployment timeline ({qf})")
    ax.grid(True, axis="x", alpha=0.25)
    ax.invert_yaxis()

    fig.tight_layout()
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _normalize_parts(parts: Iterable[str] | None) -> tuple[str, ...]:
    if parts is None:
        return DIAGNOSTIC_PARTS
    normalized: list[str] = []
    for item in parts:
        value = str(item).strip().lower()
        if not value:
            continue
        if value == "all":
            return DIAGNOSTIC_PARTS
        if value not in DIAGNOSTIC_PARTS:
            raise ValueError(f"Unsupported diagnostics part: {item}")
        if value not in normalized:
            normalized.append(value)
    return tuple(normalized) or DIAGNOSTIC_PARTS


def _summary_suffix(*, adjusted: bool) -> str:
    return "adj" if adjusted else "raw"


def _deployment_summary_cache_path(config: MeopConfig, deployment: str, *, qf: str, adjusted: bool) -> Path:
    suffix = _summary_suffix(adjusted=adjusted)
    return config.plots_by_deployment_dir / deployment / f"{deployment}_{qf}_deployment_summary_{suffix}.json"


def _write_deployment_summary_cache(config: MeopConfig, summary: DeploymentDiagnosticSummary, *, qf: str) -> Path:
    path = _deployment_summary_cache_path(config, summary.deployment, qf=qf, adjusted=summary.adjusted)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary.as_dict(), indent=2, sort_keys=True), encoding="utf-8")
    return path


def _read_deployment_summary_cache(config: MeopConfig, deployment: str, *, qf: str, adjusted: bool) -> DeploymentDiagnosticSummary | None:
    path = _deployment_summary_cache_path(config, deployment, qf=qf, adjusted=adjusted)
    if not path.exists():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        return None
    summary = DeploymentDiagnosticSummary.from_dict(payload)
    if summary.deployment != deployment:
        return None
    return summary


def _selected_deployments(config: MeopConfig, selection: Selection, *, qf: str) -> tuple[str, ...]:
    selection = selection.normalized()
    if selection.deployment:
        return (selection.deployment,)
    if selection.smru_name:
        return (selection.smru_name.split("-")[0],)
    if not config.final_dataset_dir.exists():
        return ()
    discovered: list[str] = []
    for deployment_dir in sorted(path for path in config.final_dataset_dir.iterdir() if path.is_dir()):
        if any(deployment_dir.glob(f"*_{qf}_prof.nc")):
            discovered.append(deployment_dir.name)
    return tuple(discovered)


def _source_paths_for_deployments(config: MeopConfig, deployments: Iterable[str], *, qf: str) -> tuple[Path, ...]:
    paths: list[Path] = []
    for deployment in deployments:
        paths.extend(list_fname_prof(deployment=deployment, qf=qf, config=config))
    return tuple(paths)


def _source_file_paths(config: MeopConfig, selection: Selection, *, qf: str) -> tuple[Path, ...]:
    selection = selection.normalized()
    if selection.smru_name:
        path = list_fname_prof(smru_name=selection.smru_name, deployment=selection.deployment, qf=qf, config=config)
        return tuple(path)
    if not selection.deployment:
        if not config.final_dataset_dir.exists():
            return ()
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
    parts: Iterable[str] | None = None,
    use_cached_summaries: bool = True,
) -> DiagnosticResult:
    normalized_parts = _normalize_parts(parts)
    part_set = set(normalized_parts)
    written: list[Path] = []
    processed: list[str] = []
    deployment_tags: dict[str, list[TagDiagnosticData]] = {}
    selected_deployments = _selected_deployments(config, selection, qf=qf)

    need_source_data = bool(part_set.intersection({"tag", "deployment"}))
    overview_summaries: dict[str, DeploymentDiagnosticSummary] = {}
    missing_overview: list[str] = []

    if "overview" in part_set:
        for deployment in selected_deployments:
            summary = _read_deployment_summary_cache(config, deployment, qf=qf, adjusted=adjusted) if use_cached_summaries else None
            if summary is None:
                missing_overview.append(deployment)
            else:
                overview_summaries[deployment] = summary
        if missing_overview:
            need_source_data = True

    if need_source_data:
        if selection.smru_name:
            source_paths = _source_file_paths(config, selection, qf=qf)
        elif selection.deployment:
            source_paths = _source_file_paths(config, selection, qf=qf)
        elif missing_overview and not part_set.intersection({"tag", "deployment"}):
            source_paths = _source_paths_for_deployments(config, missing_overview, qf=qf)
        else:
            source_paths = _source_file_paths(config, selection, qf=qf)
    else:
        source_paths = ()

    suffix = _summary_suffix(adjusted=adjusted)

    for source_path in source_paths:
        dataset = open_meop_netcdf(source_path)
        try:
            smru_name = source_path.name.split("_")[0]
            tag = _tag_diagnostic_data(dataset, smru_name, adjusted=adjusted, config=config)
            deployment_tags.setdefault(tag.deployment, []).append(tag)
            if "tag" in part_set:
                section_path = fname_plots(smru_name, deployment=tag.deployment, qf=qf, suffix=f"section_{suffix}", config=config)
                ts_path = fname_plots(smru_name, deployment=tag.deployment, qf=qf, suffix=f"TS_{suffix}", config=config)
                map_path = fname_plots(smru_name, deployment=tag.deployment, qf=qf, suffix=f"map_{suffix}", config=config)
                profiles_path = fname_plots(smru_name, deployment=tag.deployment, qf=qf, suffix=f"profiles_{suffix}", config=config)
                flags_path = fname_plots(smru_name, deployment=tag.deployment, qf=qf, suffix=f"flags_{suffix}", config=config)
                info_path = _info_text_path(smru_name, tag.deployment, qf, suffix=f"info_{suffix}", config=config)
                _save_section_figure(dataset, section_path, adjusted=adjusted, pmax=pmax)
                _save_ts_figure(dataset, ts_path, adjusted=adjusted)
                _save_map_figure(dataset, map_path, adjusted=adjusted)
                _save_profiles_figure(dataset, profiles_path, adjusted=adjusted, pmax=pmax)
                _save_flags_figure(dataset, flags_path, adjusted=adjusted)
                _save_info_text(dataset, info_path, adjusted=adjusted, has_hr=tag.has_hr)
                written.extend([section_path, ts_path, map_path, profiles_path, flags_path, info_path])
            processed.append(smru_name)
        finally:
            dataset.close()

    if deployment_tags:
        for deployment, tags in sorted(deployment_tags.items()):
            tag_tuple = tuple(tags)
            summary = _deployment_summary(tag_tuple)
            overview_summaries[deployment] = summary
            cache_path = _write_deployment_summary_cache(config, summary, qf=qf)
            written.append(cache_path)
            if "deployment" in part_set:
                dep_map_path = config.plots_by_deployment_dir / deployment / f"{deployment}_{qf}_map_{suffix}.png"
                dep_ts_path = config.plots_by_deployment_dir / deployment / f"{deployment}_{qf}_TS_{suffix}.png"
                dep_hist_path = config.plots_by_deployment_dir / deployment / f"{deployment}_{qf}_histograms_{suffix}.png"
                dep_timing_path = config.plots_by_deployment_dir / deployment / f"{deployment}_{qf}_timing_{suffix}.png"
                _save_deployment_map_figure(tag_tuple, dep_map_path, qf=qf)
                _save_deployment_ts_figure(tag_tuple, dep_ts_path, qf=qf)
                _save_deployment_histograms_figure(tag_tuple, dep_hist_path, qf=qf)
                _save_deployment_timing_figure(tag_tuple, dep_timing_path, qf=qf)
                written.extend([dep_map_path, dep_ts_path, dep_hist_path, dep_timing_path])

    if "overview" in part_set and len(overview_summaries) >= 2:
        summaries = tuple(summary for _, summary in sorted(overview_summaries.items()))
        overview_map = config.plots_overview_dir / f"all_deployments_{qf}_map_{suffix}.png"
        overview_hist = config.plots_overview_dir / f"all_deployments_{qf}_histograms_{suffix}.png"
        overview_timing = config.plots_overview_dir / f"all_deployments_{qf}_timing_{suffix}.png"
        _save_overview_map_figure(summaries, overview_map, qf=qf)
        _save_overview_histograms_figure(summaries, overview_hist, qf=qf)
        _save_overview_timing_figure(summaries, overview_timing, qf=qf)
        written.extend([overview_map, overview_hist, overview_timing])

    return DiagnosticResult(written_files=tuple(written), processed_tags=tuple(processed))
