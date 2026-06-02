"""Overview map generation for MEOP published datasets.

Generates cartopy-based track maps grouped by region, deployment,
or country, saved as PNG files.  Requires ``cartopy`` (optional dep).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .regions import label_regions

try:
    import cartopy.crs as ccrs  # type: ignore
    import cartopy.feature as cfeature  # type: ignore
    _CARTOPY = True
except Exception:
    ccrs = None
    cfeature = None
    _CARTOPY = False


_COLORS_TAB20 = [
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a",
    "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
    "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d",
    "#17becf", "#9edae5",
]


def enrich_profiles_dataframe(
    df: pd.DataFrame,
    *,
    catalog_path: Path | None = None,
) -> pd.DataFrame:
    """Add REGION and (optionally) COUNTRY columns to a list-profiles DataFrame.

    Parameters
    ----------
    df:
        DataFrame with at least SMRU_PLATFORM_CODE, DEPLOYMENT_CODE,
        LATITUDE, LONGITUDE columns.
    catalog_path:
        Path to ``list_deployment.csv``.  When provided, COUNTRY is merged
        from this catalog.  When absent, COUNTRY is set to ``"Unknown"``.

    Returns
    -------
    Enriched DataFrame (copy).
    """
    df = df.copy()

    if "REGION" not in df.columns:
        # Compute per-tag median positions then broadcast back to per-profile.
        tag_col = "SMRU_PLATFORM_CODE"
        if tag_col in df.columns and "LATITUDE" in df.columns and "LONGITUDE" in df.columns:
            medians = (
                df[[tag_col, "LATITUDE", "LONGITUDE"]]
                .groupby(tag_col)
                .median()
                .reset_index()
            )
            medians["REGION"] = label_regions(
                medians["LONGITUDE"].to_numpy(),
                medians["LATITUDE"].to_numpy(),
            ).astype(str)
            df = df.merge(medians[[tag_col, "REGION"]], on=tag_col, how="left")
            df["REGION"] = df["REGION"].fillna("Unknown")
        else:
            df["REGION"] = "Unknown"

    if "COUNTRY" not in df.columns:
        df["COUNTRY"] = "Unknown"
        if catalog_path is not None and catalog_path.exists():
            try:
                cat = pd.read_csv(catalog_path, usecols=["deployment_code", "country"])
                cat = cat.rename(columns={"deployment_code": "DEPLOYMENT_CODE", "country": "COUNTRY"})
                df = df.drop(columns=["COUNTRY"]).merge(cat, on="DEPLOYMENT_CODE", how="left")
                df["COUNTRY"] = df["COUNTRY"].fillna("Unknown")
            except Exception:
                df["COUNTRY"] = "Unknown"

    return df


def _central_longitude(lons: np.ndarray) -> float:
    valid = lons[np.isfinite(lons)]
    if not valid.size:
        return 0.0
    return 180.0 if float(np.nanmax(valid) - np.nanmin(valid)) > 180 else 0.0


def _group_colors(groups: list[str]) -> dict[str, Any]:
    n = len(groups)
    try:
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap("tab20", max(n, 1))
        return {g: cmap(i % max(n, 1)) for i, g in enumerate(groups)}
    except Exception:
        return {g: _COLORS_TAB20[i % len(_COLORS_TAB20)] for i, g in enumerate(groups)}


def _make_track_map(
    df: pd.DataFrame,
    groupby: str,
    *,
    title: str,
    legend: bool = True,
    legend_horiz: bool = False,
) -> Any:
    """Create a cartopy track map; returns (fig, ax) or raises if cartopy absent."""
    import matplotlib.pyplot as plt

    lons = df["LONGITUDE"].to_numpy(dtype=float)
    lats = df["LATITUDE"].to_numpy(dtype=float)
    central_lon = _central_longitude(lons)

    groups = sorted(df[groupby].dropna().unique().tolist())
    color_map = _group_colors(groups)

    if _CARTOPY:
        projection = ccrs.PlateCarree(central_longitude=central_lon)
        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot(1, 1, 1, projection=projection)
        transform = ccrs.PlateCarree()
    else:
        fig, ax = plt.subplots(figsize=(15, 10))
        transform = None

    # Group by (groupby, SMRU_PLATFORM_CODE) to get individual tag tracks
    tag_col = "SMRU_PLATFORM_CODE" if "SMRU_PLATFORM_CODE" in df.columns else groupby
    plot_cols = [groupby, tag_col, "LONGITUDE", "LATITUDE"]
    avail_cols = [c for c in plot_cols if c in df.columns]

    labeled: set[str] = set()
    for group_val, group_df in df.groupby(groupby):
        group_val = str(group_val)
        color = color_map.get(group_val, "steelblue")
        label = group_val if group_val not in labeled else ""
        labeled.add(group_val)

        # Plot each tag track within the group
        if tag_col in group_df.columns and tag_col != groupby:
            for _, tag_df in group_df.groupby(tag_col):
                # Sort by JULD (time) to connect profiles in temporal order, not arbitrary order
                if "JULD" in tag_df.columns:
                    tag_df = tag_df.sort_values("JULD")
                
                lon = tag_df["LONGITUDE"].to_numpy(dtype=float)
                lat = tag_df["LATITUDE"].to_numpy(dtype=float)
                valid = np.isfinite(lon) & np.isfinite(lat)
                if not valid.any():
                    continue
                lon_v = lon[valid]
                lat_v = lat[valid]
                # Break tracks at large jumps (>100° gap)
                lon_plot = np.where(np.abs(np.diff(lon_v, prepend=lon_v[:1])) > 100, np.nan, lon_v)
                if _CARTOPY:
                    ax.plot(lon_plot, lat_v, color=color, linewidth=0.7, alpha=0.7,
                            transform=transform, label=label)
                else:
                    ax.plot(lon_plot, lat_v, color=color, linewidth=0.7, alpha=0.7, label=label)
                label = ""  # only first track in group gets legend label
        else:
            # Sort by JULD (time) to connect profiles in temporal order, not arbitrary order
            if "JULD" in group_df.columns:
                group_df = group_df.sort_values("JULD")
            
            lon_v = group_df["LONGITUDE"].to_numpy(dtype=float)
            lat_v = group_df["LATITUDE"].to_numpy(dtype=float)
            valid = np.isfinite(lon_v) & np.isfinite(lat_v)
            if valid.any():
                if _CARTOPY:
                    ax.plot(lon_v[valid], lat_v[valid], color=color, linewidth=0.7,
                            alpha=0.7, transform=transform, label=label)
                else:
                    ax.plot(lon_v[valid], lat_v[valid], color=color, linewidth=0.7,
                            alpha=0.7, label=label)

    if _CARTOPY:
        ax.coastlines(resolution="110m", linewidth=0.6)
        ax.add_feature(cfeature.LAND, facecolor="0.88")
        gl = ax.gridlines(draw_labels=True, linewidth=0.4, alpha=0.4)
        gl.top_labels = False
        gl.right_labels = False
    else:
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.grid(True, alpha=0.25)
        valid_mask = np.isfinite(lons) & np.isfinite(lats)
        if valid_mask.any():
            lon_pad = max(2.0, 0.05 * float(np.nanmax(lons[valid_mask]) - np.nanmin(lons[valid_mask])))
            lat_pad = max(2.0, 0.05 * float(np.nanmax(lats[valid_mask]) - np.nanmin(lats[valid_mask])))
            ax.set_xlim(float(np.nanmin(lons[valid_mask])) - lon_pad, float(np.nanmax(lons[valid_mask])) + lon_pad)
            ax.set_ylim(float(np.nanmin(lats[valid_mask])) - lat_pad, float(np.nanmax(lats[valid_mask])) + lat_pad)

    ax.set_title(title, fontsize=14)

    if legend and len(groups) > 0:
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            n_groups = len(groups)
            if legend_horiz:
                ncol = min(10, max(1, n_groups))
                offset = -0.10 if n_groups > 30 else -0.07
                ax.legend(handles, labels, bbox_to_anchor=(0, offset, 1, 0),
                          loc="upper left", mode="expand", ncol=ncol, fontsize=9,
                          frameon=False)
            else:
                ax.legend(handles, labels, bbox_to_anchor=(1.02, 1.0), loc="upper left",
                          fontsize=9, frameon=False)

    plt.tight_layout()
    return fig, ax


def _save_map(
    df: pd.DataFrame,
    groupby: str,
    title: str,
    dest: Path,
    *,
    legend: bool = True,
    legend_horiz: bool = False,
    rebuild: bool = False,
) -> Path | None:
    """Generate and save a single track map PNG; returns path on success."""
    if dest.exists() and not rebuild:
        return dest
    if df.empty or groupby not in df.columns:
        return None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, _ = _make_track_map(df, groupby, title=title, legend=legend, legend_horiz=legend_horiz)
        dest.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(dest, dpi=180, bbox_inches="tight")
        plt.close(fig)
        return dest
    except Exception:
        return None


def build_overview_maps(
    df: pd.DataFrame,
    output_dir: Path,
    *,
    rebuild: bool = False,
    verbose: bool = False,
) -> list[Path]:
    """Generate the suite of MEOP overview map PNGs.

    Parameters
    ----------
    df:
        Profile-level DataFrame already enriched with REGION and COUNTRY columns
        (use :func:`enrich_profiles_dataframe` first).
    output_dir:
        Directory where map PNGs are written.
    rebuild:
        Regenerate even if output already exists.
    verbose:
        Print each generated file.

    Returns
    -------
    List of Path objects for every map file successfully written.
    """
    written: list[Path] = []

    def _emit(path: Path | None) -> None:
        if path:
            written.append(path)
            if verbose:
                print(f"  map: {path.name}")

    # Global maps
    _emit(_save_map(df, "REGION", "Distribution of profiles by ocean region",
                    output_dir / "Global_distribution_by_region.png",
                    legend=True, rebuild=rebuild))

    _emit(_save_map(df, "DEPLOYMENT_CODE", "Distribution of profiles by deployment code",
                    output_dir / "Global_distribution_by_deployment.png",
                    legend=False, rebuild=rebuild))

    _emit(_save_map(df, "DEPLOYMENT_CODE", "Distribution of profiles by deployment code",
                    output_dir / "Global_distribution_by_deployment_legend.png",
                    legend=True, legend_horiz=True, rebuild=rebuild))

    if "COUNTRY" in df.columns and df["COUNTRY"].nunique() > 1:
        _emit(_save_map(df, "COUNTRY", "Distribution of profiles by country",
                        output_dir / "Global_distribution_by_country.png",
                        legend=True, rebuild=rebuild))

    # Per-region maps
    if "REGION" in df.columns:
        for region in sorted(df["REGION"].dropna().unique()):
            if region in ("Unknown", ""):
                continue
            region_df = df[df["REGION"] == region]
            safe = region.replace(" ", "_").replace("/", "-")
            _emit(_save_map(region_df, "DEPLOYMENT_CODE",
                            f"Distribution of {region} profiles",
                            output_dir / f"Regional_distribution_{safe}.png",
                            legend=True, legend_horiz=True, rebuild=rebuild))

    # Per-country maps
    if "COUNTRY" in df.columns:
        for country in sorted(df["COUNTRY"].dropna().unique()):
            if country in ("Unknown", ""):
                continue
            country_df = df[df["COUNTRY"] == country]
            safe = country.replace(" ", "_").replace("/", "-")
            _emit(_save_map(country_df, "DEPLOYMENT_CODE",
                            f"Distribution of profiles for {country}",
                            output_dir / f"National_distribution_{safe}.png",
                            legend=True, legend_horiz=True, rebuild=rebuild))

    return written
