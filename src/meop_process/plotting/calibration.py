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
``config.plotdir / deployment / {smru_name}_lr1_calibration_N.png``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np


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
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (pres_grid, temp, psal, juld) from a MEOP lr0/lr1 prof file.

    All output arrays are float64.  ``pres_grid`` shape is (N_LEVELS,),
    ``temp`` / ``psal`` shape is (N_PROF, N_LEVELS), ``juld`` shape is (N_PROF,).
    """
    with _open_any(path) as ds:
        pres = np.asarray(
            ds["PRES_ADJUSTED"].values
            if "PRES_ADJUSTED" in ds
            else ds["PRES"].values,
            dtype=np.float64,
        )
        temp = np.asarray(
            ds["TEMP_ADJUSTED"].values
            if "TEMP_ADJUSTED" in ds
            else ds["TEMP"].values,
            dtype=np.float64,
        )
        psal = np.asarray(
            ds["PSAL_ADJUSTED"].values
            if "PSAL_ADJUSTED" in ds
            else ds["PSAL"].values,
            dtype=np.float64,
        )
        juld = np.asarray(ds["JULD"].values, dtype=np.float64)
        # Replace fill values (≥9999) with NaN
        for arr in (pres, temp, psal):
            arr[arr >= 9999.0] = np.nan

        # pres may be (N_PROF, N_LEVELS) → take median as shared grid
        if pres.ndim == 2:
            pres_grid = np.nanmedian(pres, axis=0)
        else:
            pres_grid = pres
        return pres_grid, temp, psal, juld


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def plot_ts_calibration(
    smru_name: str,
    *,
    cora_data: dict[str, np.ndarray],
    target_path: Path,
    other_paths: Sequence[Path] = (),
    output_dir: Path,
    chunk_size: int = 200,
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
    other_paths:
        Paths to other tags in the same deployment (used for context cloud).
    output_dir:
        Directory where PNG files are written.
    chunk_size:
        Maximum number of target profiles per figure.

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

    pres_grid, tgt_temp, tgt_psal, tgt_juld = _extract_profiles(target_path)
    n_tgt = tgt_temp.shape[0]

    # Build CORA median salinity per pressure level for anomaly panel
    cora_pres = cora_data.get("pres", np.array([]))
    cora_temp = cora_data.get("temp", np.empty((0, len(pres_grid))))
    cora_psal = cora_data.get("psal", np.empty((0, len(pres_grid))))

    # Interpolate CORA to the target pressure grid if needed
    if len(cora_pres) > 0 and not np.allclose(cora_pres, pres_grid, equal_nan=True):
        cora_temp_interp = np.full((cora_temp.shape[0], len(pres_grid)), np.nan)
        cora_psal_interp = np.full((cora_psal.shape[0], len(pres_grid)), np.nan)
        for i in range(cora_temp.shape[0]):
            valid = ~np.isnan(cora_pres) & ~np.isnan(cora_temp[i])
            if valid.sum() >= 2:
                cora_temp_interp[i] = np.interp(
                    pres_grid, cora_pres[valid], cora_temp[i, valid], left=np.nan, right=np.nan
                )
            valid = ~np.isnan(cora_pres) & ~np.isnan(cora_psal[i])
            if valid.sum() >= 2:
                cora_psal_interp[i] = np.interp(
                    pres_grid, cora_pres[valid], cora_psal[i, valid], left=np.nan, right=np.nan
                )
        cora_temp = cora_temp_interp
        cora_psal = cora_psal_interp
        cora_pres = pres_grid

    cora_psal_median = np.nanmedian(cora_psal, axis=0) if cora_psal.shape[0] > 0 else np.full_like(pres_grid, np.nan)

    # Build "other" profiles cloud
    other_temp_list: list[np.ndarray] = []
    other_psal_list: list[np.ndarray] = []
    for op in other_paths:
        if not op.exists():
            continue
        try:
            _, ot, op_psal, _ = _extract_profiles(op)
            other_temp_list.append(ot)
            other_psal_list.append(op_psal)
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
        for i in range(c_psal.shape[0]):
            anom = c_psal[i] - cora_psal_median
            ax_anom.plot(anom, pres_grid, color=cmap(c_norm[i]), lw=0.8, alpha=0.8)

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

    return written
