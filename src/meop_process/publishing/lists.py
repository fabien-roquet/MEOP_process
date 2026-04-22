"""Build list_profiles.csv from published NC files."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from ..io.netcdf import juld_to_datetime


_FILL_VALUE = 99999.0
_LIST_PROFILES_COLUMNS = (
    "SMRU_PLATFORM_CODE",
    "DEPLOYMENT_CODE",
    "PROFILE_INDEX",
    "JULD",
    "LATITUDE",
    "LONGITUDE",
    "N_LEVELS",
    "FILENAME",
)


def _decode_attr(ds: xr.Dataset, key: str, default: str = "") -> str:
    val = ds.attrs.get(key, default)
    if isinstance(val, (bytes, bytearray)):
        return val.decode("utf-8", errors="replace").strip()
    return str(val).strip()


def _valid_float(val: object) -> float | None:
    try:
        v = float(val)  # type: ignore[arg-type]
        if abs(v) >= _FILL_VALUE:
            return None
        return v
    except (TypeError, ValueError):
        return None


def _count_valid_levels(ds: xr.Dataset, prof_idx: int) -> int:
    for var in ("PRES", "PRES_ADJUSTED", "PRES_INTERP"):
        if var in ds:
            arr = ds[var].values[prof_idx]
            mask = np.isfinite(arr) & (arr < _FILL_VALUE)
            return int(mask.sum())
    return 0


def _profiles_from_file(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    try:
        with xr.open_dataset(path, decode_times=False) as ds:
            smru = _decode_attr(ds, "smru_platform_code", path.name.split("_")[0])
            deployment = _decode_attr(ds, "deployment_code", "")
            if not deployment:
                from ..catalog.filenames import deployment_from_smru_name

                deployment = deployment_from_smru_name(smru)

            n_prof = ds.sizes.get("N_PROF", 0)
            juld_vals: np.ndarray
            if "JULD" in ds:
                raw = np.asarray(ds["JULD"].values, dtype=float)
                juld_vals = raw
            else:
                juld_vals = np.full(n_prof, np.nan)

            lat_vals = np.asarray(ds["LATITUDE"].values, dtype=float) if "LATITUDE" in ds else np.full(n_prof, np.nan)
            lon_vals = np.asarray(ds["LONGITUDE"].values, dtype=float) if "LONGITUDE" in ds else np.full(n_prof, np.nan)

            datetimes = juld_to_datetime(juld_vals)

            for i in range(n_prof):
                dt = datetimes[i]
                juld_str = dt.strftime("%Y-%m-%d %H:%M:%S") if dt is not None and not (isinstance(dt, float) and np.isnan(dt)) else ""
                rows.append({
                    "SMRU_PLATFORM_CODE": smru,
                    "DEPLOYMENT_CODE": deployment,
                    "PROFILE_INDEX": i,
                    "JULD": juld_str,
                    "LATITUDE": _valid_float(lat_vals[i]),
                    "LONGITUDE": _valid_float(lon_vals[i]),
                    "N_LEVELS": _count_valid_levels(ds, i),
                    "FILENAME": path.name,
                })
    except Exception:
        pass
    return rows


def build_list_profiles(
    output_dir: Path,
    *,
    catalog_path: Path | None = None,
    verbose: bool = False,
) -> Path | None:
    """Scan output_dir for *_all_prof.nc files and write list_profiles.csv.

    Parameters
    ----------
    output_dir:
        Directory containing ``*_all_prof.nc`` files.
    catalog_path:
        Optional path to ``list_deployment.csv`` to enrich the output with a
        COUNTRY column.  When provided, REGION is also added via regionmask.
    verbose:
        Print progress to stdout.

    Returns
    -------
    Path to the written CSV, or None if no files were found.
    """
    nc_files = sorted(output_dir.glob("*_all_prof.nc"))
    if not nc_files:
        return None

    all_rows: list[dict[str, object]] = []
    for path in nc_files:
        rows = _profiles_from_file(path)
        all_rows.extend(rows)
        if verbose:
            print(f"  {path.name}: {len(rows)} profiles")

    if not all_rows:
        return None

    frame = pd.DataFrame(all_rows, columns=list(_LIST_PROFILES_COLUMNS))

    # Optionally enrich with REGION and COUNTRY
    try:
        from ..plotting.maps import enrich_profiles_dataframe
        frame = enrich_profiles_dataframe(frame, catalog_path=catalog_path)
    except Exception:
        pass

    dest = output_dir / "list_profiles.csv"
    frame.to_csv(dest, index=False)
    if verbose:
        print(f"  wrote {len(frame)} rows to {dest}")
    return dest
