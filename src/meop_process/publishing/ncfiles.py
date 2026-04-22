"""Create published NC profile files (_all_prof.nc)."""
from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import xarray as xr

from ..catalog.filenames import deployment_from_smru_name, fname_prof, list_smru_name
from ..models import MeopConfig


# Quality-flag preference order for the base product to publish
_QF_PREFERENCE = ("hr2", "lr1", "fr1", "hr1", "fr0", "lr0")
# Candidate qf files that supply per-standard-level interpolated data
_INTERP_QF = ("hr1",)


@dataclass(frozen=True)
class PublishResult:
    written_files: tuple[Path, ...] = field(default_factory=tuple)
    skipped_files: tuple[Path, ...] = field(default_factory=tuple)
    processed_tags: tuple[str, ...] = field(default_factory=tuple)

    @property
    def written(self) -> bool:
        return bool(self.written_files)


def _best_source(config: MeopConfig, smru_name: str, deployment: str) -> tuple[str, Path] | None:
    for qf in _QF_PREFERENCE:
        path = fname_prof(smru_name, deployment=deployment, qf=qf, config=config)
        if path.exists():
            return qf, path
    return None


def _interp_source(config: MeopConfig, smru_name: str, deployment: str, base_qf: str) -> Path | None:
    if base_qf == "hr2":
        # hr2 already contains the best data; INTERP == ADJUSTED
        return None
    for qf in _INTERP_QF:
        path = fname_prof(smru_name, deployment=deployment, qf=qf, config=config)
        if path.exists():
            return path
    return None


def _add_interp_variables(dataset: xr.Dataset, interp_path: Path) -> xr.Dataset:
    """Add PRES/TEMP/PSAL _INTERP variables from an interpolated (hr1) source file."""
    try:
        with xr.open_dataset(interp_path) as src:
            src_loaded = src.load()
            result = dataset.copy(deep=True)
            for base in ("PRES", "TEMP", "PSAL"):
                src_var = f"{base}_ADJUSTED" if f"{base}_ADJUSTED" in src_loaded else base
                interp_name = f"{base}_INTERP"
                if src_var in src_loaded and interp_name not in result:
                    interp_data = src_loaded[src_var]
                    result[interp_name] = interp_data.rename(
                        {dim: dim.replace("N_LEVELS", "N_INTERP") for dim in interp_data.dims}
                    )
            return result
    except Exception:
        return dataset


def _add_base_interp_variables(dataset: xr.Dataset) -> xr.Dataset:
    """When base is hr2, add INTERP aliases pointing at ADJUSTED data (backward compat)."""
    result = dataset.copy(deep=True)
    for base in ("PRES", "TEMP", "PSAL"):
        src_var = f"{base}_ADJUSTED"
        interp_name = f"{base}_INTERP"
        if src_var in result and interp_name not in result:
            adjusted = result[src_var]
            result[interp_name] = adjusted.rename(
                {dim: dim.replace("N_LEVELS", "N_INTERP") for dim in adjusted.dims}
            )
    return result


def create_ncfile_all(
    config: MeopConfig,
    smru_name: str,
    *,
    output_dir: Path,
    rebuild: bool = False,
) -> Path | None:
    """Create <smru_name>_all_prof.nc in output_dir from the best available processed product.

    Returns the output path on success, or None if no source file is found.
    """
    deployment = deployment_from_smru_name(smru_name)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{smru_name}_all_prof.nc"

    if output_path.exists() and not rebuild:
        return output_path

    found = _best_source(config, smru_name, deployment)
    if found is None:
        return None
    base_qf, source_path = found

    from ..processing.netcdf import save_dataset_with_compression  # local import avoids cycle

    if base_qf == "hr2":
        with xr.open_dataset(source_path) as ds:
            result = _add_base_interp_variables(ds.load())
        save_dataset_with_compression(result, output_path)
    else:
        interp_path = _interp_source(config, smru_name, deployment, base_qf)
        if interp_path is not None:
            with xr.open_dataset(source_path) as ds:
                result = _add_interp_variables(ds.load(), interp_path)
            save_dataset_with_compression(result, output_path)
        else:
            shutil.copy2(source_path, output_path)

    return output_path


def publish_ncfiles(
    config: MeopConfig,
    *,
    output_dir: Path,
    deployment: str = "",
    smru_name: str = "",
    rebuild: bool = False,
) -> PublishResult:
    """Create _all_prof.nc files for all matching tags and return a summary."""
    if smru_name:
        candidates = [smru_name]
    else:
        dep = deployment or ""
        candidates = list_smru_name(smru_name="", deployment=dep, qf="*", config=config)
        if not candidates and dep:
            # fall back: discover from any qf
            from ..metadata.summaries import _discover_processed_products

            inventory = _discover_processed_products(config)
            candidates = list(inventory.get(dep, {}).keys())
        elif not candidates and not dep:
            from ..metadata.summaries import _discover_processed_products

            inventory = _discover_processed_products(config)
            candidates = [tag for tags in inventory.values() for tag in tags]

    written: list[Path] = []
    skipped: list[Path] = []
    processed: list[str] = []

    for tag in sorted(set(candidates)):
        target = create_ncfile_all(config, tag, output_dir=output_dir, rebuild=rebuild)
        if target is not None:
            if target in written or not rebuild:
                pass
            written.append(target)
            processed.append(tag)

    return PublishResult(
        written_files=tuple(written),
        skipped_files=tuple(skipped),
        processed_tags=tuple(processed),
    )
