from __future__ import annotations

import filecmp
import shutil
from dataclasses import dataclass
from importlib.resources import as_file, files
from pathlib import Path
from typing import Any

from ..config.paths import ensure_runtime_directories
from ..models import MeopConfig


PACKAGE_TABLES: dict[str, str] = {
    "table_coeff.csv": "Calibration coefficients and algorithm constants.",
    "table_config_plots.csv": "Plot/report configuration.",
    "table_filter.csv": "Filtering controls applied during processing.",
    "table_meta.csv": "Metadata patches applied to processed netCDF files.",
    "table_param.csv": "Processing parameters and thresholds.",
    "table_salinity_offsets.csv": "Salinity offset adjustments.",
    "table_split_tags.csv": "Split-tag handling rules.",
}


PACKAGE_CATALOGS: dict[str, str] = {
    "list_deployment.csv": "Packaged sample deployment catalog used by tests and examples.",
    "list_deployment_hr.csv": "Packaged sample high-resolution catalog used by tests and examples.",
    "list_wmo_numbers.csv": "Optional packaged WMO mapping sample used by tests and examples.",
}

PROCESS_TABLES: dict[str, str] = {
    "list_deployment.csv": "Deployment registry used to resolve deployments and tags.",
    "list_deployment_hr.csv": "High-resolution deployment registry keyed by smru_platform_code.",
    "list_wmo_numbers.csv": "Optional/legacy WMO mapping table.",
}

MIRRORED_CONFIG_JSON: dict[str, str] = {
    "deployment2.json": "Legacy deployment metadata used during catalog completion.",
    "deployment3.json": "Current deployment metadata used during catalog completion.",
    "platform2.json": "Legacy platform metadata used during catalog completion.",
    "platform3.json": "Current platform metadata used during catalog completion.",
    "deployment2_patch.json": "Deployment patch metadata.",
    "platform2_patch.json": "Platform patch metadata.",
}

REFERENCE_SUBTREES: dict[str, str] = {
    "WOD_data": "World Ocean Database reference profiles.",
    "CORA_data/CORA_ncfiles": "CORA climatology/reference profiles.",
    "MEOP_last_stable_version": "Reference MEOP release used in comparisons.",
}

OPTIONAL_DATA_DIRECTORIES: dict[str, str] = {
    "crawl_locations": "Updated Argos crawl locations.",
    "smooth_cls_locations": "Smoothed CLS locations.",
    "raw_smru_data_odv": "Low-resolution ODV zip/text inputs staged directly under data/.",
    "raw_smru_hr_data": "High-resolution raw CTD text inputs staged directly under data/.",
}


@dataclass(frozen=True)
class ValidationRecord:
    category: str
    name: str
    path: Path
    required: bool
    exists: bool
    description: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "category": self.category,
            "name": self.name,
            "path": str(self.path),
            "required": self.required,
            "exists": self.exists,
            "description": self.description,
        }


@dataclass(frozen=True)
class DataLayout:
    config: MeopConfig

    @property
    def package_root(self):
        return files("meop_process.resources")

    @property
    def packaged_tables_root(self):
        return self.package_root.joinpath("tables")

    @property
    def packaged_catalog_root(self):
        return self.package_root.joinpath("catalog")

    def canonical_table_path(self, name: str) -> Path:
        return self.config.tablesdir / name

    def canonical_catalog_path(self, name: str) -> Path:
        return self.config.catalogdir / name

    def raw_odv_pattern(self) -> str:
        return str(self.config.raw_odv_dir / "<deployment>_ODV.zip")

    def raw_odv_text_pattern(self) -> str:
        return str(self.config.raw_odv_dir / "<deployment>_ODV.txt")

    def raw_hr_pattern(self) -> str:
        return str(self.config.raw_hr_dir / "<year>" / "<prefix><instr_id>_ctd.txt")

    def public_release_root(self) -> Path:
        return self.config.publicdir_ctd


def _layout(config: MeopConfig) -> DataLayout:
    return DataLayout(config)


def _copy_if_different(source: Path, destination: Path) -> bool:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and filecmp.cmp(source, destination, shallow=False):
        return False
    shutil.copy2(source, destination)
    return True


def bootstrap_packaged_tables(config: MeopConfig, *, names: list[str] | tuple[str, ...] | None = None) -> list[Path]:
    """Ensure shipped CSV tables are available under ``data/tables`` only."""

    ensure_runtime_directories(config)
    layout = _layout(config)
    requested = list(names or PACKAGE_TABLES.keys())
    changed: list[Path] = []

    for name in requested:
        canonical = layout.canonical_table_path(name)
        resource = layout.packaged_tables_root.joinpath(name)
        if canonical.exists():
            continue
        if not resource.is_file():
            continue
        with as_file(resource) as resource_path:
            if _copy_if_different(resource_path, canonical):
                changed.append(canonical)
    return changed


def bootstrap_packaged_catalogs(config: MeopConfig, *, names: list[str] | tuple[str, ...] | None = None) -> list[Path]:
    """Optionally seed ``data/catalog`` from packaged sample catalog CSV files."""

    ensure_runtime_directories(config)
    layout = _layout(config)
    requested = list(names or PACKAGE_CATALOGS.keys())
    changed: list[Path] = []

    for name in requested:
        canonical = layout.canonical_catalog_path(name)
        resource = layout.packaged_catalog_root.joinpath(name)
        if canonical.exists():
            continue
        if not resource.is_file():
            continue
        with as_file(resource) as resource_path:
            if _copy_if_different(resource_path, canonical):
                changed.append(canonical)
    return changed


def bootstrap_catalog_files(config: MeopConfig, *, names: list[str] | tuple[str, ...] | None = None) -> list[Path]:
    """Catalog files are operator-managed under ``data/catalog`` and are never mirrored elsewhere."""

    ensure_runtime_directories(config)
    changed: list[Path] = []
    for name in names or tuple(PROCESS_TABLES.keys()):
        path = config.catalogdir / name
        if path.exists():
            changed.append(path)
    return changed


def prepare_runtime_environment(config: MeopConfig) -> list[Path]:
    created = ensure_runtime_directories(config)
    changed = bootstrap_packaged_tables(config)
    return created + changed


def resolve_table_path(config: MeopConfig, name: str, *, required: bool = True) -> Path:
    layout = _layout(config)
    canonical = layout.canonical_table_path(name)

    if name in PACKAGE_TABLES:
        bootstrap_packaged_tables(config, names=[name])
        if canonical.exists():
            return canonical
    elif canonical.exists():
        return canonical

    if required:
        raise FileNotFoundError(f"MEOP table not found: {name}")
    return canonical


def resolve_catalog_path(config: MeopConfig, name: str, *, required: bool = True) -> Path:
    path = _layout(config).canonical_catalog_path(name)
    if path.exists() or not required:
        return path
    raise FileNotFoundError(f"MEOP catalog table not found: {name}")


def validate_data_layout(config: MeopConfig) -> list[ValidationRecord]:
    layout = _layout(config)
    prepare_runtime_environment(config)

    records: list[ValidationRecord] = []
    for name, description in PACKAGE_TABLES.items():
        path = layout.canonical_table_path(name)
        records.append(
            ValidationRecord(
                category="packaged_table",
                name=name,
                path=path,
                required=True,
                exists=path.exists(),
                description=description,
            )
        )

    for name, description in PROCESS_TABLES.items():
        path = layout.canonical_catalog_path(name)
        records.append(
            ValidationRecord(
                category="catalog",
                name=name,
                path=path,
                required=name in {"list_deployment.csv", "list_deployment_hr.csv"},
                exists=path.exists(),
                description=description,
            )
        )

    for name, description in MIRRORED_CONFIG_JSON.items():
        path = config.config_files_dir / name
        records.append(
            ValidationRecord(
                category="config_json",
                name=name,
                path=path,
                required=False,
                exists=path.exists(),
                description=description,
            )
        )

    for relative, description in REFERENCE_SUBTREES.items():
        path = config.refdir / relative
        records.append(
            ValidationRecord(
                category="reference_dir",
                name=relative,
                path=path,
                required=False,
                exists=path.exists(),
                description=description,
            )
        )

    for relative, description in OPTIONAL_DATA_DIRECTORIES.items():
        path = config.datadir / relative
        records.append(
            ValidationRecord(
                category="data_dir",
                name=relative,
                path=path,
                required=False,
                exists=path.exists(),
                description=description,
            )
        )

    return records


def describe_data_layout(config: MeopConfig) -> dict[str, Any]:
    prepare_runtime_environment(config)
    return {
        "roots": {
            "process_root": str(config.processdir),
            "data_root": str(config.datadir),
            "tables_root": str(config.tablesdir),
            "catalog_root": str(config.catalogdir),
            "config_files_root": str(config.config_files_dir),
            "raw_odv_root": str(config.raw_odv_dir),
            "raw_hr_root": str(config.raw_hr_dir),
            "references_root": str(config.refdir),
            "final_dataset_root": str(config.final_dataset_dir),
            "plots_root": str(config.plotdir),
            "maps_root": str(config.mapsdir),
            "doc_latex_root": str(config.texdir),
            "public_release_root": str(config.publicdir_ctd),
        },
        "patterns": {
            "raw_odv_zip": _layout(config).raw_odv_pattern(),
            "raw_odv_text": _layout(config).raw_odv_text_pattern(),
            "raw_hr_text": _layout(config).raw_hr_pattern(),
        },
        "packaged_tables": [
            {"name": name, "path": str(config.tablesdir / name), "description": description}
            for name, description in PACKAGE_TABLES.items()
        ],
        "packaged_catalogs": [
            {"name": name, "path": str(_layout(config).packaged_catalog_root.joinpath(name)), "description": description}
            for name, description in PACKAGE_CATALOGS.items()
        ],
        "catalog_tables": [
            {"name": name, "path": str(config.catalogdir / name), "description": description}
            for name, description in PROCESS_TABLES.items()
        ],
    }


def format_data_layout(config: MeopConfig) -> str:
    info = describe_data_layout(config)
    roots = info["roots"]
    patterns = info["patterns"]
    lines = [
        "MEOP runtime data layout",
        "",
        f"- process root: {roots['process_root']}",
        f"- data root: {roots['data_root']}",
        f"- tables: {roots['tables_root']}",
        f"- catalog: {roots['catalog_root']}",
        f"- config json: {roots['config_files_root']}",
        f"- raw ODV: {roots['raw_odv_root']}",
        f"- raw HR: {roots['raw_hr_root']}",
        f"- references: {roots['references_root']}",
        f"- final products: {roots['final_dataset_root']}",
        f"- plots: {roots['plots_root']}",
        f"- maps: {roots['maps_root']}",
        f"- LaTeX outputs: {roots['doc_latex_root']}",
        f"- public release: {roots['public_release_root']}",
        "",
        "Expected raw inputs:",
        f"- low-resolution archive: {patterns['raw_odv_zip']}",
        f"- low-resolution text: {patterns['raw_odv_text']}",
        f"- high-resolution text: {patterns['raw_hr_text']}",
        "",
        "Notes:",
        "- Packaged CSV defaults are seeded into data/tables only.",
        "- Catalog CSV files are read and written from data/catalog only.",
        "- No root-level configuration or table mirrors are required.",
        "- Raw ODV and HR files are staged directly under data/.",
    ]
    return "\n".join(lines)
