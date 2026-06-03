from __future__ import annotations

import filecmp
import json
import shutil
import zipfile
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

from ..catalog.deployments import DEPLOYMENT_FIELD_ORDER
from ..catalog.tables import read_csv_rows, write_csv_rows
from ..models import MeopConfig


JSON_SPECS: dict[str, str] = {
    "deployment2.json": "deployment_code",
    "deployment3.json": "deployment_code",
    "platform2.json": "smru_platform_code",
    "platform3.json": "smru_platform_code",
}


@dataclass(frozen=True)
class SourceDeployment:
    deployment: str
    path: Path
    odv_zip: Path | None
    generic_zip: Path | None
    generic_members: tuple[str, ...] = ()

    @property
    def has_odv_zip(self) -> bool:
        return self.odv_zip is not None and self.odv_zip.exists()

    @property
    def has_generic_zip(self) -> bool:
        return self.generic_zip is not None and self.generic_zip.exists()

    @property
    def generic_has_mdb(self) -> bool:
        return any(name.lower().endswith(".mdb") for name in self.generic_members)

    @property
    def generic_has_odv(self) -> bool:
        return any("odv" in name.lower() for name in self.generic_members)

    @property
    def metadata_only_mdb(self) -> bool:
        return not self.has_odv_zip and self.generic_has_mdb and not self.generic_has_odv


@dataclass(frozen=True)
class SmruSyncResult:
    source_dir: Path
    resolved_source_dir: Path
    all_dir: Path
    apply: bool
    scanned_deployments: tuple[str, ...]
    new_processable_deployments: tuple[str, ...] = ()
    metadata_only_deployments: tuple[str, ...] = ()
    copied_odv: tuple[Path, ...] = ()
    replaced_odv: tuple[Path, ...] = ()
    unchanged_odv: tuple[str, ...] = ()
    backed_up_files: tuple[Path, ...] = ()
    written_files: tuple[Path, ...] = ()
    impacted_deployments: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    json_source_records: dict[str, int] = field(default_factory=dict)
    json_added_records: dict[str, int] = field(default_factory=dict)
    json_updated_records: dict[str, int] = field(default_factory=dict)

    def as_dict(self) -> dict[str, object]:
        return {
            "source_dir": str(self.source_dir),
            "resolved_source_dir": str(self.resolved_source_dir),
            "all_dir": str(self.all_dir),
            "apply": self.apply,
            "scanned_deployments": list(self.scanned_deployments),
            "new_processable_deployments": list(self.new_processable_deployments),
            "metadata_only_deployments": list(self.metadata_only_deployments),
            "copied_odv": [str(path) for path in self.copied_odv],
            "replaced_odv": [str(path) for path in self.replaced_odv],
            "unchanged_odv": list(self.unchanged_odv),
            "backed_up_files": [str(path) for path in self.backed_up_files],
            "written_files": [str(path) for path in self.written_files],
            "impacted_deployments": list(self.impacted_deployments),
            "warnings": list(self.warnings),
            "json_source_records": dict(self.json_source_records),
            "json_added_records": dict(self.json_added_records),
            "json_updated_records": dict(self.json_updated_records),
        }

    def format_markdown(self) -> str:
        mode = "APPLY" if self.apply else "DRY RUN"
        lines = [
            f"# SMRU data sync report ({mode})",
            "",
            f"- source: `{self.source_dir}`",
            f"- resolved source: `{self.resolved_source_dir}`",
            f"- deployments scanned: {len(self.scanned_deployments)}",
            f"- new processable deployments: {len(self.new_processable_deployments)}",
            f"- metadata-only MDB deployments: {len(self.metadata_only_deployments)}",
            f"- ODV copied: {len(self.copied_odv)}",
            f"- ODV replaced: {len(self.replaced_odv)}",
            f"- ODV unchanged: {len(self.unchanged_odv)}",
            f"- impacted deployments: {len(self.impacted_deployments)}",
            "",
        ]
        if not self.apply:
            lines.extend(["No files were changed.", ""])
        if self.new_processable_deployments:
            lines.extend(["## New Processable Deployments", ""])
            lines.extend(f"- `{deployment}`" for deployment in self.new_processable_deployments)
            lines.append("")
        if self.metadata_only_deployments:
            lines.extend(["## Metadata-Only MDB Deployments", ""])
            lines.extend(f"- `{deployment}`" for deployment in self.metadata_only_deployments)
            lines.append("")
        if self.replaced_odv:
            lines.extend(["## Replaced ODV Inputs", ""])
            lines.extend(f"- `{path.name}`" for path in self.replaced_odv)
            lines.append("")
        if self.copied_odv:
            lines.extend(["## Copied ODV Inputs", ""])
            lines.extend(f"- `{path.name}`" for path in self.copied_odv)
            lines.append("")
        if self.json_source_records:
            lines.extend(["## JSON Merge Summary", ""])
            for name in sorted(self.json_source_records):
                lines.append(
                    f"- `{name}`: source={self.json_source_records.get(name, 0)}, "
                    f"added={self.json_added_records.get(name, 0)}, "
                    f"updated={self.json_updated_records.get(name, 0)}"
                )
            lines.append("")
        if self.warnings:
            lines.extend(["## Warnings", ""])
            lines.extend(f"- {warning}" for warning in self.warnings)
            lines.append("")
        return "\n".join(lines).rstrip() + "\n"


def _timestamp() -> str:
    return datetime.now().strftime("%Y%m%dT%H%M%S")


def _read_json_records(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if isinstance(payload, list):
        return [record for record in payload if isinstance(record, dict)]
    if isinstance(payload, dict):
        return [payload]
    return []


def _write_json_records(path: Path, records: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(list(records), handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def _zip_members(path: Path | None) -> tuple[str, ...]:
    if path is None or not path.exists():
        return ()
    try:
        with zipfile.ZipFile(path) as archive:
            return tuple(archive.namelist())
    except zipfile.BadZipFile:
        return ()


def _source_deployments(all_dir: Path) -> tuple[SourceDeployment, ...]:
    deployments: list[SourceDeployment] = []
    for directory in sorted(path for path in all_dir.iterdir() if path.is_dir()):
        deployment = directory.name
        odv_zip = directory / f"{deployment}_ODV.zip"
        generic_zip = directory / f"{deployment}.zip"
        deployments.append(
            SourceDeployment(
                deployment=deployment,
                path=directory,
                odv_zip=odv_zip if odv_zip.exists() else None,
                generic_zip=generic_zip if generic_zip.exists() else None,
                generic_members=_zip_members(generic_zip if generic_zip.exists() else None),
            )
        )
    return tuple(deployments)


def _merge_records(
    existing: list[dict[str, Any]],
    source: list[dict[str, Any]],
    *,
    key: str,
) -> tuple[list[dict[str, Any]], int, int]:
    positions: dict[str, int] = {}
    merged: list[dict[str, Any]] = []
    for record in existing:
        value = str(record.get(key, "")).strip()
        if value and value not in positions:
            positions[value] = len(merged)
        merged.append(dict(record))

    added = 0
    updated = 0
    for record in source:
        value = str(record.get(key, "")).strip()
        if not value:
            continue
        current = dict(record)
        if value in positions:
            index = positions[value]
            if merged[index] != current:
                merged[index] = current
                updated += 1
            continue
        positions[value] = len(merged)
        merged.append(current)
        added += 1
    return merged, added, updated


def _source_json_records(deployments: Iterable[SourceDeployment], name: str) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for deployment in deployments:
        records.extend(_read_json_records(deployment.path / name))
    return records


def _backup_file(path: Path, *, timestamp: str) -> Path:
    backup = path.with_name(f"{path.stem}_{timestamp}{path.suffix}")
    shutil.copy2(path, backup)
    return backup


def _deployment_dates(platform_records: list[dict[str, Any]]) -> tuple[str, str]:
    starts = [str(record.get("time_coverage_start", ""))[:10] for record in platform_records if str(record.get("time_coverage_start", "")).strip()]
    ends = [str(record.get("time_coverage_end", ""))[:10] for record in platform_records if str(record.get("time_coverage_end", "")).strip()]
    return (min(starts) if starts else "", max(ends) if ends else "")


def _deployment_catalog_row(deployment: SourceDeployment, *, process: str) -> dict[str, str]:
    records = _read_json_records(deployment.path / "deployment3.json") or _read_json_records(deployment.path / "deployment2.json")
    record = records[0] if records else {}
    platforms = _read_json_records(deployment.path / "platform3.json") or _read_json_records(deployment.path / "platform2.json")
    start_date, end_date = _deployment_dates(platforms)
    return {
        "deployment_code": deployment.deployment,
        "pi_code": str(record.get("pi_code", "")).strip(),
        "process": process,
        "public": "0",
        "country": "UNKNOWN",
        "task_done": "",
        "first_version": "",
        "last_version": "",
        "start_date": start_date,
        "end_date": end_date,
        "start_date_jul": "",
        "description": str(record.get("description", "")).strip(),
        "gts": str(record.get("gts", "")).strip(),
    }


def _write_catalog_if_changed(config: MeopConfig, rows: list[dict[str, str]], *, apply: bool, timestamp: str) -> tuple[Path | None, Path | None]:
    path = config.catalogdir / "list_deployment.csv"
    existing = read_csv_rows(path)
    if existing == rows:
        return None, None
    if not apply:
        return path, None
    backup = _backup_file(path, timestamp=timestamp) if path.exists() else None
    written = write_csv_rows(path, rows, field_order=DEPLOYMENT_FIELD_ORDER + ("description", "gts"))
    return written, backup


def sync_smru_data(
    config: MeopConfig,
    *,
    source_dir: str | Path,
    apply: bool = False,
    timestamp: str | None = None,
) -> SmruSyncResult:
    source_path = Path(source_dir).expanduser()
    resolved_source = source_path.resolve()
    all_dir = resolved_source / "all"
    if not all_dir.exists():
        return SmruSyncResult(
            source_dir=source_path,
            resolved_source_dir=resolved_source,
            all_dir=all_dir,
            apply=apply,
            scanned_deployments=(),
            warnings=(f"source deployment directory does not exist: {all_dir}",),
        )

    stamp = timestamp or _timestamp()
    deployments = _source_deployments(all_dir)
    scanned = tuple(deployment.deployment for deployment in deployments)

    catalog_path = config.catalogdir / "list_deployment.csv"
    catalog_rows = read_csv_rows(catalog_path)
    catalog_by_deployment = {
        str(row.get("deployment_code") or row.get("row_name") or "").strip(): dict(row)
        for row in catalog_rows
        if str(row.get("deployment_code") or row.get("row_name") or "").strip()
    }

    new_processable: list[str] = []
    metadata_only: list[str] = []
    warnings: list[str] = []
    copied_odv: list[Path] = []
    replaced_odv: list[Path] = []
    unchanged_odv: list[str] = []
    backed_up: list[Path] = []
    written_files: list[Path] = []
    impacted: set[str] = set()

    updated_catalog_rows = [dict(row) for row in catalog_rows]
    for deployment in deployments:
        in_catalog = deployment.deployment in catalog_by_deployment
        if not in_catalog:
            if deployment.has_odv_zip:
                new_processable.append(deployment.deployment)
                updated_catalog_rows.append(_deployment_catalog_row(deployment, process="1"))
                impacted.add(deployment.deployment)
            elif deployment.metadata_only_mdb:
                metadata_only.append(deployment.deployment)
                updated_catalog_rows.append(_deployment_catalog_row(deployment, process="0"))
            else:
                warnings.append(f"{deployment.deployment}: no ODV zip and no MDB-only generic zip classification")

        if not deployment.has_odv_zip:
            if deployment.has_generic_zip and not deployment.generic_has_odv:
                warnings.append(f"{deployment.deployment}: generic zip contains no ODV files; MDB conversion is not implemented")
            continue

        assert deployment.odv_zip is not None
        destination = config.raw_odv_dir / deployment.odv_zip.name
        if not destination.exists():
            copied_odv.append(destination)
            impacted.add(deployment.deployment)
            if apply:
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(deployment.odv_zip, destination)
                written_files.append(destination)
            continue
        if filecmp.cmp(deployment.odv_zip, destination, shallow=False):
            unchanged_odv.append(deployment.deployment)
            continue
        replaced_odv.append(destination)
        impacted.add(deployment.deployment)
        if apply:
            backup_dir = config.raw_odv_dir / "backups" / stamp
            backup_dir.mkdir(parents=True, exist_ok=True)
            backup = backup_dir / destination.name
            shutil.copy2(destination, backup)
            shutil.copy2(deployment.odv_zip, destination)
            backed_up.append(backup)
            written_files.append(destination)

    catalog_written, catalog_backup = _write_catalog_if_changed(config, updated_catalog_rows, apply=apply, timestamp=stamp)
    if catalog_written is not None:
        if apply:
            written_files.append(catalog_written)
            if catalog_backup is not None:
                backed_up.append(catalog_backup)
        else:
            written_files.append(catalog_written)

    json_source_counts: dict[str, int] = {}
    json_added_counts: dict[str, int] = {}
    json_updated_counts: dict[str, int] = {}
    for name, key in JSON_SPECS.items():
        source_records = _source_json_records(deployments, name)
        existing_path = config.config_files_dir / name
        existing_records = _read_json_records(existing_path)
        merged, added, updated = _merge_records(existing_records, source_records, key=key)
        json_source_counts[name] = len(source_records)
        json_added_counts[name] = added
        json_updated_counts[name] = updated
        if added == 0 and updated == 0:
            continue
        if apply:
            if existing_path.exists():
                backed_up.append(_backup_file(existing_path, timestamp=stamp))
            _write_json_records(existing_path, merged)
        written_files.append(existing_path)

    return SmruSyncResult(
        source_dir=source_path,
        resolved_source_dir=resolved_source,
        all_dir=all_dir,
        apply=apply,
        scanned_deployments=scanned,
        new_processable_deployments=tuple(sorted(new_processable)),
        metadata_only_deployments=tuple(sorted(metadata_only)),
        copied_odv=tuple(copied_odv),
        replaced_odv=tuple(replaced_odv),
        unchanged_odv=tuple(sorted(unchanged_odv)),
        backed_up_files=tuple(backed_up),
        written_files=tuple(written_files),
        impacted_deployments=tuple(sorted(impacted)),
        warnings=tuple(warnings),
        json_source_records=json_source_counts,
        json_added_records=json_added_counts,
        json_updated_records=json_updated_counts,
    )
