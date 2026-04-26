from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


REGISTRY_FILE = "versions.json"
VALID_STATUSES = {"development", "published"}


@dataclass(frozen=True)
class VersionRecord:
    version: str
    status: str
    first_seen: str
    last_updated: str
    output_dir: str



def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()



def _registry_path(public_root: Path) -> Path:
    return public_root / REGISTRY_FILE



def load_version_registry(public_root: Path) -> list[VersionRecord]:
    path = _registry_path(public_root)
    if not path.exists():
        return []
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return []
    if not isinstance(payload, dict):
        return []
    rows = payload.get("versions", [])
    records: list[VersionRecord] = []
    if not isinstance(rows, list):
        return records
    for row in rows:
        if not isinstance(row, dict):
            continue
        version = str(row.get("version", "")).strip()
        status = str(row.get("status", "")).strip().lower()
        if not version or status not in VALID_STATUSES:
            continue
        records.append(
            VersionRecord(
                version=version,
                status=status,
                first_seen=str(row.get("first_seen", "")),
                last_updated=str(row.get("last_updated", "")),
                output_dir=str(row.get("output_dir", "")),
            )
        )
    records.sort(key=lambda item: item.version)
    return records



def update_version_registry(
    public_root: Path,
    *,
    version: str,
    status: str,
    output_dir: Path,
) -> Path:
    normalized_version = str(version).strip()
    normalized_status = str(status).strip().lower()
    if not normalized_version:
        raise ValueError("version must be a non-empty string")
    if normalized_status not in VALID_STATUSES:
        raise ValueError(f"status must be one of {sorted(VALID_STATUSES)}")

    path = _registry_path(public_root)
    now = _utc_now_iso()
    existing = load_version_registry(public_root)
    by_version = {record.version: record for record in existing}
    previous = by_version.get(normalized_version)
    first_seen = previous.first_seen if previous and previous.first_seen else now
    by_version[normalized_version] = VersionRecord(
        version=normalized_version,
        status=normalized_status,
        first_seen=first_seen,
        last_updated=now,
        output_dir=str(output_dir),
    )

    rows = [
        {
            "version": record.version,
            "status": record.status,
            "first_seen": record.first_seen,
            "last_updated": record.last_updated,
            "output_dir": record.output_dir,
        }
        for record in sorted(by_version.values(), key=lambda item: item.version)
    ]
    payload = {"versions": rows}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path
