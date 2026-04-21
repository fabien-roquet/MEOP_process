from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable


def _normalize_header(name: str | None) -> str:
    if name is None:
        return ""
    normalized = name.lstrip("\ufeff").strip()
    if normalized.lower() == "unnamed: 0":
        return ""
    return normalized


def read_csv_rows(path: str | Path) -> list[dict[str, str]]:
    csv_path = Path(path)
    if not csv_path.exists():
        return []
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return []
        normalized_fieldnames = [_normalize_header(name) for name in reader.fieldnames]
        rows: list[dict[str, str]] = []
        for raw_row in reader:
            row: dict[str, str] = {}
            for original_key, normalized_key in zip(reader.fieldnames, normalized_fieldnames, strict=False):
                value = raw_row.get(original_key, "")
                row[normalized_key] = "" if value is None else str(value).strip()
            rows.append(row)
        return rows


def read_indexed_csv_rows(path: str | Path, *, index_name: str = "row_name") -> list[dict[str, str]]:
    rows = read_csv_rows(path)
    indexed_rows: list[dict[str, str]] = []
    for row in rows:
        current = dict(row)
        row_name = current.pop("", "")
        if row_name:
            current.setdefault(index_name, row_name)
        indexed_rows.append(current)
    return indexed_rows


def write_csv_rows(
    path: str | Path,
    rows: Iterable[dict[str, str]],
    *,
    field_order: Iterable[str] | None = None,
) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    materialized = [dict(row) for row in rows]

    if field_order is None:
        ordered_fields: list[str] = []
        seen: set[str] = set()
        for row in materialized:
            for key in row:
                if not key or key in seen:
                    continue
                ordered_fields.append(key)
                seen.add(key)
    else:
        ordered_fields = [field for field in field_order if field]
        for row in materialized:
            for key in row:
                if not key or key in ordered_fields:
                    continue
                ordered_fields.append(key)

    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered_fields)
        writer.writeheader()
        writer.writerows([{field: row.get(field, "") for field in ordered_fields} for row in materialized])

    return destination


def write_indexed_csv_rows(
    path: str | Path,
    rows: Iterable[dict[str, str]],
    *,
    index_name: str = "row_name",
    field_order: Iterable[str] | None = None,
) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    materialized = [dict(row) for row in rows]
    if field_order is None:
        ordered_fields: list[str] = []
        seen: set[str] = set()
        for row in materialized:
            for key in row:
                if key == index_name or key in seen:
                    continue
                ordered_fields.append(key)
                seen.add(key)
    else:
        ordered_fields = [field for field in field_order if field != index_name]
        for row in materialized:
            for key in row:
                if key == index_name or key in ordered_fields:
                    continue
                ordered_fields.append(key)

    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([""] + ordered_fields)
        for row in materialized:
            index_value = row.get(index_name, "")
            writer.writerow([index_value] + [row.get(field, "") for field in ordered_fields])
    return destination


def read_json_records(path: str | Path) -> list[dict[str, object]]:
    json_path = Path(path)
    if not json_path.exists():
        return []
    with json_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if isinstance(payload, list):
        return [record for record in payload if isinstance(record, dict)]
    if isinstance(payload, dict):
        return [payload]
    return []


def non_empty_attributes(row: dict[str, str], *, exclude: Iterable[str] = ()) -> dict[str, str]:
    excluded = set(exclude)
    return {
        key: value
        for key, value in row.items()
        if key not in excluded and value not in (None, "")
    }
