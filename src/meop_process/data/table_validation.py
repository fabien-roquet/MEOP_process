from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from ..models import MeopConfig


@dataclass(frozen=True)
class TableValidationIssue:
    level: str
    table: str
    message: str
    row_number: int | None = None
    column: str = ""

    def as_dict(self) -> dict[str, object]:
        return {
            "level": self.level,
            "table": self.table,
            "message": self.message,
            "row_number": self.row_number,
            "column": self.column,
        }


@dataclass(frozen=True)
class TableValidationResult:
    checked_tables: tuple[str, ...]
    issues: tuple[TableValidationIssue, ...]

    @property
    def errors(self) -> tuple[TableValidationIssue, ...]:
        return tuple(issue for issue in self.issues if issue.level == "error")

    @property
    def warnings(self) -> tuple[TableValidationIssue, ...]:
        return tuple(issue for issue in self.issues if issue.level == "warning")

    @property
    def ok(self) -> bool:
        return not self.errors

    def as_dict(self) -> dict[str, object]:
        return {
            "ok": self.ok,
            "checked_tables": list(self.checked_tables),
            "errors": [issue.as_dict() for issue in self.errors],
            "warnings": [issue.as_dict() for issue in self.warnings],
        }

    def format_text(self) -> str:
        lines = [
            "MEOP runtime table validation",
            f"checked tables: {len(self.checked_tables)}",
            f"errors: {len(self.errors)}",
            f"warnings: {len(self.warnings)}",
        ]
        for issue in self.issues:
            location = issue.table
            if issue.row_number is not None:
                location += f": row {issue.row_number}"
            if issue.column:
                location += f": {issue.column}"
            lines.append(f"{issue.level.upper():7s} {location}: {issue.message}")
        if self.ok:
            lines.append("OK: runtime tables passed validation")
        return "\n".join(lines)


@dataclass(frozen=True)
class _TableSchema:
    required_columns: tuple[str, ...]
    key_column: str
    numeric_columns: tuple[str, ...] = ()
    blank_columns_error: bool = True
    allow_blank_numeric: bool = False
    duplicate_key_level: str = "error"


TABLE_COEFF_COMMENT_CODES = frozenset(
    {
        "OK",
        "RM",
        "SRM",
        "NO_POS",
        "SMRU_WAIT",
        "MDB_MANY",
        "MDB_SEV",
        "MDB_SOME",
        "MDB_NONE",
        "DIVE_POS",
        "GPS_LAND",
        "INVIS",
        "S_BAD",
        "S_STOP",
        "S_SURF_HI",
        "T_WARM",
        "T_DRIFT",
        "T_BAD",
        "TS_BAD",
        "TS_DRIFT",
        "SPLIT_CHECK",
        "OFFSET_VAR",
        "OFFSET_HIGH",
        "DATE_FILTER",
        "SHALLOW",
        "PRES_MISS",
        "FREEZE_NOREF",
        "PART_WEIRD",
        "RM_CHECK",
        "ALL_BAD",
    }
)


TABLE_SCHEMAS: dict[str, _TableSchema] = {
    "table_coeff.csv": _TableSchema(
        required_columns=("smru_platform_code", "T1", "T2", "S1", "S2", "remove", "Sremove", "comment"),
        key_column="smru_platform_code",
        numeric_columns=("T1", "T2", "S1", "S2", "remove", "Sremove"),
        blank_columns_error=True,
        allow_blank_numeric=False,
    ),
    "table_coeff_comment_codes.csv": _TableSchema(
        required_columns=("code", "short_label", "details"),
        key_column="code",
        blank_columns_error=True,
    ),
    "table_config_plots.csv": _TableSchema(
        required_columns=("smru_platform_code", "Tmin", "Tmax", "Smin", "Smax", "lat_min", "lat_max", "lon_min", "lon_max"),
        key_column="smru_platform_code",
        numeric_columns=("Tmin", "Tmax", "Smin", "Smax", "lat_min", "lat_max", "lon_min", "lon_max"),
        blank_columns_error=True,
        allow_blank_numeric=True,
    ),
    "table_filter.csv": _TableSchema(
        required_columns=("smru_platform_name", "Sonly", "filter", "x1", "x2"),
        key_column="smru_platform_name",
        numeric_columns=("Sonly", "x1", "x2"),
        blank_columns_error=False,
        allow_blank_numeric=True,
        duplicate_key_level="allow",
    ),
    "table_meta.csv": _TableSchema(
        required_columns=("smru_platform_code", "location", "reference_doi", "species"),
        key_column="smru_platform_code",
        blank_columns_error=True,
    ),
    "table_param.csv": _TableSchema(
        required_columns=(
            "deployment_code",
            "temp_error",
            "psal_error",
            "minT",
            "maxT",
            "minS",
            "maxS",
            "min_Nprof",
            "pmax",
            "pmax_fluo",
            "is_lon_center_180",
        ),
        key_column="deployment_code",
        numeric_columns=("temp_error", "psal_error", "minT", "maxT", "minS", "maxS", "min_Nprof", "pmax", "pmax_fluo", "is_lon_center_180"),
        blank_columns_error=True,
        allow_blank_numeric=True,
    ),
    "table_salinity_offsets.csv": _TableSchema(
        required_columns=("smru_platform_code", "index_1", "index_2", "index_3", "index_4", "offset_1", "offset_2", "offset_3", "offset_4"),
        key_column="smru_platform_code",
        numeric_columns=("index_1", "index_2", "index_3", "index_4", "offset_1", "offset_2", "offset_3", "offset_4"),
        blank_columns_error=True,
        allow_blank_numeric=True,
        duplicate_key_level="warning",
    ),
    "table_split_tags.csv": _TableSchema(
        required_columns=("smru_platform_name", "nsplit"),
        key_column="smru_platform_name",
        numeric_columns=("nsplit",),
        blank_columns_error=True,
        allow_blank_numeric=False,
    ),
}


def _normalize_header(name: str | None) -> str:
    if name is None:
        return ""
    value = name.lstrip("\ufeff").strip()
    if value.lower() == "unnamed: 0":
        return ""
    return value


def _is_number(value: str) -> bool:
    try:
        float(value)
    except (TypeError, ValueError):
        return False
    return True


def _validate_table_coeff_comment(row: dict[str, str], *, table: str, row_index: int) -> list[TableValidationIssue]:
    issues: list[TableValidationIssue] = []
    comment = row.get("comment", "").strip()
    if not comment:
        return [TableValidationIssue("error", table, "blank standardized comment", row_number=row_index, column="comment")]
    if len(comment) > 64:
        issues.append(TableValidationIssue("error", table, "standardized comment is too long", row_number=row_index, column="comment"))
    codes = [part.strip() for part in comment.split(";") if part.strip()]
    if not codes:
        issues.append(TableValidationIssue("error", table, "blank standardized comment", row_number=row_index, column="comment"))
        return issues

    for code in codes:
        if code not in TABLE_COEFF_COMMENT_CODES:
            issues.append(TableValidationIssue("error", table, f"unknown standardized comment code {code!r}", row_number=row_index, column="comment"))
    if len(set(codes)) != len(codes):
        issues.append(TableValidationIssue("warning", table, "duplicate standardized comment code", row_number=row_index, column="comment"))
    if "OK" in codes and len(codes) > 1:
        issues.append(TableValidationIssue("error", table, "OK cannot be combined with other comment codes", row_number=row_index, column="comment"))

    remove = row.get("remove", "").strip()
    sremove = row.get("Sremove", "").strip()
    if remove == "1" and "RM" not in codes:
        issues.append(TableValidationIssue("warning", table, "remove=1 without RM comment code", row_number=row_index, column="comment"))
    if remove != "1" and "RM" in codes:
        issues.append(TableValidationIssue("warning", table, "RM comment code without remove=1", row_number=row_index, column="comment"))
    if sremove == "1" and "SRM" not in codes:
        issues.append(TableValidationIssue("warning", table, "Sremove=1 without SRM comment code", row_number=row_index, column="comment"))
    if sremove != "1" and "SRM" in codes:
        issues.append(TableValidationIssue("warning", table, "SRM comment code without Sremove=1", row_number=row_index, column="comment"))
    return issues


def _validate_csv_table(path: Path, table: str, schema: _TableSchema) -> list[TableValidationIssue]:
    issues: list[TableValidationIssue] = []
    if not path.exists():
        return [TableValidationIssue("error", table, f"missing runtime table: {path}")]

    try:
        handle = path.open("r", encoding="utf-8-sig", newline="")
    except OSError as exc:
        return [TableValidationIssue("error", table, f"cannot open table: {exc}")]

    with handle:
        try:
            reader = csv.DictReader(handle)
        except csv.Error as exc:
            return [TableValidationIssue("error", table, f"cannot parse CSV header: {exc}")]
        if reader.fieldnames is None:
            return [TableValidationIssue("error", table, "missing CSV header")]

        raw_headers = list(reader.fieldnames)
        headers = [_normalize_header(name) for name in raw_headers]
        blank_indices = [index + 1 for index, name in enumerate(headers) if not name]
        for index in blank_indices:
            level = "error" if schema.blank_columns_error else "warning"
            issues.append(TableValidationIssue(level, table, f"blank or unnamed column in header position {index}"))

        blocking_header_error = False
        duplicate_headers = sorted({name for name in headers if name and headers.count(name) > 1})
        for name in duplicate_headers:
            issues.append(TableValidationIssue("error", table, "duplicate CSV header", column=name))
            blocking_header_error = True

        header_set = set(headers)
        for column in schema.required_columns:
            if column not in header_set:
                issues.append(TableValidationIssue("error", table, "missing required column", column=column))
                blocking_header_error = True

        if blocking_header_error:
            return issues

        seen_keys: dict[str, int] = {}
        try:
            for row_index, raw_row in enumerate(reader, start=2):
                if None in raw_row:
                    issues.append(TableValidationIssue("error", table, f"row has {len(raw_row[None] or [])} extra field(s)", row_number=row_index))
                row = {
                    _normalize_header(key): "" if value is None else str(value).strip()
                    for key, value in raw_row.items()
                    if key is not None
                }
                key_value = row.get(schema.key_column, "").strip()
                if not key_value:
                    issues.append(TableValidationIssue("error", table, "blank key value", row_number=row_index, column=schema.key_column))
                elif key_value in seen_keys and schema.duplicate_key_level != "allow":
                    issues.append(
                        TableValidationIssue(
                            schema.duplicate_key_level,
                            table,
                            f"duplicate key also seen on row {seen_keys[key_value]}",
                            row_number=row_index,
                            column=schema.key_column,
                        )
                    )
                else:
                    seen_keys[key_value] = row_index

                for column in schema.numeric_columns:
                    value = row.get(column, "").strip()
                    if not value and schema.allow_blank_numeric:
                        continue
                    if not _is_number(value):
                        issues.append(TableValidationIssue("error", table, f"expected numeric value, got {value!r}", row_number=row_index, column=column))
                if table == "table_coeff.csv":
                    issues.extend(_validate_table_coeff_comment(row, table=table, row_index=row_index))
        except csv.Error as exc:
            issues.append(TableValidationIssue("error", table, f"CSV parser error: {exc}"))
    return issues


def validate_runtime_tables(config: MeopConfig, *, tables: Iterable[str] | None = None) -> TableValidationResult:
    requested = tuple(tables or TABLE_SCHEMAS.keys())
    issues: list[TableValidationIssue] = []
    for table in requested:
        schema = TABLE_SCHEMAS.get(table)
        if schema is None:
            issues.append(TableValidationIssue("error", table, "no validation schema is registered for this table"))
            continue
        issues.extend(_validate_csv_table(config.tablesdir / table, table, schema))
    return TableValidationResult(checked_tables=requested, issues=tuple(issues))
