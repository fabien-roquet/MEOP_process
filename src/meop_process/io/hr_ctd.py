from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..catalog.deployments import load_hr_catalog
from ..catalog.tables import read_json_records
from ..models import MeopConfig


PLATFORM_JSON_FILES = ("platform3.json", "platform2.json", "platform2_patch.json")


@dataclass(frozen=True)
class HrCtdPath:
    smru_platform_code: str
    year: str
    instr_id: str
    prefix: str
    continuous: bool
    expected_path: Path
    exists: bool

    def as_dict(self) -> dict[str, object]:
        return {
            "smru_platform_code": self.smru_platform_code,
            "year": self.year,
            "instr_id": self.instr_id,
            "prefix": self.prefix,
            "expected_path": str(self.expected_path),
            "exists": self.exists,
        }



def _clean_prefix(raw: str) -> str:
    value = (raw or "").strip()
    if not value:
        return ""
    return value if value.endswith("_") else f"{value}_"



def _append_unique(items: list[str], value: str) -> None:
    cleaned = (value or "").strip()
    if cleaned and cleaned not in items:
        items.append(cleaned)



def _matching_platform_records(config: MeopConfig, smru_platform_code: str) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for name in PLATFORM_JSON_FILES:
        for record in read_json_records(config.config_files_dir / name):
            if str(record.get("smru_platform_code", "")).strip() == smru_platform_code:
                records.append(record)
    return records



def _candidate_years(row: dict[str, str] | None, platform_records: list[dict[str, object]]) -> list[str]:
    years: list[str] = []
    if row is not None:
        _append_unique(years, str(row.get("year", ""))[:4])
    for record in platform_records:
        _append_unique(years, str(record.get("time_coverage_start", ""))[:4])
        _append_unique(years, str(record.get("time_coverage_end", ""))[:4])
        _append_unique(years, str(record.get("dt_created", ""))[:4])
    return years



def _candidate_instr_ids(row: dict[str, str] | None, platform_records: list[dict[str, object]]) -> list[str]:
    instr_ids: list[str] = []
    if row is not None:
        _append_unique(instr_ids, str(row.get("instr_id", "")))
    for record in platform_records:
        _append_unique(instr_ids, str(record.get("instr_id", "")))
    return instr_ids



def _candidate_prefixes(row: dict[str, str] | None) -> list[str]:
    prefixes: list[str] = []
    if row is not None:
        _append_unique(prefixes, _clean_prefix(str(row.get("prefix", ""))))
    if "" not in prefixes:
        prefixes.append("")
    return prefixes



def _iter_candidate_paths(
    config: MeopConfig,
    *,
    years: list[str],
    instr_ids: list[str],
    prefixes: list[str],
) -> list[tuple[str, str, str, Path]]:
    candidates: list[tuple[str, str, str, Path]] = []
    seen: set[str] = set()

    def _add(year: str, instr_id: str, prefix: str, path: Path) -> None:
        key = str(path)
        if key not in seen:
            candidates.append((year, instr_id, prefix, path))
            seen.add(key)

    for year in years:
        if not year:
            continue
        year_dir = config.raw_hr_dir / year
        for prefix in prefixes:
            for instr_id in instr_ids:
                if not instr_id:
                    continue
                _add(year, instr_id, prefix, year_dir / f"{prefix}{instr_id}_ctd.txt")
        for instr_id in instr_ids:
            if not instr_id or not year_dir.exists():
                continue
            for path in sorted(year_dir.glob(f"*{instr_id}*_ctd.txt")):
                _add(year, instr_id, "", path)
    return candidates



def resolve_hr_ctd_path(config: MeopConfig, smru_platform_code: str) -> HrCtdPath | None:
    catalog = load_hr_catalog(config)
    row = catalog.get(smru_platform_code)
    platform_records = _matching_platform_records(config, smru_platform_code)
    if row is None and not platform_records:
        return None

    continuous = False if row is None else str(row.get("continuous", "")).strip() in {"1", "true", "True", "yes", "Y"}
    years = _candidate_years(row, platform_records)
    instr_ids = _candidate_instr_ids(row, platform_records)
    prefixes = _candidate_prefixes(row)

    requested_year = years[0] if years else ""
    requested_instr_id = instr_ids[0] if instr_ids else ""
    requested_prefix = prefixes[0] if prefixes else ""

    selected_year = requested_year
    selected_instr_id = requested_instr_id
    selected_prefix = requested_prefix
    selected_path = config.raw_hr_dir / requested_year / f"{requested_prefix}{requested_instr_id}_ctd.txt"

    for year, instr_id, prefix, path in _iter_candidate_paths(
        config,
        years=years,
        instr_ids=instr_ids,
        prefixes=prefixes,
    ):
        if path.exists():
            selected_year = year
            selected_instr_id = instr_id
            selected_prefix = prefix
            selected_path = path
            break

    return HrCtdPath(
        smru_platform_code=smru_platform_code,
        year=selected_year,
        instr_id=selected_instr_id,
        prefix=selected_prefix,
        continuous=continuous,
        expected_path=selected_path,
        exists=selected_path.exists(),
    )



def discover_hr_ctd_paths(config: MeopConfig, deployment: str = "") -> tuple[HrCtdPath, ...]:
    catalog = load_hr_catalog(config)
    rows = []
    for smru_platform_code in sorted(catalog):
        if deployment and smru_platform_code.split("-")[0] != deployment:
            continue
        resolved = resolve_hr_ctd_path(config, smru_platform_code)
        if resolved is not None:
            rows.append(resolved)
    return tuple(rows)
