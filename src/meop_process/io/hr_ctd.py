from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..catalog.deployments import load_hr_catalog
from ..models import MeopConfig


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
    if value.lower() in {"nan", "none"}:
        return ""
    return value if value.endswith("_") else f"{value}_"


def resolve_hr_ctd_path(config: MeopConfig, smru_platform_code: str) -> HrCtdPath | None:
    """Resolve the HR raw-file path strictly from ``list_deployment_hr.csv``.

    No fallback to platform JSON metadata or filename globbing is performed here. The HR filename
    is defined solely by the catalog row: ``data/raw_smru_hr_data/<year>/<prefix><instr_id>_ctd.txt``.
    """

    catalog = load_hr_catalog(config)
    row = catalog.get(smru_platform_code)
    if row is None:
        return None

    year = str(row.get("year", "")).strip()[:4]
    instr_id = str(row.get("instr_id", "")).strip()
    prefix = _clean_prefix(str(row.get("prefix", "")))
    continuous = str(row.get("continuous", "")).strip() in {"1", "true", "True", "yes", "Y"}
    candidates = [root / year / f"{prefix}{instr_id}_ctd.txt" for root in config.raw_hr_search_dirs]
    expected_path = next((path for path in candidates if path.exists()), candidates[0])

    return HrCtdPath(
        smru_platform_code=smru_platform_code,
        year=year,
        instr_id=instr_id,
        prefix=prefix,
        continuous=continuous,
        expected_path=expected_path,
        exists=expected_path.exists(),
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
