from __future__ import annotations

from pathlib import Path

from ..data.layout import resolve_catalog_path
from ..io.raw_odv import build_odv_profile_index, discover_raw_odv_files
from ..models import DeploymentInfo, DeploymentRecord, MeopConfig, Selection
from .filenames import list_fname_prof
from .tables import read_indexed_csv_rows, read_json_records, write_indexed_csv_rows


DEPLOYMENT_FIELD_ORDER = (
    "deployment_code",
    "pi_code",
    "process",
    "public",
    "country",
    "task_done",
    "first_version",
    "last_version",
    "start_date",
    "end_date",
    "start_date_jul",
)

HR_FIELD_ORDER = (
    "smru_platform_code",
    "instr_id",
    "year",
    "prefix",
    "continuous",
)

DEPLOYMENT_JSON_FILES = ("deployment3.json", "deployment2.json", "deployment2_patch.json")
PLATFORM_JSON_FILES = ("platform3.json", "platform2.json", "platform2_patch.json")


def deployment_from_smru_name(smru_name: str) -> str:
    return smru_name.split("-")[0]


def normalize_selection(deployment: str = "", smru_name: str = "") -> Selection:
    return Selection(deployment=deployment, smru_name=smru_name).normalized()


def _clean_string(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _read_deployment_rows(config: MeopConfig) -> list[dict[str, str]]:
    path = resolve_catalog_path(config, "list_deployment.csv", required=False)
    if not path.exists():
        return []
    return read_indexed_csv_rows(path)


def _read_hr_rows(config: MeopConfig) -> list[dict[str, str]]:
    path = resolve_catalog_path(config, "list_deployment_hr.csv", required=False)
    if not path.exists():
        return []
    return read_indexed_csv_rows(path)


def _load_deployment_json_records(config: MeopConfig) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for name in DEPLOYMENT_JSON_FILES:
        records.extend(read_json_records(config.config_files_dir / name))
    return records


def _load_platform_json_records(config: MeopConfig) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for name in PLATFORM_JSON_FILES:
        records.extend(read_json_records(config.config_files_dir / name))
    return records


def load_deployment_catalog(config: MeopConfig, *, persist: bool = True) -> dict[str, DeploymentRecord]:
    """Resolve the deployment registry in Python from ``data/catalog`` and ``data/config_files``."""

    rows = _read_deployment_rows(config)

    catalog: dict[str, DeploymentRecord] = {}
    ordered_rows: list[dict[str, str]] = []
    for row in rows:
        deployment_code = _clean_string(row.get("deployment_code") or row.get("row_name"))
        if not deployment_code:
            continue
        normalized = {key: _clean_string(value) for key, value in row.items()}
        normalized["deployment_code"] = deployment_code
        normalized.setdefault("row_name", _clean_string(row.get("row_name") or deployment_code))
        ordered_rows.append(normalized)
        catalog[deployment_code] = DeploymentRecord(
            deployment_code=deployment_code,
            pi_code=normalized.get("pi_code", ""),
            country=normalized.get("country", ""),
            process=normalized.get("process", ""),
            public=normalized.get("public", ""),
            description=normalized.get("description", ""),
            gts=normalized.get("gts", ""),
            start_date=normalized.get("start_date", ""),
            end_date=normalized.get("end_date", ""),
            task_done=normalized.get("task_done", ""),
            source="catalog",
            row_name=normalized.get("row_name", deployment_code),
            extra={
                key: value
                for key, value in normalized.items()
                if key not in {"row_name", "deployment_code", "pi_code", "country", "process", "public", "description", "gts", "start_date", "end_date", "task_done"}
            },
        )

    platform_records = _load_platform_json_records(config)
    start_by_deployment: dict[str, str] = {}
    end_by_deployment: dict[str, str] = {}
    for record in platform_records:
        deployment_code = _clean_string(record.get("deployment_code"))
        if not deployment_code:
            continue
        start = _clean_string(record.get("time_coverage_start"))[:10]
        end = _clean_string(record.get("time_coverage_end"))[:10]
        if start and (deployment_code not in start_by_deployment or start < start_by_deployment[deployment_code]):
            start_by_deployment[deployment_code] = start
        if end and (deployment_code not in end_by_deployment or end > end_by_deployment[deployment_code]):
            end_by_deployment[deployment_code] = end

    changed = False
    for record in _load_deployment_json_records(config):
        deployment_code = _clean_string(record.get("deployment_code"))
        if not deployment_code:
            continue
        row = {
            "row_name": deployment_code,
            "deployment_code": deployment_code,
            "pi_code": _clean_string(record.get("pi_code")),
            "process": "1",
            "public": "0",
            "country": "UNKNOWN",
            "task_done": "",
            "first_version": "",
            "last_version": "",
            "start_date": start_by_deployment.get(deployment_code, ""),
            "end_date": end_by_deployment.get(deployment_code, ""),
            "start_date_jul": "",
            "description": _clean_string(record.get("description")),
            "gts": _clean_string(record.get("gts")),
        }
        existing = catalog.get(deployment_code)
        if existing is None:
            ordered_rows.append(row)
            catalog[deployment_code] = DeploymentRecord(
                deployment_code=deployment_code,
                pi_code=row["pi_code"],
                country=row["country"],
                process=row["process"],
                public=row["public"],
                description=row["description"],
                gts=row["gts"],
                start_date=row["start_date"],
                end_date=row["end_date"],
                task_done=row["task_done"],
                source="json",
                row_name=deployment_code,
                extra={
                    "first_version": row["first_version"],
                    "last_version": row["last_version"],
                    "start_date_jul": row["start_date_jul"],
                },
            )
            changed = True
            continue

        updated_row = None
        for candidate in ordered_rows:
            if _clean_string(candidate.get("row_name") or candidate.get("deployment_code")) == deployment_code:
                updated_row = candidate
                break
        if updated_row is None:
            continue
        row_updated = False
        if not updated_row.get("pi_code") and row["pi_code"]:
            updated_row["pi_code"] = row["pi_code"]
            changed = True
            row_updated = True
        if not updated_row.get("description") and row["description"]:
            updated_row["description"] = row["description"]
            changed = True
            row_updated = True
        if not updated_row.get("gts") and row["gts"]:
            updated_row["gts"] = row["gts"]
            changed = True
            row_updated = True
        if not updated_row.get("start_date") and row["start_date"]:
            updated_row["start_date"] = row["start_date"]
            changed = True
            row_updated = True
        if not updated_row.get("end_date") and row["end_date"]:
            updated_row["end_date"] = row["end_date"]
            changed = True
            row_updated = True
        if row_updated:
            catalog[deployment_code] = DeploymentRecord(
                deployment_code=deployment_code,
                pi_code=updated_row.get("pi_code", ""),
                country=updated_row.get("country", ""),
                process=updated_row.get("process", ""),
                public=updated_row.get("public", ""),
                description=updated_row.get("description", ""),
                gts=updated_row.get("gts", ""),
                start_date=updated_row.get("start_date", ""),
                end_date=updated_row.get("end_date", ""),
                task_done=updated_row.get("task_done", ""),
                source="catalog+json",
                row_name=updated_row.get("row_name", deployment_code),
                extra={
                    key: value
                    for key, value in updated_row.items()
                    if key not in {"row_name", "deployment_code", "pi_code", "country", "process", "public", "description", "gts", "start_date", "end_date", "task_done"}
                },
            )

    if persist and ordered_rows:
        canonical = config.catalogdir / "list_deployment.csv"
        write_indexed_csv_rows(canonical, ordered_rows, field_order=DEPLOYMENT_FIELD_ORDER + ("description", "gts"))

    return catalog


def load_hr_catalog(config: MeopConfig, *, persist: bool = True) -> dict[str, dict[str, str]]:
    """Load the high-resolution catalog from ``list_deployment_hr.csv``."""

    rows = _read_hr_rows(config)
    catalog: dict[str, dict[str, str]] = {}
    for row in rows:
        smru_platform_code = _clean_string(row.get("smru_platform_code") or row.get("row_name"))
        if not smru_platform_code:
            continue
        normalized = {key: _clean_string(value) for key, value in row.items()}
        normalized["row_name"] = normalized.get("row_name") or smru_platform_code
        normalized["smru_platform_code"] = smru_platform_code
        catalog[smru_platform_code] = normalized

    if persist and rows:
        canonical = config.catalogdir / "list_deployment_hr.csv"
        write_indexed_csv_rows(canonical, rows, field_order=HR_FIELD_ORDER)

    return catalog


def load_platform_catalog(config: MeopConfig) -> dict[str, tuple[str, ...]]:
    """Return platform codes grouped by deployment using mirrored JSON metadata."""

    grouped: dict[str, set[str]] = {}
    for record in _load_platform_json_records(config):
        smru_platform_code = _clean_string(record.get("smru_platform_code"))
        deployment_code = _clean_string(record.get("deployment_code"))
        if not smru_platform_code or not deployment_code:
            continue
        grouped.setdefault(deployment_code, set()).add(smru_platform_code)
    return {deployment: tuple(sorted(values)) for deployment, values in grouped.items()}


def _list_existing_outputs(config: MeopConfig, selection: Selection, qf: str) -> tuple[Path, ...]:
    return tuple(list_fname_prof(smru_name=selection.smru_name, deployment=selection.deployment, qf=qf, config=config))


def load_info_deployment(config: MeopConfig, deployment: str = "", smru_name: str = "") -> DeploymentInfo:
    """Resolve deployment and tag information from package-managed data.

        It validates the deployment code, resolves where package-managed raw data and outputs are
    expected to live, and inventories any already-produced MEOP netCDF files.
    """

    selection = normalize_selection(deployment=deployment, smru_name=smru_name)
    deployment_catalog = load_deployment_catalog(config)
    hr_catalog = load_hr_catalog(config)
    platform_catalog = load_platform_catalog(config)
    record = deployment_catalog.get(selection.deployment)
    directory = config.final_dataset_dir / selection.deployment

    list_tag_lr0 = _list_existing_outputs(config, selection, "lr0")
    list_tag_lr1 = _list_existing_outputs(config, selection, "lr1")
    list_tag_hr0 = _list_existing_outputs(config, selection, "hr0")
    list_tag_hr1 = _list_existing_outputs(config, selection, "hr1")
    list_tag_hr2 = _list_existing_outputs(config, selection, "hr2")
    list_tag_fr0 = _list_existing_outputs(config, selection, "fr0")
    list_tag_fr1 = _list_existing_outputs(config, selection, "fr1")

    known_platform_codes = platform_catalog.get(selection.deployment, ())
    hr_platform_codes = tuple(sorted(code for code in hr_catalog if deployment_from_smru_name(code) == selection.deployment))

    raw_files = discover_raw_odv_files(config, selection.deployment)
    raw_index = build_odv_profile_index(raw_files)

    if list_tag_lr0:
        list_smru_name = tuple(path.name.split("_")[0] for path in list_tag_lr0)
    elif selection.smru_name:
        list_smru_name = (selection.smru_name,)
    elif raw_index.smru_names:
        list_smru_name = raw_index.smru_names
    elif known_platform_codes:
        list_smru_name = known_platform_codes
    else:
        list_smru_name = ()

    return DeploymentInfo(
        selection=selection,
        record=record,
        invalid_code=record is None or not selection.deployment,
        directory=directory,
        raw_input_dir=config.raw_odv_dir,
        raw_input_zip=config.raw_odv_dir / f"{selection.deployment}_ODV.zip",
        raw_working_text=raw_files.preferred_ctd_text or (config.raw_odv_dir / f"{selection.deployment}_ODV.txt"),
        raw_working_ctd_text=raw_files.ctd_text,
        raw_working_fl_text=raw_files.fl_text,
        raw_working_fcell=config.raw_odv_dir / f"{selection.deployment}_fcell.mat",
        raw_smru_names=raw_index.smru_names,
        raw_profile_count_by_smru=raw_index.profile_count_by_smru,
        list_smru_name=list_smru_name,
        list_tag_lr0=list_tag_lr0,
        list_tag_lr1=list_tag_lr1,
        list_tag_hr0=list_tag_hr0,
        list_tag_hr1=list_tag_hr1,
        list_tag_hr2=list_tag_hr2,
        list_tag_fr0=list_tag_fr0,
        list_tag_fr1=list_tag_fr1,
        known_platform_codes=known_platform_codes,
        hr_platform_codes=hr_platform_codes,
    )
