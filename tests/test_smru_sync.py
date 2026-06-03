from __future__ import annotations

import csv
import json
import zipfile
from pathlib import Path

from meop_process.data.smru_sync import sync_smru_data


def _write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def _write_zip(path: Path, members: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w") as archive:
        for name, content in members.items():
            archive.writestr(name, content)


def _seed_source_deployment(source: Path, deployment: str, *, odv: bool = True, mdb: bool = False) -> Path:
    root = source / "all" / deployment
    root.mkdir(parents=True, exist_ok=True)
    _write_json(
        root / "deployment3.json",
        [{"deployment_code": deployment, "pi_code": "PI", "description": f"{deployment} desc", "gts": "Y"}],
    )
    _write_json(
        root / "deployment2.json",
        [{"deployment_code": deployment, "pi_code": "PI", "description": f"{deployment} desc", "gts": "Y"}],
    )
    _write_json(
        root / "platform3.json",
        [
            {
                "deployment_code": deployment,
                "smru_platform_code": f"{deployment}-AAA",
                "instr_id": "42",
                "time_coverage_start": "2026-01-01T00:00:00Z",
                "time_coverage_end": "2026-02-01T00:00:00Z",
            }
        ],
    )
    _write_json(
        root / "platform2.json",
        [
            {
                "deployment_code": deployment,
                "smru_platform_code": f"{deployment}-AAA",
                "instr_id": "42",
                "time_coverage_start": "2026-01-01T00:00:00Z",
                "time_coverage_end": "2026-02-01T00:00:00Z",
            }
        ],
    )
    if odv:
        _write_zip(root / f"{deployment}_ODV.zip", {f"{deployment}_ODV.txt": "odv"})
    if mdb:
        _write_zip(root / f"{deployment}.zip", {f"{deployment}.mdb": "mdb"})
    return root


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def test_sync_smru_data_dry_run_does_not_mutate_runtime(meop_config, tmp_path: Path) -> None:
    source = tmp_path / "source"
    _seed_source_deployment(source, "DEPNEW", odv=True)

    result = sync_smru_data(meop_config, source_dir=source, apply=False)

    assert result.new_processable_deployments == ("DEPNEW",)
    assert result.copied_odv == (meop_config.raw_odv_dir / "DEPNEW_ODV.zip",)
    assert not (meop_config.raw_odv_dir / "DEPNEW_ODV.zip").exists()
    assert not (meop_config.catalogdir / "list_deployment.csv").exists()
    assert "No files were changed." in result.format_markdown()


def test_sync_smru_data_apply_adds_processable_and_disabled_mdb_deployments(meop_config, tmp_path: Path) -> None:
    source = tmp_path / "source"
    _seed_source_deployment(source, "DEPNEW", odv=True)
    _seed_source_deployment(source, "DEPMDB", odv=False, mdb=True)

    meop_config.catalogdir.mkdir(parents=True, exist_ok=True)
    (meop_config.catalogdir / "list_deployment.csv").write_text(
        "deployment_code,pi_code,process,public,country,task_done,first_version,last_version,start_date,end_date,start_date_jul\n"
        "OLD,OLDPI,1,1,SE,,,,2020-01-01,2020-02-01,\n",
        encoding="utf-8",
    )
    meop_config.config_files_dir.mkdir(parents=True, exist_ok=True)
    _write_json(meop_config.config_files_dir / "deployment3.json", [{"deployment_code": "OLD", "description": "keep"}])
    _write_json(meop_config.config_files_dir / "deployment2.json", [{"deployment_code": "OLD", "description": "keep"}])
    _write_json(meop_config.config_files_dir / "platform3.json", [{"smru_platform_code": "OLD-AAA", "deployment_code": "OLD"}])
    _write_json(meop_config.config_files_dir / "platform2.json", [{"smru_platform_code": "OLD-AAA", "deployment_code": "OLD"}])

    result = sync_smru_data(meop_config, source_dir=source, apply=True, timestamp="20260101T000000")

    assert result.new_processable_deployments == ("DEPNEW",)
    assert result.metadata_only_deployments == ("DEPMDB",)
    assert (meop_config.raw_odv_dir / "DEPNEW_ODV.zip").exists()

    rows = {row["deployment_code"]: row for row in _read_rows(meop_config.catalogdir / "list_deployment.csv")}
    assert rows["OLD"]["process"] == "1"
    assert rows["DEPNEW"]["process"] == "1"
    assert rows["DEPNEW"]["public"] == "0"
    assert rows["DEPMDB"]["process"] == "0"
    assert rows["DEPMDB"]["public"] == "0"

    deployment_records = json.loads((meop_config.config_files_dir / "deployment3.json").read_text(encoding="utf-8"))
    assert {record["deployment_code"] for record in deployment_records} == {"OLD", "DEPNEW", "DEPMDB"}


def test_sync_smru_data_replaces_changed_odv_with_backup(meop_config, tmp_path: Path) -> None:
    source = tmp_path / "source"
    _seed_source_deployment(source, "DEP001", odv=True)

    meop_config.raw_odv_dir.mkdir(parents=True, exist_ok=True)
    _write_zip(meop_config.raw_odv_dir / "DEP001_ODV.zip", {"DEP001_ODV.txt": "old"})
    meop_config.catalogdir.mkdir(parents=True, exist_ok=True)
    (meop_config.catalogdir / "list_deployment.csv").write_text(
        "deployment_code,pi_code,process,public,country,task_done,first_version,last_version,start_date,end_date,start_date_jul\n"
        "DEP001,PI,1,1,SE,,,,2020-01-01,2020-02-01,\n",
        encoding="utf-8",
    )

    result = sync_smru_data(meop_config, source_dir=source, apply=True, timestamp="20260101T000000")

    backup = meop_config.raw_odv_dir / "backups" / "20260101T000000" / "DEP001_ODV.zip"
    assert result.replaced_odv == (meop_config.raw_odv_dir / "DEP001_ODV.zip",)
    assert backup.exists()
    with zipfile.ZipFile(meop_config.raw_odv_dir / "DEP001_ODV.zip") as archive:
        assert archive.read("DEP001_ODV.txt").decode() == "odv"
    with zipfile.ZipFile(backup) as archive:
        assert archive.read("DEP001_ODV.txt").decode() == "old"


def test_sync_smru_data_leaves_unchanged_existing_odv_alone(meop_config, tmp_path: Path) -> None:
    source = tmp_path / "source"
    _seed_source_deployment(source, "DEP001", odv=True)
    meop_config.raw_odv_dir.mkdir(parents=True, exist_ok=True)
    _write_zip(meop_config.raw_odv_dir / "DEP001_ODV.zip", {"DEP001_ODV.txt": "odv"})

    result = sync_smru_data(meop_config, source_dir=source, apply=False)

    assert result.unchanged_odv == ("DEP001",)
    assert result.replaced_odv == ()
