from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from meop_process.catalog.tables import write_indexed_csv_rows
from meop_process.data.layout import bootstrap_packaged_catalogs
from meop_process.models import MeopConfig


@pytest.fixture()
def meop_config(tmp_path: Path) -> MeopConfig:
    processdir = tmp_path / "repo"
    datadir = processdir / "data"
    refdir = processdir / "references"
    publicdir = processdir / "public"
    for path in (processdir, datadir, refdir, publicdir):
        path.mkdir(parents=True, exist_ok=True)
    for path in (datadir / "tables", datadir / "catalog", datadir / "config_files", datadir / "raw_smru_data_odv", datadir / "raw_smru_hr_data"):
        path.mkdir(parents=True, exist_ok=True)
    return MeopConfig(
        processdir=processdir,
        datadir=datadir,
        refdir=refdir,
        publicdir=publicdir,
        version="MEOP-CTD_2099-01-01",
        machine="test_machine",
    )


@pytest.fixture()
def seed_catalog(meop_config: MeopConfig):
    def _seed(*, deployment: str = "DEP001", smru_name: str = "DEP001-AAA") -> None:
        meop_config.catalogdir.mkdir(parents=True, exist_ok=True)
        write_indexed_csv_rows(
            meop_config.catalogdir / "list_deployment.csv",
            [
                {
                    "row_name": deployment,
                    "deployment_code": deployment,
                    "pi_code": "PI001",
                    "process": "1",
                    "public": "1",
                    "country": "SE",
                    "task_done": "",
                    "first_version": "",
                    "last_version": "",
                    "start_date": "2020-01-01",
                    "end_date": "2020-12-31",
                    "start_date_jul": "",
                }
            ],
        )
        write_indexed_csv_rows(
            meop_config.catalogdir / "list_deployment_hr.csv",
            [
                {
                    "row_name": smru_name,
                    "smru_platform_code": smru_name,
                    "instr_id": "42",
                    "year": "2020",
                    "prefix": "",
                    "continuous": "0",
                }
            ],
        )
    return _seed


def _seed_config_jsons(meop_config: MeopConfig, *, deployment_payload: list[dict], platform_payload: list[dict]) -> None:
    meop_config.config_files_dir.mkdir(parents=True, exist_ok=True)
    for name in ("deployment3.json", "platform3.json", "deployment2_patch.json", "platform2_patch.json"):
        (meop_config.config_files_dir / name).write_text("[]", encoding="utf-8")
    (meop_config.config_files_dir / "deployment2.json").write_text(json.dumps(deployment_payload), encoding="utf-8")
    (meop_config.config_files_dir / "platform2.json").write_text(json.dumps(platform_payload), encoding="utf-8")


def _seed_table_meta(meop_config: MeopConfig, rows: str) -> None:
    meop_config.tablesdir.mkdir(parents=True, exist_ok=True)
    (meop_config.tablesdir / "table_meta.csv").write_text(rows, encoding="utf-8")


@pytest.fixture()
def stage_ct88_example(meop_config: MeopConfig):
    fixtures_root = Path(__file__).resolve().parent / "fixtures" / "ct88"

    def _stage() -> dict[str, Path]:
        deployment = "ct88"
        meop_config.raw_odv_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(fixtures_root / "ct88_ODV.zip", meop_config.raw_odv_dir / "ct88_ODV.zip")

        bootstrap_packaged_catalogs(meop_config)

        _seed_config_jsons(
            meop_config,
            deployment_payload=[
                {
                    "deployment_code": "ct88",
                    "description": "Weddell Seal CTD Jan 2012",
                    "pi_code": "COSTA   ",
                    "gts": "Y",
                    "dt_created": "2012-04-20T10:31:17Z",
                    "dt_modified": "2022-01-13T20:53:27Z",
                }
            ],
            platform_payload=[
                {
                    "platform_code": "26172",
                    "wmo_platform_code": "454",
                    "smru_platform_code": "ct88-225-12",
                    "deployment_code": "ct88",
                    "species": "Weddell seal",
                    "time_coverage_start": "2012-02-05T00:00:00Z",
                    "time_coverage_end": "2012-07-01T00:00:00Z",
                    "location": "Ross Sea McMurdo",
                    "firmware_version": "94",
                    "firmware_parameters": "CTD_GEN_11A",
                    "instr_id": "12225",
                    "ptt": "113863",
                    "loc_algorithm": "L",
                    "dt_created": "2012-04-20T10:34:52Z",
                    "dt_modified": "2020-05-27T14:13:20Z",
                }
            ],
        )

        _seed_table_meta(meop_config, "smru_platform_code,location\nct88-225-12,ct88\n")

        return {
            "deployment": Path(deployment),
            "reference_lr0": fixtures_root / "ct88-225-12_lr0_prof.nc",
            "reference_lr1": fixtures_root / "ct88-225-12_lr1_prof.nc",
            "zip": meop_config.raw_odv_dir / "ct88_ODV.zip",
        }

    return _stage


@pytest.fixture()
def stage_ct78_example(meop_config: MeopConfig):
    fixtures_root = Path(__file__).resolve().parent / "fixtures" / "ct78"

    def _stage() -> dict[str, Path]:
        deployment = "ct78"
        meop_config.raw_odv_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(fixtures_root / "ct78_ODV.zip", meop_config.raw_odv_dir / "ct78_ODV.zip")

        bootstrap_packaged_catalogs(meop_config)

        _seed_config_jsons(
            meop_config,
            deployment_payload=[
                {
                    "deployment_code": "ct78",
                    "description": "Southern elephant seal CTD 2012",
                    "pi_code": "IMOS    ",
                    "gts": "Y",
                    "dt_created": "2012-04-20T10:31:17Z",
                    "dt_modified": "2024-03-05T18:06:07Z",
                }
            ],
            platform_payload=[
                {
                    "platform_code": "23880",
                    "wmo_platform_code": "Q9900512",
                    "smru_platform_code": "ct78-465-12",
                    "deployment_code": "ct78",
                    "species": "Southern ellie",
                    "time_coverage_start": "2012-01-01T00:00:00Z",
                    "time_coverage_end": "2012-11-01T00:00:00Z",
                    "location": "Heard Island",
                    "firmware_version": "82",
                    "firmware_parameters": "CTD_GEN_10A",
                    "instr_id": "11465",
                    "ptt": "54979",
                    "loc_algorithm": "L",
                    "dt_created": "2012-04-20T10:34:52Z",
                    "dt_modified": "2024-03-05T18:06:07Z",
                }
            ],
        )

        _seed_table_meta(meop_config, "smru_platform_code,location\nct78-465-12,ct78\n")

        return {
            "deployment": Path(deployment),
            "reference_lr0": fixtures_root / "ct78-465-12_lr0_prof.nc",
            "reference_lr1": fixtures_root / "ct78-465-12_lr1_prof.nc",
            "reference_hr2": fixtures_root / "ct78-465-12_hr2_prof_shortened.nc",
            "zip": meop_config.raw_odv_dir / "ct78_ODV.zip",
        }

    return _stage


@pytest.fixture()
def stage_ct96_example(meop_config: MeopConfig):
    fixtures_root = Path(__file__).resolve().parent / "fixtures" / "ct96"

    def _stage() -> dict[str, Path]:
        deployment = "ct96"
        meop_config.raw_odv_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(fixtures_root / "ct96_ODV.zip", meop_config.raw_odv_dir / "ct96_ODV.zip")
        hr_dir = meop_config.raw_hr_dir / "2013"
        hr_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(fixtures_root / "12664_ctd_shortened.txt", hr_dir / "12664_ctd.txt")

        bootstrap_packaged_catalogs(meop_config)

        _seed_config_jsons(
            meop_config,
            deployment_payload=[
                {
                    "deployment_code": "ct96",
                    "description": "Southern elephant seal CTD 2013",
                    "pi_code": "IMOS    ",
                    "gts": "Y",
                    "dt_created": "2012-04-20T10:31:17Z",
                    "dt_modified": "2024-03-05T18:32:44Z",
                }
            ],
            platform_payload=[
                {
                    "platform_code": "34650",
                    "wmo_platform_code": "Q9900566",
                    "smru_platform_code": "ct96-24-13",
                    "deployment_code": "ct96",
                    "species": "Southern ellie",
                    "time_coverage_start": "2013-02-05T00:00:00Z",
                    "time_coverage_end": "2014-02-01T00:00:00Z",
                    "location": "ct96",
                    "firmware_version": "101",
                    "firmware_parameters": "CTD_GEN_11A",
                    "instr_id": "12664",
                    "ptt": "122378",
                    "loc_algorithm": "K",
                    "dt_created": "2012-04-20T10:34:52Z",
                    "dt_modified": "2024-03-05T18:32:44Z",
                }
            ],
        )

        _seed_table_meta(meop_config, "smru_platform_code,location\nct96-24-13,ct96\n")

        return {
            "deployment": Path(deployment),
            "zip": meop_config.raw_odv_dir / "ct96_ODV.zip",
            "hr_file": meop_config.raw_hr_dir / "2013" / "12664_ctd.txt",
            "reference_fr0": fixtures_root / "ct96-24-13_fr0_prof_shortened.nc",
            "reference_fr1": fixtures_root / "ct96-24-13_fr1_prof_shortened.nc",
            "reference_hr2": fixtures_root / "ct96-24-13_hr2_prof_shortened.nc",
            "reference_traj": fixtures_root / "ct96-24-13_traj_shortened.nc",
        }

    return _stage
