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
    publicdir = processdir / "public"
    for path in (processdir, datadir, publicdir):
        path.mkdir(parents=True, exist_ok=True)
    for path in (
        datadir / "tables",
        datadir / "catalog",
        datadir / "data_raw" / "config_files",
        datadir / "data_raw" / "raw_smru_data_odv",
        datadir / "data_raw" / "raw_smru_hr_data",
        datadir / "data_raw" / "crawl_locations",
        datadir / "data_raw" / "smooth_cls_locations",
    ):
        path.mkdir(parents=True, exist_ok=True)
    return MeopConfig(
        processdir=processdir,
        datadir=datadir,
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


def _seed_reference_catalogs(
    meop_config: MeopConfig,
    *,
    deployment_rows: list[dict[str, str]],
    hr_rows: list[dict[str, str]],
) -> None:
    """Seed the exact catalog rows needed by a reference fixture.

    Processing regressions must not depend on the packaged sample catalogs. Those packaged CSVs
    are useful for data-layout/bootstrap tests, but the FR/HR tests should stage the authoritative
    `list_deployment.csv` and `list_deployment_hr.csv` rows explicitly.
    """

    meop_config.catalogdir.mkdir(parents=True, exist_ok=True)
    write_indexed_csv_rows(meop_config.catalogdir / "list_deployment.csv", deployment_rows)
    write_indexed_csv_rows(meop_config.catalogdir / "list_deployment_hr.csv", hr_rows)


@pytest.fixture()
def stage_ct88_example(meop_config: MeopConfig):
    fixtures_root = Path(__file__).resolve().parent / "fixtures" / "ct88"

    def _stage() -> dict[str, Path]:
        deployment = "ct88"
        meop_config.raw_odv_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(fixtures_root / "ct88_ODV.zip", meop_config.raw_odv_dir / "ct88_ODV.zip")

        _seed_reference_catalogs(
            meop_config,
            deployment_rows=[
                {
                    "row_name": deployment,
                    "deployment_code": deployment,
                    "pi_code": "COSTA",
                    "process": "1",
                    "public": "1",
                    "country": "USA",
                    "task_done": "",
                    "first_version": "",
                    "last_version": "",
                    "start_date": "2012-02-05",
                    "end_date": "2012-07-01",
                    "start_date_jul": "",
                    "description": "Weddell Seal CTD Jan 2012",
                    "gts": "Y",
                }
            ],
            hr_rows=[
                {
                    "row_name": "ct88-225-12",
                    "smru_platform_code": "ct88-225-12",
                    "instr_id": "12225",
                    "year": "2012",
                    "prefix": "",
                    "continuous": "1",
                }
            ],
        )

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

        _seed_reference_catalogs(
            meop_config,
            deployment_rows=[
                {
                    "row_name": deployment,
                    "deployment_code": deployment,
                    "pi_code": "IMOS",
                    "process": "1",
                    "public": "1",
                    "country": "AUSTRALIA",
                    "task_done": "",
                    "first_version": "",
                    "last_version": "",
                    "start_date": "2012-01-01",
                    "end_date": "2012-11-01",
                    "start_date_jul": "",
                    "description": "Southern elephant seal CTD 2012",
                    "gts": "Y",
                }
            ],
            hr_rows=[
                {
                    "row_name": "ct78-465-12",
                    "smru_platform_code": "ct78-465-12",
                    "instr_id": "11465",
                    "year": "2012",
                    "prefix": "",
                    "continuous": "1",
                }
            ],
        )

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

        _seed_reference_catalogs(
            meop_config,
            deployment_rows=[
                {
                    "row_name": deployment,
                    "deployment_code": deployment,
                    "pi_code": "IMOS",
                    "process": "1",
                    "public": "1",
                    "country": "AUSTRALIA",
                    "task_done": "",
                    "first_version": "",
                    "last_version": "",
                    "start_date": "2013-02-05",
                    "end_date": "2014-02-01",
                    "start_date_jul": "",
                    "description": "Southern elephant seal CTD 2013",
                    "gts": "Y",
                }
            ],
            hr_rows=[
                {
                    "row_name": "ct96-24-13",
                    "smru_platform_code": "ct96-24-13",
                    "instr_id": "12664",
                    "year": "2013",
                    "prefix": "",
                    "continuous": "1",
                }
            ],
        )

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


@pytest.fixture()
def stage_ct160_oxy_example(meop_config: MeopConfig):
    """ct160-Oxy726-20: CTD + fluorescence + oxygen (real CHLA and DOXY values)."""
    fixtures_root = Path(__file__).resolve().parent / "fixtures" / "ct160"

    def _stage() -> dict[str, Path]:
        deployment = "ct160"
        smru_name = "ct160-Oxy726-20"
        meop_config.raw_odv_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(fixtures_root / "ct160_ODV.zip", meop_config.raw_odv_dir / "ct160_ODV.zip")

        _seed_reference_catalogs(
            meop_config,
            deployment_rows=[
                {
                    "row_name": deployment,
                    "deployment_code": deployment,
                    "pi_code": "IMOS",
                    "process": "1",
                    "public": "1",
                    "country": "AUSTRALIA",
                    "task_done": "",
                    "first_version": "",
                    "last_version": "",
                    "start_date": "2020-12-20",
                    "end_date": "2023-05-01",
                    "start_date_jul": "",
                    "description": "IMOS Kerguelen late 2020",
                    "gts": "Y",
                }
            ],
            hr_rows=[
                {
                    "row_name": smru_name,
                    "smru_platform_code": smru_name,
                    "instr_id": "14726",
                    "year": "2021",
                    "prefix": "",
                    "continuous": "0",
                }
            ],
        )

        _seed_config_jsons(
            meop_config,
            deployment_payload=[
                {
                    "deployment_code": "ct160",
                    "description": "IMOS Kerguelen late 2020",
                    "pi_code": "IMOS    ",
                    "gts": "Y",
                    "dt_created": "2020-06-22T11:55:59Z",
                    "dt_modified": "2025-06-26T11:55:39Z",
                }
            ],
            platform_payload=[
                {
                    "platform_code": "93486",
                    "wmo_platform_code": "1446",
                    "smru_platform_code": smru_name,
                    "deployment_code": "ct160",
                    "species": "Southern elephant seal",
                    "time_coverage_start": "2021-01-26T00:00:00Z",
                    "time_coverage_end": "2021-04-30T00:00:00Z",
                    "location": "Kerguelen",
                    "firmware_version": "222",
                    "firmware_parameters": "OXY_CONT_20A",
                    "instr_id": "14726",
                    "ptt": "120372",
                    "loc_algorithm": "K",
                    "dt_created": "2020-08-14T00:38:05Z",
                    "dt_modified": "2024-03-14T18:32:33Z",
                }
            ],
        )

        _seed_table_meta(meop_config, f"smru_platform_code,location\n{smru_name},ct160\n")

        return {
            "deployment": Path(deployment),
            "smru_name": smru_name,
            "zip": meop_config.raw_odv_dir / "ct160_ODV.zip",
            "reference_lr0": fixtures_root / f"{smru_name}_lr0_prof_shortened.nc",
            "reference_lr1": fixtures_root / f"{smru_name}_lr1_prof_shortened.nc",
            "reference_hr2": fixtures_root / f"{smru_name}_hr2_prof_shortened.nc",
        }

    return _stage


@pytest.fixture()
def stage_ft11_chla_example(meop_config: MeopConfig):
    """ft11-Cy07b-12: CTD + fluorescence (real CHLA values, no DOXY)."""
    fixtures_root = Path(__file__).resolve().parent / "fixtures" / "ft11"

    def _stage() -> dict[str, Path]:
        deployment = "ft11"
        smru_name = "ft11-Cy07b-12"
        meop_config.raw_odv_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(fixtures_root / "ft11_ODV.zip", meop_config.raw_odv_dir / "ft11_ODV.zip")

        _seed_reference_catalogs(
            meop_config,
            deployment_rows=[
                {
                    "row_name": deployment,
                    "deployment_code": deployment,
                    "pi_code": "GUINET",
                    "process": "1",
                    "public": "1",
                    "country": "FRANCE",
                    "task_done": "",
                    "first_version": "",
                    "last_version": "",
                    "start_date": "2012-02-05",
                    "end_date": "2012-10-01",
                    "start_date_jul": "",
                    "description": "Fluoro 2012 rebatteried",
                    "gts": "N",
                }
            ],
            hr_rows=[
                {
                    "row_name": smru_name,
                    "smru_platform_code": smru_name,
                    "instr_id": "11404",
                    "year": "2012",
                    "prefix": "",
                    "continuous": "0",
                }
            ],
        )

        _seed_config_jsons(
            meop_config,
            deployment_payload=[
                {
                    "deployment_code": "ft11",
                    "description": "Fluoro 2012 rebatteried",
                    "pi_code": "XTOPHE  ",
                    "gts": "N",
                    "dt_created": "2012-04-20T10:31:17Z",
                    "dt_modified": "2020-07-15T10:28:56Z",
                }
            ],
            platform_payload=[
                {
                    "platform_code": "26658",
                    "wmo_platform_code": "530",
                    "smru_platform_code": smru_name,
                    "deployment_code": "ft11",
                    "species": "Southern elephant seal",
                    "time_coverage_start": "2012-02-05T00:00:00Z",
                    "time_coverage_end": "2012-10-01T00:00:00Z",
                    "location": "Kerguelen",
                    "firmware_version": "94",
                    "firmware_parameters": "FTD_10A",
                    "instr_id": "11404",
                    "ptt": "49771",
                    "loc_algorithm": "K",
                    "dt_created": "2012-04-20T10:34:52Z",
                    "dt_modified": "2020-05-27T14:13:20Z",
                }
            ],
        )

        _seed_table_meta(meop_config, f"smru_platform_code,location\n{smru_name},ft11\n")

        return {
            "deployment": Path(deployment),
            "smru_name": smru_name,
            "zip": meop_config.raw_odv_dir / "ft11_ODV.zip",
            "reference_lr0": fixtures_root / f"{smru_name}_lr0_prof_shortened.nc",
            "reference_lr1": fixtures_root / f"{smru_name}_lr1_prof_shortened.nc",
            "reference_hr2": fixtures_root / f"{smru_name}_hr2_prof_shortened.nc",
        }

    return _stage


@pytest.fixture()
def stage_ct153_lr_example(meop_config: MeopConfig):
    """ct153-Pendragon-19: LR-only tag (no hr2), Weddell seal, density removal."""
    fixtures_root = Path(__file__).resolve().parent / "fixtures" / "ct153"

    def _stage() -> dict[str, Path]:
        deployment = "ct153"
        smru_name = "ct153-Pendragon-19"
        meop_config.raw_odv_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(fixtures_root / "ct153_ODV.zip", meop_config.raw_odv_dir / "ct153_ODV.zip")

        _seed_reference_catalogs(
            meop_config,
            deployment_rows=[
                {
                    "row_name": deployment,
                    "deployment_code": deployment,
                    "pi_code": "LARS",
                    "process": "1",
                    "public": "1",
                    "country": "UK",
                    "task_done": "",
                    "first_version": "",
                    "last_version": "",
                    "start_date": "2020-02-01",
                    "end_date": "2020-11-01",
                    "start_date_jul": "",
                    "description": "TARSAN CTD 2019",
                    "gts": "Y",
                }
            ],
            hr_rows=[
                {
                    "row_name": smru_name,
                    "smru_platform_code": smru_name,
                    "instr_id": "15088",
                    "year": "2020",
                    "prefix": "",
                    "continuous": "0",
                }
            ],
        )

        _seed_config_jsons(
            meop_config,
            deployment_payload=[
                {
                    "deployment_code": "ct153",
                    "description": "TARSAN CTD 2019",
                    "pi_code": "LARS    ",
                    "gts": "Y",
                    "dt_created": "2019-07-23T15:59:34Z",
                    "dt_modified": "2021-09-28T17:31:19Z",
                }
            ],
            platform_payload=[
                {
                    "platform_code": "87048",
                    "wmo_platform_code": "1349",
                    "smru_platform_code": smru_name,
                    "deployment_code": "ct153",
                    "species": "Weddell seal",
                    "time_coverage_start": "2020-02-01T00:00:00Z",
                    "time_coverage_end": "2020-11-01T00:00:00Z",
                    "location": "Edwards Islands",
                    "firmware_version": "216",
                    "firmware_parameters": "CTD_GEN_18C",
                    "instr_id": "15088",
                    "ptt": "120359",
                    "loc_algorithm": "K",
                    "dt_created": "2019-07-23T16:04:10Z",
                    "dt_modified": "2021-07-27T15:39:59Z",
                }
            ],
        )

        _seed_table_meta(meop_config, f"smru_platform_code,location\n{smru_name},ct153\n")

        return {
            "deployment": Path(deployment),
            "smru_name": smru_name,
            "zip": meop_config.raw_odv_dir / "ct153_ODV.zip",
            "reference_lr0": fixtures_root / f"{smru_name}_lr0_prof_shortened.nc",
            "reference_lr1": fixtures_root / f"{smru_name}_lr1_prof_shortened.nc",
        }

    return _stage


@pytest.fixture()
def stage_ct107_split_example(meop_config: MeopConfig):
    """ct107-938-13-N2: split tag (-N2 suffix), northern elephant seal."""
    fixtures_root = Path(__file__).resolve().parent / "fixtures" / "ct107"

    def _stage() -> dict[str, Path]:
        deployment = "ct107"
        smru_name = "ct107-938-13-N2"
        base_name = "ct107-938-13"
        meop_config.raw_odv_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(fixtures_root / "ct107_ODV.zip", meop_config.raw_odv_dir / "ct107_ODV.zip")

        _seed_reference_catalogs(
            meop_config,
            deployment_rows=[
                {
                    "row_name": deployment,
                    "deployment_code": deployment,
                    "pi_code": "COSTA",
                    "process": "1",
                    "public": "1",
                    "country": "USA",
                    "task_done": "",
                    "first_version": "",
                    "last_version": "",
                    "start_date": "2014-01-28",
                    "end_date": "2015-01-01",
                    "start_date_jul": "",
                    "description": "Costa northern ellie CTD 2014",
                    "gts": "Y",
                }
            ],
            hr_rows=[
                {
                    "row_name": smru_name,
                    "smru_platform_code": smru_name,
                    "instr_id": "12938",
                    "year": "2014",
                    "prefix": "",
                    "continuous": "0",
                }
            ],
        )

        _seed_config_jsons(
            meop_config,
            deployment_payload=[
                {
                    "deployment_code": "ct107",
                    "description": "Costa northern ellie CTD 2014",
                    "pi_code": "COSTA   ",
                    "gts": "Y",
                    "dt_created": "2013-12-03T10:00:21Z",
                    "dt_modified": "2023-04-21T11:15:18Z",
                }
            ],
            platform_payload=[
                {
                    "platform_code": "44898",
                    "wmo_platform_code": "657",
                    "smru_platform_code": base_name,
                    "deployment_code": "ct107",
                    "species": "Northern elephant seal",
                    "time_coverage_start": "2014-01-28T00:00:00Z",
                    "time_coverage_end": "2015-01-01T00:00:00Z",
                    "location": "California",
                    "firmware_version": "119",
                    "firmware_parameters": "CTD_GEN_13B",
                    "instr_id": "12938",
                    "ptt": "133773",
                    "loc_algorithm": "L",
                    "dt_created": "2013-12-03T10:03:03Z",
                    "dt_modified": "2020-05-27T14:13:20Z",
                }
            ],
        )

        # Seed the split-tags table so the processor knows ct107-938-13 is split into 2
        meop_config.tablesdir.mkdir(parents=True, exist_ok=True)
        (meop_config.tablesdir / "table_split_tags.csv").write_text(
            "smru_platform_name,nsplit\nct107-938-13,2\n", encoding="utf-8"
        )
        _seed_table_meta(meop_config, f"smru_platform_code,location\n{smru_name},ct107\n")

        return {
            "deployment": Path(deployment),
            "smru_name": smru_name,
            "base_name": base_name,
            "zip": meop_config.raw_odv_dir / "ct107_ODV.zip",
            "reference_lr0": fixtures_root / f"{smru_name}_lr0_prof_shortened.nc",
            "reference_lr1": fixtures_root / f"{smru_name}_lr1_prof_shortened.nc",
            "reference_hr2": fixtures_root / f"{smru_name}_hr2_prof_shortened.nc",
        }

    return _stage
