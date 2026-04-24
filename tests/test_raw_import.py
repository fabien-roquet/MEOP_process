from __future__ import annotations

import zipfile

import numpy as np
import pytest

from meop_process.catalog.deployments import load_info_deployment
from meop_process.io.raw_odv import import_raw_data_zip
from meop_process.io.raw_odv import OdvProfile, discover_raw_odv_files, load_raw_odv_profiles
from meop_process.processing.ncargo import _numeric_matrix, create_ncargo_python, prepare_ncargo_inputs
from meop_process.models import Selection


def test_import_raw_data_extracts_archive_in_place(meop_config) -> None:
    deployment = "DEP001"
    meop_config.raw_odv_dir.mkdir(parents=True, exist_ok=True)
    zip_path = meop_config.raw_odv_dir / f"{deployment}_ODV.zip"

    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("profiles/example.txt", "payload")

    assert import_raw_data_zip(meop_config, deployment)

    extracted_file = meop_config.raw_odv_dir / "profiles" / "example.txt"
    assert zip_path.exists()
    assert extracted_file.read_text(encoding="utf-8") == "payload"


def test_prepare_ncargo_inputs_is_a_noop_in_pure_python(meop_config, seed_catalog) -> None:
    deployment = "DEP001"
    seed_catalog(deployment=deployment, smru_name="DEP001-AAA")
    raw_root = meop_config.raw_odv_dir
    raw_root.mkdir(parents=True, exist_ok=True)
    (raw_root / f"{deployment}_ODV.txt").write_text("ctd", encoding="utf-8")
    (raw_root / f"{deployment}_FL_ODV.txt").write_text("fl", encoding="utf-8")

    info = load_info_deployment(meop_config, deployment=deployment)

    assert prepare_ncargo_inputs(meop_config, info) is True


def test_load_raw_odv_profiles_splits_tag_on_largest_time_gap(meop_config, seed_catalog) -> None:
    deployment = "DEP001"
    smru_name = "DEP001-AAA"
    seed_catalog(deployment=deployment, smru_name=smru_name)
    meop_config.raw_odv_dir.mkdir(parents=True, exist_ok=True)
    (meop_config.tablesdir / "table_split_tags.csv").write_text(
        "smru_platform_name,nsplit\nDEP001-AAA,2\n",
        encoding="utf-8",
    )
    (meop_config.raw_odv_dir / f"{deployment}_ODV.txt").write_text(
        "// comment\n"
        "Cruise;Station;Type;yyyy-mm-dd hh:mm;Longitude;Latitude;Quality;Pressure;Temperature;Salinity\n"
        "DEP001-AAA;1;;2020-01-01 00:00;10;20;1;5;1.0;33.1\n"
        ";;;;;;;10;1.1;33.2\n"
        "DEP001-AAA;2;;2020-01-02 00:00;10;20;1;5;1.2;33.3\n"
        ";;;;;;;10;1.3;33.4\n"
        "DEP001-AAA;3;;2020-03-01 00:00;10;20;1;5;1.4;33.5\n"
        ";;;;;;;10;1.5;33.6\n"
        "DEP001-AAA;4;;2020-03-02 00:00;10;20;1;5;1.6;33.7\n"
        ";;;;;;;10;1.7;33.8\n",
        encoding="utf-8",
    )

    files = discover_raw_odv_files(meop_config, deployment)
    profiles = load_raw_odv_profiles(files, config=meop_config)

    assert [profile.smru_name for profile in profiles] == [
        "DEP001-AAA-N1",
        "DEP001-AAA-N1",
        "DEP001-AAA-N2",
        "DEP001-AAA-N2",
    ]


def test_create_ncargo_python_writes_split_tag_outputs(meop_config, seed_catalog) -> None:
    deployment = "DEP001"
    smru_name = "DEP001-AAA"
    seed_catalog(deployment=deployment, smru_name=smru_name)
    meop_config.raw_odv_dir.mkdir(parents=True, exist_ok=True)
    (meop_config.tablesdir / "table_split_tags.csv").write_text(
        "smru_platform_name,nsplit\nDEP001-AAA,2\n",
        encoding="utf-8",
    )
    (meop_config.raw_odv_dir / f"{deployment}_ODV.txt").write_text(
        "// comment\n"
        "Cruise;Station;Type;yyyy-mm-dd hh:mm;Longitude;Latitude;Quality;Pressure;Temperature;Salinity\n"
        "DEP001-AAA;1;;2020-01-01 00:00;10;20;1;5;1.0;33.1\n"
        ";;;;;;;10;1.1;33.2\n"
        "DEP001-AAA;2;;2020-01-02 00:00;10;20;1;5;1.2;33.3\n"
        ";;;;;;;10;1.3;33.4\n"
        "DEP001-AAA;3;;2020-03-01 00:00;10;20;1;5;1.4;33.5\n"
        ";;;;;;;10;1.5;33.6\n"
        "DEP001-AAA;4;;2020-03-02 00:00;10;20;1;5;1.6;33.7\n"
        ";;;;;;;10;1.7;33.8\n",
        encoding="utf-8",
    )

    result = create_ncargo_python(meop_config, Selection(deployment=deployment))
    info = load_info_deployment(meop_config, deployment=deployment)

    assert result.processed_tags == ("DEP001-AAA-N1", "DEP001-AAA-N2")
    assert info.raw_smru_names == ("DEP001-AAA-N1", "DEP001-AAA-N2")


def test_load_raw_odv_profiles_drops_profiles_with_blank_header_lon_lat(meop_config, seed_catalog) -> None:
    deployment = "DEP001"
    smru_name = "DEP001-AAA"
    seed_catalog(deployment=deployment, smru_name=smru_name)
    meop_config.raw_odv_dir.mkdir(parents=True, exist_ok=True)
    (meop_config.tablesdir / "table_split_tags.csv").write_text(
        "smru_platform_name,nsplit\nDEP001-AAA,2\n",
        encoding="utf-8",
    )
    (meop_config.raw_odv_dir / f"{deployment}_ODV.txt").write_text(
        "// comment\n"
        "Cruise;Station;Type;yyyy-mm-dd hh:mm;Longitude;Latitude;Quality;Pressure;Temperature;Salinity\n"
        "DEP001-AAA;1;;2020-01-01 00:00;;;1;5;1.0;33.1\n"
        ";;;;;;;10;1.1;33.2\n"
        "DEP001-AAA;2;;2020-01-02 00:00;10;20;1;5;1.2;33.3\n"
        ";;;;;;;10;1.3;33.4\n",
        encoding="utf-8",
    )

    files = discover_raw_odv_files(meop_config, deployment)
    profiles = load_raw_odv_profiles(files, config=meop_config)

    assert len(profiles) == 1
    assert profiles[0].station == 2


def test_numeric_matrix_uses_field_specific_max_length() -> None:
    profiles = [
        OdvProfile(
            smru_name="DEP001-AAA",
            station=1,
            timestamp="2020-01-01 00:00",
            longitude=10.0,
            latitude=20.0,
            pressure=tuple(float(i) for i in range(16)),
            temperature=tuple(1.0 + i for i in range(16)),
            salinity=tuple(33.0 + i for i in range(16)),
            fluorescence=tuple(0.1 + i for i in range(16)),
        ),
        OdvProfile(
            smru_name="DEP001-AAA",
            station=2,
            timestamp="2020-01-02 00:00",
            longitude=10.0,
            latitude=20.0,
            pressure=tuple(float(i) for i in range(16)),
            temperature=tuple(1.0 + i for i in range(16)),
            salinity=tuple(33.0 + i for i in range(16)),
            fluorescence=tuple(0.1 + i for i in range(23)),
        ),
    ]

    matrix = _numeric_matrix(profiles, "fluorescence")

    assert matrix.shape == (2, 23)
    assert float(matrix[1, 22]) == pytest.approx(22.1)
    assert np.all(np.isnan(matrix[0, 16:]))
