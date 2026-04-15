from __future__ import annotations

import zipfile

from meop_process.catalog.deployments import load_info_deployment
from meop_process.io.raw_odv import import_raw_data_zip
from meop_process.processing.ncargo import prepare_ncargo_inputs


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
