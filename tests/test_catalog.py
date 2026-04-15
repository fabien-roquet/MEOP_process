from __future__ import annotations

import json

from meop_process.api import load_info_deployment


def test_load_info_deployment_can_bootstrap_from_mirrored_json_when_catalog_csv_is_absent(meop_config) -> None:
    meop_config.config_files_dir.mkdir(parents=True, exist_ok=True)
    (meop_config.config_files_dir / "deployment3.json").write_text(
        json.dumps(
            [
                {
                    "deployment_code": "DEPJSON",
                    "description": "JSON deployment",
                    "pi_code": "PIJSON",
                    "gts": "Y",
                }
            ]
        ),
        encoding="utf-8",
    )
    (meop_config.config_files_dir / "deployment2_patch.json").write_text("[]", encoding="utf-8")
    (meop_config.config_files_dir / "platform3.json").write_text(
        json.dumps(
            [
                {
                    "deployment_code": "DEPJSON",
                    "smru_platform_code": "DEPJSON-AAA",
                    "time_coverage_start": "2021-01-05T00:00:00Z",
                    "time_coverage_end": "2021-02-10T00:00:00Z",
                }
            ]
        ),
        encoding="utf-8",
    )
    (meop_config.config_files_dir / "platform2_patch.json").write_text("[]", encoding="utf-8")

    info = load_info_deployment(deployment="DEPJSON", config=meop_config)

    assert info.invalid_code is False
    assert info.PI == "PIJSON"
    assert info.known_platform_codes == ("DEPJSON-AAA",)
    assert (meop_config.catalogdir / "list_deployment.csv").exists()
    assert (meop_config.config_files_dir / "deployment3.json").exists()
