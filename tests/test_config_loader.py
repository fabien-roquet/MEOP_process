from __future__ import annotations

import json
from pathlib import Path

from meop_process.config.loader import load_config


def test_load_config_falls_back_to_repo_local_defaults(tmp_path: Path) -> None:
    config = load_config(processdir=tmp_path)
    assert config.processdir == tmp_path
    assert config.datadir == tmp_path / "data"
    assert config.refdir == tmp_path / "references"
    assert config.publicdir == tmp_path / "public"


def test_load_config_reads_selected_machine_entry(tmp_path: Path) -> None:
    config_file = tmp_path / "configs.json"
    payload = {
        "version": {"CTDnew": "MEOP-CTD_2030-01-01"},
        "configs": {
            "test_machine": {
                "processdir": str(tmp_path / "checkout"),
                "datadir": str(tmp_path / "external_data"),
                "refdir": str(tmp_path / "external_refs"),
                "public": str(tmp_path / "external_public"),
                "pdflatex": "xelatex",
            }
        },
    }
    config_file.write_text(json.dumps(payload), encoding="utf-8")

    config = load_config(processdir=tmp_path, config_file=config_file, machine="test_machine")

    assert config.version == "MEOP-CTD_2030-01-01"
    assert config.processdir == tmp_path / "checkout"
    assert config.datadir == tmp_path / "external_data"
    assert config.refdir == tmp_path / "external_refs"
    assert config.publicdir == tmp_path / "external_public"
    assert config.pdflatex == "xelatex"
