from __future__ import annotations

import json
from pathlib import Path

from meop_process.config.loader import load_config


def test_load_config_falls_back_to_repo_local_defaults(tmp_path: Path) -> None:
    config = load_config(processdir=tmp_path)
    assert config.processdir == tmp_path
    assert config.datadir == tmp_path / "data"
    assert config.publicdir == tmp_path / "public"


def test_load_config_reads_selected_machine_entry(tmp_path: Path) -> None:
    config_file = tmp_path / "configs.json"
    payload = {
        "version": {"CTDnew": "MEOP-CTD_2030-01-01"},
        "configs": {
            "test_machine": {
                "processdir": str(tmp_path / "checkout"),
                "datadir": str(tmp_path / "external_data"),
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
    assert config.publicdir == tmp_path / "external_public"
    assert config.pdflatex == "xelatex"


def test_load_config_reads_runtime_defaults_and_notifications(tmp_path: Path) -> None:
    config_file = tmp_path / "configs.json"
    payload = {
        "defaults": {
            "diagnostics": {"qf": "hr1", "adjusted": False, "parts": ["overview"]},
            "batch": {"jobs": 4, "verbose": True},
            "notifications": {
                "email": {
                    "enabled": True,
                    "when": "failure",
                    "to": ["ops@example.org"],
                    "attach": ["summary_md", "summary_csv"],
                    "subject_prefix": "[MEOP TEST]",
                    "smtp": {
                        "host": "smtp.example.org",
                        "port": 2525,
                        "starttls": False,
                        "use_ssl": True,
                        "username_env": "SMTP_USER",
                        "password_env": "SMTP_PASS",
                        "from": "meop@example.org",
                    },
                }
            },
        },
        "configs": {
            "test_machine": {
                "processdir": str(tmp_path / "checkout"),
                "datadir": str(tmp_path / "data_runtime"),
                "public": str(tmp_path / "public_runtime"),
            }
        },
    }
    config_file.write_text(json.dumps(payload), encoding="utf-8")

    config = load_config(processdir=tmp_path, config_file=config_file, machine="test_machine")

    assert config.diagnostics_defaults.qf == "hr1"
    assert config.diagnostics_defaults.adjusted is False
    assert config.diagnostics_defaults.parts == ("overview",)
    assert config.batch_defaults.jobs == 4
    assert config.batch_defaults.verbose is True
    assert config.email_notifications.enabled is True
    assert config.email_notifications.when == "failure"
    assert config.email_notifications.to == ("ops@example.org",)
    assert config.email_notifications.attach == ("summary_md", "summary_csv")
    assert config.email_notifications.subject_prefix == "[MEOP TEST]"
    assert config.email_notifications.transport.host == "smtp.example.org"
    assert config.email_notifications.transport.port == 2525
    assert config.email_notifications.transport.use_ssl is True
    assert config.email_notifications.transport.from_address == "meop@example.org"
