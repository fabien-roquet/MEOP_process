from __future__ import annotations

from meop_process.models import Selection
from meop_process import cli as cli_module


def test_cli_runs_selected_actions(meop_config, seed_catalog, monkeypatch) -> None:
    seed_catalog(deployment="DEP001", smru_name="DEP001-AAA")

    monkeypatch.setattr(cli_module, "process_tags", lambda **kwargs: True)
    monkeypatch.setattr(cli_module, "create_hr2", lambda **kwargs: True)

    exit_code = cli_module.main(
        ["--deployment", "DEP001", "--process_data", "--create_hr2"],
        config=meop_config,
    )

    assert exit_code == 0


def test_cli_utility_mode_does_not_require_deployment_selection(meop_config, capsys) -> None:
    exit_code = cli_module.main(["--show-data-layout", "--bootstrap-data"], config=meop_config)

    captured = capsys.readouterr()
    assert exit_code == 0
    assert "MEOP runtime data layout" in captured.out


def test_cli_allows_overview_only_diagnostics_without_selection(meop_config, monkeypatch) -> None:
    calls: list[tuple[str, str, tuple[str, ...], bool]] = []

    def fake_generate_diagnostics(**kwargs):
        calls.append(
            (
                kwargs.get("deployment", ""),
                kwargs.get("smru_name", ""),
                    tuple(kwargs.get("parts", ())),
                    bool(kwargs.get("adjusted", True)),
            )
        )
        return {"written_files": [], "processed_tags": []}

    monkeypatch.setattr(cli_module, "generate_diagnostics", fake_generate_diagnostics)

    exit_code = cli_module.main(["--diagnostics", "--diagnostics-part", "overview"], config=meop_config)

    assert exit_code == 0
    assert calls == [("", "", ("overview",), True)]


def test_cli_forwards_notification_overrides_to_batch_runner(meop_config, monkeypatch) -> None:
    captured: list[dict[str, object]] = []

    def fake_run_all_deployments(**kwargs):
        captured.append(kwargs)
        return {"summary_markdown": "summary.md", "failed_count": 0}

    monkeypatch.setattr(cli_module, "run_all_deployments", fake_run_all_deployments)

    exit_code = cli_module.main(
        [
            "--run-all-deployments",
            "--diagnostics",
            "--notify-email",
            "ops@example.org",
            "--notify-when",
            "failure",
            "--notify-attach",
            "summary_md",
        ],
        config=meop_config,
    )

    assert exit_code == 0
    assert captured[0]["notify_email"] == ["ops@example.org"]
    assert captured[0]["notify_when"] == "failure"
    assert captured[0]["notify_attach"] == ["summary_md"]


def test_cli_validate_runtime_tables_utility(meop_config, monkeypatch, capsys) -> None:
    monkeypatch.setattr(
        cli_module,
        "validate_runtime_tables",
        lambda config: {
            "ok": True,
            "checked_tables": ["table_coeff.csv"],
            "errors": [],
            "warnings": [],
        },
    )

    exit_code = cli_module.main(["--validate-runtime-tables"], config=meop_config)

    captured = capsys.readouterr()
    assert exit_code == 0
    assert "MEOP runtime table validation" in captured.out
    assert "OK: runtime tables passed validation" in captured.out


def test_cli_sync_smru_data_dry_run_requires_source_and_prints_report(meop_config, monkeypatch, tmp_path, capsys) -> None:
    class FakeSyncResult:
        def format_markdown(self):
            return "# fake sync\n"

    calls: list[tuple[str, bool]] = []

    def fake_sync(config, *, source_dir, apply=False):
        calls.append((source_dir, apply))
        return FakeSyncResult()

    monkeypatch.setattr("meop_process.data.smru_sync.sync_smru_data", fake_sync)

    exit_code = cli_module.main(["--sync-smru-data", "--source-dir", str(tmp_path)], config=meop_config)

    captured = capsys.readouterr()
    assert exit_code == 0
    assert calls == [(str(tmp_path), False)]
    assert "# fake sync" in captured.out
