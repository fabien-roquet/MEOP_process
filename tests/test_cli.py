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
