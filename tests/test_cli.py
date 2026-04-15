from __future__ import annotations

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
