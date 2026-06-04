from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path

from meop_process.catalog.tables import write_indexed_csv_rows
from meop_process.compare_batch import discover_calibration_tags, main, run_calibration_plots_batch


def _seed_deployment_catalog(meop_config, rows: list[dict[str, str]]) -> None:
    write_indexed_csv_rows(meop_config.catalogdir / "list_deployment.csv", rows)
    write_indexed_csv_rows(meop_config.catalogdir / "list_deployment_hr.csv", [])


def _touch_profile(meop_config, deployment: str, smru_name: str, qf: str = "hr1") -> Path:
    path = meop_config.final_dataset_dir / deployment / f"{smru_name}_{qf}_prof.nc"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("placeholder", encoding="utf-8")
    return path


def test_discover_calibration_tags_uses_enabled_deployments_by_default(meop_config) -> None:
    _seed_deployment_catalog(
        meop_config,
        [
            {
                "row_name": "DEP001",
                "deployment_code": "DEP001",
                "process": "1",
            },
            {
                "row_name": "DEP002",
                "deployment_code": "DEP002",
                "process": "0",
            },
        ],
    )
    _touch_profile(meop_config, "DEP001", "DEP001-AAA", "hr1")
    _touch_profile(meop_config, "DEP001", "DEP001-AAA", "lr1")
    _touch_profile(meop_config, "DEP002", "DEP002-BBB", "hr1")

    assert discover_calibration_tags(meop_config) == ["DEP001-AAA"]
    assert discover_calibration_tags(meop_config, include_disabled=True) == ["DEP001-AAA", "DEP002-BBB"]
    assert discover_calibration_tags(meop_config, tags=["DEP003-CCC"]) == ["DEP003-CCC"]


def test_run_calibration_plots_batch_writes_state_and_skips_success(meop_config, monkeypatch, tmp_path) -> None:
    cfg = replace(meop_config, cora_dir=tmp_path / "cora")
    _seed_deployment_catalog(
        cfg,
        [
            {
                "row_name": "DEP001",
                "deployment_code": "DEP001",
                "process": "1",
            },
        ],
    )
    _touch_profile(cfg, "DEP001", "DEP001-AAA", "hr1")
    _touch_profile(cfg, "DEP001", "DEP001-BBB", "lr1")
    calls: list[str] = []

    def fake_generate(smru_name: str, *, config):
        calls.append(smru_name)
        outdir = config.plotdir / smru_name.split("-")[0]
        outdir.mkdir(parents=True, exist_ok=True)
        written = [
            outdir / f"{smru_name}_calibration.png",
            outdir / f"{smru_name}_calibration_offsets.csv",
            outdir / f"{smru_name}_calibration_offsets.png",
        ]
        for path in written:
            path.write_text("ok", encoding="utf-8")
        return written

    monkeypatch.setattr("meop_process.compare_batch.generate_calibration_plots", fake_generate)

    state_dir = cfg.datadir / "compare_state"
    first = run_calibration_plots_batch(config=cfg, state_dir=state_dir)

    assert first.success_count == 2
    assert first.failed_count == 0
    assert first.summary_markdown.exists()
    assert first.summary_csv.exists()
    assert first.state_path.exists()
    assert calls == ["DEP001-AAA", "DEP001-BBB"]

    calls.clear()
    second = run_calibration_plots_batch(config=cfg, state_dir=state_dir)

    assert second.skipped_count == 2
    assert calls == []

    forced = run_calibration_plots_batch(config=cfg, state_dir=state_dir, force=True)

    assert forced.success_count == 2
    assert calls == ["DEP001-AAA", "DEP001-BBB"]


def test_compare_batch_main_dry_run_lists_selected_tags(meop_config, tmp_path, capsys) -> None:
    cora_dir = tmp_path / "cora"
    cora_dir.mkdir()
    config_json = meop_config.processdir / "configs.json"
    config_json.write_text(
        json.dumps(
            {
                "defaults": {
                    "processdir": str(meop_config.processdir),
                    "datadir": str(meop_config.datadir),
                    "public": str(meop_config.publicdir),
                    "references": {"cora_dir": str(cora_dir)},
                }
            }
        ),
        encoding="utf-8",
    )
    _seed_deployment_catalog(
        meop_config,
        [
            {
                "row_name": "DEP001",
                "deployment_code": "DEP001",
                "process": "1",
            },
        ],
    )
    _touch_profile(meop_config, "DEP001", "DEP001-AAA", "hr1")

    exit_code = main(["--processdir", str(meop_config.processdir), "--dry-run"])

    output = capsys.readouterr().out
    assert exit_code == 0
    assert "Selected tags: 1" in output
    assert "DEP001-AAA" in output
