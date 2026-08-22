from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path

from meop_process.catalog.tables import write_indexed_csv_rows
from meop_process.compare_batch import (
    CALIBRATION_METHOD_VERSION,
    _should_skip,
    discover_calibration_tags,
    main,
    run_calibration_plots_batch,
)
from meop_process.compare_cli import (
    InsufficientReferenceDataError,
    InsufficientTargetDataError,
    InvalidTargetDataError,
    NoReferenceDataError,
)


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
    state = json.loads(first.state_path.read_text(encoding="utf-8"))
    assert {entry["method_version"] for entry in state.values()} == {CALIBRATION_METHOD_VERSION}

    calls.clear()
    second = run_calibration_plots_batch(config=cfg, state_dir=state_dir)

    assert second.skipped_count == 2
    assert calls == []

    forced = run_calibration_plots_batch(config=cfg, state_dir=state_dir, force=True)

    assert forced.success_count == 2
    assert calls == ["DEP001-AAA", "DEP001-BBB"]


def test_old_compare_success_state_is_not_reused(meop_config, tmp_path) -> None:
    cfg = replace(meop_config, cora_dir=tmp_path / "cora")
    _touch_profile(cfg, "DEP001", "DEP001-AAA", "hr1")
    output_dir = cfg.plotdir / "DEP001"
    output_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("calibration.png", "calibration_offsets.csv", "calibration_offsets.png"):
        (output_dir / f"DEP001-AAA_{suffix}").write_text("old", encoding="utf-8")

    skip, reason = _should_skip(
        "DEP001-AAA",
        state={"DEP001-AAA": {"status": "success", "method_version": 3}},
        force=False,
        config=cfg,
    )

    assert skip is False
    assert CALIBRATION_METHOD_VERSION == 4
    assert "older calibration method version" in reason


def test_compare_batch_records_data_availability_statuses(meop_config, monkeypatch, tmp_path) -> None:
    cfg = replace(meop_config, cora_dir=tmp_path / "cora")
    _seed_deployment_catalog(
        cfg,
        [{"row_name": "DEP001", "deployment_code": "DEP001", "process": "1"}],
    )
    for smru_name in ("DEP001-NONE", "DEP001-THIN", "DEP001-SPARSE", "DEP001-BAD"):
        _touch_profile(cfg, "DEP001", smru_name, "hr1")

    thin_output = cfg.plotdir / "DEP001" / "DEP001-THIN_calibration_offsets.csv"

    def fake_generate(smru_name: str, *, config):
        if smru_name.endswith("NONE"):
            raise NoReferenceDataError("no profiles in the requested CORA cells")
        if smru_name.endswith("THIN"):
            thin_output.parent.mkdir(parents=True, exist_ok=True)
            thin_output.write_text("diagnostic", encoding="utf-8")
            raise InsufficientReferenceDataError("no usable deep overlap", written_files=[thin_output])
        if smru_name.endswith("SPARSE"):
            raise InsufficientTargetDataError("only four usable target profiles")
        raise InvalidTargetDataError("target has no valid position")

    monkeypatch.setattr("meop_process.compare_batch.generate_calibration_plots", fake_generate)

    result = run_calibration_plots_batch(config=cfg, state_dir=cfg.datadir / "compare_state")

    statuses = {item.smru_name: item.status for item in result.tag_results}
    assert statuses == {
        "DEP001-BAD": "invalid_target",
        "DEP001-NONE": "no_reference",
        "DEP001-SPARSE": "insufficient_target",
        "DEP001-THIN": "insufficient_reference",
    }
    assert result.success_count == 0
    assert result.failed_count == 0
    assert result.no_reference_count == 1
    assert result.insufficient_reference_count == 1
    assert result.insufficient_target_count == 1
    assert result.invalid_target_count == 1
    assert result.tag_results[-1].written_files == (thin_output,)
    summary = result.summary_markdown.read_text(encoding="utf-8")
    assert "No reference: **1**" in summary
    assert "Insufficient reference: **1**" in summary
    assert "Insufficient target: **1**" in summary
    assert "Invalid target: **1**" in summary


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
