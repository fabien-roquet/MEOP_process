from __future__ import annotations

from pathlib import Path
import zipfile

from meop_process.processing.adjustments import AdjustmentResult
from meop_process.processing.hr import HrResult
from meop_process.workflows import process as process_module



def test_process_tags_runs_pure_python_tlc_sequence(meop_config, monkeypatch, seed_catalog) -> None:
    seed_catalog(deployment="DEP001", smru_name="DEP001-AAA")
    metadata_calls = []
    adjustment_calls = []
    tlc_calls = []
    fr0_calls = []
    fr1_calls = []
    hr2_calls = []
    monkeypatch.setattr(process_module, "sync_external_config_files", lambda config: [config.config_files_dir])
    monkeypatch.setattr(process_module, "import_raw_data_zip", lambda config, deployment: True)
    monkeypatch.setattr(process_module, "remove_deployment_outputs", lambda config, info: [])
    monkeypatch.setattr(
        process_module,
        "create_ncargo_python",
        lambda config, selection: HrResult(written_files=(Path("lr0.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(
        process_module,
        "create_hr0_python",
        lambda config, selection: HrResult(written_files=(Path("hr0.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(
        process_module,
        "create_fr0_python",
        lambda config, selection: fr0_calls.append((selection.deployment, selection.smru_name)) or HrResult(written_files=(Path("fr0.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(process_module, "apply_location_adjustment_placeholder", lambda config, selection: None)
    monkeypatch.setattr(
        process_module,
        "update_metadata_from_table",
        lambda config, deployment, smru_name: metadata_calls.append((deployment, smru_name)),
    )
    monkeypatch.setattr(
        process_module,
        "apply_adjustments",
        lambda config, selection: adjustment_calls.append((selection.deployment, selection.smru_name)) or AdjustmentResult(written_files=(Path("lr0.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(
        process_module,
        "apply_tlc",
        lambda config, selection: tlc_calls.append((selection.deployment, selection.smru_name)) or HrResult(written_files=(Path("hr1.nc"), Path("lr1.nc")), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(
        process_module,
        "apply_tlc_fr",
        lambda config, selection: fr1_calls.append((selection.deployment, selection.smru_name)) or HrResult(written_files=(Path("fr1.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(
        process_module,
        "create_hr2_python",
        lambda config, selection: hr2_calls.append((selection.deployment, selection.smru_name)) or HrResult(written_files=(Path("hr2.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )

    result = process_module.process_tags(meop_config, deployment="DEP001")

    assert result.success is True
    assert bool(result) is True
    assert fr0_calls == [("DEP001", "")]
    assert metadata_calls == [("DEP001", "")]
    assert adjustment_calls == [("DEP001", "")]
    assert tlc_calls == [("DEP001", "")]
    assert fr1_calls == [("DEP001", "")]
    assert hr2_calls == [("DEP001", "")]



def test_process_tags_uses_python_notlc_branch(meop_config, monkeypatch, seed_catalog) -> None:
    seed_catalog(deployment="DEP001", smru_name="DEP001-AAA")
    notlc_calls = []
    adjustment_calls = []
    fr1_calls = []
    hr2_calls = []
    monkeypatch.setattr(process_module, "sync_external_config_files", lambda config: [])
    monkeypatch.setattr(process_module, "import_raw_data_zip", lambda config, deployment: True)
    monkeypatch.setattr(process_module, "remove_deployment_outputs", lambda config, info: [])
    monkeypatch.setattr(
        process_module,
        "create_ncargo_python",
        lambda config, selection: HrResult(written_files=(Path("lr0.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(
        process_module,
        "create_hr0_python",
        lambda config, selection: HrResult(written_files=(Path("hr0.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(
        process_module,
        "create_fr0_python",
        lambda config, selection: HrResult(written_files=(Path("fr0.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(process_module, "apply_location_adjustment_placeholder", lambda config, selection: None)
    monkeypatch.setattr(process_module, "update_metadata_from_table", lambda config, deployment, smru_name: [])
    monkeypatch.setattr(
        process_module,
        "apply_adjustments",
        lambda config, selection: adjustment_calls.append((selection.deployment, selection.smru_name)) or AdjustmentResult(written_files=(Path("lr0.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(
        process_module,
        "apply_notlc",
        lambda config, selection: notlc_calls.append((selection.deployment, selection.smru_name)) or HrResult(written_files=(Path("hr1.nc"), Path("lr1.nc")), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(
        process_module,
        "apply_notlc_fr",
        lambda config, selection: fr1_calls.append((selection.deployment, selection.smru_name)) or HrResult(written_files=(Path("fr1.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(
        process_module,
        "create_hr2_python",
        lambda config, selection: hr2_calls.append((selection.deployment, selection.smru_name)) or HrResult(written_files=(Path("hr2.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )

    result = process_module.process_tags(meop_config, smru_name="DEP001-AAA", notlc=True)

    assert result.success is True
    assert bool(result) is True
    assert notlc_calls == [("DEP001", "DEP001-AAA")]
    assert adjustment_calls == [("DEP001", "DEP001-AAA")]
    assert fr1_calls == [("DEP001", "DEP001-AAA")]
    assert hr2_calls == [("DEP001", "DEP001-AAA")]



def test_process_tags_stops_early_for_invalid_selection(meop_config, monkeypatch) -> None:
    monkeypatch.setattr(process_module, "sync_external_config_files", lambda config: [])
    monkeypatch.setattr(process_module, "import_raw_data_zip", lambda config, deployment: True)
    result = process_module.process_tags(meop_config, deployment="DEP001")
    assert result.success is False
    assert result.reason == "invalid deployment code"
    assert bool(result) is False



def test_process_tags_notlc_prunes_intermediate_products_by_default(meop_config, stage_ct88_example) -> None:
    stage_ct88_example()

    result = process_module.process_tags(
        meop_config,
        deployment="ct88",
        smru_name="ct88-225-12",
        notlc=True,
    )

    assert result.success is True
    assert bool(result) is True
    assert not (meop_config.final_dataset_dir / "ct88" / "ct88-225-12_lr0_prof.nc").exists()
    assert not (meop_config.final_dataset_dir / "ct88" / "ct88-225-12_hr0_prof.nc").exists()
    assert (meop_config.final_dataset_dir / "ct88" / "ct88-225-12_hr1_prof.nc").exists()
    assert (meop_config.final_dataset_dir / "ct88" / "ct88-225-12_lr1_prof.nc").exists()
    assert (meop_config.final_dataset_dir / "ct88" / "ct88-225-12_hr2_prof.nc").exists()
    assert any(path.endswith("_lr0_prof.nc") for path in result.pruned_files)
    assert any(path.endswith("_hr0_prof.nc") for path in result.pruned_files)


def test_process_tags_notlc_can_keep_intermediate_products(meop_config, stage_ct88_example) -> None:
    stage_ct88_example()

    result = process_module.process_tags(
        meop_config,
        deployment="ct88",
        smru_name="ct88-225-12",
        notlc=True,
        keep_intermediate_products=True,
    )

    assert result.success is True
    assert bool(result) is True
    assert (meop_config.final_dataset_dir / "ct88" / "ct88-225-12_lr0_prof.nc").exists()
    assert (meop_config.final_dataset_dir / "ct88" / "ct88-225-12_hr0_prof.nc").exists()
    assert (meop_config.final_dataset_dir / "ct88" / "ct88-225-12_hr1_prof.nc").exists()
    assert (meop_config.final_dataset_dir / "ct88" / "ct88-225-12_lr1_prof.nc").exists()
    assert (meop_config.final_dataset_dir / "ct88" / "ct88-225-12_hr2_prof.nc").exists()
    assert result.pruned_files == ()


def test_process_tags_reports_missing_hr2_files(meop_config, monkeypatch, seed_catalog) -> None:
    seed_catalog(deployment="DEP001", smru_name="DEP001-AAA")
    monkeypatch.setattr(process_module, "sync_external_config_files", lambda config: [])
    monkeypatch.setattr(process_module, "import_raw_data_zip", lambda config, deployment: True)
    monkeypatch.setattr(process_module, "remove_deployment_outputs", lambda config, info: [])
    monkeypatch.setattr(
        process_module,
        "create_ncargo_python",
        lambda config, selection: HrResult(written_files=(Path("lr0.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(process_module, "apply_location_adjustment_placeholder", lambda config, selection: None)
    monkeypatch.setattr(
        process_module,
        "create_hr0_python",
        lambda config, selection: HrResult(written_files=(Path("hr0.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(
        process_module,
        "create_fr0_python",
        lambda config, selection: HrResult(written_files=(Path("fr0.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(process_module, "update_metadata_from_table", lambda config, deployment, smru_name: [])
    monkeypatch.setattr(
        process_module,
        "apply_adjustments",
        lambda config, selection: AdjustmentResult(written_files=(Path("lr0.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(
        process_module,
        "apply_tlc",
        lambda config, selection: HrResult(written_files=(Path("hr1.nc"), Path("lr1.nc")), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(
        process_module,
        "apply_tlc_fr",
        lambda config, selection: HrResult(written_files=(Path("fr1.nc"),), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(
        process_module,
        "create_hr2_python",
        lambda config, selection: HrResult(written_files=(), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )

    result = process_module.process_tags(meop_config, deployment="DEP001")

    assert result.success is False
    assert result.reason == "no HR2 files written"


def test_process_tags_reports_empty_raw_odv_archive(meop_config, seed_catalog) -> None:
    seed_catalog(deployment="DEP001", smru_name="DEP001-AAA")
    archive = meop_config.raw_odv_dir / "DEP001_ODV.zip"
    archive.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive, "w"):
        pass

    result = process_module.process_tags(meop_config, deployment="DEP001")

    assert result.success is False
    assert result.reason == "raw ODV archive is empty"


def test_process_tags_reports_hr0_failure_with_lr0_input_count(meop_config, monkeypatch, seed_catalog) -> None:
    seed_catalog(deployment="DEP001", smru_name="DEP001-AAA")
    monkeypatch.setattr(process_module, "sync_external_config_files", lambda config: [])
    monkeypatch.setattr(process_module, "import_raw_data_zip", lambda config, deployment: True)
    monkeypatch.setattr(process_module, "remove_deployment_outputs", lambda config, info: [])
    monkeypatch.setattr(
        process_module,
        "create_ncargo_python",
        lambda config, selection: HrResult(written_files=(meop_config.final_dataset_dir / "DEP001" / "DEP001-AAA_lr0_prof.nc",), processed_tags=(selection.smru_name or "DEP001-AAA",)),
    )
    monkeypatch.setattr(process_module, "apply_location_adjustment_placeholder", lambda config, selection: None)
    monkeypatch.setattr(
        process_module,
        "create_hr0_python",
        lambda config, selection: HrResult(written_files=(), processed_tags=()),
    )

    lr0_path = meop_config.final_dataset_dir / "DEP001" / "DEP001-AAA_lr0_prof.nc"
    lr0_path.parent.mkdir(parents=True, exist_ok=True)
    lr0_path.write_text("placeholder", encoding="utf-8")

    result = process_module.process_tags(meop_config, deployment="DEP001")

    assert result.success is False
    assert result.reason == "no HR0 files written from 1 LR0 inputs"
