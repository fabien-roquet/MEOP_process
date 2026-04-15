from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

from meop_process.metadata.summaries import update_metadata_summaries


def _copy_output(reference: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(reference, destination)


def test_update_metadata_summaries_incremental(meop_config, stage_ct88_example, stage_ct78_example):
    ct88 = stage_ct88_example()
    ct78 = stage_ct78_example()

    _copy_output(
        ct88["reference_lr1"],
        meop_config.final_dataset_dir / "ct88" / "ct88-225-12_lr1_prof.nc",
    )
    first = update_metadata_summaries(meop_config, processed_deployments=["ct88"])
    assert first.written is True
    assert first.list_tags_path.exists()
    assert first.list_deployments_path.exists()

    tags = pd.read_csv(first.list_tags_path)
    deployments = pd.read_csv(first.list_deployments_path)
    assert set(tags["SMRU_PLATFORM_CODE"]) == {"ct88-225-12"}
    assert set(deployments["DEPLOYMENT_CODE"]) == {"ct88"}

    _copy_output(
        ct78["reference_hr2"],
        meop_config.final_dataset_dir / "ct78" / "ct78-465-12_hr2_prof.nc",
    )
    second = update_metadata_summaries(meop_config, processed_deployments=["ct78"])
    assert second.written is True
    tags = pd.read_csv(second.list_tags_path)
    deployments = pd.read_csv(second.list_deployments_path)
    assert set(tags["SMRU_PLATFORM_CODE"]) == {"ct88-225-12", "ct78-465-12"}
    assert set(deployments["DEPLOYMENT_CODE"]) == {"ct88", "ct78"}

    third = update_metadata_summaries(meop_config)
    assert third.written is False
