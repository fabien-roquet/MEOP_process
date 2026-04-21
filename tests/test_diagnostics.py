from __future__ import annotations

import shutil
from pathlib import Path
import json

import matplotlib.image as mpimg
import numpy as np
import xarray as xr

from meop_process.api import generate_diagnostics
from meop_process.catalog.filenames import fname_plots, fname_prof
from meop_process.io.netcdf import open_meop_netcdf


def _stage_reference_product(meop_config, deployment: str, smru_name: str, reference_file: Path) -> Path:
    target = fname_prof(smru_name, deployment=deployment, qf="lr1", config=meop_config)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(reference_file, target)
    return target


def _assert_nonempty_png(path: Path) -> None:
    assert path.is_file()
    assert path.stat().st_size > 20_000
    image = mpimg.imread(path)
    assert image.ndim in (2, 3)
    assert float(image.std()) > 0.001


def _relative_inventory(root: Path) -> set[str]:
    return {
        str(path.relative_to(root))
        for path in root.rglob("*")
        if path.is_file()
    }


def test_open_meop_netcdf_reads_reference_lr1(stage_ct88_example, meop_config) -> None:
    example = stage_ct88_example()
    staged = _stage_reference_product(meop_config, "ct88", "ct88-225-12", example["reference_lr1"])

    dataset = open_meop_netcdf(staged)
    try:
        assert int(dataset.sizes["N_PROF"]) == 306
        assert dataset.attrs["deployment_code"] == "ct88"
        assert dataset.attrs["smru_platform_code"] == "ct88-225-12"
    finally:
        dataset.close()


def test_generate_diagnostics_ct88_writes_overview_and_section_pngs(stage_ct88_example, meop_config) -> None:
    example = stage_ct88_example()
    _stage_reference_product(meop_config, "ct88", "ct88-225-12", example["reference_lr1"])

    result = generate_diagnostics(smru_name="ct88-225-12", qf="lr1", config=meop_config)

    overview = fname_plots("ct88-225-12", deployment="ct88", qf="lr1", suffix="diags_TS_adj", config=meop_config)
    section = fname_plots("ct88-225-12", deployment="ct88", qf="lr1", suffix="transect_adj", config=meop_config)
    deployment = meop_config.plots_by_deployment_dir / "ct88" / "ct88_lr1_deployment_overview_adj.png"
    summary = meop_config.plots_by_deployment_dir / "ct88" / "ct88_lr1_deployment_summary_adj.json"
    assert result.processed_tags == ("ct88-225-12",)
    _assert_nonempty_png(overview)
    _assert_nonempty_png(section)
    _assert_nonempty_png(deployment)
    assert summary.is_file()
    payload = json.loads(summary.read_text(encoding="utf-8"))
    assert payload["deployment"] == "ct88"
    assert payload["smru_names"] == ["ct88-225-12"]


def test_generate_diagnostics_ct78_writes_overview_and_section_pngs(stage_ct78_example, meop_config) -> None:
    example = stage_ct78_example()
    _stage_reference_product(meop_config, "ct78", "ct78-465-12", example["reference_lr1"])

    result = generate_diagnostics(smru_name="ct78-465-12", qf="lr1", config=meop_config)

    overview = fname_plots("ct78-465-12", deployment="ct78", qf="lr1", suffix="diags_TS_adj", config=meop_config)
    section = fname_plots("ct78-465-12", deployment="ct78", qf="lr1", suffix="transect_adj", config=meop_config)
    deployment = meop_config.plots_by_deployment_dir / "ct78" / "ct78_lr1_deployment_overview_adj.png"
    assert result.processed_tags == ("ct78-465-12",)
    _assert_nonempty_png(overview)
    _assert_nonempty_png(section)
    _assert_nonempty_png(deployment)


def test_generate_diagnostics_by_deployment_writes_deployment_overview(stage_ct88_example, meop_config) -> None:
    example = stage_ct88_example()
    _stage_reference_product(meop_config, "ct88", "ct88-225-12", example["reference_lr1"])

    result = generate_diagnostics(deployment="ct88", qf="lr1", config=meop_config)

    overview = fname_plots("ct88-225-12", deployment="ct88", qf="lr1", suffix="diags_TS_adj", config=meop_config)
    section = fname_plots("ct88-225-12", deployment="ct88", qf="lr1", suffix="transect_adj", config=meop_config)
    deployment = meop_config.plots_by_deployment_dir / "ct88" / "ct88_lr1_deployment_overview_adj.png"
    assert result.processed_tags == ("ct88-225-12",)
    _assert_nonempty_png(overview)
    _assert_nonempty_png(section)
    _assert_nonempty_png(deployment)


def test_generate_diagnostics_all_nan_track_still_writes_pngs(stage_ct88_example, meop_config) -> None:
    example = stage_ct88_example()
    staged = _stage_reference_product(meop_config, "ct88", "ct88-225-12", example["reference_lr1"])

    with xr.open_dataset(staged, decode_times=False) as dataset:
        updated = dataset.load()
    updated["LATITUDE"] = (updated["LATITUDE"].dims, np.full(updated["LATITUDE"].shape, np.nan, dtype=np.float64))
    updated["LONGITUDE"] = (updated["LONGITUDE"].dims, np.full(updated["LONGITUDE"].shape, np.nan, dtype=np.float64))
    updated.to_netcdf(staged)

    result = generate_diagnostics(smru_name="ct88-225-12", qf="lr1", config=meop_config)

    overview = fname_plots("ct88-225-12", deployment="ct88", qf="lr1", suffix="diags_TS_adj", config=meop_config)
    section = fname_plots("ct88-225-12", deployment="ct88", qf="lr1", suffix="transect_adj", config=meop_config)
    deployment = meop_config.plots_by_deployment_dir / "ct88" / "ct88_lr1_deployment_overview_adj.png"
    assert result.processed_tags == ("ct88-225-12",)
    _assert_nonempty_png(overview)
    _assert_nonempty_png(section)
    _assert_nonempty_png(deployment)


def test_generate_diagnostics_all_deployments_writes_global_overview(stage_ct88_example, stage_ct78_example, meop_config) -> None:
    ct88 = stage_ct88_example()
    ct78 = stage_ct78_example()
    _stage_reference_product(meop_config, "ct88", "ct88-225-12", ct88["reference_lr1"])
    _stage_reference_product(meop_config, "ct78", "ct78-465-12", ct78["reference_lr1"])

    result = generate_diagnostics(qf="lr1", config=meop_config)

    ct88_deployment = meop_config.plots_by_deployment_dir / "ct88" / "ct88_lr1_deployment_overview_adj.png"
    ct78_deployment = meop_config.plots_by_deployment_dir / "ct78" / "ct78_lr1_deployment_overview_adj.png"
    overview = meop_config.plots_overview_dir / "all_deployments_lr1_overview_adj.png"
    assert set(result.processed_tags) == {"ct88-225-12", "ct78-465-12"}
    _assert_nonempty_png(ct88_deployment)
    _assert_nonempty_png(ct78_deployment)
    _assert_nonempty_png(overview)


def test_generate_diagnostics_matches_expected_fixture_inventory(stage_ct88_example, stage_ct78_example, meop_config) -> None:
    ct88 = stage_ct88_example()
    ct78 = stage_ct78_example()
    _stage_reference_product(meop_config, "ct88", "ct88-225-12", ct88["reference_lr1"])
    _stage_reference_product(meop_config, "ct78", "ct78-465-12", ct78["reference_lr1"])

    generate_diagnostics(qf="lr1", config=meop_config)

    expected_root = Path(__file__).resolve().parent / "fixtures" / "diagnostics_expected"
    actual_root = meop_config.datadir
    actual_inventory = _relative_inventory(actual_root / "plots_by_tags")
    actual_inventory |= {f"deployments/{item}" for item in _relative_inventory(actual_root / "plots_by_deployments")}
    actual_inventory |= {f"overview/{item}" for item in _relative_inventory(actual_root / "plots_overview")}

    expected_inventory = _relative_inventory(expected_root / "plots_by_tags")
    expected_inventory |= {f"deployments/{item}" for item in _relative_inventory(expected_root / "plots_by_deployments")}
    expected_inventory |= {f"overview/{item}" for item in _relative_inventory(expected_root / "plots_overview")}

    assert actual_inventory == expected_inventory


def test_generate_diagnostics_overview_only_uses_cached_summaries(stage_ct88_example, stage_ct78_example, meop_config) -> None:
    ct88 = stage_ct88_example()
    ct78 = stage_ct78_example()
    _stage_reference_product(meop_config, "ct88", "ct88-225-12", ct88["reference_lr1"])
    _stage_reference_product(meop_config, "ct78", "ct78-465-12", ct78["reference_lr1"])

    generate_diagnostics(qf="lr1", config=meop_config)
    shutil.rmtree(meop_config.final_dataset_dir)

    result = generate_diagnostics(qf="lr1", config=meop_config, parts=("overview",))

    overview = meop_config.plots_overview_dir / "all_deployments_lr1_overview_adj.png"
    assert result.processed_tags == ()
    _assert_nonempty_png(overview)
