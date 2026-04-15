from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib.image as mpimg

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
    assert result.processed_tags == ("ct88-225-12",)
    _assert_nonempty_png(overview)
    _assert_nonempty_png(section)


def test_generate_diagnostics_ct78_writes_overview_and_section_pngs(stage_ct78_example, meop_config) -> None:
    example = stage_ct78_example()
    _stage_reference_product(meop_config, "ct78", "ct78-465-12", example["reference_lr1"])

    result = generate_diagnostics(smru_name="ct78-465-12", qf="lr1", config=meop_config)

    overview = fname_plots("ct78-465-12", deployment="ct78", qf="lr1", suffix="diags_TS_adj", config=meop_config)
    section = fname_plots("ct78-465-12", deployment="ct78", qf="lr1", suffix="transect_adj", config=meop_config)
    assert result.processed_tags == ("ct78-465-12",)
    _assert_nonempty_png(overview)
    _assert_nonempty_png(section)
