"""Tests for the meop-publish workflow."""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest
import xarray as xr

from meop_process.api import publish_data
from meop_process.catalog.filenames import fname_prof
from meop_process.models import MeopConfig
from meop_process.publishing.attributes import update_global_attributes
from meop_process.publishing.lists import build_list_profiles
from meop_process.publishing.ncfiles import create_ncfile_all
from meop_process.publishing.site import build_site, build_overview_page, build_deployment_pages, build_tag_pages
from meop_process.publish_cli import main as publish_main


FIXTURES_CT96 = Path(__file__).resolve().parent / "fixtures" / "ct96"
SMRU_NAME = "ct96-24-13"
DEPLOYMENT = "ct96"


def _plant_hr2(config: MeopConfig) -> Path:
    """Copy reference hr2 fixture into the expected final_dataset_dir location."""
    src = FIXTURES_CT96 / "ct96-24-13_hr2_prof_shortened.nc"
    dest = fname_prof(SMRU_NAME, deployment=DEPLOYMENT, qf="hr2", config=config)
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    return dest


class TestCreateNcfileAll:
    def test_creates_all_prof_nc_from_hr2(self, meop_config: MeopConfig, tmp_path: Path) -> None:
        _plant_hr2(meop_config)
        output = create_ncfile_all(meop_config, SMRU_NAME, output_dir=tmp_path)
        assert output is not None
        assert output.exists()
        assert output.name == f"{SMRU_NAME}_all_prof.nc"

    def test_returns_existing_file_without_rebuild(self, meop_config: MeopConfig, tmp_path: Path) -> None:
        _plant_hr2(meop_config)
        first = create_ncfile_all(meop_config, SMRU_NAME, output_dir=tmp_path)
        assert first is not None
        mtime = first.stat().st_mtime
        second = create_ncfile_all(meop_config, SMRU_NAME, output_dir=tmp_path)
        assert second is not None
        # file not recreated
        assert second.stat().st_mtime == mtime

    def test_rebuild_recreates_file(self, meop_config: MeopConfig, tmp_path: Path) -> None:
        _plant_hr2(meop_config)
        first = create_ncfile_all(meop_config, SMRU_NAME, output_dir=tmp_path)
        assert first is not None
        mtime = first.stat().st_mtime
        second = create_ncfile_all(meop_config, SMRU_NAME, output_dir=tmp_path, rebuild=True)
        assert second is not None
        assert second.stat().st_mtime >= mtime

    def test_returns_none_when_no_source(self, meop_config: MeopConfig, tmp_path: Path) -> None:
        result = create_ncfile_all(meop_config, "unknown-00-00", output_dir=tmp_path)
        assert result is None


class TestUpdateGlobalAttributes:
    def test_patches_data_type_and_date_update(self, meop_config: MeopConfig, tmp_path: Path) -> None:
        _plant_hr2(meop_config)
        nc_path = create_ncfile_all(meop_config, SMRU_NAME, output_dir=tmp_path)
        assert nc_path is not None

        patched = update_global_attributes(tmp_path, verbose=False)
        assert nc_path in patched

        with xr.open_dataset(nc_path) as ds:
            assert ds.attrs.get("data_type") == "Argo profile"
            assert "date_update" in ds.attrs

    def test_sets_version_when_provided(self, meop_config: MeopConfig, tmp_path: Path) -> None:
        _plant_hr2(meop_config)
        nc_path = create_ncfile_all(meop_config, SMRU_NAME, output_dir=tmp_path)
        assert nc_path is not None

        update_global_attributes(tmp_path, version="TEST-v1")
        with xr.open_dataset(nc_path) as ds:
            assert ds.attrs.get("version_database") == "TEST-v1"

    def test_empty_dir_returns_empty_list(self, tmp_path: Path) -> None:
        assert update_global_attributes(tmp_path) == []


class TestBuildListProfiles:
    def test_writes_list_profiles_csv(self, meop_config: MeopConfig, tmp_path: Path) -> None:
        _plant_hr2(meop_config)
        create_ncfile_all(meop_config, SMRU_NAME, output_dir=tmp_path)

        dest = build_list_profiles(tmp_path)
        assert dest is not None
        assert dest.exists()
        assert dest.name == "list_profiles.csv"

    def test_csv_has_expected_columns(self, meop_config: MeopConfig, tmp_path: Path) -> None:
        import pandas as pd

        _plant_hr2(meop_config)
        create_ncfile_all(meop_config, SMRU_NAME, output_dir=tmp_path)
        build_list_profiles(tmp_path)

        frame = pd.read_csv(tmp_path / "list_profiles.csv")
        for col in ("SMRU_PLATFORM_CODE", "DEPLOYMENT_CODE", "PROFILE_INDEX", "JULD", "LATITUDE", "LONGITUDE"):
            assert col in frame.columns, f"missing column {col}"

    def test_returns_none_when_no_files(self, tmp_path: Path) -> None:
        assert build_list_profiles(tmp_path) is None


class TestPublishData:
    def test_creates_all_files(self, meop_config: MeopConfig, seed_catalog: None, tmp_path: Path) -> None:
        _plant_hr2(meop_config)
        output_dir = tmp_path / "public"

        result = publish_data(
            config=meop_config,
            smru_name=SMRU_NAME,
            output_dir=output_dir,
            verbose=False,
        )

        assert Path(result["output_dir"]) == output_dir
        assert len(result["published_files"]) == 1
        assert result["list_profiles_path"] is not None
        assert Path(result["list_profiles_path"]).exists()

    def test_skip_flags_respected(self, meop_config: MeopConfig, seed_catalog: None, tmp_path: Path) -> None:
        _plant_hr2(meop_config)
        output_dir = tmp_path / "public"

        result = publish_data(
            config=meop_config,
            smru_name=SMRU_NAME,
            output_dir=output_dir,
            create_files=True,
            update_attrs=False,
            list_profiles=False,
            list_tags=False,
            verbose=False,
        )

        assert len(result["published_files"]) == 1
        assert result["list_profiles_path"] is None
        assert result["list_tags_path"] is None


class TestPublishCli:
    def test_cli_help_exits_zero(self) -> None:
        with pytest.raises(SystemExit) as exc_info:
            publish_main(["--help"])
        assert exc_info.value.code == 0

    def test_cli_runs_with_output_dir(self, meop_config: MeopConfig, seed_catalog: None, tmp_path: Path) -> None:
        _plant_hr2(meop_config)
        output_dir = tmp_path / "pub"

        rc = publish_main(
            [
                "--smru_name", SMRU_NAME,
                "--output-dir", str(output_dir),
                "--no-list-tags",
            ],
            config=meop_config,
        )
        assert rc == 0
        assert (output_dir / f"{SMRU_NAME}_all_prof.nc").exists()
        assert (output_dir / "list_profiles.csv").exists()


# ---------------------------------------------------------------------------
# Site generation tests
# ---------------------------------------------------------------------------

def _plant_fake_plots(base: Path) -> dict[str, Path]:
    """Create minimal fake PNG files mimicking the diagnostics layout.

    Returns a dict with keys 'tags_dir', 'dep_dir', 'overview_dir'.
    """
    dep = "DEP01"
    smru = "DEP01-001-24"

    tags_dir = base / "plots_by_tags" / dep
    tags_dir.mkdir(parents=True)
    dep_dir = base / "plots_by_deployments" / dep
    dep_dir.mkdir(parents=True)
    overview_dir = base / "plots_overview"
    overview_dir.mkdir(parents=True)

    # Per-tag plots (1x1 pixel PNG is fine for existence tests)
    _tiny_png = b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01\x08\x02\x00\x00\x00\x90wS\xde\x00\x00\x00\x0cIDATx\x9cc\xf8\x0f\x00\x00\x01\x01\x00\x05\x18\xd8N\x00\x00\x00\x00IEND\xaeB`\x82"
    for plot_type in ("section", "TS", "map", "profiles", "flags"):
        (tags_dir / f"{smru}_lr1_{plot_type}_adj.png").write_bytes(_tiny_png)
    (tags_dir / f"{smru}_lr1_info_adj.txt").write_text(
        f"{smru}\ndeployment: {dep}\nperiod: 2024-01-01 to 2024-02-01\n", encoding="utf-8"
    )

    # Per-deployment plots
    for plot_type in ("map", "TS", "histograms", "timing"):
        (dep_dir / f"{dep}_lr1_{plot_type}_adj.png").write_bytes(_tiny_png)

    # Overview plots
    for plot_type in ("map", "histograms", "timing"):
        (overview_dir / f"all_deployments_lr1_{plot_type}_adj.png").write_bytes(_tiny_png)

    return {
        "tags_dir": base / "plots_by_tags",
        "dep_dir": base / "plots_by_deployments",
        "overview_dir": overview_dir,
        "dep": dep,
        "smru": smru,
    }


class TestBuildSite:
    def test_build_site_returns_html_paths(self, tmp_path: Path) -> None:
        dirs = _plant_fake_plots(tmp_path)
        paths = build_site(dirs["tags_dir"], dirs["dep_dir"], dirs["overview_dir"])
        assert len(paths) > 0
        assert all(p.suffix == ".html" for p in paths)

    def test_overview_index_created(self, tmp_path: Path) -> None:
        dirs = _plant_fake_plots(tmp_path)
        build_site(dirs["tags_dir"], dirs["dep_dir"], dirs["overview_dir"])
        overview_html = dirs["overview_dir"] / "index.html"
        assert overview_html.exists()

    def test_overview_lists_deployment_link(self, tmp_path: Path) -> None:
        dirs = _plant_fake_plots(tmp_path)
        build_site(dirs["tags_dir"], dirs["dep_dir"], dirs["overview_dir"])
        content = (dirs["overview_dir"] / "index.html").read_text(encoding="utf-8")
        assert dirs["dep"] in content

    def test_deployment_index_created(self, tmp_path: Path) -> None:
        dirs = _plant_fake_plots(tmp_path)
        build_site(dirs["tags_dir"], dirs["dep_dir"], dirs["overview_dir"])
        dep_html = dirs["dep_dir"] / dirs["dep"] / "index.html"
        assert dep_html.exists()

    def test_deployment_index_links_to_tag(self, tmp_path: Path) -> None:
        dirs = _plant_fake_plots(tmp_path)
        build_site(dirs["tags_dir"], dirs["dep_dir"], dirs["overview_dir"])
        content = (dirs["dep_dir"] / dirs["dep"] / "index.html").read_text(encoding="utf-8")
        assert dirs["smru"] in content

    def test_tag_page_created(self, tmp_path: Path) -> None:
        dirs = _plant_fake_plots(tmp_path)
        build_site(dirs["tags_dir"], dirs["dep_dir"], dirs["overview_dir"])
        tag_html = dirs["tags_dir"] / dirs["dep"] / f"{dirs['smru']}.html"
        assert tag_html.exists()

    def test_tag_page_contains_info_text(self, tmp_path: Path) -> None:
        dirs = _plant_fake_plots(tmp_path)
        build_site(dirs["tags_dir"], dirs["dep_dir"], dirs["overview_dir"])
        tag_html = dirs["tags_dir"] / dirs["dep"] / f"{dirs['smru']}.html"
        content = tag_html.read_text(encoding="utf-8")
        assert "2024-01-01" in content  # from info.txt

    def test_rebuild_false_skips_existing(self, tmp_path: Path) -> None:
        dirs = _plant_fake_plots(tmp_path)
        build_site(dirs["tags_dir"], dirs["dep_dir"], dirs["overview_dir"])
        overview_html = dirs["overview_dir"] / "index.html"
        mtime = overview_html.stat().st_mtime
        build_site(dirs["tags_dir"], dirs["dep_dir"], dirs["overview_dir"], rebuild=False)
        assert overview_html.stat().st_mtime == mtime

    def test_rebuild_true_overwrites_existing(self, tmp_path: Path) -> None:
        dirs = _plant_fake_plots(tmp_path)
        build_site(dirs["tags_dir"], dirs["dep_dir"], dirs["overview_dir"])
        overview_html = dirs["overview_dir"] / "index.html"
        mtime = overview_html.stat().st_mtime
        build_site(dirs["tags_dir"], dirs["dep_dir"], dirs["overview_dir"], rebuild=True)
        assert overview_html.stat().st_mtime >= mtime

    def test_empty_dirs_returns_empty_list(self, tmp_path: Path) -> None:
        result = build_site(
            tmp_path / "plots_by_tags",
            tmp_path / "plots_by_deployments",
            tmp_path / "plots_overview",
        )
        assert result == []

    def test_build_site_cli_flag(self, meop_config: MeopConfig, tmp_path: Path) -> None:
        """--build-site flag should not raise even when plots dirs are empty."""
        dirs = _plant_fake_plots(tmp_path)
        # Point config to our fake plots dirs via patching
        import meop_process.models as _models
        original_plotdir = _models.MeopConfig.plotdir.fget  # type: ignore[attr-defined]
        original_dep_dir = _models.MeopConfig.plots_by_deployment_dir.fget  # type: ignore[attr-defined]
        original_ov_dir = _models.MeopConfig.plots_overview_dir.fget  # type: ignore[attr-defined]
        try:
            _models.MeopConfig.plotdir = property(lambda self: dirs["tags_dir"])  # type: ignore[attr-defined]
            _models.MeopConfig.plots_by_deployment_dir = property(lambda self: dirs["dep_dir"])  # type: ignore[attr-defined]
            _models.MeopConfig.plots_overview_dir = property(lambda self: dirs["overview_dir"])  # type: ignore[attr-defined]
            from meop_process.workflows.publish import publish
            result = publish(meop_config, build_site=True, create_files=False, update_attrs=False,
                             list_profiles=False, list_tags=False, output_dir=tmp_path / "pub")
            assert len(result.site_paths) > 0
        finally:
            _models.MeopConfig.plotdir = property(original_plotdir)  # type: ignore[attr-defined]
            _models.MeopConfig.plots_by_deployment_dir = property(original_dep_dir)  # type: ignore[attr-defined]
            _models.MeopConfig.plots_overview_dir = property(original_ov_dir)  # type: ignore[attr-defined]
