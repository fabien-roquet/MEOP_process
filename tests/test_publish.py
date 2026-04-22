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
