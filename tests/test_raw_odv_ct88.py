from __future__ import annotations

from meop_process.catalog.deployments import load_info_deployment
from meop_process.io.raw_odv import build_odv_profile_index, discover_raw_odv_files, import_raw_data_zip, load_raw_odv_profiles


def test_ct88_raw_index_is_available_from_runtime_store(meop_config, stage_ct88_example) -> None:
    stage_ct88_example()

    assert import_raw_data_zip(meop_config, "ct88") is True

    files = discover_raw_odv_files(meop_config, "ct88")
    index = build_odv_profile_index(files)

    assert files.preferred_ctd_text is not None
    assert index.total_profiles == 6497
    assert len(index.smru_names) == 19
    assert index.profile_count_by_smru["ct88-225-12"] == 306
    assert index.max_levels == 16


def test_load_info_deployment_uses_raw_odv_index_when_catalog_lacks_tag_rows(meop_config, stage_ct88_example) -> None:
    stage_ct88_example()
    assert import_raw_data_zip(meop_config, "ct88") is True

    info = load_info_deployment(meop_config, deployment="ct88")

    assert info.invalid_code is False
    assert len(info.raw_smru_names) == 19
    assert info.raw_profile_count_by_smru["ct88-225-12"] == 306
    assert "ct88-225-12" in info.list_smru_name




def test_load_info_deployment_accepts_deployment2_json_when_catalog_csv_is_missing(meop_config, stage_ct88_example) -> None:
    stage_ct88_example()
    catalog_path = meop_config.catalogdir / "list_deployment.csv"
    if catalog_path.exists():
        catalog_path.unlink()

    assert import_raw_data_zip(meop_config, "ct88") is True
    info = load_info_deployment(meop_config, deployment="ct88")

    assert info.invalid_code is False
    assert info.record is not None
    assert info.record.description == "Weddell Seal CTD Jan 2012"
    assert info.record.pi_code.strip() == "COSTA"

def test_load_raw_odv_profiles_keeps_ct88_profile_structure(meop_config, stage_ct88_example) -> None:
    stage_ct88_example()
    assert import_raw_data_zip(meop_config, "ct88") is True

    files = discover_raw_odv_files(meop_config, "ct88")
    profiles = load_raw_odv_profiles(files)
    first = next(profile for profile in profiles if profile.smru_name == "ct88-225-12")

    assert len(profiles) == 6497
    assert first.station == 1
    assert first.n_levels == 16
    assert first.timestamp == "2012-02-06 06:00"
    assert first.pressure[:4] == (4.0, 10.0, 20.0, 26.0)
    assert first.temperature[:4] == (-1.5893, -1.5733, -1.5532, -1.5322)
    assert first.salinity[:4] == (33.302, 33.3199, 33.3491, 33.3662)
