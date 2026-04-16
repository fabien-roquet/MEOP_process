from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from meop_process.catalog.filenames import fname_prof, fname_traj
from meop_process.io.raw_odv import import_raw_data_zip
from meop_process.metadata.patch import update_metadata_from_table
from meop_process.models import Selection
from meop_process.processing.adjustments import apply_adjustments
from meop_process.processing.fr0 import create_fr0_python
from meop_process.processing.hr import create_fr1_python, create_hr0_python, create_hr1_python
from meop_process.processing.hr2 import create_hr2_python
from meop_process.processing.ncargo import create_ncargo_python


TIMESTAMP = datetime(2024, 3, 6, 2, 15, 16, tzinfo=timezone.utc)
TIMESTAMP_CT78 = datetime(2024, 3, 5, 18, 6, 7, tzinfo=timezone.utc)
CORE_DIMS = {"N_PROF", "N_LEVELS", "N_PARAM", "N_CALIB"}


def _open_any(path: Path):
    try:
        return xr.open_dataset(path, decode_times=False)
    except Exception:
        return xr.open_dataset(path, engine="scipy", decode_times=False)


def _seed_table_param_min_profiles_one(meop_config) -> None:
    table_param = meop_config.tablesdir / "table_param.csv"
    table_param.write_text(
        "row_name,temp_error,psal_error,minT,maxT,minS,maxS,min_Nprof,pmax,pmax_fluo,is_lon_centre_180\n"
        "ct96,0.1,0.2,-3,32,4,40,1,1000,200,0\n",
        encoding="utf-8",
    )


def _stage_ct96_with_fr0(meop_config, stage_ct96_example):
    staged = stage_ct96_example()
    assert import_raw_data_zip(meop_config, "ct96") is True
    create_ncargo_python(
        meop_config,
        Selection(deployment="ct96", smru_name="ct96-24-13"),
        now=TIMESTAMP,
    )
    _seed_table_param_min_profiles_one(meop_config)
    create_fr0_python(
        meop_config,
        Selection(deployment="ct96", smru_name="ct96-24-13"),
        now=TIMESTAMP,
    )
    apply_adjustments(
        meop_config,
        Selection(deployment="ct96", smru_name="ct96-24-13"),
    )
    return staged


def _stage_ct78_hr1(meop_config, stage_ct78_example):
    staged = stage_ct78_example()
    assert import_raw_data_zip(meop_config, "ct78") is True
    create_ncargo_python(
        meop_config,
        Selection(deployment="ct78", smru_name="ct78-465-12"),
        now=TIMESTAMP_CT78,
    )
    update_metadata_from_table(meop_config, smru_name="ct78-465-12", modes=("lr0",))
    create_hr0_python(
        meop_config,
        Selection(deployment="ct78", smru_name="ct78-465-12"),
        now=TIMESTAMP_CT78,
    )
    apply_adjustments(
        meop_config,
        Selection(deployment="ct78", smru_name="ct78-465-12"),
    )
    create_hr1_python(
        meop_config,
        Selection(deployment="ct78", smru_name="ct78-465-12"),
        thermal_lag=True,
        now=TIMESTAMP_CT78,
    )
    return staged


def test_create_fr0_python_matches_ct96_traj_reference_core_fields(meop_config, stage_ct96_example) -> None:
    staged = _stage_ct96_with_fr0(meop_config, stage_ct96_example)

    candidate_path = fname_traj("ct96-24-13", config=meop_config)
    assert candidate_path.exists()

    with _open_any(candidate_path) as candidate, _open_any(staged["reference_traj"]) as reference:
        assert 0 <= int(reference.sizes["TIME"]) - int(candidate.sizes["TIME"]) <= 1
        shared = min(int(candidate.sizes["TIME"]), int(reference.sizes["TIME"]))
        np.testing.assert_allclose(candidate["TIME"].values[:shared], reference["TIME"].values[:shared], atol=1e-6)
        np.testing.assert_allclose(candidate["LATITUDE"].values[:shared], reference["LATITUDE"].values[:shared], atol=1e-3)
        np.testing.assert_allclose(candidate["LONGITUDE"].values[:shared], reference["LONGITUDE"].values[:shared], atol=1e-3)
        np.testing.assert_allclose(candidate["PRES"].values[:shared], reference["PRES"].values[:shared], atol=1e-6)
        np.testing.assert_allclose(candidate["TEMP"].values[:shared], reference["TEMP"].values[:shared], atol=1e-4)
        np.testing.assert_allclose(candidate["PSAL"].values[:shared], reference["PSAL"].values[:shared], atol=1e-4)


def test_create_fr0_python_writes_classic_compressed_outputs(meop_config, stage_ct96_example) -> None:
    _stage_ct96_with_fr0(meop_config, stage_ct96_example)

    netcdf4 = pytest.importorskip("netCDF4")
    prof_path = fname_prof("ct96-24-13", qf="fr0", config=meop_config)
    traj_path = fname_traj("ct96-24-13", config=meop_config)

    with netcdf4.Dataset(prof_path, mode="r") as ds:
        assert ds.file_format == "NETCDF4_CLASSIC"
        assert ds.variables["TEMP"].filters().get("zlib") is True

    with netcdf4.Dataset(traj_path, mode="r") as ds:
        assert ds.file_format == "NETCDF4_CLASSIC"
        assert ds.variables["TEMP"].filters().get("zlib") is True

    with _open_any(prof_path) as candidate:
        assert CORE_DIMS.issubset(set(candidate.sizes))


def test_create_fr1_python_matches_ct96_shortened_reference_core_fields(meop_config, stage_ct96_example) -> None:
    staged = _stage_ct96_with_fr0(meop_config, stage_ct96_example)

    result = create_fr1_python(
        meop_config,
        Selection(deployment="ct96", smru_name="ct96-24-13"),
        thermal_lag=True,
        now=TIMESTAMP,
    )

    assert result.processed_tags == ("ct96-24-13",)
    candidate_path = fname_prof("ct96-24-13", qf="fr1", config=meop_config)
    assert candidate_path.exists()

    with _open_any(candidate_path) as candidate, _open_any(staged["reference_fr1"]) as reference:
        assert candidate.attrs["thermal_lag_adjustment"] == "yes"
        for dim in CORE_DIMS:
            assert int(candidate.sizes[dim]) == int(reference.sizes[dim])
        np.testing.assert_allclose(candidate["JULD"].values, reference["JULD"].values)
        np.testing.assert_allclose(candidate["LATITUDE"].values, reference["LATITUDE"].values, atol=2e-5)
        np.testing.assert_allclose(candidate["LONGITUDE"].values, reference["LONGITUDE"].values, atol=2e-5)
        np.testing.assert_allclose(candidate["PRES"].values, reference["PRES"].values, atol=1e-6, equal_nan=True)

        temp_diff = np.abs(candidate["TEMP_ADJUSTED"].values - reference["TEMP_ADJUSTED"].values)
        psal_diff = np.abs(candidate["PSAL_ADJUSTED"].values - reference["PSAL_ADJUSTED"].values)
        finite_temp = np.isfinite(temp_diff)
        finite_psal = np.isfinite(psal_diff)
        assert finite_temp.any()
        assert finite_psal.any()
        assert float(np.nanmean(temp_diff[finite_temp])) < 0.02
        assert float(np.nanmean(psal_diff[finite_psal])) < 0.01
        assert float(np.nanpercentile(temp_diff[finite_temp], 99)) < 0.10
        assert float(np.nanpercentile(psal_diff[finite_psal], 99)) < 0.05
        assert float(np.nanmax(temp_diff[finite_temp])) < 0.20
        assert float(np.nanmax(psal_diff[finite_psal])) < 0.10


def test_create_hr2_python_uses_fr1_when_available_matches_ct96_reference(meop_config, stage_ct96_example) -> None:
    staged = _stage_ct96_with_fr0(meop_config, stage_ct96_example)
    create_fr1_python(
        meop_config,
        Selection(deployment="ct96", smru_name="ct96-24-13"),
        thermal_lag=True,
        now=TIMESTAMP,
    )

    result = create_hr2_python(
        meop_config,
        Selection(deployment="ct96", smru_name="ct96-24-13"),
        now=TIMESTAMP,
    )

    assert result.processed_tags == ("ct96-24-13",)
    candidate_path = fname_prof("ct96-24-13", qf="hr2", config=meop_config)
    assert candidate_path.exists()

    with _open_any(candidate_path) as candidate, _open_any(staged["reference_hr2"]) as reference:
        for dim in CORE_DIMS:
            assert int(candidate.sizes[dim]) == int(reference.sizes[dim])
        np.testing.assert_allclose(candidate["JULD"].values, reference["JULD"].values)
        np.testing.assert_allclose(candidate["LATITUDE"].values, reference["LATITUDE"].values, atol=2e-5)
        np.testing.assert_allclose(candidate["LONGITUDE"].values, reference["LONGITUDE"].values, atol=2e-5)
        np.testing.assert_allclose(candidate["PRES"].values, reference["PRES"].values, atol=1e-6, equal_nan=True)

        temp_diff = np.abs(candidate["TEMP_ADJUSTED"].values - reference["TEMP_ADJUSTED"].values)
        psal_diff = np.abs(candidate["PSAL_ADJUSTED"].values - reference["PSAL_ADJUSTED"].values)
        finite_temp = np.isfinite(temp_diff)
        finite_psal = np.isfinite(psal_diff)
        assert finite_temp.any()
        assert finite_psal.any()
        assert float(np.nanmean(temp_diff[finite_temp])) < 0.02
        assert float(np.nanmean(psal_diff[finite_psal])) < 0.01
        assert float(np.nanpercentile(temp_diff[finite_temp], 99)) < 0.10
        assert float(np.nanpercentile(psal_diff[finite_psal], 99)) < 0.05


def test_create_hr2_python_uses_hr1_when_fr_data_is_absent(meop_config, stage_ct78_example) -> None:
    staged = _stage_ct78_hr1(meop_config, stage_ct78_example)

    result = create_hr2_python(
        meop_config,
        Selection(deployment="ct78", smru_name="ct78-465-12"),
        now=TIMESTAMP_CT78,
    )

    assert result.processed_tags == ("ct78-465-12",)
    candidate_path = fname_prof("ct78-465-12", qf="hr2", config=meop_config)
    assert candidate_path.exists()

    with _open_any(candidate_path) as candidate, _open_any(staged["reference_hr2"]) as reference:
        subset = candidate.isel(N_PROF=slice(0, int(reference.sizes["N_PROF"])))
        assert int(subset.sizes["N_LEVELS"]) == int(reference.sizes["N_LEVELS"])
        np.testing.assert_allclose(subset["JULD"].values, reference["JULD"].values)
        np.testing.assert_allclose(subset["LATITUDE"].values, reference["LATITUDE"].values)
        np.testing.assert_allclose(subset["LONGITUDE"].values, reference["LONGITUDE"].values)
        np.testing.assert_allclose(subset["PRES"].values, reference["PRES"].values, atol=1e-6, equal_nan=True)

        temp_diff = np.abs(subset["TEMP_ADJUSTED"].values - reference["TEMP_ADJUSTED"].values)
        psal_diff = np.abs(subset["PSAL_ADJUSTED"].values - reference["PSAL_ADJUSTED"].values)
        finite_temp = np.isfinite(temp_diff)
        finite_psal = np.isfinite(psal_diff)
        assert finite_temp.any()
        assert finite_psal.any()
        assert float(np.nanmean(temp_diff[finite_temp])) < 0.03
        assert float(np.nanmean(psal_diff[finite_psal])) < 0.03
        assert float(np.nanpercentile(temp_diff[finite_temp], 99)) < 0.10
        assert float(np.nanpercentile(psal_diff[finite_psal], 99)) < 0.10
