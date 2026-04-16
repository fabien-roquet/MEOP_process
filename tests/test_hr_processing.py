from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pytest
import xarray as xr

from meop_process.catalog.filenames import fname_prof
from meop_process.io.raw_odv import import_raw_data_zip
from meop_process.metadata.patch import update_metadata_from_table
from meop_process.models import Selection
from meop_process.processing.adjustments import apply_adjustments
from meop_process.processing.hr import create_hr0_python, create_hr1_python
from meop_process.processing.stabilise_sa_const_ct import SolverInfo
from meop_process.processing.ncargo import create_ncargo_python


TIMESTAMP = datetime(2024, 3, 6, 2, 15, 16, tzinfo=timezone.utc)
CORE_DIMS = {"N_PROF", "N_LEVELS", "N_PARAM", "N_CALIB"}



def _open_any(path):
    try:
        return xr.open_dataset(path, decode_times=False)
    except Exception:
        return xr.open_dataset(path, engine="scipy", decode_times=False)



def _stage_lr0(meop_config, stage_ct88_example) -> dict[str, str]:
    staged = stage_ct88_example()
    assert import_raw_data_zip(meop_config, "ct88") is True
    create_ncargo_python(
        meop_config,
        Selection(deployment="ct88", smru_name="ct88-225-12"),
        now=TIMESTAMP,
    )
    update_metadata_from_table(meop_config, smru_name="ct88-225-12", modes=("lr0",))
    return staged



def _stage_hr0_adjusted(meop_config, stage_ct88_example):
    staged = _stage_lr0(meop_config, stage_ct88_example)
    create_hr0_python(
        meop_config,
        Selection(deployment="ct88", smru_name="ct88-225-12"),
        now=TIMESTAMP,
    )
    apply_adjustments(
        meop_config,
        Selection(deployment="ct88", smru_name="ct88-225-12"),
    )
    return staged



def test_create_hr0_python_builds_1000_level_interpolated_grid(meop_config, stage_ct88_example) -> None:
    staged = _stage_lr0(meop_config, stage_ct88_example)
    result = create_hr0_python(
        meop_config,
        Selection(deployment="ct88", smru_name="ct88-225-12"),
        now=TIMESTAMP,
    )

    assert result.processed_tags == ("ct88-225-12",)
    candidate = fname_prof("ct88-225-12", qf="hr0", config=meop_config)
    assert candidate.exists()

    with _open_any(staged["reference_lr0"]) as lr0, _open_any(candidate) as hr0:
        assert int(hr0.sizes["N_PROF"]) == int(lr0.sizes["N_PROF"])
        assert int(hr0.sizes["N_LEVELS"]) == 1000
        np.testing.assert_allclose(hr0["PRES"].values[0, :5], np.array([1, 2, 3, 4, 5], dtype=np.float32))
        np.testing.assert_allclose(hr0["TEMP"].values[0, 3], lr0["TEMP"].values[0, 0], atol=1e-4)
        np.testing.assert_allclose(hr0["TEMP"].values[0, 9], lr0["TEMP"].values[0, 1], atol=1e-4)
        np.testing.assert_allclose(hr0["PSAL"].values[0, 3], lr0["PSAL"].values[0, 0], atol=1e-4)
        np.testing.assert_allclose(hr0["PSAL"].values[0, 9], lr0["PSAL"].values[0, 1], atol=1e-4)
        assert np.isnan(hr0["TEMP"].values[0, 156])
        assert np.isnan(hr0["PSAL"].values[0, 156])
        assert hr0.attrs["thermal_lag_adjustment"] == "no"


def test_create_hr0_python_writes_classic_compressed_file(meop_config, stage_ct88_example) -> None:
    _stage_lr0(meop_config, stage_ct88_example)
    create_hr0_python(
        meop_config,
        Selection(deployment="ct88", smru_name="ct88-225-12"),
        now=TIMESTAMP,
    )

    candidate = fname_prof("ct88-225-12", qf="hr0", config=meop_config)
    netcdf4 = pytest.importorskip("netCDF4")
    with netcdf4.Dataset(candidate, mode="r") as ds:
        assert ds.file_format == "NETCDF4_CLASSIC"
        filters = ds.variables["TEMP"].filters()
        assert filters.get("zlib") is True
        assert ds.variables["TEMP"].dtype == "float32"

    with _open_any(candidate) as hr0:
        assert CORE_DIMS.issubset(set(hr0.sizes))


def test_create_hr1_python_notlc_writes_hr1_and_lr1(meop_config, stage_ct88_example) -> None:
    _stage_hr0_adjusted(meop_config, stage_ct88_example)
    result = create_hr1_python(
        meop_config,
        Selection(deployment="ct88", smru_name="ct88-225-12"),
        thermal_lag=False,
        now=TIMESTAMP,
    )

    assert result.processed_tags == ("ct88-225-12",)
    hr0_path = fname_prof("ct88-225-12", qf="hr0", config=meop_config)
    hr1_path = fname_prof("ct88-225-12", qf="hr1", config=meop_config)
    lr1_path = fname_prof("ct88-225-12", qf="lr1", config=meop_config)
    assert hr1_path.exists()
    assert lr1_path.exists()

    with _open_any(hr0_path) as hr0, _open_any(hr1_path) as hr1, _open_any(lr1_path) as lr1:
        assert hr1.sizes == hr0.sizes
        np.testing.assert_allclose(hr1["TEMP"].values, hr0["TEMP"].values, equal_nan=True)
        np.testing.assert_allclose(hr1["PSAL"].values, hr0["PSAL"].values, equal_nan=True)
        assert hr1.attrs["thermal_lag_adjustment"] == "no"
        assert hr1.attrs["density_inversion_adjustment"] in {"python-gsw", "python-fallback"}
        finite = np.isfinite(hr1["TEMP_ADJUSTED"].values) & np.isfinite(hr0["TEMP_ADJUSTED"].values)
        assert np.any(np.abs(hr1["TEMP_ADJUSTED"].values[finite] - hr0["TEMP_ADJUSTED"].values[finite]) > 0)
        assert lr1.attrs["thermal_lag_adjustment"] == "no"
        assert int(lr1.sizes["N_PROF"]) == int(hr0.sizes["N_PROF"])



def test_create_hr1_python_tlc_matches_ct88_lr1_reference_core_fields(meop_config, stage_ct88_example) -> None:
    staged = _stage_hr0_adjusted(meop_config, stage_ct88_example)
    result = create_hr1_python(
        meop_config,
        Selection(deployment="ct88", smru_name="ct88-225-12"),
        thermal_lag=True,
        now=TIMESTAMP,
    )

    assert result.processed_tags == ("ct88-225-12",)
    hr1_path = fname_prof("ct88-225-12", qf="hr1", config=meop_config)
    lr1_path = fname_prof("ct88-225-12", qf="lr1", config=meop_config)
    assert hr1_path.exists()
    assert lr1_path.exists()

    with _open_any(staged["reference_lr1"]) as reference, _open_any(lr1_path) as candidate:
        assert candidate.attrs["thermal_lag_adjustment"] == "yes"
        for dim in CORE_DIMS:
            assert int(candidate.sizes[dim]) == int(reference.sizes[dim])
        np.testing.assert_allclose(candidate["JULD"].values, reference["JULD"].values)
        np.testing.assert_allclose(candidate["LATITUDE"].values, reference["LATITUDE"].values)
        np.testing.assert_allclose(candidate["LONGITUDE"].values, reference["LONGITUDE"].values)
        np.testing.assert_allclose(candidate["PRES"].values, reference["PRES"].values, atol=1e-6, equal_nan=True)

        temp_diff = np.abs(candidate["TEMP_ADJUSTED"].values - reference["TEMP_ADJUSTED"].values)
        psal_diff = np.abs(candidate["PSAL_ADJUSTED"].values - reference["PSAL_ADJUSTED"].values)
        finite_temp = np.isfinite(temp_diff)
        finite_psal = np.isfinite(psal_diff)
        assert finite_temp.any()
        assert finite_psal.any()
        assert float(np.nanmean(temp_diff[finite_temp])) < 0.01
        assert float(np.nanmean(psal_diff[finite_psal])) < 0.02
        assert float(np.nanmax(temp_diff[finite_temp])) < 0.12
        assert float(np.nanmax(psal_diff[finite_psal])) < 1.0



def test_create_hr1_python_passes_level_major_pressure_to_stabiliser(meop_config, stage_ct88_example, monkeypatch) -> None:
    _stage_hr0_adjusted(meop_config, stage_ct88_example)

    import meop_process.processing.hr as hr_module

    class _FakeGsw:
        @staticmethod
        def CT_from_t(sa, temp, pres):
            _ = sa, pres
            return np.asarray(temp, dtype=np.float64)

    captured = {}

    def _fake_stabilise(sp, ct, p, **kwargs):
        _ = ct, kwargs
        captured['shape'] = np.asarray(sp).shape
        captured['p'] = np.asarray(p, dtype=np.float64)
        metadata = [None] * captured['shape'][1]
        return np.asarray(sp, dtype=np.float64), metadata

    monkeypatch.setattr(hr_module, '_gsw', _FakeGsw())
    monkeypatch.setattr(hr_module, 'stabilise_SP_const_CT', _fake_stabilise)

    create_hr1_python(
        meop_config,
        Selection(deployment='ct88', smru_name='ct88-225-12'),
        thermal_lag=False,
        now=TIMESTAMP,
    )

    assert captured['shape'][0] == 1000
    assert captured['shape'][1] == 306
    finite = np.isfinite(captured['p'])
    first_profile = captured['p'][:, 0]
    dp = np.diff(first_profile[np.isfinite(first_profile)])
    assert dp.size > 0
    assert np.all(dp > 0)


def test_create_hr1_python_records_skipped_stabilisation_profile_indices(meop_config, stage_ct88_example, monkeypatch) -> None:
    _stage_hr0_adjusted(meop_config, stage_ct88_example)

    import meop_process.processing.hr as hr_module

    def _fake_stabilise(sp, ct, p, **kwargs):
        _ = ct, p, kwargs
        metadata = [None] * sp.shape[1]
        metadata[7] = SolverInfo(
            solver="isotonic-gsw",
            status="skipped: numeric overflow while pooling density profile",
            success=False,
            profile_index=7,
        )
        for index in range(sp.shape[1]):
            if metadata[index] is None:
                metadata[index] = SolverInfo(
                    solver="isotonic-gsw",
                    status="already_stable",
                    success=True,
                    profile_index=index,
                )
        return np.asarray(sp, dtype=np.float64), metadata

    monkeypatch.setattr(hr_module, 'stabilise_SP_const_CT', _fake_stabilise)

    create_hr1_python(
        meop_config,
        Selection(deployment='ct88', smru_name='ct88-225-12'),
        thermal_lag=False,
        now=TIMESTAMP,
    )

    hr1_path = fname_prof("ct88-225-12", qf="hr1", config=meop_config)
    with _open_any(hr1_path) as hr1:
        assert hr1.attrs["stabilisation_skipped_profile_indices"] == "7"
        assert int(hr1.attrs["stabilisation_skipped_profile_count"]) == 1
