from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from meop_process.catalog.filenames import fname_prof, fname_traj
from meop_process.models import Selection
from meop_process.io.raw_odv import import_raw_data_zip
from meop_process.processing.fr0 import HrRawData, _pick_column, create_fr0_python, detect_profiles, load_hr_data
from meop_process.processing.ncargo import create_ncargo_python


TIMESTAMP = datetime(2024, 3, 6, 2, 15, 16, tzinfo=timezone.utc)
JULD_REF = datetime(1950, 1, 1, tzinfo=timezone.utc)


def _open_any(path: Path):
    try:
        return xr.open_dataset(path, decode_times=False)
    except Exception:
        return xr.open_dataset(path, engine="scipy", decode_times=False)


def _make_synthetic_dive_series(start: datetime, *, dives: int = 2) -> tuple[pd.DatetimeIndex, np.ndarray, np.ndarray, np.ndarray]:
    timestamps: list[datetime] = []
    pressure: list[float] = []
    temp: list[float] = []
    sal: list[float] = []

    current = start
    for _ in range(dives):
        def _append_block(values: list[float]) -> None:
            nonlocal current
            for value in values:
                timestamps.append(current)
                pressure.append(float(value))
                temp.append(6.0 - 0.02 * float(value))
                sal.append(33.4 + 0.01 * float(value))
                current += timedelta(seconds=1)

        _append_block([0.0] * 30)
        _append_block(np.linspace(0.0, 120.0, 180).tolist())
        _append_block((120.0 + 0.2 * np.sin(np.linspace(0.0, 4.0 * np.pi, 40))).tolist())
        _append_block(np.linspace(120.0, 0.0, 180).tolist())
        _append_block([0.0] * 45)

    return pd.DatetimeIndex(timestamps, tz="UTC"), np.asarray(pressure), np.asarray(temp), np.asarray(sal)


def _seed_table_param_min_profiles_one(meop_config) -> None:
    table_param = meop_config.tablesdir / "table_param.csv"
    table_param.write_text(
        "row_name,temp_error,psal_error,minT,maxT,minS,maxS,min_Nprof,pmax,pmax_fluo,is_lon_centre_180\n"
        "ct96,0.1,0.2,-3,32,4,40,1,1000,200,0\n",
        encoding="utf-8",
    )


def _stage_ct96_lr0(meop_config, stage_ct96_example):
    staged = stage_ct96_example()
    assert import_raw_data_zip(meop_config, "ct96") is True
    create_ncargo_python(
        meop_config,
        Selection(deployment="ct96", smru_name="ct96-24-13"),
        now=TIMESTAMP,
    )
    _seed_table_param_min_profiles_one(meop_config)
    return staged


def test_pick_column_pressure_re_surface_variants() -> None:
    # Both spaced and no-space variants of PRESSURE_RE_SURFACE must be detected
    for col in ("PRESSURE_RE_SURFACE (dbar)", "PRESSURE_RE_SURFACE(dbar)"):
        assert _pick_column([col], "pressure") == col, f"Expected {col!r} to be picked"
    # ABS_PRESSURE_AT_SURFACE must still be excluded
    assert _pick_column(["ABS_PRESSURE_AT_SURFACE (dbar)"], "pressure") is None


def test_load_hr_data_skips_malformed_rows(tmp_path) -> None:
    path = tmp_path / "bad_hr.txt"
    path.write_text(
        "datetime\tpressure\ttemperature\tsalinity\n"
        "2020-01-01 00:00:00\t1\t2.0\t33.1\n"
        "2020-01-01 00:00:01\t2\t2.1\t33.2\textra\n"
        "2020-01-01 00:00:02\t3\t2.2\t33.3\n",
        encoding="utf-8",
    )

    with pytest.warns(UserWarning, match="Malformed HR rows skipped"):
        data = load_hr_data(path, continuous=True)

    assert data.timestamp.size == 2
    np.testing.assert_allclose(data.pressure, np.asarray([1.0, 3.0]))
    # Standard variants must still work
    assert _pick_column(["PRESSURE (dbar)"], "pressure") == "PRESSURE (dbar)"
    assert _pick_column(["corrected depth"], "pressure") == "corrected depth"


def test_detect_profiles_matches_continuous_dive_logic() -> None:
    timestamps, pressure, temp, sal = _make_synthetic_dive_series(datetime(2013, 2, 13, 16, 0, 0, tzinfo=timezone.utc))
    juld = np.asarray([(ts.to_pydatetime() - JULD_REF).total_seconds() / 86400.0 for ts in timestamps], dtype=np.float64)
    raw = HrRawData(
        timestamp=timestamps,
        juld=juld,
        pressure=pressure,
        temperature=temp,
        salinity=sal,
        conductivity=np.full(pressure.shape, np.nan, dtype=np.float64),
        fluorescence=np.full(pressure.shape, np.nan, dtype=np.float64),
        oxygen=np.full(pressure.shape, np.nan, dtype=np.float64),
        light=np.full(pressure.shape, np.nan, dtype=np.float64),
        continuous=True,
        isfluo=False,
        isoxy=False,
        islight=False,
    )

    detected = detect_profiles(raw)

    assert detected.n_profiles == 2
    assert np.all(detected.profile_end_index > detected.ascent_start_index)
    assert int(detected.ascent_start_index[0]) > 150
    assert int(detected.profile_end_index[0] - detected.ascent_start_index[0]) > 100


def test_create_fr0_python_matches_ct96_shortened_reference_core_fields(meop_config, stage_ct96_example) -> None:
    staged = _stage_ct96_lr0(meop_config, stage_ct96_example)

    result = create_fr0_python(
        meop_config,
        Selection(deployment="ct96", smru_name="ct96-24-13"),
        now=TIMESTAMP,
    )

    assert result.processed_tags == ("ct96-24-13",)
    fr0_path = fname_prof("ct96-24-13", qf="fr0", config=meop_config)
    assert fr0_path.exists()

    with _open_any(fr0_path) as dataset, _open_any(staged["reference_fr0"]) as reference:
        assert int(dataset.sizes["N_PROF"]) == 20
        assert int(dataset.sizes["N_LEVELS"]) == 1000
        assert dataset.attrs["thermal_lag_adjustment"] == "no"
        assert dataset.attrs["profile_source"] == "full-resolution"

        np.testing.assert_allclose(dataset["JULD"].values, reference["JULD"].values, atol=1e-4)
        np.testing.assert_allclose(dataset["LATITUDE"].values, reference["LATITUDE"].values, atol=1e-3)
        np.testing.assert_allclose(dataset["LONGITUDE"].values, reference["LONGITUDE"].values, atol=1e-3)
        temp_finite = np.isfinite(dataset["TEMP"].values)
        ref_temp_finite = np.isfinite(reference["TEMP"].values)
        psal_finite = np.isfinite(dataset["PSAL"].values)
        ref_psal_finite = np.isfinite(reference["PSAL"].values)
        assert int(np.count_nonzero(temp_finite != ref_temp_finite)) <= 20
        assert int(np.count_nonzero(psal_finite != ref_psal_finite)) <= 20
        np.testing.assert_allclose(temp_finite.sum(axis=1), ref_temp_finite.sum(axis=1), atol=2)
        np.testing.assert_allclose(psal_finite.sum(axis=1), ref_psal_finite.sum(axis=1), atol=2)

        temp_mask = temp_finite & ref_temp_finite
        psal_mask = psal_finite & ref_psal_finite
        np.testing.assert_allclose(dataset["TEMP"].values[temp_mask], reference["TEMP"].values[temp_mask], atol=6e-2)
        np.testing.assert_allclose(dataset["PSAL"].values[psal_mask], reference["PSAL"].values[psal_mask], atol=6e-2)
        temp_qc = dataset["TEMP_QC"].values.astype(str)
        ref_temp_qc = reference["TEMP_QC"].values.astype(str)
        psal_qc = dataset["PSAL_QC"].values.astype(str)
        ref_psal_qc = reference["PSAL_QC"].values.astype(str)
        assert int(np.count_nonzero(temp_qc != ref_temp_qc)) <= 20
        assert int(np.count_nonzero(psal_qc != ref_psal_qc)) <= 20
        np.testing.assert_array_equal(dataset["PROFILE_TEMP_QC"].values.astype(str), reference["PROFILE_TEMP_QC"].values.astype(str))
        np.testing.assert_array_equal(dataset["PROFILE_PSAL_QC"].values.astype(str), reference["PROFILE_PSAL_QC"].values.astype(str))


def test_create_fr0_python_uses_only_hr_catalog_filename_and_skips_missing_raw_file(meop_config, stage_ct96_example) -> None:
    stage_ct96_example()

    # Strict behavior: the HR filename comes only from list_deployment_hr.csv. If that catalog row
    # points to a missing raw file, the FR pipeline must skip processing instead of guessing.
    (meop_config.catalogdir / "list_deployment_hr.csv").write_text(
        ",smru_platform_code,instr_id,year,prefix,continuous\n"
        "ct96-24-13,ct96-24-13,12664,2013,,1\n",
        encoding="utf-8",
    )

    # Remove the HR raw file to verify that missing files cause skipping
    hr_file = meop_config.raw_hr_dir / "2013" / "12664_ctd.txt"
    if hr_file.exists():
        hr_file.unlink()

    assert import_raw_data_zip(meop_config, "ct96") is True
    create_ncargo_python(
        meop_config,
        Selection(deployment="ct96", smru_name="ct96-24-13"),
        now=TIMESTAMP,
    )
    _seed_table_param_min_profiles_one(meop_config)

    result = create_fr0_python(
        meop_config,
        Selection(deployment="ct96", smru_name="ct96-24-13"),
        now=TIMESTAMP,
    )

    assert result.processed_tags == ()
    fr0_path = fname_prof("ct96-24-13", qf="fr0", config=meop_config)
    assert not fr0_path.exists()


def test_create_fr0_python_writes_ct96_traj_file(meop_config, stage_ct96_example) -> None:
    staged = _stage_ct96_lr0(meop_config, stage_ct96_example)
    create_fr0_python(
        meop_config,
        Selection(deployment="ct96", smru_name="ct96-24-13"),
        now=TIMESTAMP,
    )

    traj_path = fname_traj("ct96-24-13", config=meop_config)
    assert traj_path.exists()

    expected_samples = sum(1 for _ in staged["hr_file"].open(encoding="utf-8")) - 1
    with _open_any(traj_path) as dataset:
        assert int(dataset.sizes["TIME"]) == expected_samples
        assert {"TIME", "LATITUDE", "LONGITUDE", "PRES", "TEMP", "PSAL"}.issubset(dataset.variables)
        assert np.isfinite(dataset["LATITUDE"].values).all()
        assert np.isfinite(dataset["LONGITUDE"].values).all()
        assert dataset.attrs["profile_source"] == "full-resolution-traj"
