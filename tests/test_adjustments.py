from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pytest
import xarray as xr

from meop_process.catalog.tables import read_csv_rows, write_indexed_csv_rows
from meop_process.catalog.filenames import fname_prof
from meop_process.io.raw_odv import import_raw_data_zip
from meop_process.metadata.patch import update_metadata_from_table
from meop_process.models import Selection
from meop_process.processing.adjustments import _ensure_default_coefficients, _load_salinity_offset_series, apply_adjustments
from meop_process.processing.hr import create_hr0_python
from meop_process.processing.ncargo import create_ncargo_python
from meop_process.workflows.compare import compare_netcdf_outputs


TIMESTAMP = datetime(2024, 3, 6, 2, 15, 16, tzinfo=timezone.utc)


ADJUSTMENT_VARIABLES = [
    "PRES_ADJUSTED",
    "TEMP_ADJUSTED",
    "PSAL_ADJUSTED",
    "TEMP_ADJUSTED_ERROR",
    "PSAL_ADJUSTED_ERROR",
    "SCIENTIFIC_CALIB_COEFFICIENT",
]



def _open_any(path):
    try:
        return xr.open_dataset(path, decode_times=False)
    except Exception:
        return xr.open_dataset(path, engine="scipy", decode_times=False)



def _stage_lr0(meop_config, stage_ct88_example):
    staged = stage_ct88_example()
    assert import_raw_data_zip(meop_config, "ct88") is True
    create_ncargo_python(
        meop_config,
        Selection(deployment="ct88", smru_name="ct88-225-12"),
        now=TIMESTAMP,
    )
    update_metadata_from_table(meop_config, smru_name="ct88-225-12", modes=("lr0",))
    return staged



def test_apply_adjustments_ct88_matches_reference_lr0_adjusted_fields(meop_config, stage_ct88_example) -> None:
    staged = _stage_lr0(meop_config, stage_ct88_example)

    result = apply_adjustments(
        meop_config,
        Selection(deployment="ct88", smru_name="ct88-225-12"),
    )

    assert result.processed_tags == ("ct88-225-12",)
    candidate = fname_prof("ct88-225-12", qf="lr0", config=meop_config)
    report = compare_netcdf_outputs(
        staged["reference_lr0"],
        candidate,
        variables=ADJUSTMENT_VARIABLES,
        attributes=(),
        atol=1e-5,
    )
    assert report.is_equal, report.variable_differences + report.attribute_differences + report.dimension_differences



def test_apply_adjustments_updates_hr0_from_lr0_coefficients_and_errors(meop_config, stage_ct88_example) -> None:
    _stage_lr0(meop_config, stage_ct88_example)
    create_hr0_python(
        meop_config,
        Selection(deployment="ct88", smru_name="ct88-225-12"),
        now=TIMESTAMP,
    )

    result = apply_adjustments(
        meop_config,
        Selection(deployment="ct88", smru_name="ct88-225-12"),
    )

    assert fname_prof("ct88-225-12", qf="lr0", config=meop_config) in result.written_files
    assert fname_prof("ct88-225-12", qf="hr0", config=meop_config) in result.written_files

    hr0_path = fname_prof("ct88-225-12", qf="hr0", config=meop_config)
    with _open_any(hr0_path) as hr0:
        psal = hr0["PSAL"].values
        psal_adjusted = hr0["PSAL_ADJUSTED"].values
        mask = np.isfinite(psal) & np.isfinite(psal_adjusted)
        assert mask.any()
        np.testing.assert_allclose(psal_adjusted[mask] - psal[mask], 0.4, atol=1e-5)

        temp_errors = hr0["TEMP_ADJUSTED_ERROR"].values
        psal_errors = hr0["PSAL_ADJUSTED_ERROR"].values
        assert np.isfinite(temp_errors).any()
        assert np.isfinite(psal_errors).any()
        unique_temp = np.unique(temp_errors[np.isfinite(temp_errors)])
        unique_psal = np.unique(psal_errors[np.isfinite(psal_errors)])
        assert all(np.isclose(value, 0.1) or np.isclose(value, 0.2) for value in unique_temp)
        assert all(np.isclose(value, 0.2) or np.isclose(value, 0.4) for value in unique_psal)


def test_ensure_default_coefficients_rewrites_plain_csv_columns(meop_config) -> None:
    meop_config.tablesdir.mkdir(parents=True, exist_ok=True)
    write_indexed_csv_rows(
        meop_config.tablesdir / "table_coeff.csv",
        [
            {
                "row_name": "existing-tag",
                "smru_platform_code": "existing-tag",
                "T1": "0",
                "T2": "0",
                "S1": "0",
                "S2": "0",
                "remove": "0",
                "Sremove": "0",
                "comment": "no comment",
            }
        ],
    )

    updated = _ensure_default_coefficients(meop_config, ["new-tag"])

    assert updated == meop_config.tablesdir / "table_coeff.csv"
    rows = read_csv_rows(updated)
    assert "" not in rows[0]
    assert any(row.get("smru_platform_code") == "new-tag" for row in rows)


def test_load_salinity_offset_series_interpolates_using_table_order(meop_config) -> None:
    meop_config.tablesdir.mkdir(parents=True, exist_ok=True)
    (meop_config.tablesdir / "table_salinity_offsets.csv").write_text(
        "smru_platform_code,index_1,index_2,index_3,index_4,offset_1,offset_2,offset_3,offset_4\n"
        "tag-1,1,3,0,0,0.0,0.2,0.5,0.5\n",
        encoding="utf-8",
    )

    offsets = _load_salinity_offset_series(meop_config, "tag-1", 5)

    np.testing.assert_allclose(offsets, np.array([0.0, 0.1, 0.2, 0.35, 0.5]))


def test_load_salinity_offset_series_uses_first_duplicate_row(meop_config) -> None:
    meop_config.tablesdir.mkdir(parents=True, exist_ok=True)
    (meop_config.tablesdir / "table_salinity_offsets.csv").write_text(
        "smru_platform_code,index_1,index_2,index_3,index_4,offset_1,offset_2,offset_3,offset_4\n"
        "tag-1,1,3,0,0,0.0,0.2,0.5,0.5\n"
        "tag-1,1,4,0,0,0.0,-0.1,-0.2,-0.2\n",
        encoding="utf-8",
    )

    offsets = _load_salinity_offset_series(meop_config, "tag-1", 5)

    np.testing.assert_allclose(offsets, np.array([0.0, 0.1, 0.2, 0.35, 0.5]))


def test_load_salinity_offset_series_sorts_non_monotonic_indices(meop_config) -> None:
    meop_config.tablesdir.mkdir(parents=True, exist_ok=True)
    (meop_config.tablesdir / "table_salinity_offsets.csv").write_text(
        "smru_platform_code,index_1,index_2,index_3,index_4,offset_1,offset_2,offset_3,offset_4\n"
        "tag-1,1,5,4,0,0.0,0.2,0.5,0.5\n",
        encoding="utf-8",
    )

    offsets = _load_salinity_offset_series(meop_config, "tag-1", 5)

    np.testing.assert_allclose(offsets, np.array([0.0, 0.16666667, 0.33333333, 0.5, 0.5]))
