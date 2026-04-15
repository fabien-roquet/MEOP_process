from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import xarray as xr

from meop_process.catalog.filenames import fname_prof
from meop_process.io.raw_odv import import_raw_data_zip
from meop_process.metadata.patch import update_metadata_from_table
from meop_process.models import Selection
from meop_process.processing.adjustments import apply_adjustments
from meop_process.processing.hr import create_hr0_python, create_hr1_python
from meop_process.processing.ncargo import create_ncargo_python


TIMESTAMP = datetime(2024, 3, 5, 18, 6, 7, tzinfo=timezone.utc)


def _open_any(path):
    try:
        return xr.open_dataset(path, decode_times=False)
    except Exception:
        return xr.open_dataset(path, engine="scipy", decode_times=False)


def _stage_hr0_adjusted(meop_config, stage_ct78_example):
    staged = stage_ct78_example()
    assert import_raw_data_zip(meop_config, "ct78") is True
    create_ncargo_python(
        meop_config,
        Selection(deployment="ct78", smru_name="ct78-465-12"),
        now=TIMESTAMP,
    )
    update_metadata_from_table(meop_config, smru_name="ct78-465-12", modes=("lr0",))
    create_hr0_python(
        meop_config,
        Selection(deployment="ct78", smru_name="ct78-465-12"),
        now=TIMESTAMP,
    )
    apply_adjustments(
        meop_config,
        Selection(deployment="ct78", smru_name="ct78-465-12"),
    )
    return staged


def test_create_hr1_python_tlc_matches_ct78_lr1_reference_core_fields(meop_config, stage_ct78_example) -> None:
    staged = _stage_hr0_adjusted(meop_config, stage_ct78_example)
    result = create_hr1_python(
        meop_config,
        Selection(deployment="ct78", smru_name="ct78-465-12"),
        thermal_lag=True,
        now=TIMESTAMP,
    )

    assert result.processed_tags == ("ct78-465-12",)
    hr1_path = fname_prof("ct78-465-12", qf="hr1", config=meop_config)
    lr1_path = fname_prof("ct78-465-12", qf="lr1", config=meop_config)
    assert hr1_path.exists()
    assert lr1_path.exists()

    with _open_any(staged["reference_lr1"]) as reference, _open_any(lr1_path) as candidate:
        assert candidate.attrs["thermal_lag_adjustment"] == "yes"
        assert candidate.sizes == reference.sizes
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
        assert float(np.nanmean(temp_diff[finite_temp])) < 0.03
        assert float(np.nanmean(psal_diff[finite_psal])) < 0.03
        assert float(np.nanpercentile(temp_diff[finite_temp], 99)) < 0.10
        assert float(np.nanpercentile(psal_diff[finite_psal], 99)) < 0.10
        assert float(np.nanmax(temp_diff[finite_temp])) < 0.20
        assert float(np.nanmax(psal_diff[finite_psal])) < 1.0
