from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import xarray as xr

from meop_process.catalog.filenames import fname_prof
from meop_process.io.raw_odv import import_raw_data_zip
from meop_process.metadata.patch import update_metadata_from_table
from meop_process.models import Selection
from meop_process.processing.ncargo import create_ncargo_python


TIMESTAMP = datetime(2024, 3, 6, 2, 15, 16, tzinfo=timezone.utc)


def _open_any(path):
    try:
        return xr.open_dataset(path, decode_times=False)
    except Exception:
        return xr.open_dataset(path, engine="scipy", decode_times=False)


def _stage_lr0(meop_config, stage_ct160_oxy_example):
    staged = stage_ct160_oxy_example()
    assert import_raw_data_zip(meop_config, "ct160") is True
    create_ncargo_python(
        meop_config,
        Selection(deployment="ct160", smru_name=staged["smru_name"]),
        now=TIMESTAMP,
    )
    update_metadata_from_table(meop_config, smru_name=staged["smru_name"], modes=("lr0",))
    return staged


def test_ct160_oxy_raw_odv_import_succeeds(meop_config, stage_ct160_oxy_example) -> None:
    stage_ct160_oxy_example()
    assert import_raw_data_zip(meop_config, "ct160") is True


def test_ct160_oxy_lr0_has_doxy_variable(meop_config, stage_ct160_oxy_example) -> None:
    staged = _stage_lr0(meop_config, stage_ct160_oxy_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    assert candidate.exists()
    with _open_any(candidate) as ds:
        assert "DOXY" in ds, "DOXY variable missing from lr0"
        assert "DOXY_QC" in ds
        assert "DOXY_ADJUSTED" in ds
        assert "DOXY_ADJUSTED_ERROR" in ds


def test_ct160_oxy_lr0_doxy_has_valid_values(meop_config, stage_ct160_oxy_example) -> None:
    staged = _stage_lr0(meop_config, stage_ct160_oxy_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(candidate) as ds:
        doxy = ds["DOXY"].values.astype(float)
        n_valid = int(np.sum(~np.all(np.isnan(doxy), axis=1)))
        assert n_valid > 0, "All DOXY profiles are NaN — expected real oxygen values"


def test_ct160_oxy_lr0_chla_not_present_when_fluorescence_all_fill(meop_config, stage_ct160_oxy_example) -> None:
    """CHLA should not be added when all fluorescence levels are 999 (fill value)."""
    staged = _stage_lr0(meop_config, stage_ct160_oxy_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(candidate) as ds:
        assert "CHLA" not in ds, "CHLA should be absent when all fluorescence values are fill"


def test_ct160_oxy_lr0_profile_count_matches_reference(meop_config, stage_ct160_oxy_example) -> None:
    staged = _stage_lr0(meop_config, stage_ct160_oxy_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(staged["reference_lr0"]) as ref, _open_any(candidate) as ds:
        assert int(ds.sizes["N_PROF"]) == int(ref.sizes["N_PROF"])


def test_ct160_oxy_lr0_juld_matches_reference(meop_config, stage_ct160_oxy_example) -> None:
    staged = _stage_lr0(meop_config, stage_ct160_oxy_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(staged["reference_lr0"]) as ref, _open_any(candidate) as ds:
        np.testing.assert_allclose(ds["JULD"].values[:5], ref["JULD"].values[:5])
