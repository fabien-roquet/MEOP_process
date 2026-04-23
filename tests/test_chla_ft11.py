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


def _stage_lr0(meop_config, stage_ft11_chla_example):
    staged = stage_ft11_chla_example()
    assert import_raw_data_zip(meop_config, "ft11") is True
    create_ncargo_python(
        meop_config,
        Selection(deployment="ft11", smru_name=staged["smru_name"]),
        now=TIMESTAMP,
    )
    update_metadata_from_table(meop_config, smru_name=staged["smru_name"], modes=("lr0",))
    return staged


def test_ft11_raw_odv_import_succeeds(meop_config, stage_ft11_chla_example) -> None:
    stage_ft11_chla_example()
    assert import_raw_data_zip(meop_config, "ft11") is True


def test_ft11_lr0_has_chla_but_no_doxy(meop_config, stage_ft11_chla_example) -> None:
    staged = _stage_lr0(meop_config, stage_ft11_chla_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    assert candidate.exists()
    with _open_any(candidate) as ds:
        assert "CHLA" in ds, "CHLA variable missing from lr0"
        assert "DOXY" not in ds, "DOXY should not be present — ft11 oxygen values are all fill values"


def test_ft11_lr0_chla_has_valid_values(meop_config, stage_ft11_chla_example) -> None:
    """Fluorescence values from the ODV should be decoded into CHLA with real (non-NaN) values."""
    staged = _stage_lr0(meop_config, stage_ft11_chla_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(candidate) as ds:
        chla = ds["CHLA"].values.astype(float)
        n_valid = int(np.sum(~np.all(np.isnan(chla), axis=1)))
        assert n_valid == int(ds.sizes["N_PROF"]), (
            f"Expected all {ds.sizes['N_PROF']} profiles to have CHLA data, got {n_valid}"
        )


def test_ft11_lr0_chla_adjusted_present(meop_config, stage_ft11_chla_example) -> None:
    staged = _stage_lr0(meop_config, stage_ft11_chla_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(candidate) as ds:
        assert "CHLA_ADJUSTED" in ds
        assert "CHLA_ADJUSTED_QC" in ds
        assert "CHLA_ADJUSTED_ERROR" in ds


def test_ft11_lr0_profile_count_matches_reference(meop_config, stage_ft11_chla_example) -> None:
    staged = _stage_lr0(meop_config, stage_ft11_chla_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(staged["reference_lr0"]) as ref, _open_any(candidate) as ds:
        assert int(ds.sizes["N_PROF"]) == int(ref.sizes["N_PROF"])


def test_ft11_lr0_juld_matches_reference(meop_config, stage_ft11_chla_example) -> None:
    staged = _stage_lr0(meop_config, stage_ft11_chla_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(staged["reference_lr0"]) as ref, _open_any(candidate) as ds:
        np.testing.assert_allclose(ds["JULD"].values[:5], ref["JULD"].values[:5])
