from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import xarray as xr

from meop_process.catalog.filenames import fname_prof
from meop_process.io.raw_odv import import_raw_data_zip
from meop_process.metadata.patch import update_metadata_from_table
from meop_process.models import Selection
from meop_process.processing.adjustments import apply_adjustments
from meop_process.processing.ncargo import create_ncargo_python


TIMESTAMP = datetime(2024, 3, 6, 2, 15, 16, tzinfo=timezone.utc)


def _open_any(path):
    try:
        return xr.open_dataset(path, decode_times=False)
    except Exception:
        return xr.open_dataset(path, engine="scipy", decode_times=False)


def _stage_lr0(meop_config, stage_ct153_lr_example):
    staged = stage_ct153_lr_example()
    assert import_raw_data_zip(meop_config, "ct153") is True
    create_ncargo_python(
        meop_config,
        Selection(deployment="ct153", smru_name=staged["smru_name"]),
        now=TIMESTAMP,
    )
    update_metadata_from_table(meop_config, smru_name=staged["smru_name"], modes=("lr0",))
    # Seed calibration coefficients matching the reference lr1 (T2=-0.08, S2=+0.20)
    meop_config.tablesdir.mkdir(parents=True, exist_ok=True)
    coeff_row = (
        "smru_platform_code,T1,T2,S1,S2,remove,Sremove,comment\n"
        f"{staged['smru_name']},0,-0.08,0,0.20,0,0,reference coefficients\n"
    )
    (meop_config.tablesdir / "table_coeff.csv").write_text(coeff_row, encoding="utf-8")
    return staged


def test_ct153_lr_raw_odv_import_succeeds(meop_config, stage_ct153_lr_example) -> None:
    stage_ct153_lr_example()
    assert import_raw_data_zip(meop_config, "ct153") is True


def test_ct153_lr_ncargo_produces_lr0_file(meop_config, stage_ct153_lr_example) -> None:
    staged = _stage_lr0(meop_config, stage_ct153_lr_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    assert candidate.exists()


def test_ct153_lr0_has_only_ctd_variables(meop_config, stage_ct153_lr_example) -> None:
    """LR-only tag should not have CHLA or DOXY in lr0."""
    staged = _stage_lr0(meop_config, stage_ct153_lr_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(candidate) as ds:
        assert "TEMP" in ds
        assert "PSAL" in ds
        assert "PRES" in ds
        assert "CHLA" not in ds
        assert "DOXY" not in ds


def test_ct153_lr0_profile_count_matches_reference(meop_config, stage_ct153_lr_example) -> None:
    staged = _stage_lr0(meop_config, stage_ct153_lr_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(staged["reference_lr0"]) as ref, _open_any(candidate) as ds:
        assert int(ds.sizes["N_PROF"]) == int(ref.sizes["N_PROF"])


def test_ct153_lr0_juld_matches_reference(meop_config, stage_ct153_lr_example) -> None:
    staged = _stage_lr0(meop_config, stage_ct153_lr_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(staged["reference_lr0"]) as ref, _open_any(candidate) as ds:
        np.testing.assert_allclose(ds["JULD"].values[:5], ref["JULD"].values[:5])


def test_ct153_lr0_temp_and_psal_match_reference(meop_config, stage_ct153_lr_example) -> None:
    staged = _stage_lr0(meop_config, stage_ct153_lr_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(staged["reference_lr0"]) as ref, _open_any(candidate) as ds:
        np.testing.assert_allclose(
            ds["TEMP"].values[0, :], ref["TEMP"].values[0, :], atol=1e-4,
        )
        np.testing.assert_allclose(
            ds["PSAL"].values[0, :], ref["PSAL"].values[0, :], atol=1e-4,
        )


def test_ct153_lr_apply_adjustments_produces_lr1(meop_config, stage_ct153_lr_example) -> None:
    """apply_adjustments should work on an LR-only tag (no hr2 required)."""
    staged = _stage_lr0(meop_config, stage_ct153_lr_example)
    result = apply_adjustments(
        meop_config,
        Selection(deployment="ct153", smru_name=staged["smru_name"]),
    )
    assert result.processed_tags == (staged["smru_name"],)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(candidate) as ds:
        assert "TEMP_ADJUSTED" in ds
        assert "PSAL_ADJUSTED" in ds
        assert "PRES_ADJUSTED" in ds


def _apply_and_open_lr0(meop_config, staged):
    """Run apply_adjustments on a staged lr0 and return the opened dataset."""
    apply_adjustments(meop_config, Selection(deployment="ct153", smru_name=staged["smru_name"]))
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    return _open_any(candidate)


def _psal_median_at_pressure(ds: xr.Dataset, target_p: float) -> float:
    """Return the median PSAL_ADJUSTED across all profiles with finite data near target_p."""
    pres = ds["PRES_ADJUSTED"].values.astype(float)
    psal = ds["PSAL_ADJUSTED"].values.astype(float)
    pres_1d = np.nanmedian(pres, axis=0)
    col = int(np.argmin(np.abs(pres_1d - target_p)))
    col_vals = psal[:, col]
    finite_vals = col_vals[np.isfinite(col_vals)]
    assert len(finite_vals) > 0, f"No finite PSAL_ADJUSTED values at ~{target_p} dbar"
    return float(np.median(finite_vals))


def test_ct153_psal_adjusted_median_matches_reference_at_50dbar(meop_config, stage_ct153_lr_example) -> None:
    """Median PSAL_ADJUSTED at ~50 dbar should match the reference lr1 within 0.05 PSU."""
    staged = _stage_lr0(meop_config, stage_ct153_lr_example)
    with _apply_and_open_lr0(meop_config, staged) as ds, _open_any(staged["reference_lr1"]) as ref:
        candidate_median = _psal_median_at_pressure(ds, 50.0)
        reference_median = _psal_median_at_pressure(ref, 50.0)
        np.testing.assert_allclose(
            candidate_median, reference_median, atol=0.05,
            err_msg="Median PSAL_ADJUSTED at ~50 dbar diverges from reference lr1",
        )


def test_ct153_psal_adjusted_median_matches_reference_at_100dbar(meop_config, stage_ct153_lr_example) -> None:
    """Median PSAL_ADJUSTED at ~100 dbar should match the reference lr1 within 0.05 PSU."""
    staged = _stage_lr0(meop_config, stage_ct153_lr_example)
    with _apply_and_open_lr0(meop_config, staged) as ds, _open_any(staged["reference_lr1"]) as ref:
        candidate_median = _psal_median_at_pressure(ds, 100.0)
        reference_median = _psal_median_at_pressure(ref, 100.0)
        np.testing.assert_allclose(
            candidate_median, reference_median, atol=0.05,
            err_msg="Median PSAL_ADJUSTED at ~100 dbar diverges from reference lr1",
        )


def test_ct153_psal_adjusted_median_matches_reference_at_200dbar(meop_config, stage_ct153_lr_example) -> None:
    """Median PSAL_ADJUSTED at ~200 dbar should match the reference lr1 within 0.05 PSU."""
    staged = _stage_lr0(meop_config, stage_ct153_lr_example)
    with _apply_and_open_lr0(meop_config, staged) as ds, _open_any(staged["reference_lr1"]) as ref:
        candidate_median = _psal_median_at_pressure(ds, 200.0)
        reference_median = _psal_median_at_pressure(ref, 200.0)
        np.testing.assert_allclose(
            candidate_median, reference_median, atol=0.05,
            err_msg="Median PSAL_ADJUSTED at ~200 dbar diverges from reference lr1",
        )
