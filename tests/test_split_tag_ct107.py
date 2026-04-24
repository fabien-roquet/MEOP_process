from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import xarray as xr

from meop_process.catalog.filenames import fname_prof
from meop_process.io.raw_odv import import_raw_data_zip, load_raw_odv_profiles, discover_raw_odv_files
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


def _stage_lr0_n2(meop_config, stage_ct107_split_example):
    staged = stage_ct107_split_example()
    assert import_raw_data_zip(meop_config, "ct107") is True
    create_ncargo_python(
        meop_config,
        Selection(deployment="ct107", smru_name=staged["smru_name"]),
        now=TIMESTAMP,
    )
    update_metadata_from_table(meop_config, smru_name=staged["smru_name"], modes=("lr0",))
    # Seed calibration coefficients matching the reference lr1 (S2=-0.02)
    meop_config.tablesdir.mkdir(parents=True, exist_ok=True)
    coeff_row = (
        "smru_platform_code,T1,T2,S1,S2,remove,Sremove,comment\n"
        f"{staged['smru_name']},0,0,0,-0.02,0,0,reference coefficients\n"
    )
    (meop_config.tablesdir / "table_coeff.csv").write_text(coeff_row, encoding="utf-8")
    return staged


def test_ct107_split_raw_odv_import_succeeds(meop_config, stage_ct107_split_example) -> None:
    stage_ct107_split_example()
    assert import_raw_data_zip(meop_config, "ct107") is True


def test_ct107_split_load_profiles_contains_n1_and_n2_suffixes(meop_config, stage_ct107_split_example) -> None:
    """load_raw_odv_profiles (which applies split rules) should produce -N1 and -N2 variants."""
    stage_ct107_split_example()
    assert import_raw_data_zip(meop_config, "ct107") is True
    files = discover_raw_odv_files(meop_config, "ct107")
    profiles = load_raw_odv_profiles(files, config=meop_config)
    smru_names = {p.smru_name for p in profiles}
    assert "ct107-938-13-N1" in smru_names, f"Expected -N1 split tag, got: {smru_names}"
    assert "ct107-938-13-N2" in smru_names, f"Expected -N2 split tag, got: {smru_names}"


def test_ct107_split_n1_plus_n2_profile_count_equals_raw_total(meop_config, stage_ct107_split_example) -> None:
    """N1 + N2 profile counts should equal the 40 raw profiles in the fixture ODV."""
    stage_ct107_split_example()
    assert import_raw_data_zip(meop_config, "ct107") is True
    files = discover_raw_odv_files(meop_config, "ct107")
    profiles = load_raw_odv_profiles(files, config=meop_config)
    n1 = sum(1 for p in profiles if p.smru_name == "ct107-938-13-N1")
    n2 = sum(1 for p in profiles if p.smru_name == "ct107-938-13-N2")
    assert n1 + n2 == 40, f"Expected 40 total profiles (N1={n1} + N2={n2})"
    assert n2 > 0, "N2 half is empty"


def test_ct107_split_ncargo_produces_n2_lr0(meop_config, stage_ct107_split_example) -> None:
    staged = _stage_lr0_n2(meop_config, stage_ct107_split_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    assert candidate.exists(), f"Expected lr0 file to be created at {candidate}"


def test_ct107_split_n2_lr0_has_nonzero_profiles(meop_config, stage_ct107_split_example) -> None:
    staged = _stage_lr0_n2(meop_config, stage_ct107_split_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(candidate) as ds:
        assert int(ds.sizes["N_PROF"]) > 0


def test_ct107_split_n2_lr0_smru_name_attribute(meop_config, stage_ct107_split_example) -> None:
    """The created lr0 should carry the -N2 smru_platform_code, not the base name."""
    staged = _stage_lr0_n2(meop_config, stage_ct107_split_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(candidate) as ds:
        assert ds.attrs.get("smru_platform_code") == staged["smru_name"], (
            f"Expected smru_platform_code='{staged['smru_name']}', "
            f"got '{ds.attrs.get('smru_platform_code')}'"
        )


def test_ct107_split_n2_lr0_temp_psal_pres_present(meop_config, stage_ct107_split_example) -> None:
    staged = _stage_lr0_n2(meop_config, stage_ct107_split_example)
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(candidate) as ds:
        assert "TEMP" in ds
        assert "PSAL" in ds
        assert "PRES" in ds


def _apply_and_open_lr0_n2(meop_config, staged):
    """Run apply_adjustments on a staged N2 lr0 and return the opened dataset."""
    apply_adjustments(meop_config, Selection(deployment="ct107", smru_name=staged["smru_name"]))
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    return _open_any(candidate)


def test_ct107_split_n2_apply_adjustments_produces_psal_adjusted(meop_config, stage_ct107_split_example) -> None:
    """apply_adjustments should produce PSAL_ADJUSTED in the lr0 file for the N2 split tag."""
    staged = _stage_lr0_n2(meop_config, stage_ct107_split_example)
    result = apply_adjustments(
        meop_config,
        Selection(deployment="ct107", smru_name=staged["smru_name"]),
    )
    assert staged["smru_name"] in result.processed_tags
    candidate = fname_prof(staged["smru_name"], qf="lr0", config=meop_config)
    with _open_any(candidate) as ds:
        assert "TEMP_ADJUSTED" in ds
        assert "PSAL_ADJUSTED" in ds
        assert "PRES_ADJUSTED" in ds


def test_ct107_split_n2_psal_adjusted_applies_s2_offset(meop_config, stage_ct107_split_example) -> None:
    """PSAL_ADJUSTED should equal PSAL - S2 (seeded offset) within float32 precision."""
    staged = _stage_lr0_n2(meop_config, stage_ct107_split_example)
    s2 = -0.02
    with _apply_and_open_lr0_n2(meop_config, staged) as ds:
        psal_raw = ds["PSAL"].values.astype(float)
        psal_adj = ds["PSAL_ADJUSTED"].values.astype(float)
        finite = np.isfinite(psal_raw) & np.isfinite(psal_adj)
        assert finite.sum() > 0, "No finite PSAL / PSAL_ADJUSTED values in candidate"
        expected = psal_raw[finite] - s2
        np.testing.assert_allclose(
            psal_adj[finite], expected, atol=1e-4,
            err_msg="PSAL_ADJUSTED should equal PSAL minus S2 offset",
        )


def test_ct107_split_n2_psal_adjusted_in_physical_range(meop_config, stage_ct107_split_example) -> None:
    """All finite PSAL_ADJUSTED values should fall within the physical salinity range."""
    staged = _stage_lr0_n2(meop_config, stage_ct107_split_example)
    with _apply_and_open_lr0_n2(meop_config, staged) as ds:
        psal_adj = ds["PSAL_ADJUSTED"].values.astype(float)
        finite_vals = psal_adj[np.isfinite(psal_adj)]
        assert len(finite_vals) > 0, "No finite PSAL_ADJUSTED values"
        assert float(finite_vals.min()) > 30.0, f"PSAL_ADJUSTED too low: {finite_vals.min()}"
        assert float(finite_vals.max()) < 38.0, f"PSAL_ADJUSTED too high: {finite_vals.max()}"

