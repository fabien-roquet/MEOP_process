from __future__ import annotations

from datetime import datetime, timezone
import warnings

import numpy as np
import xarray as xr

from meop_process.catalog.deployments import load_info_deployment
from meop_process.io.raw_odv import OdvProfile
from meop_process.processing.ncargo import _build_ncargo_dataset
from meop_process.processing.qc import apply_lr0_qc_filters


def _qc_ready_dataset() -> xr.Dataset:
    return xr.Dataset(
        data_vars={
            "PRES": (("N_PROF", "N_LEVELS"), np.asarray([[1.0, 2.0]], dtype=np.float32)),
            "TEMP": (("N_PROF", "N_LEVELS"), np.asarray([[10.0, 9.0]], dtype=np.float32)),
            "PSAL": (("N_PROF", "N_LEVELS"), np.asarray([[34.0, 34.1]], dtype=np.float32)),
            "PRES_QC": (("N_PROF", "N_LEVELS"), np.asarray([["1", "1"]], dtype=object)),
            "TEMP_QC": (("N_PROF", "N_LEVELS"), np.asarray([["1", "1"]], dtype=object)),
            "PSAL_QC": (("N_PROF", "N_LEVELS"), np.asarray([["1", "1"]], dtype=object)),
            "PRES_ADJUSTED_QC": (("N_PROF", "N_LEVELS"), np.asarray([["1", "1"]], dtype=object)),
            "TEMP_ADJUSTED_QC": (("N_PROF", "N_LEVELS"), np.asarray([["1", "1"]], dtype=object)),
            "PSAL_ADJUSTED_QC": (("N_PROF", "N_LEVELS"), np.asarray([["1", "1"]], dtype=object)),
            "JULD": (("N_PROF",), np.asarray([12345.0], dtype=np.float64)),
            "LATITUDE": (("N_PROF",), np.asarray([np.nan], dtype=np.float64)),
            "LONGITUDE": (("N_PROF",), np.asarray([np.nan], dtype=np.float64)),
        },
        coords={"N_PROF": np.arange(1), "N_LEVELS": np.arange(2)},
    )


def test_build_ncargo_dataset_handles_missing_geolocation_without_warnings(meop_config, seed_catalog) -> None:
    seed_catalog(deployment="DEP001", smru_name="DEP001-AAA")
    info = load_info_deployment(meop_config, deployment="DEP001")
    profiles = [
        OdvProfile(
            smru_name="DEP001-AAA",
            station=1,
            timestamp="2020-01-01 00:00",
            longitude=np.nan,
            latitude=np.nan,
            pressure=(5.0, 10.0),
            temperature=(1.0, 0.5),
            salinity=(34.0, 34.1),
        )
    ]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("error", RuntimeWarning)
        dataset = _build_ncargo_dataset(meop_config, info, "DEP001-AAA", profiles, now=datetime.now(timezone.utc))

    assert caught == []
    assert np.isnan(dataset.attrs["geospatial_lat_min"])
    assert np.isnan(dataset.attrs["geospatial_lon_min"])


def test_apply_lr0_qc_filters_handles_all_nan_locations_without_warnings(meop_config, seed_catalog) -> None:
    seed_catalog(deployment="DEP001", smru_name="DEP001-AAA")
    info = load_info_deployment(meop_config, deployment="DEP001")
    dataset = _qc_ready_dataset()

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("error", RuntimeWarning)
        result = apply_lr0_qc_filters(meop_config, info, "DEP001-AAA", dataset)

    assert caught == []
    assert result.dataset.attrs["geospatial_lat_min"] == " "
    assert result.dataset.attrs["geospatial_lon_min"] == " "