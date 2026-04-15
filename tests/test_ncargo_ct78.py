from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import xarray as xr

from meop_process.catalog.filenames import fname_prof
from meop_process.io.raw_odv import import_raw_data_zip
from meop_process.metadata.patch import update_metadata_from_table
from meop_process.models import Selection
from meop_process.processing.ncargo import create_ncargo_python
from meop_process.workflows.compare import compare_netcdf_outputs


CORE_VARIABLES = [
    "PLATFORM_NUMBER",
    "JULD",
    "LATITUDE",
    "LONGITUDE",
    "PRES",
    "TEMP",
    "PSAL",
    "STATION_PARAMETERS",
    "PARAMETER",
]

CORE_ATTRIBUTES = [
    "platform_code",
    "wmo_platform_code",
    "smru_platform_code",
    "deployment_code",
    "species",
    "time_coverage_start",
    "time_coverage_end",
    "location",
    "loc_algorithm",
    "firmware_version",
    "firmware_parameters",
    "instr_id",
    "ptt",
    "nation",
    "PI",
]


def _open_any(path):
    try:
        return xr.open_dataset(path, decode_times=False)
    except Exception:
        return xr.open_dataset(path, engine="scipy", decode_times=False)


def test_create_ncargo_python_writes_ct78_lr0_core_file(meop_config, stage_ct78_example) -> None:
    staged = stage_ct78_example()
    assert import_raw_data_zip(meop_config, "ct78") is True

    result = create_ncargo_python(
        meop_config,
        Selection(deployment="ct78", smru_name="ct78-465-12"),
        now=datetime(2024, 3, 5, 18, 6, 7, tzinfo=timezone.utc),
    )
    update_metadata_from_table(meop_config, smru_name="ct78-465-12", modes=("lr0",))

    assert result.processed_tags == ("ct78-465-12",)
    candidate = fname_prof("ct78-465-12", qf="lr0", config=meop_config)
    assert candidate.exists()

    with _open_any(staged["reference_lr0"]) as reference, _open_any(candidate) as created:
        assert int(created.sizes["N_PROF"]) == 286
        assert int(created.sizes["N_LEVELS"]) == 17
        np.testing.assert_allclose(created["JULD"].values[:5], reference["JULD"].values[:5])
        np.testing.assert_allclose(created["PRES"].values[0, :], reference["PRES"].values[0, :])
        np.testing.assert_allclose(created["TEMP"].values[0, :], reference["TEMP"].values[0, :], atol=1e-4)
        np.testing.assert_allclose(created["PSAL"].values[0, :], reference["PSAL"].values[0, :], atol=1e-4)
        assert created.attrs["platform_code"] == reference.attrs["platform_code"]
        assert created.attrs["wmo_platform_code"] == reference.attrs["wmo_platform_code"]
        assert created.attrs["smru_platform_code"] == reference.attrs["smru_platform_code"]
        assert created.attrs["species"] == reference.attrs["species"]
        assert created.attrs["location"] == reference.attrs["location"]
        assert created.attrs["loc_algorithm"] == reference.attrs["loc_algorithm"]


def test_ct78_core_subset_matches_reference_dataset(meop_config, stage_ct78_example) -> None:
    staged = stage_ct78_example()
    assert import_raw_data_zip(meop_config, "ct78") is True
    create_ncargo_python(
        meop_config,
        Selection(deployment="ct78", smru_name="ct78-465-12"),
        now=datetime(2024, 3, 5, 18, 6, 7, tzinfo=timezone.utc),
    )
    update_metadata_from_table(meop_config, smru_name="ct78-465-12", modes=("lr0",))
    candidate = fname_prof("ct78-465-12", qf="lr0", config=meop_config)

    report = compare_netcdf_outputs(
        staged["reference_lr0"],
        candidate,
        variables=CORE_VARIABLES,
        attributes=CORE_ATTRIBUTES,
        atol=1e-4,
    )

    assert report.is_equal, report.variable_differences + report.attribute_differences + report.dimension_differences
