from __future__ import annotations

import xarray as xr

from meop_process.catalog.filenames import fname_prof
from meop_process.metadata.patch import update_metadata_from_table


def test_update_metadata_patches_all_matching_files(meop_config) -> None:
    tag = "DEP001-AAA"
    deployment_dir = meop_config.processdir / "final_dataset_prof" / "DEP001"
    deployment_dir.mkdir(parents=True, exist_ok=True)

    for qf in ("lr0", "hr0", "fr0"):
        path = fname_prof(tag, qf=qf, config=meop_config)
        xr.Dataset({"TEMP": ("N_LEVELS", [1.0, 2.0])}).to_netcdf(path, engine="h5netcdf")

    table_meta = meop_config.tablesdir / "table_meta.csv"
    table_meta.write_text(
        "smru_platform_code,location,principal_investigator\n"
        "DEP001-AAA,South Ocean,Dr Example\n",
        encoding="utf-8",
    )

    updated = update_metadata_from_table(meop_config, deployment="DEP001")

    assert len(updated) == 3
    for qf in ("lr0", "hr0", "fr0"):
        path = fname_prof(tag, qf=qf, config=meop_config)
        with xr.open_dataset(path, engine="h5netcdf") as dataset:
            assert dataset.attrs["location"] == "South Ocean"
            assert dataset.attrs["principal_investigator"] == "Dr Example"
