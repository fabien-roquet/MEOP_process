from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import xarray as xr

from meop_process.reference.cora import load_cora_tiles
from meop_process.reference.cora_download import (
    CORA_DATASET_ID,
    EASYCORA_DATASET_ID,
    RemoteFile,
    create_download_plan,
    download_plan,
    select_remote_files,
)
from meop_process.reference.cora_prepare import prepare_cora_tiles


def test_select_remote_files_requires_global_and_filters_year() -> None:
    rows = [
        RemoteFile("s3://bucket/dataset_202511/Global/CO_X_20180101_PR_CT.nc", 10),
        RemoteFile("s3://bucket/dataset_202511/Global/CO_X_20240101_PR_OC.nc", 20),
        RemoteFile("s3://bucket/dataset_202511/Arctic/CO_X_20240101_PR_CT.nc", 30),
        RemoteFile("s3://bucket/dataset_202511/Global/CO_X_20240101_PR_PF.nc", 40),
    ]

    selected = select_remote_files(
        rows,
        file_classes=("CT", "OC", "TE"),
        start_year=2020,
        end_year=2025,
    )

    assert [Path(row.filename).name for row in selected] == ["CO_X_20240101_PR_OC.nc"]


class _FakeCopernicusMarine:
    def __init__(self) -> None:
        self.downloads: list[dict[str, object]] = []

    def get(self, **kwargs: object) -> None:
        if "create_file_list" in kwargs:
            output = Path(kwargs["output_directory"]) / str(kwargs["create_file_list"])
            dataset_id = str(kwargs["dataset_id"])
            if dataset_id == CORA_DATASET_ID:
                names = [
                    "s3://bucket/product/"
                    f"{CORA_DATASET_ID}_202511/Global/CO_X_20240101_PR_CT.nc",
                    "s3://bucket/product/"
                    f"{CORA_DATASET_ID}_202511/Global/CO_X_20240101_PR_PF.nc",
                ]
            else:
                names = [
                    "s3://bucket/product/"
                    f"{EASYCORA_DATASET_ID}_202511/Global/ECO_X_20240101_PR_PF.nc"
                ]
            with output.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=("filename", "size", "last_modified_datetime", "etag"),
                )
                writer.writeheader()
                for index, name in enumerate(names):
                    writer.writerow(
                        {
                            "filename": name,
                            "size": 100 + index,
                            "last_modified_datetime": "2026-01-01T00:00:00Z",
                            "etag": str(index),
                        }
                    )
            return
        self.downloads.append(kwargs)
        destination = Path(kwargs["output_directory"])
        destination.mkdir(parents=True, exist_ok=True)
        for line in Path(kwargs["file_list"]).read_text(encoding="utf-8").splitlines():
            if line:
                (destination / Path(line).name).write_bytes(b"x" * 100)


def test_download_plan_pins_inventory_version_and_exact_lists(tmp_path: Path) -> None:
    fake = _FakeCopernicusMarine()

    plan = create_download_plan(tmp_path, toolbox=fake)

    assert plan["argo_source"] == "easycora"
    assert [item["name"] for item in plan["selections"]] == ["ctd", "argo"]
    assert [item["file_count"] for item in plan["selections"]] == [1, 1]
    assert {item["dataset_version"] for item in plan["selections"]} == {"202511"}
    assert (tmp_path / "download_manifest.json").is_file()

    download_plan(tmp_path, plan, toolbox=fake)
    assert len(fake.downloads) == 2
    assert all(item["dataset_version"] == "202511" for item in fake.downloads)
    assert all(item["skip_existing"] is True for item in fake.downloads)


def _write_source_file(path: Path, *, probe_types: list[int], adjusted_temp: bool) -> None:
    n_prof = len(probe_types)
    pressure = np.tile(np.asarray([0.0, 20.0, 40.0]), (n_prof, 1))
    temperature = np.tile(np.asarray([1.0, 2.0, 3.0]), (n_prof, 1))
    salinity = np.tile(np.asarray([34.0, 34.2, 34.4]), (n_prof, 1))
    data_vars: dict[str, object] = {
        "LATITUDE": ("N_PROF", np.full(n_prof, -65.0)),
        "LONGITUDE": ("N_PROF", np.full(n_prof, -35.0)),
        "TIME": ("N_PROF", np.arange(n_prof, dtype=float) + 27000.0),
        "POSITION_QC": ("N_PROF", np.ones(n_prof, dtype=np.int8)),
        "TIME_QC": ("N_PROF", np.ones(n_prof, dtype=np.int8)),
        "PROBE_TYPE": ("N_PROF", np.asarray(probe_types, dtype=np.int8)),
        "PRES": (("N_PROF", "N_LEVELS"), pressure),
        "PRES_QC": (("N_PROF", "N_LEVELS"), np.ones_like(pressure, dtype=np.int8)),
        "TEMP": (("N_PROF", "N_LEVELS"), temperature),
        "TEMP_QC": (("N_PROF", "N_LEVELS"), np.ones_like(pressure, dtype=np.int8)),
        "PSAL": (("N_PROF", "N_LEVELS"), salinity),
        "PSAL_QC": (("N_PROF", "N_LEVELS"), np.ones_like(pressure, dtype=np.int8)),
        "PLATFORM_CODE": ("N_PROF", np.asarray([f"P{index}" for index in range(n_prof)])),
        "DC_REFERENCE": ("N_PROF", np.asarray([f"D{index}" for index in range(n_prof)])),
        "DATA_MODE": ("N_PROF", np.asarray(["D"] * n_prof)),
        "WMO_INST_TYPE": ("N_PROF", np.asarray(["830"] * n_prof)),
    }
    if adjusted_temp:
        data_vars["TEMP_ADJUSTED"] = (
            ("N_PROF", "N_LEVELS"),
            temperature + 10.0,
        )
        data_vars["TEMP_ADJUSTED_QC"] = (
            ("N_PROF", "N_LEVELS"),
            np.ones_like(pressure, dtype=np.int8),
        )
    dataset = xr.Dataset(data_vars)
    dataset["TIME"].attrs.update(units="days since 1950-01-01 00:00:00", calendar="standard")
    path.parent.mkdir(parents=True, exist_ok=True)
    dataset.to_netcdf(path)


def test_prepare_cora_tiles_keeps_only_ctd_and_argo_with_provenance(tmp_path: Path) -> None:
    source = tmp_path / "source_archive"
    _write_source_file(
        source / "source" / "cora_ctd_candidates" / "CO_X_20240101_PR_CT.nc",
        probe_types=[2, 12],
        adjusted_temp=True,
    )
    _write_source_file(
        source / "source" / "easycora_argo" / "ECO_X_20240102_PR_PF.nc",
        probe_types=[5, 2],
        adjusted_temp=False,
    )
    output = tmp_path / "tiles"

    manifest = prepare_cora_tiles(
        source,
        output,
        pressure_start=10.0,
        pressure_stop=30.0,
        pressure_step=10.0,
    )

    assert manifest["counters"]["profiles_written_ctd"] == 1
    assert manifest["counters"]["profiles_written_argo"] == 1
    tile = output / "CORA_lon40W_lat70S.nc"
    assert tile.is_file()
    with xr.open_dataset(tile, decode_times=False) as dataset:
        assert dataset.sizes["N_PROF"] == 2
        assert set(dataset["REFERENCE_KIND"].values.astype(str)) == {"ARGO", "CTD"}
        assert set(dataset["PROBE_TYPE"].values.tolist()) == {2, 5}
        ctd_index = int(np.flatnonzero(dataset["REFERENCE_KIND"].values.astype(str) == "CTD")[0])
        np.testing.assert_allclose(dataset["TEMP"].values[ctd_index], [11.5, 12.0, 12.5])
        assert dataset["JULD"].attrs["units"].startswith("days since 1950-01-01")
        assert str(dataset["SOURCE_FILE"].values[ctd_index]).endswith("PR_CT.nc")

    loaded = load_cora_tiles(
        output,
        lon_min=-40.0,
        lon_max=-30.0,
        lat_min=-70.0,
        lat_max=-60.0,
    )
    assert loaded["temp"].shape == (2, 3)
    assert loaded["psal"].shape == (2, 3)
