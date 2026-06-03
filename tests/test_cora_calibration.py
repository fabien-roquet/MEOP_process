"""Tests for the CORA reference tile loader and CORA-based calibration plots."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Helpers to create minimal mock CORA tiles
# ---------------------------------------------------------------------------


def _write_mock_cora_tile(
    path: Path,
    n_prof: int = 10,
    pres_grid: list[float] | None = None,
    lat_range: tuple[float, float] = (-70.0, -60.0),
    lon_range: tuple[float, float] = (-40.0, -30.0),
) -> None:
    """Write a minimal CORA-style NetCDF tile to *path*."""
    import xarray as xr

    if pres_grid is None:
        pres_grid = [10.0 * i for i in range(1, 11)]  # 10, 20, ..., 100 dbar
    n_pres = len(pres_grid)

    rng = np.random.default_rng(42)
    lats = rng.uniform(lat_range[0], lat_range[1], n_prof)
    lons = rng.uniform(lon_range[0], lon_range[1], n_prof)
    temp = rng.uniform(0.0, 5.0, (n_prof, n_pres))
    psal = rng.uniform(33.5, 34.5, (n_prof, n_pres))

    ds = xr.Dataset(
        {
            "LATITUDE": ("N_PROF", lats),
            "LONGITUDE": ("N_PROF", lons),
            "JULD": ("N_PROF", np.arange(n_prof, dtype=np.float64) * 1.0),
            "TEMP": (("N_PROF", "PRES_GRID"), temp),
            "PSAL": (("N_PROF", "PRES_GRID"), psal),
        },
        coords={"PRES_GRID": pres_grid},
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(path)


# ---------------------------------------------------------------------------
# CORA tile loader tests
# ---------------------------------------------------------------------------


def test_load_cora_tiles_basic(tmp_path: Path) -> None:
    """load_cora_tiles returns correct keys and shapes from a single tile."""
    from meop_process.reference.cora import load_cora_tiles

    tile_path = tmp_path / "CORA_lon040W_lat70S.nc"
    _write_mock_cora_tile(
        tile_path,
        n_prof=15,
        lat_range=(-70.0, -61.0),
        lon_range=(-40.0, -31.0),
    )

    result = load_cora_tiles(
        tmp_path,
        lon_min=-45.0,
        lon_max=-25.0,
        lat_min=-75.0,
        lat_max=-55.0,
    )
    assert set(result.keys()) == {"lat", "lon", "juld", "temp", "psal", "pres"}
    assert result["lat"].ndim == 1
    assert result["temp"].ndim == 2
    assert result["psal"].shape == result["temp"].shape
    assert result["temp"].shape[0] == result["lat"].shape[0]
    assert result["pres"].ndim == 1
    assert result["temp"].shape[1] == result["pres"].shape[0]


def test_load_cora_tiles_filters_by_bbox(tmp_path: Path) -> None:
    """Profiles outside the bounding box are excluded."""
    from meop_process.reference.cora import load_cora_tiles

    tile_path = tmp_path / "CORA_lon040W_lat70S.nc"
    _write_mock_cora_tile(
        tile_path,
        n_prof=20,
        lat_range=(-75.0, -65.0),
        lon_range=(-45.0, -35.0),
    )

    # Tight bbox that excludes many profiles
    result = load_cora_tiles(
        tmp_path,
        lon_min=-42.0,
        lon_max=-38.0,
        lat_min=-72.0,
        lat_max=-68.0,
    )
    # All returned profiles must be within the bbox
    assert np.all(result["lat"] >= -72.0)
    assert np.all(result["lat"] <= -68.0)
    assert np.all(result["lon"] >= -42.0)
    assert np.all(result["lon"] <= -38.0)


def test_load_cora_tiles_missing_tiles_returns_empty(tmp_path: Path) -> None:
    """Returns zero-length arrays when no tiles exist for the bbox."""
    from meop_process.reference.cora import load_cora_tiles

    result = load_cora_tiles(
        tmp_path,
        lon_min=-45.0,
        lon_max=-35.0,
        lat_min=-75.0,
        lat_max=-65.0,
    )
    assert result["lat"].shape[0] == 0
    assert result["temp"].shape[0] == 0
    assert result["psal"].shape[0] == 0


def test_load_cora_tiles_multiple_tiles(tmp_path: Path) -> None:
    """Profiles from multiple tiles are concatenated."""
    from meop_process.reference.cora import load_cora_tiles

    # Two adjacent tiles: lon 040W lat 70S and lon 040W lat 60S
    _write_mock_cora_tile(
        tmp_path / "CORA_lon040W_lat70S.nc",
        n_prof=8,
        lat_range=(-70.0, -61.0),
        lon_range=(-40.0, -31.0),
    )
    _write_mock_cora_tile(
        tmp_path / "CORA_lon040W_lat60S.nc",
        n_prof=6,
        lat_range=(-60.0, -51.0),
        lon_range=(-40.0, -31.0),
    )

    result = load_cora_tiles(
        tmp_path,
        lon_min=-42.0,
        lon_max=-29.0,
        lat_min=-75.0,
        lat_max=-50.0,
    )
    assert result["lat"].shape[0] > 0
    assert result["temp"].shape[0] == result["lat"].shape[0]


def test_load_cora_tiles_skips_malformed_existing_tile(tmp_path: Path) -> None:
    import xarray as xr
    from meop_process.reference.cora import load_cora_tiles

    xr.Dataset().to_netcdf(tmp_path / "CORA_lon050W_lat80S.nc")
    _write_mock_cora_tile(
        tmp_path / "CORA_lon040W_lat70S.nc",
        n_prof=8,
        lat_range=(-70.0, -61.0),
        lon_range=(-40.0, -31.0),
    )

    result = load_cora_tiles(
        tmp_path,
        lon_min=-45.0,
        lon_max=-25.0,
        lat_min=-75.0,
        lat_max=-55.0,
    )

    assert result["lat"].shape[0] > 0
    assert result["temp"].shape[1] == result["pres"].shape[0]


# ---------------------------------------------------------------------------
# Tile filename helper tests
# ---------------------------------------------------------------------------


def test_tile_filename_east_north() -> None:
    from meop_process.reference.cora import _tile_filename

    assert _tile_filename(0, 0) == "CORA_lon000E_lat00N.nc"
    assert _tile_filename(10, 60) == "CORA_lon010E_lat60N.nc"
    assert _tile_filename(170, 0) == "CORA_lon170E_lat00N.nc"


def test_tile_filename_west_south() -> None:
    from meop_process.reference.cora import _tile_filename

    assert _tile_filename(-10, -70) == "CORA_lon010W_lat70S.nc"
    assert _tile_filename(-40, -60) == "CORA_lon040W_lat60S.nc"
    assert _tile_filename(-180, -80) == "CORA_lon180W_lat80S.nc"


# ---------------------------------------------------------------------------
# Calibration plot tests
# ---------------------------------------------------------------------------


def _write_mock_meop_prof(
    path: Path,
    n_prof: int = 5,
    n_levels: int = 10,
    lat: float = -65.0,
    lon: float = -35.0,
) -> None:
    """Write a minimal MEOP-style lr1 prof NetCDF file."""
    import xarray as xr

    rng = np.random.default_rng(0)
    pres_vals = np.tile(np.arange(1, n_levels + 1, dtype=np.float64) * 10.0, (n_prof, 1))
    temp_vals = rng.uniform(0.0, 5.0, (n_prof, n_levels))
    psal_vals = rng.uniform(33.8, 34.5, (n_prof, n_levels))
    ds = xr.Dataset(
        {
            "PRES_ADJUSTED": (("N_PROF", "N_LEVELS"), pres_vals),
            "TEMP_ADJUSTED": (("N_PROF", "N_LEVELS"), temp_vals),
            "PSAL_ADJUSTED": (("N_PROF", "N_LEVELS"), psal_vals),
            "LATITUDE": ("N_PROF", np.full(n_prof, lat)),
            "LONGITUDE": ("N_PROF", np.full(n_prof, lon)),
            "JULD": ("N_PROF", np.arange(n_prof, dtype=np.float64)),
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(path)


def test_plot_ts_calibration_produces_file(tmp_path: Path) -> None:
    """plot_ts_calibration writes at least one PNG file."""
    from meop_process.plotting.calibration import plot_ts_calibration

    target_path = tmp_path / "mock_lr1.nc"
    _write_mock_meop_prof(target_path, n_prof=5)

    # Empty CORA (no tiles) — still should produce a plot
    cora_data = {
        "lat": np.empty(0),
        "lon": np.empty(0),
        "juld": np.empty(0),
        "temp": np.empty((0, 10)),
        "psal": np.empty((0, 10)),
        "pres": np.arange(1, 11, dtype=np.float64) * 10.0,
    }

    written = plot_ts_calibration(
        "test-tag-00",
        cora_data=cora_data,
        target_path=target_path,
        output_dir=tmp_path / "plots",
    )
    assert len(written) >= 1
    for p in written:
        assert p.exists()
        assert p.suffix == ".png"


def test_plot_ts_calibration_with_cora_data(tmp_path: Path) -> None:
    """plot_ts_calibration runs successfully with real CORA mock data."""
    from meop_process.plotting.calibration import plot_ts_calibration
    from meop_process.reference.cora import load_cora_tiles

    # Write mock CORA tile
    _write_mock_cora_tile(
        tmp_path / "cora" / "CORA_lon040W_lat70S.nc",
        n_prof=20,
        lat_range=(-70.0, -61.0),
        lon_range=(-40.0, -31.0),
    )
    cora_data = load_cora_tiles(
        tmp_path / "cora",
        lon_min=-45.0,
        lon_max=-25.0,
        lat_min=-75.0,
        lat_max=-55.0,
    )
    assert cora_data["lat"].shape[0] > 0

    target_path = tmp_path / "mock_lr1.nc"
    _write_mock_meop_prof(target_path, n_prof=10, lat=-65.0, lon=-35.0)

    written = plot_ts_calibration(
        "mock-tag-01",
        cora_data=cora_data,
        target_path=target_path,
        output_dir=tmp_path / "plots",
    )
    assert len(written) == 1
    assert written[0].name == "mock-tag-01_calibration.png"


def test_plot_ts_calibration_chunking(tmp_path: Path) -> None:
    """With 250 profiles and chunk_size=200, two files should be written."""
    from meop_process.plotting.calibration import plot_ts_calibration

    target_path = tmp_path / "mock_lr1.nc"
    _write_mock_meop_prof(target_path, n_prof=250)

    cora_data = {
        "lat": np.empty(0),
        "lon": np.empty(0),
        "juld": np.empty(0),
        "temp": np.empty((0, 10)),
        "psal": np.empty((0, 10)),
        "pres": np.arange(1, 11, dtype=np.float64) * 10.0,
    }

    written = plot_ts_calibration(
        "test-tag-chunk",
        cora_data=cora_data,
        target_path=target_path,
        output_dir=tmp_path / "plots",
        chunk_size=200,
    )
    assert len(written) == 2
    assert written[0].name == "test-tag-chunk_calibration_part0.png"
    assert written[1].name == "test-tag-chunk_calibration_part1.png"


def test_plot_ts_calibration_interpolates_cora_and_context_grids(tmp_path: Path) -> None:
    from meop_process.plotting.calibration import plot_ts_calibration

    target_path = tmp_path / "mock_hr1.nc"
    other_path = tmp_path / "other_hr1.nc"
    _write_mock_meop_prof(target_path, n_prof=4, n_levels=16)
    _write_mock_meop_prof(other_path, n_prof=3, n_levels=1000)

    cora_pres = np.linspace(0.0, 990.0, 100)
    cora_data = {
        "lat": np.asarray([-65.0, -66.0]),
        "lon": np.asarray([-35.0, -36.0]),
        "juld": np.asarray([0.0, 1.0]),
        "temp": np.tile(np.linspace(0.0, 4.0, cora_pres.size), (2, 1)),
        "psal": np.tile(np.linspace(33.8, 34.4, cora_pres.size), (2, 1)),
        "pres": cora_pres,
    }

    written = plot_ts_calibration(
        "test-tag-grid",
        cora_data=cora_data,
        target_path=target_path,
        other_paths=[other_path],
        output_dir=tmp_path / "plots",
    )

    assert len(written) == 1
    assert written[0].exists()


# ---------------------------------------------------------------------------
# CLI --plot1 tests
# ---------------------------------------------------------------------------


def test_compare_cli_plot1_no_cora_dir_error(tmp_path: Path, meop_config) -> None:
    """--plot1 exits with code 2 when cora_dir is not configured."""
    import json
    from meop_process.compare_cli import main

    # meop_config has no cora_dir, build a config.json that points at it
    config_json = tmp_path / "configs.json"
    config_json.write_text(
        json.dumps(
            {
                "defaults": {
                    "datadir": str(meop_config.datadir),
                    "public": str(meop_config.publicdir),
                }
            }
        )
    )
    ret = main(["--plot1", "some-tag", "--config", str(config_json)])
    assert ret == 2


def test_compare_cli_plot1_prefers_hr1_product(meop_config) -> None:
    from meop_process.compare_cli import _select_calibration_target

    smru_name = "DEP001-AAA"
    deployment_dir = meop_config.final_dataset_dir / "DEP001"
    deployment_dir.mkdir(parents=True, exist_ok=True)
    for qf in ("lr1", "hr1", "hr2"):
        (deployment_dir / f"{smru_name}_{qf}_prof.nc").write_text("placeholder", encoding="utf-8")

    target, qf = _select_calibration_target(smru_name, config=meop_config)

    assert qf == "hr1"
    assert target is not None
    assert target.name == f"{smru_name}_hr1_prof.nc"


def test_compare_cli_requires_reference_and_candidate() -> None:
    """Calling meop-compare without --plot1 and without positional args fails."""
    from meop_process.compare_cli import main

    with pytest.raises(SystemExit) as exc_info:
        main([])
    assert exc_info.value.code != 0
