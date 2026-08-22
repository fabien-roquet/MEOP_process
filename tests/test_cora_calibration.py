"""Tests for the CORA reference tile loader and CORA-based calibration plots."""

from __future__ import annotations

import csv
from dataclasses import replace
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


def test_load_cora_tiles_discovers_unpadded_archive_names(tmp_path: Path) -> None:
    """The real local archive uses lon40W/lon00E rather than lon040W/lon000E."""
    from meop_process.reference.cora import discover_cora_tiles, load_cora_tiles

    tile_path = tmp_path / "CORA_lon40W_lat70S.nc"
    _write_mock_cora_tile(
        tile_path,
        n_prof=12,
        lat_range=(-70.0, -61.0),
        lon_range=(-40.0, -31.0),
    )

    inventory = discover_cora_tiles(tmp_path)
    assert inventory[(-40, -70)] == tile_path

    result = load_cora_tiles(
        tmp_path,
        lon_min=-45.0,
        lon_max=-25.0,
        lat_min=-75.0,
        lat_max=-55.0,
    )

    assert result["lat"].shape[0] == 12
    assert result["temp"].shape == (12, 10)


def test_discover_cora_tiles_deduplicates_padding_variants(tmp_path: Path, caplog) -> None:
    from meop_process.reference.cora import discover_cora_tiles

    short = tmp_path / "CORA_lon40W_lat70S.nc"
    padded = tmp_path / "CORA_lon040W_lat70S.nc"
    short.touch()
    padded.touch()

    inventory = discover_cora_tiles(tmp_path)

    assert inventory == {(-40, -70): short}
    assert "Multiple CORA tiles resolve to cell" in caplog.text


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
    lat: float | np.ndarray = -65.0,
    lon: float | np.ndarray = -35.0,
    juld: np.ndarray | None = None,
) -> None:
    """Write a minimal MEOP-style lr1 prof NetCDF file."""
    import xarray as xr

    rng = np.random.default_rng(0)
    pres_vals = np.tile(np.arange(1, n_levels + 1, dtype=np.float64) * 10.0, (n_prof, 1))
    temp_vals = rng.uniform(0.0, 5.0, (n_prof, n_levels))
    psal_vals = rng.uniform(33.8, 34.5, (n_prof, n_levels))
    latitude = np.broadcast_to(np.asarray(lat, dtype=np.float64), (n_prof,)).copy()
    longitude = np.broadcast_to(np.asarray(lon, dtype=np.float64), (n_prof,)).copy()
    time = (
        np.arange(n_prof, dtype=np.float64)
        if juld is None
        else np.asarray(juld, dtype=np.float64)
    )
    ds = xr.Dataset(
        {
            "PRES_ADJUSTED": (("N_PROF", "N_LEVELS"), pres_vals),
            "TEMP_ADJUSTED": (("N_PROF", "N_LEVELS"), temp_vals),
            "PSAL_ADJUSTED": (("N_PROF", "N_LEVELS"), psal_vals),
            "LATITUDE": ("N_PROF", latitude),
            "LONGITUDE": ("N_PROF", longitude),
            "JULD": ("N_PROF", time),
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(path)


def _write_offset_meop_prof(
    path: Path,
    *,
    psal_offset: float,
    psal_slope: float,
    temp_offset: float,
    n_prof: int = 12,
) -> np.ndarray:
    """Write a small profile file with controlled tag-minus-CORA anomalies."""
    import xarray as xr

    pres_grid = np.asarray([100.0, 300.0, 400.0, 500.0, 600.0, 800.0], dtype=np.float64)
    pres_vals = np.tile(pres_grid, (n_prof, 1))
    cora_temp = np.asarray([1.4, 1.2, 1.0, 0.8, 0.6, 0.4], dtype=np.float64)
    cora_psal = np.asarray([34.0, 34.1, 34.2, 34.3, 34.4, 34.5], dtype=np.float64)
    temp_vals = np.tile(cora_temp + temp_offset, (n_prof, 1))
    psal_vals = np.tile(cora_psal + psal_offset + psal_slope * pres_grid, (n_prof, 1))
    ds = xr.Dataset(
        {
            "PRES_ADJUSTED": (("N_PROF", "N_LEVELS"), pres_vals),
            "TEMP_ADJUSTED": (("N_PROF", "N_LEVELS"), temp_vals),
            "PSAL_ADJUSTED": (("N_PROF", "N_LEVELS"), psal_vals),
            "LATITUDE": ("N_PROF", np.full(n_prof, -65.0)),
            "LONGITUDE": ("N_PROF", np.full(n_prof, -35.0)),
            "JULD": ("N_PROF", np.arange(n_prof, dtype=np.float64)),
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(path)
    return pres_grid


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


def test_plot_ts_calibration_writes_offset_diagnostics(tmp_path: Path) -> None:
    from meop_process.compare_cli import _diagnostic_support
    from meop_process.plotting.calibration import plot_ts_calibration

    target_path = tmp_path / "offset_hr1.nc"
    pres_grid = _write_offset_meop_prof(target_path, psal_offset=0.01, psal_slope=2e-5, temp_offset=-0.03)
    cora_temp = np.asarray([1.4, 1.2, 1.0, 0.8, 0.6, 0.4], dtype=np.float64)
    cora_psal = np.asarray([34.0, 34.1, 34.2, 34.3, 34.4, 34.5], dtype=np.float64)
    cora_data = {
        "lat": np.asarray([-65.0, -66.0]),
        "lon": np.asarray([-35.0, -36.0]),
        "juld": np.asarray([0.0, 1.0]),
        "temp": np.tile(cora_temp, (2, 1)),
        "psal": np.tile(cora_psal, (2, 1)),
        "pres": pres_grid,
    }

    written = plot_ts_calibration(
        "test-tag-offset",
        cora_data=cora_data,
        target_path=target_path,
        output_dir=tmp_path / "plots",
        write_diagnostics=True,
    )

    names = {path.name for path in written}
    assert "test-tag-offset_calibration.png" in names
    assert "test-tag-offset_calibration_offsets.png" in names
    csv_path = tmp_path / "plots" / "test-tag-offset_calibration_offsets.csv"
    assert csv_path.name in names

    rows = list(csv.DictReader(csv_path.open(encoding="utf-8", newline="")))
    psal_band = next(row for row in rows if row["variable"] == "PSAL" and row["diagnostic"] == "band_400_600")
    psal_linear = next(row for row in rows if row["variable"] == "PSAL" and row["diagnostic"] == "linear_200_1000")
    temp_band = next(row for row in rows if row["variable"] == "TEMP" and row["diagnostic"] == "band_400_600")

    assert float(psal_band["suggested_S2"]) == pytest.approx(0.02, abs=1e-6)
    assert float(psal_linear["suggested_S1"]) == pytest.approx(0.02, abs=1e-6)
    assert float(psal_linear["suggested_S2"]) == pytest.approx(0.01, abs=1e-6)
    assert float(temp_band["suggested_T2"]) == pytest.approx(-0.03, abs=1e-6)
    assert _diagnostic_support(written) == {"TEMP": True, "PSAL": True}


def test_diagnostic_support_rejects_sparse_profiles(tmp_path: Path) -> None:
    from meop_process.compare_cli import _diagnostic_support
    from meop_process.plotting.calibration import plot_ts_calibration

    target_path = tmp_path / "offset_hr1.nc"
    pres_grid = _write_offset_meop_prof(
        target_path,
        psal_offset=0.01,
        psal_slope=2e-5,
        temp_offset=-0.03,
        n_prof=4,
    )
    cora_data = {
        "lat": np.asarray([-65.0, -66.0]),
        "lon": np.asarray([-35.0, -36.0]),
        "juld": np.asarray([0.0, 1.0]),
        "temp": np.tile(np.asarray([1.4, 1.2, 1.0, 0.8, 0.6, 0.4]), (2, 1)),
        "psal": np.tile(np.asarray([34.0, 34.1, 34.2, 34.3, 34.4, 34.5]), (2, 1)),
        "pres": pres_grid,
    }

    written = plot_ts_calibration(
        "test-tag-sparse",
        cora_data=cora_data,
        target_path=target_path,
        output_dir=tmp_path / "plots",
        write_diagnostics=True,
    )

    assert _diagnostic_support(written) == {"TEMP": False, "PSAL": False}


def test_diagnostic_support_rejects_incomplete_vertical_coverage(tmp_path: Path) -> None:
    from meop_process.compare_cli import _diagnostic_support

    csv_path = tmp_path / "test-tag_calibration_offsets.csv"
    fieldnames = (
        "variable",
        "diagnostic",
        "pressure_min",
        "pressure_max",
        "n_profiles",
        "n_levels",
        "suggested_T1",
        "suggested_T2",
        "suggested_S1",
        "suggested_S2",
    )
    rows: list[dict[str, str]] = []
    for variable in ("TEMP", "PSAL"):
        prefix = "T" if variable == "TEMP" else "S"
        coefficient_fields = {
            f"suggested_{prefix}1": "0.01",
            f"suggested_{prefix}2": "0.02",
        }
        rows.extend(
            (
                {
                    "variable": variable,
                    "diagnostic": "linear_200_1000",
                    "n_profiles": "12",
                    "n_levels": "3",
                    **coefficient_fields,
                },
                {
                    "variable": variable,
                    "diagnostic": "band_400_600",
                    "n_profiles": "12",
                    "n_levels": "3",
                    **coefficient_fields,
                },
                {
                    "variable": variable,
                    "diagnostic": "all_valid_profile_median",
                    "pressure_min": "300",
                    "pressure_max": "600",
                },
            )
        )
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    assert _diagnostic_support([csv_path]) == {"TEMP": False, "PSAL": False}


def test_cora_cells_for_track_is_antimeridian_safe_and_avoids_envelope() -> None:
    from meop_process.reference.cora import cora_cells_for_track

    antimeridian = cora_cells_for_track(
        np.asarray([-65.0, -64.0]),
        np.asarray([179.0, -179.0]),
        margin=5.0,
    )
    assert not any(-170 < lon < 170 for lon, _ in antimeridian)

    disconnected = cora_cells_for_track(
        np.asarray([-44.5, 55.0]),
        np.asarray([146.4, -1.3]),
        margin=5.0,
    )
    assert len(disconnected) < 20
    assert (70, 0) not in disconnected


def test_load_cora_track_loads_only_profiles_near_accepted_positions(tmp_path: Path) -> None:
    from meop_process.reference.cora import load_cora_track

    tile_path = tmp_path / "CORA_lon040W_lat70S.nc"
    _write_mock_cora_tile(
        tile_path,
        n_prof=20,
        lat_range=(-70.0, -61.0),
        lon_range=(-40.0, -31.0),
    )
    import xarray as xr

    with xr.open_dataset(tile_path, decode_times=False) as dataset:
        target_lat = float(dataset["LATITUDE"].values[0])
        target_lon = float(dataset["LONGITUDE"].values[0])

    result = load_cora_track(
        tmp_path,
        latitudes=np.asarray([target_lat]),
        longitudes=np.asarray([target_lon]),
        margin=1e-4,
    )

    assert result["lat"].shape == (1,)
    assert result["lat"][0] == pytest.approx(target_lat)
    assert result["lon"][0] == pytest.approx(target_lon)
    assert result["temp"].shape == (1, 10)


def test_target_preflight_removes_null_island_and_minor_test_component(tmp_path: Path) -> None:
    from meop_process.compare_cli import (
        _preflight_calibration_target,
        _target_support_is_sufficient,
    )

    target_path = tmp_path / "mixed_hr1_prof.nc"
    local_lat = np.linspace(-50.0, -54.0, 12)
    local_lon = np.linspace(140.0, 146.0, 12)
    _write_mock_meop_prof(
        target_path,
        n_prof=14,
        n_levels=80,
        lat=np.concatenate((local_lat, [55.0, 0.0])),
        lon=np.concatenate((local_lon, [-2.0, 0.0])),
        juld=np.arange(14, dtype=np.float64),
    )

    selected = _preflight_calibration_target(target_path)

    assert selected.indices.tolist() == list(range(12))
    assert selected.component_sizes == (12, 1)
    assert selected.dropped_invalid == 1
    assert selected.dropped_minor_segments == 1
    assert all(_target_support_is_sufficient(item) for item in selected.support.values())


def test_generate_calibration_plots_rejects_sparse_target_before_cora(
    meop_config, monkeypatch, tmp_path: Path
) -> None:
    from meop_process.compare_cli import (
        InsufficientTargetDataError,
        generate_calibration_plots,
    )

    target_path = tmp_path / "DEP001-SPARSE_hr1_prof.nc"
    _write_mock_meop_prof(target_path, n_prof=4, n_levels=80)
    cfg = replace(meop_config, cora_dir=tmp_path / "cora")
    monkeypatch.setattr(
        "meop_process.compare_cli._select_calibration_target",
        lambda smru_name, config: (target_path, "hr1"),
    )

    def unexpected_cora_load(*args, **kwargs):
        raise AssertionError("CORA must not be loaded for an impossible target")

    monkeypatch.setattr("meop_process.compare_cli.load_cora_track", unexpected_cora_load)

    with pytest.raises(InsufficientTargetDataError, match="before CORA loading"):
        generate_calibration_plots("DEP001-SPARSE", config=cfg)


# ---------------------------------------------------------------------------
# CLI --plot1 tests
# ---------------------------------------------------------------------------


def test_generate_calibration_plots_reports_no_reference(meop_config, monkeypatch, tmp_path: Path) -> None:
    from meop_process.compare_cli import NoReferenceDataError, generate_calibration_plots

    target_path = tmp_path / "DEP001-AAA_hr1_prof.nc"
    _write_mock_meop_prof(target_path, n_prof=12, n_levels=80)
    cfg = replace(meop_config, cora_dir=tmp_path / "cora")
    monkeypatch.setattr(
        "meop_process.compare_cli._select_calibration_target",
        lambda smru_name, config: (target_path, "hr1"),
    )
    monkeypatch.setattr(
        "meop_process.compare_cli.load_cora_track",
        lambda *args, **kwargs: {
            "lat": np.empty(0),
            "lon": np.empty(0),
            "juld": np.empty(0),
            "temp": np.empty((0, 0)),
            "psal": np.empty((0, 0)),
            "pres": np.empty(0),
        },
    )

    with pytest.raises(NoReferenceDataError, match="no CORA profiles found"):
        generate_calibration_plots("DEP001-AAA", config=cfg)


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
