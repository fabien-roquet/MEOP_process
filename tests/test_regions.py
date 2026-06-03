"""Tests for meop_process.plotting.regions — geographic region labelling."""
from __future__ import annotations

import numpy as np
import pytest

from meop_process.plotting import regions
from meop_process.plotting.regions import label_region, label_regions


def test_southern_ocean():
    assert label_region(0.0, -60.0) == "Southern Ocean"


def test_north_atlantic():
    assert label_region(-30.0, 45.0) == "North Atlantic"


def test_south_pacific():
    assert label_region(160.0, -45.0) == "South Pacific"


def test_indian_ocean():
    assert label_region(70.0, -35.0) == "Indian Ocean"


def test_north_pacific():
    assert label_region(-150.0, 40.0) == "North Pacific"


def test_south_atlantic():
    assert label_region(-20.0, -35.0) == "South Atlantic"


def test_label_regions_array():
    lons = np.array([0.0, -30.0, 160.0])
    lats = np.array([-60.0, 45.0, -45.0])
    result = label_regions(lons, lats)
    assert result[0] == "Southern Ocean"
    assert result[1] == "North Atlantic"
    assert result[2] == "South Pacific"


def test_label_regions_2d():
    lons = np.array([[0.0, -30.0], [160.0, -150.0]])
    lats = np.array([[-60.0, 45.0], [-45.0, 40.0]])
    result = label_regions(lons, lats)
    assert result.shape == (2, 2)
    assert result[0, 0] == "Southern Ocean"
    assert result[0, 1] == "North Atlantic"


def test_label_region_returns_string():
    result = label_region(0.0, -60.0)
    assert isinstance(result, str)
    assert len(result) > 0


def test_label_regions_returns_ndarray():
    result = label_regions(np.array([0.0]), np.array([-60.0]))
    assert isinstance(result, np.ndarray)


def test_regionmask_is_opt_in(monkeypatch):
    class ExplodingRegionmask:
        @property
        def defined_regions(self):  # pragma: no cover - should not be touched
            raise AssertionError("regionmask should not be used by default")

    monkeypatch.delenv("MEOP_USE_REGIONMASK_AR6", raising=False)
    monkeypatch.setattr(regions, "_AVAILABLE", True)
    monkeypatch.setattr(regions, "_interpolator", None)
    monkeypatch.setattr(regions, "_region_names", None)
    monkeypatch.setattr(regions, "regionmask", ExplodingRegionmask(), raising=False)

    assert label_region(-30.0, 45.0) == "North Atlantic"


def test_regionmask_failure_falls_back(monkeypatch):
    class ExplodingRegionmask:
        @property
        def defined_regions(self):
            raise RuntimeError("offline cache missing")

    monkeypatch.setenv("MEOP_USE_REGIONMASK_AR6", "1")
    monkeypatch.setattr(regions, "_AVAILABLE", True)
    monkeypatch.setattr(regions, "_interpolator", None)
    monkeypatch.setattr(regions, "_region_names", None)
    monkeypatch.setattr(regions, "regionmask", ExplodingRegionmask(), raising=False)

    assert label_region(160.0, -45.0) == "South Pacific"
