from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import xarray as xr

from meop_process.catalog.filenames import fname_prof
from meop_process.models import Selection
from meop_process.processing.locations import apply_location_adjustment_placeholder


JULD_REF = datetime(1950, 1, 1, tzinfo=timezone.utc)


def _juld(year: int, month: int, day: int) -> float:
    moment = datetime(year, month, day, tzinfo=timezone.utc)
    return (moment - JULD_REF).total_seconds() / 86400.0


def test_apply_location_adjustment_uses_cls_locations(meop_config, seed_catalog) -> None:
    seed_catalog(deployment="DEP001", smru_name="DEP001-AAA")
    path = fname_prof("DEP001-AAA", qf="lr0", config=meop_config)
    path.parent.mkdir(parents=True, exist_ok=True)
    xr.Dataset(
        {
            "JULD": (("N_PROF",), np.asarray([_juld(2014, 1, 2), _juld(2014, 1, 3), _juld(2014, 1, 4)], dtype=np.float64)),
            "LATITUDE": (("N_PROF",), np.asarray([-60.0, -60.1, -60.2], dtype=np.float64)),
            "LONGITUDE": (("N_PROF",), np.asarray([10.0, 10.1, 10.2], dtype=np.float64)),
        },
        coords={"N_PROF": np.arange(3)},
        attrs={"ptt": "12345", "loc_algorithm": "CLS KALMAN", "smru_platform_code": "DEP001-AAA"},
    ).to_netcdf(path)

    cls_path = meop_config.cls_locdir / "12345_2014-01-01_2014-12-31.smoothing.csv"
    cls_path.parent.mkdir(parents=True, exist_ok=True)
    cls_path.write_text(
        "Platform ID No.;Latitude;Longitude;Latitude 2;Longitude 2;Loc. quality;Loc. date=yyyy/MM/dd HH:mm:ss;Altitude;Pass;Sat.;Frequency;Msg;Error radius;Semi-major axis;Semi-minor axis;Ellipse orientation;GDOP (m/Hz);Dist. Subsat Track (deg);\n"
        "12345;-61.0;11.0;-61.0;11.0;A;2014/01/02 00:00:00;0;0;0;0;0;0;0;0;0;0;0;\n"
        "12345;-61.1;11.1;-61.1;11.1;A;2014/01/03 00:00:00;0;0;0;0;0;0;0;0;0;0;0;\n"
        "12345;-61.2;11.2;-61.2;11.2;A;2014/01/04 00:00:00;0;0;0;0;0;0;0;0;0;0;0;\n",
        encoding="utf-8",
    )

    result = apply_location_adjustment_placeholder(meop_config, Selection(deployment="DEP001", smru_name="DEP001-AAA"))

    assert result.placeholder is False
    assert path in result.written_files
    with xr.open_dataset(path, decode_times=False) as updated:
        np.testing.assert_allclose(updated["LATITUDE"].values, np.asarray([-61.0, -61.1, -61.2]))
        np.testing.assert_allclose(updated["LONGITUDE"].values, np.asarray([11.0, 11.1, 11.2]))
        assert updated.attrs["loc_algorithm"] == "CLS SMOOTH KALMAN"


def test_apply_location_adjustment_noop_without_external_sources(meop_config, seed_catalog) -> None:
    seed_catalog(deployment="DEP001", smru_name="DEP001-AAA")
    path = fname_prof("DEP001-AAA", qf="lr0", config=meop_config)
    path.parent.mkdir(parents=True, exist_ok=True)
    xr.Dataset(
        {
            "JULD": (("N_PROF",), np.asarray([_juld(2014, 1, 2), _juld(2014, 1, 3), _juld(2014, 1, 4)], dtype=np.float64)),
            "LATITUDE": (("N_PROF",), np.asarray([-60.0, -60.1, -60.2], dtype=np.float64)),
            "LONGITUDE": (("N_PROF",), np.asarray([10.0, 10.1, 10.2], dtype=np.float64)),
        },
        coords={"N_PROF": np.arange(3)},
        attrs={"ptt": "12345", "loc_algorithm": "CLS KALMAN", "smru_platform_code": "DEP001-AAA"},
    ).to_netcdf(path)

    result = apply_location_adjustment_placeholder(meop_config, Selection(deployment="DEP001", smru_name="DEP001-AAA"))

    assert result.written_files == ()


def test_apply_location_adjustment_ignores_cls_file_missing_required_columns(meop_config, seed_catalog, capsys) -> None:
    seed_catalog(deployment="DEP001", smru_name="DEP001-AAA")
    path = fname_prof("DEP001-AAA", qf="lr0", config=meop_config)
    path.parent.mkdir(parents=True, exist_ok=True)
    xr.Dataset(
        {
            "JULD": (("N_PROF",), np.asarray([_juld(2014, 1, 2), _juld(2014, 1, 3), _juld(2014, 1, 4)], dtype=np.float64)),
            "LATITUDE": (("N_PROF",), np.asarray([-60.0, -60.1, -60.2], dtype=np.float64)),
            "LONGITUDE": (("N_PROF",), np.asarray([10.0, 10.1, 10.2], dtype=np.float64)),
        },
        coords={"N_PROF": np.arange(3)},
        attrs={"ptt": "12345", "loc_algorithm": "CLS KALMAN", "smru_platform_code": "DEP001-AAA"},
    ).to_netcdf(path)

    cls_path = meop_config.cls_locdir / "12345_2014-01-01_2014-12-31.smoothing.csv"
    cls_path.parent.mkdir(parents=True, exist_ok=True)
    cls_path.write_text(
        "Platform ID No.;Latitude;Longitude;Latitude 2;Longitude 2;Loc. quality\n"
        "12345;-61.0;11.0;-61.0;11.0;A\n",
        encoding="utf-8",
    )

    result = apply_location_adjustment_placeholder(meop_config, Selection(deployment="DEP001", smru_name="DEP001-AAA"))

    assert result.written_files == ()
    assert "missing columns" in capsys.readouterr().out


def test_apply_location_adjustment_updates_positioning_system_with_char_dim(meop_config, seed_catalog) -> None:
    seed_catalog(deployment="DEP001", smru_name="DEP001-AAA")
    path = fname_prof("DEP001-AAA", qf="lr0", config=meop_config)
    path.parent.mkdir(parents=True, exist_ok=True)
    xr.Dataset(
        {
            "JULD": (("N_PROF",), np.asarray([_juld(2014, 1, 2), _juld(2014, 1, 3), _juld(2014, 1, 4)], dtype=np.float64)),
            "LATITUDE": (("N_PROF",), np.asarray([-60.0, -60.1, -60.2], dtype=np.float64)),
            "LONGITUDE": (("N_PROF",), np.asarray([10.0, 10.1, 10.2], dtype=np.float64)),
            "POSITIONING_SYSTEM": (("N_PROF", "STRING8"), np.full((3, 8), b" ", dtype="S1")),
        },
        coords={"N_PROF": np.arange(3), "STRING8": np.arange(8)},
        attrs={"ptt": "12345", "loc_algorithm": "CLS LEAST SQUARES", "smru_platform_code": "DEP001-AAA"},
    ).to_netcdf(path)

    cls_path = meop_config.cls_locdir / "12345_2014-01-01_2014-12-31.smoothing.csv"
    cls_path.parent.mkdir(parents=True, exist_ok=True)
    cls_path.write_text(
        "Platform ID No.;Latitude;Longitude;Latitude 2;Longitude 2;Loc. quality;Loc. date=yyyy/MM/dd HH:mm:ss;Altitude;Pass;Sat.;Frequency;Msg;Error radius;Semi-major axis;Semi-minor axis;Ellipse orientation;GDOP (m/Hz);Dist. Subsat Track (deg);\n"
        "12345;-61.0;11.0;-61.0;11.0;A;2014/01/02 00:00:00;0;0;0;0;0;0;0;0;0;0;0;\n"
        "12345;-61.1;11.1;-61.1;11.1;A;2014/01/03 00:00:00;0;0;0;0;0;0;0;0;0;0;0;\n"
        "12345;-61.2;11.2;-61.2;11.2;A;2014/01/04 00:00:00;0;0;0;0;0;0;0;0;0;0;0;\n",
        encoding="utf-8",
    )

    result = apply_location_adjustment_placeholder(meop_config, Selection(deployment="DEP001", smru_name="DEP001-AAA"))

    assert path in result.written_files
    with xr.open_dataset(path, decode_times=False) as updated:
        assert updated["POSITIONING_SYSTEM"].dims == ("N_PROF", "STRING8")
        values = updated["POSITIONING_SYSTEM"].values.astype("S1")
        decoded = [b"".join(row).decode("ascii").strip() for row in values]
        assert decoded == ["LS", "LS", "LS"]
