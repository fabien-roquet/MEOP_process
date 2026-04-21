from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from meop_process.compare_cli import main
from meop_process.processing.adjustments import _ensure_default_coefficients


def _write_dataset(path: Path, *, values: list[float], title: str = "same") -> None:
    dataset = xr.Dataset(
        data_vars={"TEMP": (("N_PROF",), np.asarray(values, dtype=np.float32))},
        attrs={"title": title},
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    dataset.to_netcdf(path)


def test_compare_cli_reports_equal_files(tmp_path) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    _write_dataset(reference, values=[1.0, 2.0])
    _write_dataset(candidate, values=[1.0, 2.0])

    assert main([str(reference), str(candidate)]) == 0


def test_compare_cli_reports_directory_difference(tmp_path, capsys) -> None:
    reference_dir = tmp_path / "ref"
    candidate_dir = tmp_path / "cand"
    _write_dataset(reference_dir / "a.nc", values=[1.0, 2.0])
    _write_dataset(candidate_dir / "a.nc", values=[1.0, 3.0])

    assert main([str(reference_dir), str(candidate_dir)]) == 1
    assert "variable: variable TEMP: values differ" in capsys.readouterr().out


def test_compare_cli_variable_filter_skips_attributes_by_default(tmp_path, capsys) -> None:
    reference = tmp_path / "reference.nc"
    candidate = tmp_path / "candidate.nc"
    _write_dataset(reference, values=[1.0, 2.0], title="left")
    _write_dataset(candidate, values=[1.0, 3.0], title="right")

    assert main([str(reference), str(candidate), "--variable", "TEMP"]) == 1
    output = capsys.readouterr().out
    assert "variable: variable TEMP: values differ" in output
    assert "attribute:" not in output


def test_adjustments_default_coefficients_rewrite_plain_csv(meop_config) -> None:
    path = meop_config.tablesdir / "table_coeff.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        ",smru_platform_code,T1,T2,S1,S2,remove,Sremove,comment\n"
        "OLD,OLD,0,0,0,0,0,0,no comment\n",
        encoding="utf-8",
    )

    _ensure_default_coefficients(meop_config, ["NEW-TAG"])

    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "smru_platform_code,T1,T2,S1,S2,remove,Sremove,comment"
    assert all(not line.startswith(",") for line in lines[1:])
