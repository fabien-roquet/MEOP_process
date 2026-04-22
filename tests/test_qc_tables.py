from __future__ import annotations

from importlib.resources import as_file, files

from meop_process.data.layout import resolve_table_path
from meop_process.processing.qc import ensure_processing_parameters, _set_coeff_flag


def _resource_text(name: str) -> str:
    resource = files("meop_process._bundled").joinpath("tables").joinpath(name)
    with as_file(resource) as path:
        return path.read_text(encoding="utf-8")


def test_ensure_processing_parameters_does_not_modify_packaged_table(meop_config) -> None:
    row = ensure_processing_parameters(meop_config, "DEP001")

    assert row["deployment_code"] == "DEP001"
    assert (meop_config.tablesdir / "table_param.csv").exists()
    assert (meop_config.tablesdir / "table_param.csv").read_text(encoding="utf-8") == _resource_text("table_param.csv")


def test_set_coeff_flag_does_not_modify_packaged_table(meop_config) -> None:
    resolve_table_path(meop_config, "table_coeff.csv")
    _set_coeff_flag(meop_config, "DEP001-AAA", "remove", "1")

    assert (meop_config.tablesdir / "table_coeff.csv").exists()
    assert (meop_config.tablesdir / "table_coeff.csv").read_text(encoding="utf-8") == _resource_text("table_coeff.csv")