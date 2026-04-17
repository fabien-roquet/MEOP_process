from __future__ import annotations

from meop_process.processing.qc import ensure_processing_parameters, _set_coeff_flag


def test_ensure_processing_parameters_does_not_create_root_table_copy(meop_config) -> None:
    ensure_processing_parameters(meop_config, "DEP001")

    assert (meop_config.tablesdir / "table_param.csv").exists()
    assert not (meop_config.processdir / "table_param.csv").exists()


def test_set_coeff_flag_does_not_create_root_table_copy(meop_config) -> None:
    _set_coeff_flag(meop_config, "DEP001-AAA", "remove", "1")

    assert (meop_config.tablesdir / "table_coeff.csv").exists()
    assert not (meop_config.processdir / "table_coeff.csv").exists()