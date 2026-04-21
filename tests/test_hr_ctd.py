from __future__ import annotations

from meop_process.io.hr_ctd import resolve_hr_ctd_path


def test_resolve_hr_ctd_path_treats_nan_prefix_as_empty(meop_config, seed_catalog) -> None:
    seed_catalog(deployment="ft18", smru_name="ft18-F2-12")
    year_dir = meop_config.raw_hr_dir / "2020"
    year_dir.mkdir(parents=True, exist_ok=True)
    expected = year_dir / "42_ctd.txt"
    expected.write_text("ok", encoding="utf-8")

    meop_config.catalogdir.joinpath("list_deployment_hr.csv").write_text(
        "smru_platform_code,instr_id,year,prefix,continuous\n"
        "ft18-F2-12,42,2020,NaN,0\n",
        encoding="utf-8",
    )

    resolved = resolve_hr_ctd_path(meop_config, "ft18-F2-12")

    assert resolved is not None
    assert resolved.prefix == ""
    assert resolved.expected_path == expected
    assert resolved.exists is True