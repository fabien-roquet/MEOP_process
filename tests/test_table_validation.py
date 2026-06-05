from __future__ import annotations

from meop_process.data.table_validation import validate_runtime_tables


def test_validate_table_coeff_detects_blank_column_duplicate_key_and_nonnumeric(meop_config) -> None:
    meop_config.tablesdir.mkdir(parents=True, exist_ok=True)
    (meop_config.tablesdir / "table_coeff.csv").write_text(
        ",smru_platform_code,T1,T2,S1,S2,remove,Sremove,comment\n"
        "idx1,DEP001-AAA,0,0,0,0,0,0,OK\n"
        "idx2,DEP001-AAA,bad,0,0,0,0,0,RM\n",
        encoding="utf-8",
    )

    result = validate_runtime_tables(meop_config, tables=("table_coeff.csv",))

    assert not result.ok
    messages = [issue.message for issue in result.errors]
    assert any("blank or unnamed column" in message for message in messages)
    assert any("duplicate key" in message for message in messages)
    assert any("expected numeric value" in message for message in messages)


def test_validate_table_coeff_rejects_unknown_standardized_comment(meop_config) -> None:
    meop_config.tablesdir.mkdir(parents=True, exist_ok=True)
    (meop_config.tablesdir / "table_coeff.csv").write_text(
        "smru_platform_code,T1,T2,S1,S2,remove,Sremove,comment\n"
        "DEP001-AAA,0,0,0,0,0,0,no comment\n",
        encoding="utf-8",
    )

    result = validate_runtime_tables(meop_config, tables=("table_coeff.csv",))

    assert not result.ok
    assert any("unknown standardized comment code" in issue.message for issue in result.errors)


def test_validate_table_filter_allows_legacy_blank_columns_as_warnings(meop_config) -> None:
    meop_config.tablesdir.mkdir(parents=True, exist_ok=True)
    (meop_config.tablesdir / "table_filter.csv").write_text(
        "smru_platform_name,Sonly,filter,x1,x2,,\n"
        "DEP001-AAA,0,date_min,0,NaN,,\n"
        "DEP001-AAA,0,date_max,1,NaN,,\n",
        encoding="utf-8",
    )

    result = validate_runtime_tables(meop_config, tables=("table_filter.csv",))

    assert result.ok
    assert len(result.warnings) == 2
    assert all("blank or unnamed column" in issue.message for issue in result.warnings)


def test_validate_runtime_tables_reports_missing_required_table(meop_config) -> None:
    result = validate_runtime_tables(meop_config, tables=("table_coeff.csv",))

    assert not result.ok
    assert result.errors[0].message.startswith("missing runtime table")
