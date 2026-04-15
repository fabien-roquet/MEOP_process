from __future__ import annotations

from meop_process.api import bootstrap_data_store, describe_runtime_data_layout, validate_runtime_data_layout
from meop_process.data.layout import bootstrap_packaged_catalogs, resolve_catalog_path, resolve_table_path


def test_bootstrap_data_store_seeds_packaged_tables(meop_config) -> None:
    changed = bootstrap_data_store(meop_config)

    assert (meop_config.tablesdir / "table_meta.csv").exists()
    assert changed


def test_resolve_table_path_returns_runtime_table(meop_config) -> None:
    bootstrap_data_store(meop_config)
    resolved = resolve_table_path(meop_config, "table_meta.csv")
    assert resolved == meop_config.tablesdir / "table_meta.csv"


def test_resolve_catalog_path_uses_runtime_catalog_root(meop_config) -> None:
    path = meop_config.catalogdir / "list_deployment.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        ",deployment_code,pi_code,process,public,country\nDEP001,DEP001,PI001,1,1,SE\n",
        encoding="utf-8",
    )

    resolved = resolve_catalog_path(meop_config, "list_deployment.csv")

    assert resolved == path
    assert resolved.read_text(encoding="utf-8") == path.read_text(encoding="utf-8")


def test_describe_runtime_data_layout_reports_expected_sections(meop_config) -> None:
    summary = describe_runtime_data_layout(meop_config)

    assert set(summary) == {"roots", "patterns", "packaged_tables", "packaged_catalogs", "catalog_tables"}
    assert summary["roots"]["tables_root"] == str(meop_config.tablesdir)
    assert summary["roots"]["catalog_root"] == str(meop_config.catalogdir)
    assert summary["roots"]["raw_odv_root"] == str(meop_config.raw_odv_dir)


def test_validate_runtime_data_layout_marks_missing_runtime_catalog_files(meop_config) -> None:
    records = validate_runtime_data_layout(meop_config)
    by_name = {record["name"]: record for record in records}

    assert by_name["table_meta.csv"]["exists"] is True
    assert by_name["list_deployment.csv"]["required"] is True
    assert by_name["list_deployment.csv"]["exists"] is False


def test_bootstrap_packaged_catalogs_seeds_reference_catalogs(meop_config) -> None:
    changed = bootstrap_packaged_catalogs(meop_config)

    assert (meop_config.catalogdir / "list_deployment_hr.csv").exists()
    assert (meop_config.catalogdir / "list_deployment.csv").exists()
    assert changed
