from __future__ import annotations

from meop_process.api import bootstrap_data_store, describe_runtime_data_layout, validate_runtime_data_layout
from meop_process.data.layout import bootstrap_packaged_catalogs, resolve_catalog_path, resolve_table_path
from meop_process.models import DeploymentInfo, Selection
from meop_process.processing.cleanup import remove_deployment_outputs


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
    path.write_text("deployment_code,pi_code,process,public,country\nBROKEN,PI001,1,1,SE\n", encoding="utf-8")

    resolved = resolve_catalog_path(meop_config, "list_deployment.csv")

    assert resolved == path
    assert resolved.read_text(encoding="utf-8") == path.read_text(encoding="utf-8")


def test_resolve_catalog_path_returns_indexed_format_unchanged(meop_config) -> None:
    path = meop_config.catalogdir / "list_deployment.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    indexed_content = ",deployment_code,pi_code,process,public,country\nDEP001,DEP001,PI001,1,1,SE\n"
    path.write_text(indexed_content, encoding="utf-8")

    resolved = resolve_catalog_path(meop_config, "list_deployment.csv")

    assert resolved == path
    assert resolved.read_text(encoding="utf-8") == indexed_content


def test_describe_runtime_data_layout_reports_expected_sections(meop_config) -> None:
    summary = describe_runtime_data_layout(meop_config)

    assert set(summary) == {"roots", "patterns", "packaged_tables", "packaged_catalogs", "catalog_tables"}
    assert summary["roots"]["tables_root"] == str(meop_config.tablesdir)
    assert summary["roots"]["catalog_root"] == str(meop_config.catalogdir)
    assert summary["roots"]["data_raw_root"] == str(meop_config.data_raw_dir)
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


def test_remove_deployment_outputs_cleans_plots_by_deployments(meop_config) -> None:
    dep = "ct1"
    raw_marker = meop_config.data_raw_dir / f"{dep}_ODV.txt"
    raw_marker.parent.mkdir(parents=True, exist_ok=True)
    raw_marker.write_text("marker", encoding="utf-8")

    dep_dir = meop_config.final_dataset_dir / dep
    dep_dir.mkdir(parents=True, exist_ok=True)
    plots_by_tags_dir = meop_config.plotdir / dep
    plots_by_tags_dir.mkdir(parents=True, exist_ok=True)
    plots_by_dep_dir = meop_config.plots_by_deployment_dir / dep
    plots_by_dep_dir.mkdir(parents=True, exist_ok=True)

    stale_nc = dep_dir / "ct1-001_lr0_prof.nc"
    stale_tag_plot = plots_by_tags_dir / "ct1-001_lr1_TS.png"
    stale_dep_plot = plots_by_dep_dir / "ct1_lr1_deployment_overview_adj.png"
    for f in (stale_nc, stale_tag_plot, stale_dep_plot):
        f.write_text("old", encoding="utf-8")

    info = DeploymentInfo(
        selection=Selection(deployment=dep),
        record=None,
        invalid_code=False,
        directory=dep_dir,
        raw_input_dir=meop_config.raw_odv_dir,
        raw_input_zip=meop_config.raw_odv_dir / f"{dep}_ODV.zip",
        raw_working_text=raw_marker,
    )

    removed = remove_deployment_outputs(meop_config, info)

    assert stale_nc not in removed or not stale_nc.exists()
    assert not stale_tag_plot.exists(), "stale tag plot should be removed"
    assert not stale_dep_plot.exists(), "stale deployment plot should be removed"
