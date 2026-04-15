from __future__ import annotations

from meop_process.batch.runner import SummaryUpdateResult, run_all_deployments
from meop_process.catalog.tables import write_indexed_csv_rows


def test_run_all_deployments_resumes_success_and_continues_after_failure(meop_config, monkeypatch):
    write_indexed_csv_rows(
        meop_config.catalogdir / "list_deployment.csv",
        [
            {
                "row_name": "DEP001",
                "deployment_code": "DEP001",
                "pi_code": "PI1",
                "process": "1",
                "public": "1",
                "country": "SE",
                "task_done": "",
                "first_version": "",
                "last_version": "",
                "start_date": "2020-01-01",
                "end_date": "2020-12-31",
                "start_date_jul": "",
            },
            {
                "row_name": "DEP002",
                "deployment_code": "DEP002",
                "pi_code": "PI2",
                "process": "1",
                "public": "1",
                "country": "SE",
                "task_done": "",
                "first_version": "",
                "last_version": "",
                "start_date": "2020-01-01",
                "end_date": "2020-12-31",
                "start_date_jul": "",
            },
        ],
    )
    write_indexed_csv_rows(meop_config.catalogdir / "list_deployment_hr.csv", [])

    calls: list[str] = []

    def fake_process_tags(config, *, deployment: str, smru_name: str = "", notlc: bool = False):
        calls.append(deployment)
        if deployment == "DEP002":
            raise RuntimeError("boom")
        out = meop_config.final_dataset_dir / deployment / f"{deployment}-AAA_hr2_prof.nc"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("ok", encoding="utf-8")
        return True

    def fake_summaries(config, impacted_deployments=None, force=False, output_dir=None):
        root = config.publicdir_ctd
        root.mkdir(parents=True, exist_ok=True)
        tags = root / "list_tags.csv"
        deps = root / "list_deployments.csv"
        tags.write_text("SMRU_PLATFORM_CODE,DEPLOYMENT_CODE\n", encoding="utf-8")
        deps.write_text("DEPLOYMENT_CODE\n", encoding="utf-8")
        return SummaryUpdateResult(
            output_dir=root,
            list_tags_path=tags,
            list_deployments_path=deps,
            impacted_deployments=tuple(impacted_deployments or []),
            written=True,
        )

    monkeypatch.setattr("meop_process.batch.runner.process_tags_workflow", fake_process_tags)
    monkeypatch.setattr("meop_process.batch.runner.update_metadata_summaries", fake_summaries)

    state_dir = meop_config.datadir / "batch_state"
    first = run_all_deployments(config=meop_config, state_dir=state_dir)
    assert first.success_count == 1
    assert first.failed_count == 1
    assert first.summary_markdown.exists()
    assert first.summary_csv.exists()
    assert first.state_path.exists()
    assert set(calls) == {"DEP001", "DEP002"}

    calls.clear()
    second = run_all_deployments(config=meop_config, state_dir=state_dir)
    assert second.skipped_count == 1
    assert second.failed_count == 1
    assert calls == ["DEP002"]
