from __future__ import annotations

import concurrent.futures
from dataclasses import replace
import json

from meop_process.batch.runner import SummaryUpdateResult, run_all_deployments
from meop_process.catalog.tables import write_indexed_csv_rows
from meop_process.models import EmailNotificationSettings, EmailTransportSettings


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

    def fake_summaries(config, processed_deployments=None, force=False, output_dir=None):
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
            impacted_deployments=tuple(processed_deployments or []),
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


def test_run_all_deployments_runs_diagnostics_for_already_successful_deployments(meop_config, monkeypatch):
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
        ],
    )
    write_indexed_csv_rows(meop_config.catalogdir / "list_deployment_hr.csv", [])

    calls: list[str] = []
    diagnostics_calls: list[tuple[str, tuple[str, ...], str, bool]] = []

    def fake_process_tags(config, *, deployment: str, smru_name: str = "", notlc: bool = False):
        calls.append(deployment)
        out = meop_config.final_dataset_dir / deployment / f"{deployment}-AAA_hr2_prof.nc"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("ok", encoding="utf-8")
        return True

    def fake_generate_diagnostics(config, selection, qf, adjusted, parts=None, use_cached_summaries=True):
        normalized_parts = tuple(parts or ())
        diagnostics_calls.append((selection.deployment, normalized_parts, qf, adjusted))
        return ["diag"]

    def fake_summaries(config, processed_deployments=None, force=False, output_dir=None):
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
            impacted_deployments=tuple(processed_deployments or []),
            written=True,
        )

    monkeypatch.setattr("meop_process.batch.runner.process_tags_workflow", fake_process_tags)
    monkeypatch.setattr("meop_process.batch.runner.generate_diagnostics_plotting", fake_generate_diagnostics)
    monkeypatch.setattr("meop_process.batch.runner.update_metadata_summaries", fake_summaries)

    state_dir = meop_config.datadir / "batch_state"
    first = run_all_deployments(config=meop_config, state_dir=state_dir)
    assert first.success_count == 1
    assert diagnostics_calls == []
    assert calls == ["DEP001"]

    calls.clear()
    second = run_all_deployments(config=meop_config, state_dir=state_dir, diagnostics=True)
    assert second.success_count == 1
    assert second.skipped_count == 0
    assert calls == []
    assert second.deployment_results[0].message == "already completed successfully; diagnostics regenerated (1 files)"
    assert diagnostics_calls == [
        ("DEP001", ("tag", "deployment"), "lr1", True),
        ("", ("overview",), "lr1", True),
    ]


def test_run_all_deployments_prunes_stale_success_entries_without_outputs(meop_config, monkeypatch):
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

    def fake_process_tags(config, *, deployment: str, smru_name: str = "", notlc: bool = False):
        out = meop_config.final_dataset_dir / deployment / f"{deployment}-AAA_hr2_prof.nc"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("ok", encoding="utf-8")
        return True

    def fake_summaries(config, processed_deployments=None, force=False, output_dir=None):
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
            impacted_deployments=tuple(processed_deployments or []),
            written=True,
        )

    monkeypatch.setattr("meop_process.batch.runner.process_tags_workflow", fake_process_tags)
    monkeypatch.setattr("meop_process.batch.runner.update_metadata_summaries", fake_summaries)

    state_dir = meop_config.datadir / "batch_state"
    latest = state_dir / "latest" / "deployment_status.json"
    latest.parent.mkdir(parents=True, exist_ok=True)
    latest.write_text(
        json.dumps(
            {
                "DEP001": {"deployment": "DEP001", "status": "success"},
                "DEP002": {"deployment": "DEP002", "status": "success"},
            }
        ),
        encoding="utf-8",
    )

    result = run_all_deployments(config=meop_config, state_dir=state_dir, deployments=["DEP001"])

    assert result.success_count == 1
    payload = json.loads(latest.read_text(encoding="utf-8"))
    assert payload["DEP001"]["status"] == "success"
    assert "DEP002" not in payload


def test_run_all_deployments_sends_summary_email_when_enabled(meop_config, monkeypatch):
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
        ],
    )
    write_indexed_csv_rows(meop_config.catalogdir / "list_deployment_hr.csv", [])

    cfg = replace(
        meop_config,
        email_notifications=EmailNotificationSettings(
            enabled=True,
            when="always",
            to=("ops@example.org",),
            attach=("summary_md",),
            subject_prefix="[MEOP TEST]",
            transport=EmailTransportSettings(host="smtp.example.org", from_address="meop@example.org"),
        ),
    )

    sent: list[dict[str, object]] = []

    def fake_process_tags(config, *, deployment: str, smru_name: str = "", notlc: bool = False):
        out = cfg.final_dataset_dir / deployment / f"{deployment}-AAA_hr2_prof.nc"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("ok", encoding="utf-8")
        return True

    def fake_summaries(config, processed_deployments=None, force=False, output_dir=None):
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
            impacted_deployments=tuple(processed_deployments or []),
            written=True,
        )

    def fake_send(settings, *, subject, body, attachments=()):
        sent.append(
            {
                "to": settings.to,
                "subject": subject,
                "body": body,
                "attachments": tuple(str(path.name) for path in attachments),
            }
        )

    monkeypatch.setattr("meop_process.batch.runner.process_tags_workflow", fake_process_tags)
    monkeypatch.setattr("meop_process.batch.runner.update_metadata_summaries", fake_summaries)
    monkeypatch.setattr("meop_process.batch.runner.send_email_message", fake_send)

    result = run_all_deployments(config=cfg, state_dir=cfg.datadir / "batch_state")

    assert result.success_count == 1
    assert len(sent) == 1
    assert sent[0]["to"] == ("ops@example.org",)
    assert "summary.md" in sent[0]["attachments"]


def test_run_all_deployments_verbose_prints_deployment_logs(meop_config, monkeypatch, capsys):
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
        ],
    )
    write_indexed_csv_rows(meop_config.catalogdir / "list_deployment_hr.csv", [])

    def fake_process_tags(config, *, deployment: str, smru_name: str = "", notlc: bool = False):
        print(f"processing {deployment} verbose")
        out = meop_config.final_dataset_dir / deployment / f"{deployment}-AAA_hr2_prof.nc"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("ok", encoding="utf-8")
        return True

    def fake_summaries(config, processed_deployments=None, force=False, output_dir=None):
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
            impacted_deployments=tuple(processed_deployments or []),
            written=True,
        )

    monkeypatch.setattr("meop_process.batch.runner.process_tags_workflow", fake_process_tags)
    monkeypatch.setattr("meop_process.batch.runner.update_metadata_summaries", fake_summaries)

    run_all_deployments(config=meop_config, state_dir=meop_config.datadir / "batch_state", verbose=True)
    captured = capsys.readouterr()
    assert "processing DEP001 verbose" in captured.out


def test_run_all_deployments_jobs_uses_parallel_executor(meop_config, monkeypatch):
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
    observed_workers: list[int] = []

    def fake_process_tags(config, *, deployment: str, smru_name: str = "", notlc: bool = False):
        calls.append(deployment)
        out = meop_config.final_dataset_dir / deployment / f"{deployment}-AAA_hr2_prof.nc"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("ok", encoding="utf-8")
        return True

    def fake_summaries(config, processed_deployments=None, force=False, output_dir=None):
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
            impacted_deployments=tuple(processed_deployments or []),
            written=True,
        )

    class FakeExecutor:
        def __init__(self, *, max_workers):
            observed_workers.append(max_workers)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def submit(self, fn, **kwargs):
            future = concurrent.futures.Future()
            future.set_result(fn(**kwargs))
            return future

    monkeypatch.setattr("meop_process.batch.runner.process_tags_workflow", fake_process_tags)
    monkeypatch.setattr("meop_process.batch.runner.update_metadata_summaries", fake_summaries)
    monkeypatch.setattr("meop_process.batch.runner._PARALLEL_EXECUTOR_FACTORY", FakeExecutor)

    result = run_all_deployments(config=meop_config, state_dir=meop_config.datadir / "batch_state", jobs=2)
    assert result.success_count == 2
    assert set(calls) == {"DEP001", "DEP002"}
    assert observed_workers == [2]


def test_run_all_deployments_records_explicit_workflow_failure_reason(meop_config, monkeypatch):
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
        ],
    )
    write_indexed_csv_rows(meop_config.catalogdir / "list_deployment_hr.csv", [])

    class FakeWorkflowResult:
        def __init__(self, success: bool, reason: str) -> None:
            self.success = success
            self.reason = reason

        def __bool__(self) -> bool:
            return self.success

    def fake_process_tags(config, *, deployment: str, smru_name: str = "", notlc: bool = False):
        return FakeWorkflowResult(False, "no HR2 files written")

    def fake_summaries(config, processed_deployments=None, force=False, output_dir=None):
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
            impacted_deployments=tuple(processed_deployments or []),
            written=True,
        )

    monkeypatch.setattr("meop_process.batch.runner.process_tags_workflow", fake_process_tags)
    monkeypatch.setattr("meop_process.batch.runner.update_metadata_summaries", fake_summaries)

    result = run_all_deployments(config=meop_config, state_dir=meop_config.datadir / "batch_state")

    assert result.failed_count == 1
    assert result.deployment_results[0].message == "no HR2 files written"
