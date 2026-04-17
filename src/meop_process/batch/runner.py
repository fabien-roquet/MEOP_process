from __future__ import annotations

import argparse
import concurrent.futures
import contextlib
import io
import json
import os
import sys
import time
import traceback
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from ..catalog.deployments import load_deployment_catalog
from ..config.loader import load_config
from ..metadata.summaries import SummaryUpdateResult, update_metadata_summaries
from ..models import MeopConfig, Selection
from ..plotting.diagnostics import generate_diagnostics as generate_diagnostics_plotting
from ..workflows.process import process_tags as process_tags_workflow


STATE_FILE_NAME = "deployment_status.json"
SUCCESS_STATUSES = {"success", "skipped"}
_PARALLEL_EXECUTOR_FACTORY = concurrent.futures.ProcessPoolExecutor


@dataclass(frozen=True)
class DeploymentRunResult:
    deployment: str
    status: str
    started_at: str
    finished_at: str
    duration_seconds: float
    log_path: Path
    message: str = ""

    def as_dict(self) -> dict[str, object]:
        return {
            "deployment": self.deployment,
            "status": self.status,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "duration_seconds": self.duration_seconds,
            "log_path": str(self.log_path),
            "message": self.message,
        }


@dataclass(frozen=True)
class BatchRunResult:
    run_id: str
    output_dir: Path
    log_dir: Path
    state_path: Path
    summary_markdown: Path
    summary_csv: Path
    deployment_results: tuple[DeploymentRunResult, ...]
    metadata_summary: SummaryUpdateResult

    @property
    def success_count(self) -> int:
        return sum(item.status == "success" for item in self.deployment_results)

    @property
    def skipped_count(self) -> int:
        return sum(item.status == "skipped" for item in self.deployment_results)

    @property
    def failed_count(self) -> int:
        return sum(item.status == "failed" for item in self.deployment_results)

    def as_dict(self) -> dict[str, object]:
        return {
            "run_id": self.run_id,
            "output_dir": str(self.output_dir),
            "log_dir": str(self.log_dir),
            "state_path": str(self.state_path),
            "summary_markdown": str(self.summary_markdown),
            "summary_csv": str(self.summary_csv),
            "success_count": self.success_count,
            "skipped_count": self.skipped_count,
            "failed_count": self.failed_count,
            "metadata_summary": self.metadata_summary.as_dict(),
            "deployment_results": [item.as_dict() for item in self.deployment_results],
        }


@dataclass(frozen=True)
class _PendingDeployment:
    index: int
    deployment: str
    log_path: Path
    diagnostics_only: bool = False


class _Tee(io.TextIOBase):
    def __init__(self, *targets: io.TextIOBase) -> None:
        self.targets = targets

    def write(self, s: str) -> int:  # pragma: no cover - trivial
        for target in self.targets:
            try:
                target.write(s)
                target.flush()
            except ValueError:
                continue
        return len(s)

    def flush(self) -> None:  # pragma: no cover - trivial
        for target in self.targets:
            try:
                target.flush()
            except ValueError:
                continue


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _timestamp(dt: datetime | None = None) -> str:
    value = dt or _utc_now()
    return value.strftime("%Y%m%dT%H%M%SZ")


def _state_root(config: MeopConfig, state_dir: str | Path | None = None) -> Path:
    root = Path(state_dir) if state_dir is not None else config.datadir / "batch"
    root.mkdir(parents=True, exist_ok=True)
    (root / "runs").mkdir(parents=True, exist_ok=True)
    (root / "latest").mkdir(parents=True, exist_ok=True)
    return root


def _load_state(path: Path) -> dict[str, dict[str, object]]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        return {}
    return {str(key): value for key, value in payload.items() if isinstance(value, dict)}


def _write_state(path: Path, payload: dict[str, dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def _truthy_process(value: object) -> bool:
    text = str(value).strip().lower()
    return text in {"", "1", "y", "yes", "true"}


def _eligible_deployments(config: MeopConfig, *, include_disabled: bool = False, selected: Iterable[str] | None = None) -> list[str]:
    catalog = load_deployment_catalog(config, persist=False)
    requested = {item.strip() for item in (selected or []) if str(item).strip()}
    deployments: list[str] = []
    for deployment, record in sorted(catalog.items()):
        if requested and deployment not in requested:
            continue
        if not include_disabled and not _truthy_process(getattr(record, "process", "")):
            continue
        deployments.append(deployment)
    if requested:
        missing = sorted(requested.difference(deployments))
        for deployment in missing:
            if deployment in catalog:
                deployments.append(deployment)
    return sorted(set(deployments))


def _has_any_outputs(config: MeopConfig, deployment: str) -> bool:
    deployment_dir = config.final_dataset_dir / deployment
    return deployment_dir.exists() and any(deployment_dir.glob(f"{deployment}-*_prof.nc"))


def _reconcile_state(config: MeopConfig, state: dict[str, dict[str, object]]) -> dict[str, dict[str, object]]:
    reconciled: dict[str, dict[str, object]] = {}
    for deployment, entry in state.items():
        status = str(entry.get("status", "")).strip().lower()
        if status in SUCCESS_STATUSES and not _has_any_outputs(config, deployment):
            continue
        reconciled[deployment] = entry
    return reconciled


def _should_skip(
    deployment: str,
    *,
    state: dict[str, dict[str, object]],
    force: bool,
    force_failed: bool,
    config: MeopConfig,
) -> tuple[bool, str]:
    if force:
        return False, "forced"
    entry = state.get(deployment)
    if not entry:
        return False, "not yet processed"
    if force_failed and str(entry.get("status", "")).strip().lower() == "failed":
        return False, "rerunning failed deployment"
    if str(entry.get("status", "")).strip().lower() != "success":
        return False, "previous run not successful"
    if not _has_any_outputs(config, deployment):
        return False, "success marker exists but outputs are missing"
    return True, "already completed successfully"


def _write_csv(path: Path, results: Iterable[DeploymentRunResult]) -> None:
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["deployment", "status", "started_at", "finished_at", "duration_seconds", "log_path", "message"])
        for result in results:
            writer.writerow([
                result.deployment,
                result.status,
                result.started_at,
                result.finished_at,
                f"{result.duration_seconds:.3f}",
                str(result.log_path),
                result.message,
            ])


def _write_markdown(path: Path, run_id: str, results: Iterable[DeploymentRunResult], metadata_summary: SummaryUpdateResult) -> None:
    items = list(results)
    success = [item for item in items if item.status == "success"]
    skipped = [item for item in items if item.status == "skipped"]
    failed = [item for item in items if item.status == "failed"]
    lines = [
        f"# MEOP batch processing report ({run_id})",
        "",
        f"Success: **{len(success)}**  ",
        f"Skipped: **{len(skipped)}**  ",
        f"Failed: **{len(failed)}**",
        "",
        "## Result table",
        "",
        "| Deployment | Status | Duration (s) | Message | Log |",
        "|---|---:|---:|---|---|",
    ]
    for item in items:
        lines.append(
            f"| {item.deployment} | {item.status} | {item.duration_seconds:.1f} | {item.message or ''} | {item.log_path.name} |"
        )
    lines.extend(
        [
            "",
            "## Metadata summaries",
            "",
            f"- list_tags.csv: `{metadata_summary.list_tags_path}`",
            f"- list_deployments.csv: `{metadata_summary.list_deployments_path}`",
            f"- impacted deployments: {', '.join(metadata_summary.impacted_deployments) or 'none'}",
            f"- files rewritten: {'yes' if metadata_summary.written else 'no'}",
            "",
        ]
    )
    if failed:
        lines.extend(["## Failed deployments", ""])
        for item in failed:
            lines.append(f"- **{item.deployment}**: {item.message}")
        lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def _diagnostics_count(result: object) -> int:
    written = getattr(result, "written_files", None)
    if written is not None:
        return len(written)
    return len(result)


def _workflow_failure_message(result: object) -> str:
    reason = getattr(result, "reason", "")
    if isinstance(reason, str) and reason.strip():
        return reason.strip()
    return "workflow returned False"


def _print_log_to_stdout(deployment: str, log_path: Path) -> None:
    text = log_path.read_text(encoding="utf-8") if log_path.exists() else ""
    if not text:
        return
    print(f"[{deployment}] log start")
    print(text, end="" if text.endswith("\n") else "\n")
    print(f"[{deployment}] log end")


def _execute_deployment(
    *,
    cfg: MeopConfig,
    deployment: str,
    log_path: Path,
    notlc: bool,
    diagnostics: bool,
    diagnostics_qf: str,
    diagnostics_raw: bool,
    diagnostics_only: bool,
    stream_stdout: bool,
) -> tuple[DeploymentRunResult, bool]:
    started = _utc_now()
    targets: list[io.TextIOBase] = []
    if stream_stdout:
        targets.append(sys.stdout)
    with log_path.open("w", encoding="utf-8") as log_handle:
        targets.insert(0, log_handle)
        tee = _Tee(*targets)
        try:
            with contextlib.redirect_stdout(tee), contextlib.redirect_stderr(tee):
                if diagnostics_only:
                    print(f"[{deployment}] skipped processing, running diagnostics instead")
                    generated = generate_diagnostics_plotting(
                        cfg,
                        Selection(deployment=deployment, smru_name="").normalized(),
                        qf=diagnostics_qf,
                        adjusted=not diagnostics_raw,
                    )
                    print(f"[{deployment}] diagnostics generated: {_diagnostics_count(generated)}")
                    status = "success"
                    message = "skipped processing; diagnostics generated"
                    processed_for_summary = False
                else:
                    print(f"[{deployment}] notlc={notlc} diagnostics={diagnostics}")
                    started_timer = time.perf_counter()
                    workflow_result = process_tags_workflow(cfg, deployment=deployment, smru_name="", notlc=notlc)
                    ok = bool(workflow_result)
                    if ok and diagnostics:
                        generated = generate_diagnostics_plotting(
                            cfg,
                            Selection(deployment=deployment, smru_name="").normalized(),
                            qf=diagnostics_qf,
                            adjusted=not diagnostics_raw,
                        )
                        print(f"[{deployment}] diagnostics generated: {_diagnostics_count(generated)}")
                    finished_timer = time.perf_counter()
                    if ok:
                        status = "success"
                        message = ""
                        processed_for_summary = True
                    else:
                        status = "failed"
                        message = _workflow_failure_message(workflow_result)
                        processed_for_summary = False
        except Exception as exc:  # pragma: no cover - exercised by dedicated tests
            finished_timer = time.perf_counter()
            status = "failed"
            message = f"diagnostics failed: {exc}" if diagnostics_only else str(exc)
            processed_for_summary = False
            with contextlib.redirect_stdout(tee), contextlib.redirect_stderr(tee):
                traceback.print_exc()
        finished = _utc_now()
        duration = (
            finished_timer - started_timer
            if not diagnostics_only and "started_timer" in locals()
            else (finished - started).total_seconds()
        )
        result = DeploymentRunResult(
            deployment=deployment,
            status=status,
            started_at=started.isoformat(),
            finished_at=finished.isoformat(),
            duration_seconds=duration,
            log_path=log_path,
            message=message,
        )
    return result, processed_for_summary


def _run_pending_deployments(
    *,
    cfg: MeopConfig,
    pending: list[_PendingDeployment],
    notlc: bool,
    diagnostics: bool,
    diagnostics_qf: str,
    diagnostics_raw: bool,
    jobs: int,
    verbose: bool,
):
    if jobs <= 1 or len(pending) <= 1:
        for task in pending:
            result, processed = _execute_deployment(
                cfg=cfg,
                deployment=task.deployment,
                log_path=task.log_path,
                notlc=notlc,
                diagnostics=diagnostics,
                diagnostics_qf=diagnostics_qf,
                diagnostics_raw=diagnostics_raw,
                diagnostics_only=task.diagnostics_only,
                stream_stdout=verbose,
            )
            yield task.index, result, processed
        return

    max_workers = min(jobs, len(pending))
    with _PARALLEL_EXECUTOR_FACTORY(max_workers=max_workers) as executor:
        futures: dict[concurrent.futures.Future, _PendingDeployment] = {}
        for task in pending:
            future = executor.submit(
                _execute_deployment,
                cfg=cfg,
                deployment=task.deployment,
                log_path=task.log_path,
                notlc=notlc,
                diagnostics=diagnostics,
                diagnostics_qf=diagnostics_qf,
                diagnostics_raw=diagnostics_raw,
                diagnostics_only=task.diagnostics_only,
                stream_stdout=False,
            )
            futures[future] = task
        for future in concurrent.futures.as_completed(futures):
            task = futures[future]
            try:
                result, processed = future.result()
            except Exception as exc:  # pragma: no cover - defensive
                finished = _utc_now()
                task.log_path.write_text(f"[{task.deployment}] worker failure: {exc}\n", encoding="utf-8")
                result = DeploymentRunResult(
                    deployment=task.deployment,
                    status="failed",
                    started_at=finished.isoformat(),
                    finished_at=finished.isoformat(),
                    duration_seconds=0.0,
                    log_path=task.log_path,
                    message=str(exc),
                )
                processed = False
            if verbose:
                _print_log_to_stdout(task.deployment, task.log_path)
            yield task.index, result, processed


def run_all_deployments(
    *,
    config: MeopConfig | None = None,
    processdir: str | Path | None = None,
    config_file: str | Path | None = None,
    machine: str | None = None,
    notlc: bool = False,
    diagnostics: bool = False,
    diagnostics_qf: str = "lr1",
    diagnostics_raw: bool = False,
    force: bool = False,
    force_failed: bool = False,
    include_disabled: bool = False,
    deployments: Iterable[str] | None = None,
    state_dir: str | Path | None = None,
    jobs: int = 1,
    verbose: bool = False,
) -> BatchRunResult:
    if jobs < 1:
        raise ValueError("jobs must be at least 1")
    cfg = config or load_config(processdir=processdir, config_file=config_file, machine=machine)
    root = _state_root(cfg, state_dir)
    run_id = _timestamp()
    output_dir = root / "runs" / run_id
    log_dir = output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    state_path = root / "latest" / STATE_FILE_NAME
    state = _reconcile_state(cfg, _load_state(state_path))
    _write_state(state_path, state)

    selected_deployments = _eligible_deployments(cfg, include_disabled=include_disabled, selected=deployments)
    results_by_index: list[tuple[int, DeploymentRunResult]] = []
    processed_for_summary: list[str] = []
    pending: list[_PendingDeployment] = []

    for index, deployment in enumerate(selected_deployments):
        skip, reason = _should_skip(
            deployment,
            state=state,
            force=force,
            force_failed=force_failed,
            config=cfg,
        )
        started = _utc_now()
        log_path = log_dir / f"{deployment}.log"

        if skip:
            if diagnostics:
                pending.append(_PendingDeployment(index=index, deployment=deployment, log_path=log_path, diagnostics_only=True))
                continue
            else:
                finished = _utc_now()
                result = DeploymentRunResult(
                    deployment=deployment,
                    status="skipped",
                    started_at=started.isoformat(),
                    finished_at=finished.isoformat(),
                    duration_seconds=(finished - started).total_seconds(),
                    log_path=log_path,
                    message=reason,
                )
                log_path.write_text(f"[{deployment}] skipped: {reason}\n", encoding="utf-8")
            if verbose:
                _print_log_to_stdout(deployment, log_path)
            results_by_index.append((index, result))
            state[deployment] = {**result.as_dict(), "run_id": run_id}
            continue

        pending.append(_PendingDeployment(index=index, deployment=deployment, log_path=log_path))

    for index, result, processed in _run_pending_deployments(
        cfg=cfg,
        pending=pending,
        notlc=notlc,
        diagnostics=diagnostics,
        diagnostics_qf=diagnostics_qf,
        diagnostics_raw=diagnostics_raw,
        jobs=jobs,
        verbose=verbose,
    ):
        results_by_index.append((index, result))
        state[result.deployment] = {**result.as_dict(), "run_id": run_id}
        if processed and result.status == "success":
            processed_for_summary.append(result.deployment)
        _write_state(state_path, state)

    results = [result for _, result in sorted(results_by_index, key=lambda item: item[0])]

    metadata_summary = update_metadata_summaries(cfg, processed_deployments=processed_for_summary, force=force)
    summary_csv = output_dir / "summary.csv"
    summary_markdown = output_dir / "summary.md"
    _write_csv(summary_csv, results)
    _write_markdown(summary_markdown, run_id, results, metadata_summary)
    _write_state(state_path, state)

    return BatchRunResult(
        run_id=run_id,
        output_dir=output_dir,
        log_dir=log_dir,
        state_path=state_path,
        summary_markdown=summary_markdown,
        summary_csv=summary_csv,
        deployment_results=tuple(results),
        metadata_summary=metadata_summary,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Process all MEOP deployments with resumable state and readable reports.")
    parser.add_argument("--processdir", default=None, help="Override the MEOP process directory.")
    parser.add_argument("--config-file", default=None, help="Explicit path to a runtime config JSON file.")
    parser.add_argument("--machine", default=None, help="Machine entry key from the runtime config JSON.")
    parser.add_argument("--notlc", action="store_true", help="Use the no-TLC branch.")
    parser.add_argument("--diagnostics", action="store_true", help="Generate diagnostics after successful processing.")
    parser.add_argument("--diagnostics-qf", default="lr1", help="Quality flag product to use for diagnostics (default: lr1).")
    parser.add_argument("--diagnostics-raw", action="store_true", help="Use raw rather than adjusted variables for diagnostics.")
    parser.add_argument("--force", action="store_true", help="Force reprocessing of deployments even if they previously completed successfully.")
    parser.add_argument("--force-failed", action="store_true", help="Re-run deployments whose latest status is failed.")
    parser.add_argument("--include-disabled", action="store_true", help="Include deployments whose PROCESS flag is disabled in the catalog.")
    parser.add_argument("--deployment", action="append", default=[], help="Restrict the run to one deployment; may be supplied multiple times.")
    parser.add_argument("--state-dir", default=None, help="Override the directory used for batch state and reports.")
    parser.add_argument("-j", "--jobs", type=int, default=1, help="Run up to N deployments in parallel (default: 1).")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print deployment log output to the terminal.")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    result = run_all_deployments(
        processdir=args.processdir,
        config_file=args.config_file,
        machine=args.machine,
        notlc=args.notlc,
        diagnostics=args.diagnostics,
        diagnostics_qf=args.diagnostics_qf,
        diagnostics_raw=args.diagnostics_raw,
        force=args.force,
        force_failed=args.force_failed,
        include_disabled=args.include_disabled,
        deployments=args.deployment,
        state_dir=args.state_dir,
        jobs=args.jobs,
        verbose=args.verbose,
    )
    print(result.summary_markdown)
    return 0 if result.failed_count == 0 else 1
