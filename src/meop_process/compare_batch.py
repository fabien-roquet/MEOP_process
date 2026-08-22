from __future__ import annotations

import argparse
import concurrent.futures
import contextlib
import csv
import io
import json
import sys
import time
import traceback
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from .catalog.deployments import load_deployment_catalog
from .catalog.filenames import deployment_from_smru_name, list_fname_prof, smru_name_from_fname_prof
from .compare_cli import CalibrationStatusError, generate_calibration_plots
from .config.loader import load_config
from .models import MeopConfig


STATE_FILE_NAME = "calibration_plot_status.json"
CALIBRATION_METHOD_VERSION = 4
_PARALLEL_EXECUTOR_FACTORY = concurrent.futures.ProcessPoolExecutor


@dataclass(frozen=True)
class CalibrationPlotRunResult:
    smru_name: str
    deployment: str
    status: str
    started_at: str
    finished_at: str
    duration_seconds: float
    log_path: Path
    written_files: tuple[Path, ...] = ()
    message: str = ""
    method_version: int = CALIBRATION_METHOD_VERSION

    def as_dict(self) -> dict[str, object]:
        return {
            "smru_name": self.smru_name,
            "deployment": self.deployment,
            "status": self.status,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "duration_seconds": self.duration_seconds,
            "log_path": str(self.log_path),
            "written_files": [str(path) for path in self.written_files],
            "message": self.message,
            "method_version": self.method_version,
        }


@dataclass(frozen=True)
class CompareBatchRunResult:
    run_id: str
    output_dir: Path
    log_dir: Path
    state_path: Path
    summary_markdown: Path
    summary_csv: Path
    tag_results: tuple[CalibrationPlotRunResult, ...]

    @property
    def success_count(self) -> int:
        return sum(item.status == "success" for item in self.tag_results)

    @property
    def skipped_count(self) -> int:
        return sum(item.status == "skipped" for item in self.tag_results)

    @property
    def failed_count(self) -> int:
        return sum(item.status == "failed" for item in self.tag_results)

    @property
    def no_reference_count(self) -> int:
        return sum(item.status == "no_reference" for item in self.tag_results)

    @property
    def insufficient_reference_count(self) -> int:
        return sum(item.status == "insufficient_reference" for item in self.tag_results)

    @property
    def insufficient_target_count(self) -> int:
        return sum(item.status == "insufficient_target" for item in self.tag_results)

    @property
    def invalid_target_count(self) -> int:
        return sum(item.status == "invalid_target" for item in self.tag_results)

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
            "no_reference_count": self.no_reference_count,
            "insufficient_reference_count": self.insufficient_reference_count,
            "insufficient_target_count": self.insufficient_target_count,
            "invalid_target_count": self.invalid_target_count,
            "tag_results": [item.as_dict() for item in self.tag_results],
        }


@dataclass(frozen=True)
class _PendingTag:
    index: int
    smru_name: str
    log_path: Path


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
    root = Path(state_dir) if state_dir is not None else config.datadir / "batch" / "compare"
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


def _eligible_deployments(
    config: MeopConfig,
    *,
    include_disabled: bool = False,
    selected: Iterable[str] | None = None,
) -> list[str]:
    requested = {item.strip() for item in (selected or []) if str(item).strip()}
    if requested:
        return sorted(requested)

    catalog = load_deployment_catalog(config, persist=False)
    deployments: list[str] = []
    for deployment, record in sorted(catalog.items()):
        if not include_disabled and not _truthy_process(getattr(record, "process", "")):
            continue
        deployments.append(deployment)
    if deployments:
        return sorted(set(deployments))

    if not config.final_dataset_dir.exists():
        return []
    return sorted(path.name for path in config.final_dataset_dir.iterdir() if path.is_dir())


def discover_calibration_tags(
    config: MeopConfig,
    *,
    deployments: Iterable[str] | None = None,
    tags: Iterable[str] | None = None,
    include_disabled: bool = False,
) -> list[str]:
    """Return SMRU tags with processed profile files, suitable for calibration plots."""
    explicit_tags = sorted({item.strip() for item in (tags or []) if str(item).strip()})
    if explicit_tags:
        return explicit_tags

    smru_names: set[str] = set()
    for deployment in _eligible_deployments(config, include_disabled=include_disabled, selected=deployments):
        for path in list_fname_prof(deployment=deployment, config=config):
            try:
                smru_names.add(smru_name_from_fname_prof(path.name))
            except Exception:
                continue
    return sorted(smru_names)


def _expected_outputs_exist(config: MeopConfig, smru_name: str) -> bool:
    deployment = deployment_from_smru_name(smru_name)
    output_dir = config.plotdir / deployment
    if not (output_dir / f"{smru_name}_calibration_offsets.csv").exists():
        return False
    if not (output_dir / f"{smru_name}_calibration_offsets.png").exists():
        return False
    calibration_pngs = [
        path
        for path in output_dir.glob(f"{smru_name}_calibration*.png")
        if not path.name.endswith("_calibration_offsets.png")
    ]
    return bool(calibration_pngs)


def _should_skip(
    smru_name: str,
    *,
    state: dict[str, dict[str, object]],
    force: bool,
    config: MeopConfig,
) -> tuple[bool, str]:
    if force:
        return False, "forced"
    entry = state.get(smru_name)
    if not entry:
        return False, "not yet plotted"
    if str(entry.get("status", "")).strip().lower() != "success":
        return False, "previous run not successful"
    if entry.get("method_version") != CALIBRATION_METHOD_VERSION:
        return False, "success marker belongs to an older calibration method version"
    if not _expected_outputs_exist(config, smru_name):
        return False, "success marker exists but plot outputs are missing"
    return True, "already completed successfully"


def _write_csv(path: Path, results: Iterable[CalibrationPlotRunResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "smru_name",
                "deployment",
                "status",
                "started_at",
                "finished_at",
                "duration_seconds",
                "written_count",
                "log_path",
                "message",
                "method_version",
            ]
        )
        for result in results:
            writer.writerow(
                [
                    result.smru_name,
                    result.deployment,
                    result.status,
                    result.started_at,
                    result.finished_at,
                    f"{result.duration_seconds:.3f}",
                    len(result.written_files),
                    str(result.log_path),
                    result.message,
                    result.method_version,
                ]
            )


def _write_markdown(path: Path, run_id: str, results: Iterable[CalibrationPlotRunResult]) -> None:
    items = list(results)
    success = [item for item in items if item.status == "success"]
    skipped = [item for item in items if item.status == "skipped"]
    no_reference = [item for item in items if item.status == "no_reference"]
    insufficient_reference = [item for item in items if item.status == "insufficient_reference"]
    insufficient_target = [item for item in items if item.status == "insufficient_target"]
    invalid_target = [item for item in items if item.status == "invalid_target"]
    failed = [item for item in items if item.status == "failed"]
    lines = [
        f"# MEOP CORA calibration plot batch report ({run_id})",
        "",
        f"Calibration method version: **{CALIBRATION_METHOD_VERSION}**  ",
        f"Success: **{len(success)}**  ",
        f"Skipped: **{len(skipped)}**  ",
        f"No reference: **{len(no_reference)}**  ",
        f"Insufficient reference: **{len(insufficient_reference)}**  ",
        f"Insufficient target: **{len(insufficient_target)}**  ",
        f"Invalid target: **{len(invalid_target)}**  ",
        f"Failed: **{len(failed)}**",
        "",
        "## Result table",
        "",
        "| Tag | Deployment | Status | Files | Duration (s) | Message | Log |",
        "|---|---|---:|---:|---:|---|---|",
    ]
    for item in items:
        lines.append(
            f"| {item.smru_name} | {item.deployment} | {item.status} | "
            f"{len(item.written_files)} | {item.duration_seconds:.1f} | {item.message or ''} | {item.log_path.name} |"
        )
    if failed:
        lines.extend(["", "## Failed tags", ""])
        for item in failed:
            lines.append(f"- **{item.smru_name}**: {item.message}")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def _print_log_to_stdout(smru_name: str, log_path: Path) -> None:
    text = log_path.read_text(encoding="utf-8") if log_path.exists() else ""
    if not text:
        return
    print(f"[{smru_name}] log start")
    print(text, end="" if text.endswith("\n") else "\n")
    print(f"[{smru_name}] log end")


def _execute_tag(
    *,
    config: MeopConfig,
    smru_name: str,
    log_path: Path,
    stream_stdout: bool,
) -> CalibrationPlotRunResult:
    started = _utc_now()
    started_timer = time.perf_counter()
    targets: list[io.TextIOBase] = []
    if stream_stdout:
        targets.append(sys.stdout)
    written: tuple[Path, ...] = ()
    with log_path.open("w", encoding="utf-8") as log_handle:
        targets.insert(0, log_handle)
        tee = _Tee(*targets)
        try:
            with contextlib.redirect_stdout(tee), contextlib.redirect_stderr(tee):
                print(f"[{smru_name}] generating CORA calibration plots")
                written = tuple(generate_calibration_plots(smru_name, config=config))
                for path in written:
                    print(f"wrote: {path}")
                status = "success"
                message = f"{len(written)} files"
        except CalibrationStatusError as exc:
            status = exc.status
            written = exc.written_files
            message = str(exc)
            with contextlib.redirect_stdout(tee), contextlib.redirect_stderr(tee):
                print(f"[{smru_name}] {status}: {message}")
                for path in written:
                    print(f"wrote: {path}")
        except Exception as exc:  # pragma: no cover - exercised by dedicated tests
            status = "failed"
            message = str(exc)
            with contextlib.redirect_stdout(tee), contextlib.redirect_stderr(tee):
                traceback.print_exc()
    finished = _utc_now()
    duration = time.perf_counter() - started_timer
    return CalibrationPlotRunResult(
        smru_name=smru_name,
        deployment=deployment_from_smru_name(smru_name),
        status=status,
        started_at=started.isoformat(),
        finished_at=finished.isoformat(),
        duration_seconds=duration,
        log_path=log_path,
        written_files=written,
        message=message,
    )


def _run_pending_tags(
    *,
    config: MeopConfig,
    pending: list[_PendingTag],
    jobs: int,
    verbose: bool,
):
    if jobs <= 1 or len(pending) <= 1:
        for task in pending:
            result = _execute_tag(
                config=config,
                smru_name=task.smru_name,
                log_path=task.log_path,
                stream_stdout=verbose,
            )
            yield task.index, result
        return

    max_workers = min(jobs, len(pending))
    with _PARALLEL_EXECUTOR_FACTORY(max_workers=max_workers) as executor:
        futures: dict[concurrent.futures.Future, _PendingTag] = {}
        for task in pending:
            future = executor.submit(
                _execute_tag,
                config=config,
                smru_name=task.smru_name,
                log_path=task.log_path,
                stream_stdout=False,
            )
            futures[future] = task
        for future in concurrent.futures.as_completed(futures):
            task = futures[future]
            try:
                result = future.result()
            except Exception as exc:  # pragma: no cover - defensive
                finished = _utc_now()
                task.log_path.write_text(f"[{task.smru_name}] worker failure: {exc}\n", encoding="utf-8")
                result = CalibrationPlotRunResult(
                    smru_name=task.smru_name,
                    deployment=deployment_from_smru_name(task.smru_name),
                    status="failed",
                    started_at=finished.isoformat(),
                    finished_at=finished.isoformat(),
                    duration_seconds=0.0,
                    log_path=task.log_path,
                    message=str(exc),
                )
            if verbose:
                _print_log_to_stdout(task.smru_name, task.log_path)
            yield task.index, result


def run_calibration_plots_batch(
    *,
    config: MeopConfig | None = None,
    processdir: str | Path | None = None,
    config_file: str | Path | None = None,
    machine: str | None = None,
    deployments: Iterable[str] | None = None,
    tags: Iterable[str] | None = None,
    include_disabled: bool = False,
    force: bool = False,
    state_dir: str | Path | None = None,
    jobs: int = 1,
    verbose: bool = False,
) -> CompareBatchRunResult:
    if jobs < 1:
        raise ValueError("jobs must be at least 1")
    cfg = config or load_config(processdir=processdir, config_file=config_file, machine=machine)
    if cfg.cora_dir is None:
        raise ValueError("references.cora_dir is not configured")

    selected_tags = discover_calibration_tags(
        cfg,
        deployments=deployments,
        tags=tags,
        include_disabled=include_disabled,
    )
    root = _state_root(cfg, state_dir)
    run_id = _timestamp()
    output_dir = root / "runs" / run_id
    log_dir = output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    state_path = root / "latest" / STATE_FILE_NAME
    state = _load_state(state_path)

    results_by_index: list[tuple[int, CalibrationPlotRunResult]] = []
    pending: list[_PendingTag] = []
    for index, smru_name in enumerate(selected_tags):
        log_path = log_dir / f"{smru_name}.log"
        skip, reason = _should_skip(smru_name, state=state, force=force, config=cfg)
        if skip:
            started = _utc_now()
            finished = _utc_now()
            result = CalibrationPlotRunResult(
                smru_name=smru_name,
                deployment=deployment_from_smru_name(smru_name),
                status="skipped",
                started_at=started.isoformat(),
                finished_at=finished.isoformat(),
                duration_seconds=(finished - started).total_seconds(),
                log_path=log_path,
                message=reason,
            )
            log_path.write_text(f"[{smru_name}] skipped: {reason}\n", encoding="utf-8")
            if verbose:
                _print_log_to_stdout(smru_name, log_path)
            results_by_index.append((index, result))
            state[smru_name] = {**result.as_dict(), "run_id": run_id}
            continue
        pending.append(_PendingTag(index=index, smru_name=smru_name, log_path=log_path))

    for index, result in _run_pending_tags(config=cfg, pending=pending, jobs=jobs, verbose=verbose):
        results_by_index.append((index, result))
        state[result.smru_name] = {**result.as_dict(), "run_id": run_id}
        _write_state(state_path, state)

    results = [result for _, result in sorted(results_by_index, key=lambda item: item[0])]
    summary_csv = output_dir / "summary.csv"
    summary_markdown = output_dir / "summary.md"
    _write_csv(summary_csv, results)
    _write_markdown(summary_markdown, run_id, results)
    _write_state(state_path, state)

    return CompareBatchRunResult(
        run_id=run_id,
        output_dir=output_dir,
        log_dir=log_dir,
        state_path=state_path,
        summary_markdown=summary_markdown,
        summary_csv=summary_csv,
        tag_results=tuple(results),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate CORA calibration plots for many MEOP tags.")
    parser.add_argument("--processdir", default=None, help="Override the MEOP process directory.")
    parser.add_argument("--config-file", default=None, help="Explicit path to a runtime config JSON file.")
    parser.add_argument("--machine", default=None, help="Machine entry key from the runtime config JSON.")
    parser.add_argument("--deployment", action="append", default=[], help="Restrict the run to one deployment; may be supplied multiple times.")
    parser.add_argument("--tag", action="append", default=[], help="Restrict the run to one SMRU tag; may be supplied multiple times.")
    parser.add_argument("--include-disabled", action="store_true", help="Include deployments whose PROCESS flag is disabled in the catalog.")
    parser.add_argument("--force", action="store_true", help="Regenerate plots even when the latest batch state is successful.")
    parser.add_argument("--dry-run", action="store_true", help="List selected tags without generating plots or writing batch state.")
    parser.add_argument("--state-dir", default=None, help="Override the directory used for compare batch state and reports.")
    parser.add_argument("-j", "--jobs", type=int, default=None, help="Run up to N tags in parallel (default: config or 1).")
    parser.add_argument("-v", "--verbose", action="store_true", default=None, help="Print tag log output to the terminal.")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    try:
        cfg = load_config(processdir=args.processdir, config_file=args.config_file, machine=args.machine)
    except (FileNotFoundError, ValueError) as exc:
        parser.error(str(exc))
    if cfg.config_path is None:
        parser.error(
            f"runtime config is required and was not found (expected {cfg.processdir / 'configs.json'}). "
            "Run 'meop-process --bootstrap-data' once, or pass --config-file / set MEOP_CONFIG_FILE."
        )
    if cfg.cora_dir is None:
        parser.error("references.cora_dir is not configured in configs.json")

    jobs = args.jobs if args.jobs is not None else cfg.batch_defaults.jobs
    verbose = args.verbose if args.verbose is not None else cfg.batch_defaults.verbose
    if args.dry_run:
        selected_tags = discover_calibration_tags(
            cfg,
            deployments=args.deployment,
            tags=args.tag,
            include_disabled=args.include_disabled,
        )
        print(f"Selected tags: {len(selected_tags)}")
        for smru_name in selected_tags:
            print(smru_name)
        return 0

    result = run_calibration_plots_batch(
        config=cfg,
        deployments=args.deployment,
        tags=args.tag,
        include_disabled=args.include_disabled,
        force=args.force,
        state_dir=args.state_dir,
        jobs=jobs,
        verbose=verbose,
    )
    print(result.summary_markdown)
    return 0 if result.failed_count == 0 else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
