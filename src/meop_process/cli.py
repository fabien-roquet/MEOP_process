from __future__ import annotations

import argparse
from collections.abc import Sequence

from .api import (
    apply_adjustments,
    bootstrap_data_store,
    create_fr0,
    create_hr2,
    describe_runtime_data_layout,
    generate_diagnostics,
    import_raw_data,
    process_tags,
    run_all_deployments,
    update_config_files,
    update_metadata_summaries,
    validate_runtime_data_layout,
)
from .catalog.filenames import deployment_from_smru_name
from .config.loader import load_config


ACTIONS = (
    "do_all",
    "update_config_files",
    "import_data",
    "process_data",
    "apply_adjustments",
    "create_fr0",
    "create_hr2",
    "diagnostics",
)

UTILITY_ACTIONS = (
    "bootstrap_data",
    "show_config",
    "show_data_layout",
    "validate_data_layout",
    "refresh_metadata_summaries",
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run MEOP processing workflows with the pure-Python package.")
    parser.add_argument("--smru_name", default="", help="Process only one SMRU PLATFORM CODE.")
    parser.add_argument("--deployment", default="", help="Process all tags in a deployment.")
    parser.add_argument("--processdir", default=None, help="Override the MEOP process directory.")
    parser.add_argument("--config-file", default=None, help="Explicit path to a runtime config JSON file.")
    parser.add_argument("--machine", default=None, help="Machine entry key from the runtime config JSON.")
    parser.add_argument("--bootstrap-data", action="store_true", help="Create runtime data folders and seed packaged CSV tables.")
    parser.add_argument("--show-config", action="store_true", help="Print the resolved runtime config source and key dataset/output paths.")
    parser.add_argument("--show-data-layout", action="store_true", help="Print where tables, raw inputs, references, and outputs are expected.")
    parser.add_argument("--validate-data-layout", action="store_true", help="Check presence of fixed-path data files and directories.")
    parser.add_argument("--do_all", action="store_true", help="Run the main pure-Python workflow stages.")
    parser.add_argument("--update_config_files", action="store_true", help="Optionally sync JSON config files into data/data_raw/config_files.")
    parser.add_argument("--import_data", action="store_true", help="Prepare raw ODV data staged under data/data_raw/raw_smru_data_odv.")
    parser.add_argument("--process_data", action="store_true", help="Run the main processing chain.")
    parser.add_argument("--apply-adjustments", dest="apply_adjustments", action="store_true", help="Apply delayed-mode offsets and error estimates to existing products.")
    parser.add_argument("--create_fr0", action="store_true", help="Create the FR0 product from full-resolution HR text files when available.")
    parser.add_argument("--notlc", action="store_true", help="Use the no-TLC branch for process_data.")
    parser.add_argument("--create_hr2", action="store_true", help="Create the HR2 product.")
    parser.add_argument("--diagnostics", action="store_true", help="Generate standard diagnostics plots for processed products.")
    parser.add_argument("--diagnostics-qf", default=None, help="Quality flag product to use for diagnostics (default: config or lr1).")
    parser.add_argument("--diagnostics-raw", action="store_true", default=None, help="Use raw rather than adjusted variables for diagnostics.")
    parser.add_argument(
        "--diagnostics-part",
        action="append",
        choices=("tag", "deployment", "overview", "all"),
        default=None,
        help="Restrict diagnostics to one or more parts: tag, deployment, overview, or all (default: all).",
    )
    parser.add_argument("--run-all-deployments", action="store_true", help="Process all deployments from the catalog, continuing past failures and keeping resumable state.")
    parser.add_argument("--refresh-metadata-summaries", dest="refresh_metadata_summaries", action="store_true", help="Refresh list_tags.csv and list_deployments.csv without processing deployments.")
    parser.add_argument("--force", action="store_true", help="Force reprocessing of deployments even if they previously completed successfully.")
    parser.add_argument("--force-failed", action="store_true", help="Re-run deployments whose latest status is failed.")
    parser.add_argument("--include-disabled", action="store_true", help="Include deployments whose PROCESS flag is disabled in the catalog.")
    parser.add_argument("--state-dir", default=None, help="Override the directory used for batch state and reports.")
    parser.add_argument("--notify-email", action="append", default=None, help="Send the batch summary email to this address; may be supplied multiple times.")
    parser.add_argument("--notify-when", choices=("always", "success", "failure"), default=None, help="When to send completion emails.")
    parser.add_argument("--notify-attach", action="append", choices=("summary_md", "summary_csv", "comparison_md"), default=None, help="Attachments to include in completion emails.")
    parser.add_argument("--no-notify", action="store_true", help="Disable completion email even if enabled in the runtime config.")
    parser.add_argument("-j", "--jobs", type=int, default=None, help="Run up to N deployments in parallel during batch mode (default: config or 1).")
    parser.add_argument("-v", "--verbose", action="store_true", default=None, help="Print batch deployment logs to the terminal.")
    return parser


def main(argv: Sequence[str] | None = None, *, config=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    workflow_requested = any(getattr(args, action) for action in ACTIONS) or args.run_all_deployments
    utility_requested = any(getattr(args, action) for action in UTILITY_ACTIONS)
    if not (workflow_requested or utility_requested):
        parser.error("at least one action flag is required")

    requires_selection = any(getattr(args, action) for action in ACTIONS if action != "diagnostics")
    if workflow_requested and not args.run_all_deployments and requires_selection and not (args.smru_name or args.deployment):
        parser.error("one of --smru_name or --deployment is required for workflow actions")

    deployment = args.deployment or (deployment_from_smru_name(args.smru_name) if args.smru_name else "")
    try:
        cfg = config or load_config(processdir=args.processdir, config_file=args.config_file, machine=args.machine)
    except (FileNotFoundError, ValueError) as exc:
        parser.error(str(exc))
    bootstrapped_now = False

    if args.bootstrap_data and cfg.config_path is None:
        bootstrap_data_store(cfg)
        bootstrapped_now = True
        try:
            cfg = config or load_config(processdir=args.processdir, config_file=args.config_file, machine=args.machine)
        except (FileNotFoundError, ValueError) as exc:
            parser.error(str(exc))

    if cfg.config_path is None and not args.bootstrap_data and config is None:
        expected = cfg.processdir / "configs.json"
        parser.error(
            f"runtime config is required and was not found (expected {expected}). "
            "Run --bootstrap-data once, or pass --config-file / set MEOP_CONFIG_FILE."
        )
    diagnostics_qf = args.diagnostics_qf or cfg.diagnostics_defaults.qf
    diagnostics_raw = args.diagnostics_raw if args.diagnostics_raw is not None else (not cfg.diagnostics_defaults.adjusted)
    diagnostics_parts = args.diagnostics_part or list(cfg.diagnostics_defaults.parts)
    jobs = args.jobs if args.jobs is not None else cfg.batch_defaults.jobs
    verbose = args.verbose if args.verbose is not None else cfg.batch_defaults.verbose

    success = True

    if args.bootstrap_data and not bootstrapped_now:
        bootstrap_data_store(cfg)

    if args.show_config:
        config_source = str(cfg.config_path) if cfg.config_path else "(none; defaults only)"
        print(f"config source     : {config_source}")
        print(f"machine           : {cfg.machine}")
        print(f"processdir        : {cfg.processdir}")
        print(f"datadir           : {cfg.datadir}")
        print(f"publicdir         : {cfg.publicdir}")
        print(f"raw_odv_dir       : {cfg.raw_odv_dir}")
        print(f"raw_hr_dir        : {cfg.raw_hr_dir}")
        print(f"final_dataset_dir : {cfg.final_dataset_dir}")
        print(f"plots_by_tags     : {cfg.plotdir}")
        print(f"plots_by_deploy   : {cfg.plots_by_deployment_dir}")
        print(f"plots_overview    : {cfg.plots_overview_dir}")
        print(f"public_ctd_dir    : {cfg.publicdir_ctd}")
        print(f"cora_dir          : {cfg.cora_dir or '(unset)'}")
        print(f"reference_dataset : {cfg.reference_dataset_dir or '(unset)'}")

    if args.show_data_layout:
        print(describe_runtime_data_layout(cfg, as_text=True))

    if args.validate_data_layout:
        records = validate_runtime_data_layout(cfg)
        for record in records:
            status = "OK" if record["exists"] else "MISSING"
            required = "required" if record["required"] else "optional"
            print(f"{status:7s} {required:8s} {record['category']:16s} {record['name']}: {record['path']}")
        if any(record["required"] and not record["exists"] for record in records):
            success = False

    if args.refresh_metadata_summaries:
        summary = update_metadata_summaries(config=cfg, force=args.force)
        print(summary["list_tags_path"])
        print(summary["list_deployments_path"])

    if not workflow_requested:
        return 0 if success else 1

    if args.run_all_deployments:
        result = run_all_deployments(
            config=cfg,
            notlc=args.notlc,
            diagnostics=args.diagnostics,
            diagnostics_qf=diagnostics_qf,
            diagnostics_raw=diagnostics_raw,
            diagnostics_parts=diagnostics_parts,
            notify_email=args.notify_email,
            notify_when=args.notify_when,
            notify_attach=args.notify_attach,
            notifications_enabled=False if args.no_notify else None,
            force=args.force,
            force_failed=args.force_failed,
            include_disabled=args.include_disabled,
            deployments=[args.deployment] if args.deployment else [],
            state_dir=args.state_dir,
            jobs=jobs,
            verbose=verbose,
        )
        print(result["summary_markdown"])
        return 0 if result.get("failed_count", 0) == 0 else 1

    if args.update_config_files or args.do_all:
        update_config_files(config=cfg)
    if args.import_data or args.do_all:
        success = import_raw_data(deployment=deployment, config=cfg) and success
    if args.process_data or args.do_all:
        success = process_tags(
            deployment=deployment,
            smru_name=args.smru_name,
            notlc=args.notlc,
            config=cfg,
        ) and success
    if args.apply_adjustments:
        apply_adjustments(deployment=deployment, smru_name=args.smru_name, config=cfg)
    if args.create_fr0:
        create_fr0(deployment=deployment, smru_name=args.smru_name, config=cfg)
    if args.create_hr2 or args.do_all:
        success = create_hr2(deployment=deployment, smru_name=args.smru_name, config=cfg) and success
    if args.diagnostics:
        generate_diagnostics(
            deployment=deployment,
            smru_name=args.smru_name,
            qf=diagnostics_qf,
            adjusted=not diagnostics_raw,
            parts=diagnostics_parts,
            config=cfg,
        )

    return 0 if success else 1
