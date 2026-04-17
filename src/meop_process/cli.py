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
    parser.add_argument("--diagnostics-qf", default="lr1", help="Quality flag product to use for diagnostics (default: lr1).")
    parser.add_argument("--diagnostics-raw", action="store_true", help="Use raw rather than adjusted variables for diagnostics.")
    parser.add_argument("--run-all-deployments", action="store_true", help="Process all deployments from the catalog, continuing past failures and keeping resumable state.")
    parser.add_argument("--refresh-metadata-summaries", dest="refresh_metadata_summaries", action="store_true", help="Refresh list_tags.csv and list_deployments.csv without processing deployments.")
    parser.add_argument("--force", action="store_true", help="Force reprocessing of deployments even if they previously completed successfully.")
    parser.add_argument("--force-failed", action="store_true", help="Re-run deployments whose latest status is failed.")
    parser.add_argument("--include-disabled", action="store_true", help="Include deployments whose PROCESS flag is disabled in the catalog.")
    parser.add_argument("--state-dir", default=None, help="Override the directory used for batch state and reports.")
    parser.add_argument("-j", "--jobs", type=int, default=1, help="Run up to N deployments in parallel during batch mode (default: 1).")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print batch deployment logs to the terminal.")
    return parser


def main(argv: Sequence[str] | None = None, *, config=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    workflow_requested = any(getattr(args, action) for action in ACTIONS) or args.run_all_deployments
    utility_requested = any(getattr(args, action) for action in UTILITY_ACTIONS)
    if not (workflow_requested or utility_requested):
        parser.error("at least one action flag is required")

    if workflow_requested and not args.run_all_deployments and not (args.smru_name or args.deployment):
        parser.error("one of --smru_name or --deployment is required for workflow actions")

    deployment = args.deployment or (deployment_from_smru_name(args.smru_name) if args.smru_name else "")
    cfg = config or load_config(processdir=args.processdir, config_file=args.config_file, machine=args.machine)

    success = True

    if args.bootstrap_data:
        bootstrap_data_store(cfg)

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
            diagnostics_qf=args.diagnostics_qf,
            diagnostics_raw=args.diagnostics_raw,
            force=args.force,
            force_failed=args.force_failed,
            include_disabled=args.include_disabled,
            deployments=[args.deployment] if args.deployment else [],
            state_dir=args.state_dir,
            jobs=args.jobs,
            verbose=args.verbose,
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
            qf=args.diagnostics_qf,
            adjusted=not args.diagnostics_raw,
            config=cfg,
        )

    return 0 if success else 1
