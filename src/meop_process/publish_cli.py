"""CLI entry point for ``meop-publish``."""
from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from .config.loader import load_config
from .publishing.versions import load_version_registry
from .workflows.publish import publish


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="meop-publish",
        description="Publish MEOP-CTD processed products to the public release directory.",
    )
    parser.add_argument("--smru_name", default="", help="Publish only this SMRU platform code.")
    parser.add_argument("--deployment", default="", help="Publish all tags in this deployment.")
    parser.add_argument("--processdir", default=None, help="Override the MEOP process directory.")
    parser.add_argument("--config-file", default=None, help="Explicit path to a runtime config JSON file.")
    parser.add_argument("--machine", default=None, help="Machine entry key from the runtime config JSON.")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Override output directory (default: publicdir_ctd from config).",
    )
    parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Force recreation of already-published files.",
    )
    parser.add_argument(
        "--no-create-files",
        dest="create_files",
        action="store_false",
        default=True,
        help="Skip creation of *_all_prof.nc files.",
    )
    parser.add_argument(
        "--no-update-attrs",
        dest="update_attrs",
        action="store_false",
        default=True,
        help="Skip patching global attributes on published NC files.",
    )
    parser.add_argument(
        "--no-list-profiles",
        dest="list_profiles",
        action="store_false",
        default=True,
        help="Skip writing list_profiles.csv.",
    )
    parser.add_argument(
        "--no-list-tags",
        dest="list_tags",
        action="store_false",
        default=True,
        help="Skip writing list_tags.csv and list_deployments.csv.",
    )
    parser.add_argument(
        "--build-maps",
        dest="build_maps",
        action="store_true",
        default=None,
        help="Generate overview map PNGs from list_profiles.csv.",
    )
    parser.add_argument(
        "--no-build-maps",
        dest="build_maps",
        action="store_false",
        help="Disable overview map generation.",
    )
    parser.add_argument(
        "--build-plots",
        dest="build_plots",
        action="store_true",
        default=None,
        help="Generate per-tag and deployment diagnostic figures.",
    )
    parser.add_argument(
        "--no-build-plots",
        dest="build_plots",
        action="store_false",
        help="Disable diagnostic figure generation.",
    )
    parser.add_argument(
        "--build-site",
        dest="build_site",
        action="store_true",
        default=None,
        help="Generate static HTML index pages from existing diagnostic plots.",
    )
    parser.add_argument(
        "--no-build-site",
        dest="build_site",
        action="store_false",
        help="Disable static HTML site generation.",
    )
    parser.add_argument(
        "--release-status",
        choices=("development", "published"),
        default=None,
        help="Version lifecycle status to record in public/versions.json.",
    )
    parser.add_argument(
        "--list-versions",
        action="store_true",
        default=False,
        help="List known dataset versions from public/versions.json and exit.",
    )
    parser.add_argument("-v", "--verbose", action="store_true", default=False, help="Print progress to stdout.")
    return parser


def main(argv: Sequence[str] | None = None, *, config=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        cfg = config or load_config(processdir=args.processdir, config_file=args.config_file, machine=args.machine, require_config=True)
    except (FileNotFoundError, ValueError) as exc:
        parser.error(str(exc))

    if args.list_versions:
        records = load_version_registry(cfg.publicdir)
        if not records:
            print("No versions registered yet.")
            return 0
        print("Version\tStatus\tFirst Seen\tLast Updated\tOutput Dir")
        for record in records:
            print(f"{record.version}\t{record.status}\t{record.first_seen}\t{record.last_updated}\t{record.output_dir}")
        return 0

    output_dir = Path(args.output_dir) if args.output_dir else None
    build_maps = args.build_maps if args.build_maps is not None else cfg.publish_defaults.build_maps
    build_plots = args.build_plots if args.build_plots is not None else cfg.publish_defaults.build_plots
    build_site = args.build_site if args.build_site is not None else cfg.publish_defaults.build_site
    release_status = args.release_status or cfg.publish_defaults.release_status

    result = publish(
        cfg,
        deployment=args.deployment,
        smru_name=args.smru_name,
        output_dir=output_dir,
        create_files=args.create_files,
        update_attrs=args.update_attrs,
        list_profiles=args.list_profiles,
        list_tags=args.list_tags,
        build_maps=build_maps,
        build_plots=build_plots,
        build_site=build_site,
        release_status=release_status,
        rebuild=args.rebuild,
        verbose=args.verbose,
    )

    info = result.as_dict()
    print(f"Output directory : {info['output_dir']}")
    print(f"Published files  : {len(info['published_files'])}")
    if info["list_profiles_path"]:
        print(f"list_profiles.csv: {info['list_profiles_path']}")
    if info["list_tags_path"]:
        print(f"list_tags.csv    : {info['list_tags_path']}")
    if info["list_deployments_path"]:
        print(f"list_deployments : {info['list_deployments_path']}")
    if info["map_paths"]:
        print(f"Maps written     : {len(info['map_paths'])}")
    if info["plot_paths"]:
        print(f"Plots written    : {len(info['plot_paths'])}")
    if info["site_paths"]:
        print(f"HTML pages       : {len(info['site_paths'])}")

    return 0
