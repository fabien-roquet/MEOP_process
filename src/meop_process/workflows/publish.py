"""Orchestrate the full publish workflow."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..models import MeopConfig
from ..publishing.attributes import update_global_attributes
from ..publishing.lists import build_list_profiles
from ..publishing.ncfiles import PublishResult, publish_ncfiles


@dataclass(frozen=True)
class PublishWorkflowResult:
    output_dir: Path
    nc_result: PublishResult
    patched_attrs: tuple[Path, ...]
    list_profiles_path: Path | None
    list_tags_path: Path | None
    list_deployments_path: Path | None

    def as_dict(self) -> dict[str, object]:
        return {
            "output_dir": str(self.output_dir),
            "published_files": [str(p) for p in self.nc_result.written_files],
            "processed_tags": list(self.nc_result.processed_tags),
            "patched_attrs": [str(p) for p in self.patched_attrs],
            "list_profiles_path": str(self.list_profiles_path) if self.list_profiles_path else None,
            "list_tags_path": str(self.list_tags_path) if self.list_tags_path else None,
            "list_deployments_path": str(self.list_deployments_path) if self.list_deployments_path else None,
        }


def publish(
    config: MeopConfig,
    *,
    deployment: str = "",
    smru_name: str = "",
    output_dir: Path | str | None = None,
    create_files: bool = True,
    update_attrs: bool = True,
    list_profiles: bool = True,
    list_tags: bool = True,
    rebuild: bool = False,
    verbose: bool = False,
) -> PublishWorkflowResult:
    """Run the full (or partial) publish workflow.

    Parameters
    ----------
    config:
        Runtime configuration.
    deployment:
        Restrict to this deployment code.
    smru_name:
        Restrict to this single tag.
    output_dir:
        Output directory for published files.  Defaults to ``config.publicdir_ctd``.
    create_files:
        Whether to create ``*_all_prof.nc`` files.
    update_attrs:
        Whether to patch global attributes on published NC files.
    list_profiles:
        Whether to write ``list_profiles.csv``.
    list_tags:
        Whether to write ``list_tags.csv`` and ``list_deployments.csv`` via
        the existing metadata summaries tool.
    rebuild:
        Force recreation of existing output files.
    verbose:
        Print progress to stdout.
    """
    dest = Path(output_dir) if output_dir else config.publicdir_ctd
    dest.mkdir(parents=True, exist_ok=True)

    nc_result = PublishResult()
    if create_files:
        if verbose:
            print(f"Creating _all_prof.nc files in {dest} …")
        nc_result = publish_ncfiles(
            config,
            output_dir=dest,
            deployment=deployment,
            smru_name=smru_name,
            rebuild=rebuild,
        )
        if verbose:
            print(f"  {len(nc_result.written_files)} files written, {len(nc_result.processed_tags)} tags processed")

    patched: tuple[Path, ...] = ()
    if update_attrs:
        if verbose:
            print("Updating global attributes …")
        patched_list = update_global_attributes(dest, version=config.version, verbose=verbose)
        patched = tuple(patched_list)

    profiles_path: Path | None = None
    if list_profiles:
        if verbose:
            print("Building list_profiles.csv …")
        profiles_path = build_list_profiles(dest, verbose=verbose)

    tags_path: Path | None = None
    deployments_path: Path | None = None
    if list_tags:
        if verbose:
            print("Refreshing list_tags.csv and list_deployments.csv …")
        from ..metadata.summaries import update_metadata_summaries

        summary = update_metadata_summaries(config, force=rebuild, output_dir=dest)
        tags_path = summary.list_tags_path
        deployments_path = summary.list_deployments_path

    return PublishWorkflowResult(
        output_dir=dest,
        nc_result=nc_result,
        patched_attrs=patched,
        list_profiles_path=profiles_path,
        list_tags_path=tags_path,
        list_deployments_path=deployments_path,
    )
