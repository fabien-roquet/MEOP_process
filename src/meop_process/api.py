from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from .batch.runner import run_all_deployments as run_all_deployments_batch
from .catalog.deployments import load_info_deployment as load_info_deployment_catalog
from .config.loader import load_config
from .config.sync import sync_external_config_files
from .data.layout import describe_data_layout, format_data_layout, prepare_runtime_environment, validate_data_layout
from .data.smru_sync import sync_smru_data as sync_smru_data_data
from .data.table_validation import validate_runtime_tables as validate_runtime_tables_data
from .io.raw_odv import import_raw_data_zip
from .metadata.patch import update_metadata_from_table
from .metadata.summaries import update_metadata_summaries as update_metadata_summaries_metadata
from .models import MeopConfig, Selection
from .plotting.diagnostics import generate_diagnostics as generate_diagnostics_plotting
from .processing.adjustments import apply_adjustments as apply_adjustments_processing
from .processing.fr0 import create_fr0_python as create_fr0_processing
from .workflows.compare import compare_reference_outputs as compare_reference_outputs_workflow
from .workflows.outputs import create_hr2 as create_hr2_workflow
from .workflows.process import process_tags as process_tags_workflow
from .workflows.publish import publish as publish_workflow


def _resolve_config(
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
) -> MeopConfig:
    return config or load_config(processdir=processdir, config_file=config_file, machine=machine)


def bootstrap_data_store(
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
) -> list[os.PathLike[str]]:
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    return prepare_runtime_environment(cfg)


def describe_runtime_data_layout(
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
    as_text: bool = False,
):
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    if as_text:
        return format_data_layout(cfg)
    return describe_data_layout(cfg)


def validate_runtime_data_layout(
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
) -> list[dict[str, Any]]:
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    return [record.as_dict() for record in validate_data_layout(cfg)]


def validate_runtime_tables(
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
) -> dict[str, Any]:
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    return validate_runtime_tables_data(cfg).as_dict()


def sync_smru_data(
    config: MeopConfig | None = None,
    *,
    source_dir: str | os.PathLike[str],
    apply: bool = False,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
) -> dict[str, Any]:
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    return sync_smru_data_data(cfg, source_dir=source_dir, apply=apply).as_dict()


def update_config_files(
    config: MeopConfig | None = None,
    *,
    source_dir: str | os.PathLike[str] | None = None,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
) -> list[Path]:
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    return sync_external_config_files(cfg, source_dir=source_dir)


def load_info_deployment(
    deployment: str = "",
    smru_name: str = "",
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
):
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    return load_info_deployment_catalog(cfg, deployment=deployment, smru_name=smru_name)


def import_raw_data(
    deployment: str,
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
) -> bool:
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    prepare_runtime_environment(cfg)
    return import_raw_data_zip(cfg, deployment)


def update_metadata(
    deployment: str = "",
    smru_name: str = "",
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
) -> bool:
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    return update_metadata_from_table(cfg, deployment=deployment, smru_name=smru_name)


def process_tags(
    deployment: str = "",
    smru_name: str = "",
    notlc: bool = False,
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
) -> bool:
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    return bool(process_tags_workflow(cfg, deployment=deployment, smru_name=smru_name, notlc=notlc))


def apply_adjustments(
    deployment: str = "",
    smru_name: str = "",
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
):
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    info = load_info_deployment_catalog(cfg, deployment=deployment, smru_name=smru_name)
    if info.invalid_code:
        return None
    return apply_adjustments_processing(cfg, info.selection)


def create_fr0(
    deployment: str = "",
    smru_name: str = "",
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
):
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    info = load_info_deployment_catalog(cfg, deployment=deployment, smru_name=smru_name)
    if info.invalid_code:
        return None
    return create_fr0_processing(cfg, info.selection)


def create_hr2(
    deployment: str = "",
    smru_name: str = "",
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
) -> bool:
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    return create_hr2_workflow(cfg, deployment=deployment, smru_name=smru_name)


def generate_diagnostics(
    deployment: str = "",
    smru_name: str = "",
    qf: str = "lr1",
    adjusted: bool = True,
    parts: list[str] | tuple[str, ...] | None = None,
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
    use_cached_summaries: bool = True,
):
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    selection = Selection(deployment=deployment, smru_name=smru_name).normalized()
    return generate_diagnostics_plotting(
        cfg,
        selection,
        qf=qf,
        adjusted=adjusted,
        parts=parts,
        use_cached_summaries=use_cached_summaries,
    )


def update_metadata_summaries(
    config: MeopConfig | None = None,
    *,
    output_dir: str | os.PathLike[str] | None = None,
    impacted_deployments: list[str] | tuple[str, ...] | None = None,
    force: bool = False,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
) -> dict[str, Any]:
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    return update_metadata_summaries_metadata(
        cfg,
        output_dir=output_dir,
        processed_deployments=impacted_deployments,
        force=force,
    ).as_dict()


def compare_reference_outputs(
    deployment: str = "",
    smru_name: str = "",
    qf: str = "lr0",
    reference_path: str | os.PathLike[str] | None = None,
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
) -> dict[str, Any]:
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    return compare_reference_outputs_workflow(
        cfg,
        deployment=deployment,
        smru_name=smru_name,
        qf=qf,
        reference_path=reference_path,
    )


def run_all_deployments(
    config: MeopConfig | None = None,
    *,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
    notlc: bool = False,
    diagnostics: bool = False,
    diagnostics_qf: str = "lr1",
    diagnostics_raw: bool = False,
    diagnostics_parts: list[str] | tuple[str, ...] | None = None,
    notify_email: list[str] | tuple[str, ...] | None = None,
    notify_when: str | None = None,
    notify_attach: list[str] | tuple[str, ...] | None = None,
    notifications_enabled: bool | None = None,
    force: bool = False,
    force_failed: bool = False,
    include_disabled: bool = False,
    deployments: list[str] | tuple[str, ...] | None = None,
    state_dir: str | os.PathLike[str] | None = None,
    jobs: int = 1,
    verbose: bool = False,
) -> dict[str, Any]:
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    return run_all_deployments_batch(
        config=cfg,
        notlc=notlc,
        diagnostics=diagnostics,
        diagnostics_qf=diagnostics_qf,
        diagnostics_raw=diagnostics_raw,
        diagnostics_parts=diagnostics_parts,
        notify_email=notify_email,
        notify_when=notify_when,
        notify_attach=notify_attach,
        notifications_enabled=notifications_enabled,
        force=force,
        force_failed=force_failed,
        include_disabled=include_disabled,
        deployments=deployments,
        state_dir=state_dir,
        jobs=jobs,
        verbose=verbose,
    ).as_dict()


def publish_data(
    config: MeopConfig | None = None,
    *,
    deployment: str = "",
    smru_name: str = "",
    output_dir: str | os.PathLike[str] | None = None,
    create_files: bool = True,
    update_attrs: bool = True,
    list_profiles: bool = True,
    list_tags: bool = True,
    rebuild: bool = False,
    verbose: bool = False,
    processdir: str | os.PathLike[str] | None = None,
    config_file: str | os.PathLike[str] | None = None,
    machine: str | None = None,
) -> dict[str, Any]:
    """Publish MEOP-CTD processed products to the public release directory.

    Creates ``*_all_prof.nc`` files combining the best available processed
    product with standard-level interpolated variables, patches global
    attributes, and generates ``list_profiles.csv``, ``list_tags.csv``,
    and ``list_deployments.csv``.

    Parameters
    ----------
    deployment:
        Restrict to this deployment code.
    smru_name:
        Restrict to this single tag.
    output_dir:
        Output directory.  Defaults to ``config.publicdir_ctd``.
    create_files:
        Whether to create ``*_all_prof.nc`` output files.
    update_attrs:
        Whether to patch global attributes on published NC files.
    list_profiles:
        Whether to write ``list_profiles.csv``.
    list_tags:
        Whether to write ``list_tags.csv`` and ``list_deployments.csv``.
    rebuild:
        Force recreation of existing output files.
    verbose:
        Print progress to stdout.
    """
    cfg = _resolve_config(config, processdir=processdir, config_file=config_file, machine=machine)
    result = publish_workflow(
        cfg,
        deployment=deployment,
        smru_name=smru_name,
        output_dir=Path(output_dir) if output_dir else None,
        create_files=create_files,
        update_attrs=update_attrs,
        list_profiles=list_profiles,
        list_tags=list_tags,
        rebuild=rebuild,
        verbose=verbose,
    )
    return result.as_dict()
