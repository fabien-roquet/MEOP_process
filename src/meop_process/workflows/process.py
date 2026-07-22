from __future__ import annotations

from dataclasses import dataclass
import zipfile

from ..catalog.filenames import list_fname_prof
from ..catalog.deployments import load_info_deployment
from ..config.paths import ensure_runtime_directories
from ..config.sync import sync_external_config_files
from ..io.raw_odv import discover_raw_odv_files, import_raw_data_zip, load_raw_odv_profiles
from ..metadata.patch import update_metadata_from_table
from ..models import MeopConfig
from ..processing.cleanup import prune_profile_products, remove_deployment_outputs
from ..processing.adjustments import apply_adjustments, apply_notlc, apply_notlc_fr, apply_tlc, apply_tlc_fr
from ..processing.fr0 import create_fr0_python
from ..processing.hr import create_hr0_python
from ..processing.hr2 import create_hr2_python
from ..processing.locations import apply_location_adjustment_placeholder
from ..processing.ncargo import create_ncargo_python


@dataclass(frozen=True)
class WorkflowResult:
    success: bool
    reason: str = ""
    pruned_files: tuple[str, ...] = ()

    def __bool__(self) -> bool:
        return self.success


def _describe_lr0_failure(config: MeopConfig, deployment: str) -> str:
    raw_files = discover_raw_odv_files(config, deployment)
    if raw_files.archive.exists():
        try:
            with zipfile.ZipFile(raw_files.archive) as archive:
                if len(archive.namelist()) == 0:
                    return "raw ODV archive is empty"
        except zipfile.BadZipFile:
            return "raw ODV archive is invalid"
    if raw_files.combined_text is None and raw_files.ctd_text is None and raw_files.fl_text is None:
        return "raw ODV archive contains no extracted ODV text"
    if not load_raw_odv_profiles(raw_files, config=config):
        return "raw ODV inputs produced no parsable profiles"
    return "no LR0 files written"


def _retained_products(config: MeopConfig, *, keep_intermediates: bool) -> tuple[str, ...]:
    products = list(config.processing_defaults.keep_products)
    if keep_intermediates:
        for qf in config.processing_defaults.debug_products:
            if qf not in products:
                products.append(qf)
    return tuple(products)


def _prune_after_success(config: MeopConfig, selection: Selection, *, keep_intermediates: bool) -> tuple[str, ...]:
    retained = _retained_products(config, keep_intermediates=keep_intermediates)
    return tuple(
        str(path)
        for path in prune_profile_products(
            config,
            selection,
            keep_products=retained,
        )
    )


def process_tags(
    config: MeopConfig,
    *,
    deployment: str = "",
    smru_name: str = "",
    notlc: bool = False,
    keep_intermediate_products: bool | None = None,
) -> WorkflowResult:
    """Run the full pure-Python MEOP processing chain for one deployment or tag."""

    ensure_runtime_directories(config)
    sync_external_config_files(config)
    info = load_info_deployment(config, deployment=deployment, smru_name=smru_name)
    if info.invalid_code:
        print(f"{info.EXP} is not a valid deployment code. Update data/catalog/list_deployment.csv or data/data_raw/config_files JSON metadata.")
        return WorkflowResult(False, "invalid deployment code")

    if not import_raw_data_zip(config, info.EXP):
        return WorkflowResult(False, _describe_lr0_failure(config, info.EXP))

    remove_deployment_outputs(config, info)
    selection = info.selection

    lr0 = create_ncargo_python(config, selection)
    if not lr0.written_files:
        return WorkflowResult(False, _describe_lr0_failure(config, info.EXP))

    apply_location_adjustment_placeholder(config, selection)

    hr0 = create_hr0_python(config, selection)
    if not hr0.written_files:
        lr0_inputs = list_fname_prof(deployment=selection.deployment, qf="lr0", config=config)
        if not lr0_inputs:
            return WorkflowResult(False, "no LR0 inputs available for HR0")
        return WorkflowResult(False, f"no HR0 files written from {len(lr0_inputs)} LR0 inputs")

    create_fr0_python(config, selection)

    update_metadata_from_table(
        config,
        deployment=selection.deployment,
        smru_name=selection.smru_name,
    )
    apply_adjustments(config, selection)

    keep_intermediates = (
        config.processing_defaults.keep_intermediate
        if keep_intermediate_products is None
        else keep_intermediate_products
    )

    if notlc:
        hr1 = apply_notlc(config, selection)
        if not hr1.written_files:
            return WorkflowResult(False, "no HR1 files written")
        apply_notlc_fr(config, selection)
        hr2 = create_hr2_python(config, selection)
        if not hr2.written_files:
            return WorkflowResult(False, "no HR2 files written")
        pruned = _prune_after_success(config, selection, keep_intermediates=keep_intermediates)
        return WorkflowResult(True, pruned_files=pruned)

    hr1 = apply_tlc(config, selection)
    if not hr1.written_files:
        return WorkflowResult(False, "no HR1 files written")
    apply_tlc_fr(config, selection)
    hr2 = create_hr2_python(config, selection)
    if not hr2.written_files:
        return WorkflowResult(False, "no HR2 files written")
    pruned = _prune_after_success(config, selection, keep_intermediates=keep_intermediates)
    return WorkflowResult(True, pruned_files=pruned)
