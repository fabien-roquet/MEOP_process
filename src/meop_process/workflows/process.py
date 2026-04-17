from __future__ import annotations

from dataclasses import dataclass

from ..catalog.deployments import load_info_deployment
from ..config.paths import ensure_runtime_directories
from ..config.sync import sync_external_config_files
from ..io.raw_odv import import_raw_data_zip
from ..metadata.patch import update_metadata_from_table
from ..models import MeopConfig
from ..processing.cleanup import remove_deployment_outputs
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

    def __bool__(self) -> bool:
        return self.success


def process_tags(
    config: MeopConfig,
    *,
    deployment: str = "",
    smru_name: str = "",
    notlc: bool = False,
) -> WorkflowResult:
    """Run the full pure-Python MEOP processing chain for one deployment or tag."""

    ensure_runtime_directories(config)
    sync_external_config_files(config)
    info = load_info_deployment(config, deployment=deployment, smru_name=smru_name)
    if info.invalid_code:
        print(f"{info.EXP} is not a valid deployment code. Update data/catalog/list_deployment.csv or data/config_files JSON metadata.")
        return WorkflowResult(False, "invalid deployment code")

    if not import_raw_data_zip(config, info.EXP):
        return WorkflowResult(False, "missing raw ODV input")

    remove_deployment_outputs(config, info)
    selection = info.selection

    lr0 = create_ncargo_python(config, selection)
    if not lr0.written_files:
        return WorkflowResult(False, "no LR0 files written")

    apply_location_adjustment_placeholder(config, selection)

    hr0 = create_hr0_python(config, selection)
    if not hr0.written_files:
        return WorkflowResult(False, "no HR0 files written")

    create_fr0_python(config, selection)

    update_metadata_from_table(
        config,
        deployment=selection.deployment,
        smru_name=selection.smru_name,
    )
    apply_adjustments(config, selection)

    if notlc:
        hr1 = apply_notlc(config, selection)
        if not hr1.written_files:
            return WorkflowResult(False, "no HR1 files written")
        apply_notlc_fr(config, selection)
        hr2 = create_hr2_python(config, selection)
        if not hr2.written_files:
            return WorkflowResult(False, "no HR2 files written")
        return WorkflowResult(True)

    hr1 = apply_tlc(config, selection)
    if not hr1.written_files:
        return WorkflowResult(False, "no HR1 files written")
    apply_tlc_fr(config, selection)
    hr2 = create_hr2_python(config, selection)
    if not hr2.written_files:
        return WorkflowResult(False, "no HR2 files written")
    return WorkflowResult(True)
