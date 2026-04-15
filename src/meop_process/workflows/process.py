from __future__ import annotations

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


def process_tags(
    config: MeopConfig,
    *,
    deployment: str = "",
    smru_name: str = "",
    notlc: bool = False,
) -> bool:
    """Run the full pure-Python MEOP processing chain for one deployment or tag."""

    ensure_runtime_directories(config)
    sync_external_config_files(config)
    info = load_info_deployment(config, deployment=deployment, smru_name=smru_name)
    if info.invalid_code:
        print(f"{info.EXP} is not a valid deployment code. Update data/catalog/list_deployment.csv or data/config_files JSON metadata.")
        return False

    if not import_raw_data_zip(config, info.EXP):
        return False

    remove_deployment_outputs(config, info)
    selection = info.selection

    lr0 = create_ncargo_python(config, selection)
    if not lr0.written_files:
        return False

    apply_location_adjustment_placeholder(config, selection)

    hr0 = create_hr0_python(config, selection)
    if not hr0.written_files:
        return False

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
            return False
        apply_notlc_fr(config, selection)
        hr2 = create_hr2_python(config, selection)
        return bool(hr2.written_files)

    hr1 = apply_tlc(config, selection)
    if not hr1.written_files:
        return False
    apply_tlc_fr(config, selection)
    hr2 = create_hr2_python(config, selection)
    return bool(hr2.written_files)
