from __future__ import annotations

from pathlib import Path

from ..models import DeploymentInfo, MeopConfig, Selection


OUTPUT_CLEANUP_ROOTS = (
    lambda cfg, info: cfg.plotdir / info.EXP,
    lambda cfg, info: cfg.plots_by_deployment_dir / info.EXP,
    lambda cfg, info: info.directory,
)

OUTPUT_CLEANUP_SUFFIXES = (".nc", ".txt", ".json", ".png")
PROCESS_PRODUCT_QFS = ("lr0", "hr0", "fr0", "hr1", "fr1", "lr1", "hr2")


def remove_deployment_outputs(config: MeopConfig, info: DeploymentInfo) -> list[Path]:
    """Python replacement for the legacy ``remove_deployment`` task.

    The cleanup only runs when the deployment has an imported raw input in the working store,
    which mirrors the guard used by ``remove_deployment.m``.
    """

    raw_marker = info.raw_working_fcell if info.nomfic.endswith("_fcell.mat") else info.raw_working_text
    if not raw_marker.exists():
        return []

    prefix = info.requested_smru_name or ""
    removed: list[Path] = []
    for root_builder in OUTPUT_CLEANUP_ROOTS:
        root = root_builder(config, info)
        if not root.exists():
            continue
        for suffix in OUTPUT_CLEANUP_SUFFIXES:
            pattern = f"{prefix}*{suffix}" if prefix else f"*{suffix}"
            for path in root.glob(pattern):
                if not path.is_file():
                    continue
                path.unlink()
                removed.append(path)
    return removed


def prune_profile_products(
    config: MeopConfig,
    selection: Selection,
    *,
    keep_products: tuple[str, ...],
    product_qfs: tuple[str, ...] = PROCESS_PRODUCT_QFS,
) -> list[Path]:
    """Remove rebuildable profile products not listed in ``keep_products``.

    The function is intentionally narrow: it only removes ``*_prof.nc`` files under the selected
    deployment's processed-profile directory. Raw inputs, diagnostics, metadata summaries, and
    public outputs are left untouched.
    """

    selected = selection.normalized()
    if not selected.deployment:
        return []
    root = config.final_dataset_dir / selected.deployment
    if not root.exists():
        return []

    keep = {item.lower() for item in keep_products}
    prefix = selected.smru_name or ""
    removed: list[Path] = []
    for qf in product_qfs:
        qf_name = qf.lower()
        if qf_name in keep:
            continue
        pattern = f"{prefix}*_{qf_name}_prof.nc" if prefix else f"*_{qf_name}_prof.nc"
        for path in sorted(root.glob(pattern)):
            if not path.is_file():
                continue
            path.unlink()
            removed.append(path)
    return removed
