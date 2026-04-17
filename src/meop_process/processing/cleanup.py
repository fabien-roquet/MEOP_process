from __future__ import annotations

from pathlib import Path

from ..models import DeploymentInfo, MeopConfig


OUTPUT_CLEANUP_ROOTS = (
    lambda cfg, info: cfg.plotdir / info.EXP,
    lambda cfg, info: cfg.legacy_plotdir / info.EXP,
    lambda cfg, info: info.directory,
    lambda cfg, info: cfg.legacy_final_dataset_dir / info.EXP,
)

OUTPUT_CLEANUP_SUFFIXES = (".nc", ".txt", ".json", ".png")


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
