from __future__ import annotations

import filecmp
import shutil
from datetime import datetime
from pathlib import Path

from ..models import MeopConfig
from .paths import ensure_runtime_directories


CONFIG_BASENAMES = (
    "deployment2",
    "deployment3",
    "platform2",
    "platform3",
    "deployment2_patch",
    "platform2_patch",
)


def sync_external_config_files(
    config: MeopConfig,
    *,
    source_dir: str | Path | None = None,
    timestamp: str | None = None,
) -> list[Path]:
    """Optionally mirror JSON config files into ``data/data_raw/config_files``.

    The pure-Python package does not require any root-level mirrors. When ``source_dir`` is
    omitted, the helper looks for ``data/external_config_files`` and performs a quiet no-op if
    that directory is absent.
    """

    ensure_runtime_directories(config)
    datestamp = timestamp or datetime.now().strftime("%Y%m%d")
    destination_dir = config.config_files_dir
    source_path = Path(source_dir) if source_dir is not None else config.datadir / "external_config_files"
    updated: list[Path] = []

    if not source_path.exists():
        return updated

    for basename in CONFIG_BASENAMES:
        source = source_path / f"{basename}.json"
        destination = destination_dir / f"{basename}.json"
        backup = destination_dir / f"{basename}_{datestamp}.json"
        if not source.exists():
            continue
        if not destination.exists():
            shutil.copy2(source, destination)
            updated.append(destination)
            continue
        if filecmp.cmp(source, destination, shallow=False):
            continue
        shutil.copy2(destination, backup)
        shutil.copy2(source, destination)
        updated.append(destination)

    return updated
