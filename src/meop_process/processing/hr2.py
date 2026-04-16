from __future__ import annotations

from datetime import datetime
from pathlib import Path

from ..catalog.filenames import fname_prof
from ..models import MeopConfig, Selection
from .hr import HrResult, _as_utc, _open_dataset, _selected_tags, _update_update_timestamps
from .netcdf import save_dataset_with_compression


def _candidate_tags(config: MeopConfig, selection: Selection) -> tuple[str, ...]:
    if selection.smru_name:
        return (selection.smru_name,)
    tags: list[str] = []
    for source_qf in ("fr1", "hr1"):
        for tag in _selected_tags(config, selection, source_qf=source_qf):
            if tag not in tags:
                tags.append(tag)
    return tuple(tags)


def _select_source_path(config: MeopConfig, deployment: str, smru_name: str) -> Path | None:
    for source_qf in ("fr1", "hr1"):
        path = fname_prof(smru_name, deployment=deployment, qf=source_qf, config=config)
        if path.exists():
            return path
    return None


def create_hr2_python(
    config: MeopConfig,
    selection: Selection,
    *,
    now: datetime | None = None,
) -> HrResult:
    tags = _candidate_tags(config, selection)
    written: list[Path] = []
    processed: list[str] = []
    timestamp = _as_utc(now)

    for smru_name in tags:
        source_path = _select_source_path(config, selection.deployment, smru_name)
        if source_path is None:
            continue
        dataset = _open_dataset(source_path)
        try:
            result = dataset.copy(deep=True)
            _update_update_timestamps(result, timestamp)
            target = fname_prof(smru_name, deployment=selection.deployment, qf="hr2", config=config)
            target.parent.mkdir(parents=True, exist_ok=True)
            save_dataset_with_compression(result, target)
            written.append(target)
            processed.append(smru_name)
        finally:
            dataset.close()

    return HrResult(
        written_files=tuple(written),
        processed_tags=tuple(processed),
    )



def create_hr2(config: MeopConfig, selection: Selection, *, now: datetime | None = None) -> HrResult:
    return create_hr2_python(config, selection, now=now)
