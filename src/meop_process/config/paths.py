from __future__ import annotations

from pathlib import Path

from ..models import MeopConfig


RUNTIME_DIRECTORIES = (
    lambda cfg: cfg.datadir,
    lambda cfg: cfg.tablesdir,
    lambda cfg: cfg.catalogdir,
    lambda cfg: cfg.data_raw_dir,
    lambda cfg: cfg.config_files_dir,
    lambda cfg: cfg.raw_odv_dir,
    lambda cfg: cfg.raw_hr_dir,
    lambda cfg: cfg.crawl_locdir,
    lambda cfg: cfg.cls_locdir,
    lambda cfg: cfg.publicdir,
    lambda cfg: cfg.final_dataset_dir,
    lambda cfg: cfg.trajectory_dataset_dir,
    lambda cfg: cfg.mapsdir,
    lambda cfg: cfg.plotdir,
    lambda cfg: cfg.plots_by_deployment_dir,
    lambda cfg: cfg.plots_overview_dir,
)


def ensure_runtime_directories(config: MeopConfig) -> list[Path]:
    created: list[Path] = []
    for builder in RUNTIME_DIRECTORIES:
        path = builder(config)
        path.mkdir(parents=True, exist_ok=True)
        created.append(path)
    return created
