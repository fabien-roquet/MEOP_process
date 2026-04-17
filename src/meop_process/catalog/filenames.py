from __future__ import annotations

import shutil
from pathlib import Path

from ..config.loader import load_config
from ..models import MeopConfig


def deployment_from_smru_name(smru_name: str) -> str:
    return smru_name.split("-")[0]


def smru_name_from_fname_prof(fname_prof: str | Path) -> str:
    return Path(fname_prof).name.split("_")[0]


def fname_prof(smru_name: str, deployment: str = "", qf: str = "lr0", *, config: MeopConfig | None = None) -> Path:
    cfg = config or load_config()
    deployment_code = deployment or deployment_from_smru_name(smru_name)
    return cfg.final_dataset_dir / deployment_code / f"{smru_name}_{qf}_prof.nc"


def fname_traj(smru_name: str, deployment: str = "", *, config: MeopConfig | None = None) -> Path:
    cfg = config or load_config()
    deployment_code = deployment or deployment_from_smru_name(smru_name)
    return cfg.trajectory_dataset_dir / deployment_code / f"{smru_name}_traj.nc"


def _processed_roots(cfg: MeopConfig, folder: str | Path | None = None) -> tuple[Path, ...]:
    if folder is None:
        return cfg.final_dataset_search_dirs
    root = Path(folder)
    candidates = []
    new_root = root / "data" / "data_prof"
    legacy_root = root / "final_dataset_prof"
    for candidate in (new_root, legacy_root):
        if candidate not in candidates:
            candidates.append(candidate)
    return tuple(candidates)


def list_fname_prof(
    smru_name: str = "",
    deployment: str = "",
    qf: str = "*",
    *,
    config: MeopConfig | None = None,
    folder: str | Path | None = None,
) -> list[Path]:
    cfg = config or load_config()
    deployment_code = deployment or deployment_from_smru_name(smru_name)
    prefix = smru_name if smru_name else f"{deployment_code}-*"
    matches: set[Path] = set()
    for root in _processed_roots(cfg, folder):
        directory = root / deployment_code
        matches.update(directory.glob(f"{prefix}_{qf}_prof.nc"))
    return sorted(matches)


def list_smru_name(smru_name: str = "", deployment: str = "", qf: str = "*", *, config: MeopConfig | None = None) -> list[str]:
    return sorted({path.name.split("_")[0] for path in list_fname_prof(smru_name, deployment, qf, config=config)})


def fname_plots(smru_name: str, deployment: str = "", qf: str = "lr0", suffix: str = "_plot", *, config: MeopConfig | None = None) -> Path:
    cfg = config or load_config()
    deployment_code = deployment or deployment_from_smru_name(smru_name)
    return cfg.plotdir / deployment_code / f"{smru_name}_{qf}_{suffix}.png"


def list_fname_plots(smru_name: str = "", deployment: str = "", qf: str = "*", suffix: str = "_plot", *, config: MeopConfig | None = None) -> list[Path]:
    cfg = config or load_config()
    deployment_code = deployment or deployment_from_smru_name(smru_name)
    prefix = smru_name if smru_name else f"{deployment_code}-*"
    matches: set[Path] = set()
    for root in cfg.plot_search_dirs:
        matches.update((root / deployment_code).glob(f"{prefix}_{qf}_{suffix}.png"))
    return sorted(matches)


def copy_file(file_name: str, src_dir: str | Path, dst_dir: str | Path) -> Path:
    destination = Path(dst_dir) / file_name
    shutil.copyfile(Path(src_dir) / file_name, destination)
    return destination
