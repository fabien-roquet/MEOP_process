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
    return cfg.processdir / "final_dataset_prof" / deployment_code / f"{smru_name}_{qf}_prof.nc"


def fname_traj(smru_name: str, deployment: str = "", *, config: MeopConfig | None = None) -> Path:
    cfg = config or load_config()
    deployment_code = deployment or deployment_from_smru_name(smru_name)
    return cfg.processdir / "final_dataset_prof" / deployment_code / f"{smru_name}_traj.nc"


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
    root = Path(folder) if folder is not None else cfg.processdir
    directory = root / "final_dataset_prof" / deployment_code
    prefix = smru_name if smru_name else f"{deployment_code}-*"
    return sorted(set(directory.glob(f"{prefix}_{qf}_prof.nc")))


def list_smru_name(smru_name: str = "", deployment: str = "", qf: str = "*", *, config: MeopConfig | None = None) -> list[str]:
    return sorted({path.name.split("_")[0] for path in list_fname_prof(smru_name, deployment, qf, config=config)})


def fname_plots(smru_name: str, deployment: str = "", qf: str = "lr0", suffix: str = "_plot", *, config: MeopConfig | None = None) -> Path:
    cfg = config or load_config()
    deployment_code = deployment or deployment_from_smru_name(smru_name)
    return cfg.processdir / "plots" / deployment_code / f"{smru_name}_{qf}_{suffix}.png"


def list_fname_plots(smru_name: str = "", deployment: str = "", qf: str = "*", suffix: str = "_plot", *, config: MeopConfig | None = None) -> list[Path]:
    cfg = config or load_config()
    deployment_code = deployment or deployment_from_smru_name(smru_name)
    directory = cfg.processdir / "plots" / deployment_code
    prefix = smru_name if smru_name else f"{deployment_code}-*"
    return sorted(directory.glob(f"{prefix}_{qf}_{suffix}.png"))


def copy_file(file_name: str, src_dir: str | Path, dst_dir: str | Path) -> Path:
    destination = Path(dst_dir) / file_name
    shutil.copyfile(Path(src_dir) / file_name, destination)
    return destination
