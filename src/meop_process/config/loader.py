from __future__ import annotations

import json
import os
import platform
import socket
from pathlib import Path
from typing import Any

from ..models import DEFAULT_VERSION, MeopConfig
from .schema import normalize_config_entry


def detect_machine_key() -> str:
    user = os.getenv("USER") or os.getenv("USERNAME") or "unknown"
    hostname = socket.gethostname() or "unknown"
    system = platform.system().lower() or "unknown"
    return f"{user}_{hostname}_{system}".replace(".", "_").replace("-", "_")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _candidate_config_files(explicit: Path | None, processdir: Path) -> list[Path]:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit)
    env_path = os.getenv("MEOP_CONFIG_FILE")
    if env_path:
        candidates.append(Path(env_path).expanduser())
    data_candidate = processdir / "data" / "configs.json"
    if data_candidate not in candidates:
        candidates.append(data_candidate)
    cwd_data_candidate = Path.cwd() / "data" / "configs.json"
    if cwd_data_candidate not in candidates:
        candidates.append(cwd_data_candidate)
    return candidates


def _read_json_file(path: Path | None) -> dict[str, Any]:
    if path is None or not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _default_paths(processdir: Path, version: str, machine: str, config_path: Path | None) -> dict[str, Any]:
    return {
        "processdir": processdir,
        "datadir": processdir / "data",
        "refdir": processdir / "references",
        "public": processdir / "public",
        "pdflatex": "pdflatex",
        "version": version,
        "machine": machine,
        "config_path": config_path,
    }


def load_config(
    *,
    processdir: str | Path | None = None,
    config_file: str | Path | None = None,
    machine: str | None = None,
) -> MeopConfig:
    """Load a MEOP runtime configuration.

    The loader prefers an explicit machine entry from a runtime JSON configuration when available.
    Otherwise it falls back to a repository-local, self-contained layout with all managed
    data living under ``data/``.
    """

    chosen_machine = machine or os.getenv("MEOP_MACHINE") or detect_machine_key()
    chosen_processdir = Path(processdir).expanduser() if processdir is not None else _repo_root()
    explicit_config = Path(config_file).expanduser() if config_file is not None else None

    config_path: Path | None = None
    payload: dict[str, Any] = {}
    for candidate in _candidate_config_files(explicit_config, chosen_processdir):
        if candidate.exists():
            config_path = candidate
            payload = _read_json_file(candidate)
            break

    version = payload.get("version", {}).get("CTDnew", DEFAULT_VERSION)
    defaults = _default_paths(chosen_processdir, version, chosen_machine, config_path)
    selected_entry = normalize_config_entry(payload.get("configs", {}).get(chosen_machine))

    if selected_entry.get("processdir"):
        defaults["processdir"] = Path(selected_entry["processdir"]).expanduser()
    processdir_path = Path(defaults["processdir"]).expanduser()
    defaults = _default_paths(processdir_path, version, chosen_machine, config_path)

    path_like_keys = {"processdir", "datadir", "refdir", "public"}
    resolved = {
        **defaults,
        **{
            key: Path(value).expanduser() if key in path_like_keys else value
            for key, value in selected_entry.items()
        },
    }

    return MeopConfig(
        processdir=Path(resolved["processdir"]),
        datadir=Path(resolved["datadir"]),
        refdir=Path(resolved["refdir"]),
        publicdir=Path(resolved["public"]),
        pdflatex=str(resolved.get("pdflatex", "pdflatex")),
        version=str(resolved.get("version", version)),
        machine=str(resolved.get("machine", chosen_machine)),
        config_path=config_path,
    )
