from __future__ import annotations

import json
import os
import platform
import socket
from pathlib import Path
from typing import Any

from ..models import (
    DEFAULT_VERSION,
    BatchDefaults,
    DiagnosticsDefaults,
    EmailNotificationSettings,
    EmailTransportSettings,
    MeopConfig,
)
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


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        existing = merged.get(key)
        if isinstance(existing, dict) and isinstance(value, dict):
            merged[key] = _deep_merge(existing, value)
        else:
            merged[key] = value
    return merged


def _as_bool(value: object, *, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "no", "n", "off"}:
        return False
    return default


def _as_str_tuple(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, (list, tuple)):
        return tuple(str(item).strip() for item in value if str(item).strip())
    text = str(value).strip()
    return (text,) if text else ()


def _parse_diagnostics_defaults(entry: dict[str, Any]) -> DiagnosticsDefaults:
    diagnostics = entry.get("diagnostics", {}) if isinstance(entry.get("diagnostics", {}), dict) else {}
    parts = _as_str_tuple(diagnostics.get("parts")) or DiagnosticsDefaults().parts
    return DiagnosticsDefaults(
        qf=str(diagnostics.get("qf", DiagnosticsDefaults().qf) or DiagnosticsDefaults().qf),
        adjusted=_as_bool(diagnostics.get("adjusted"), default=DiagnosticsDefaults().adjusted),
        parts=parts,
    )


def _parse_batch_defaults(entry: dict[str, Any]) -> BatchDefaults:
    batch = entry.get("batch", {}) if isinstance(entry.get("batch", {}), dict) else {}
    jobs_raw = batch.get("jobs", BatchDefaults().jobs)
    try:
        jobs = int(jobs_raw)
    except (TypeError, ValueError):
        jobs = BatchDefaults().jobs
    if jobs < 1:
        jobs = BatchDefaults().jobs
    return BatchDefaults(
        jobs=jobs,
        verbose=_as_bool(batch.get("verbose"), default=BatchDefaults().verbose),
    )


def _parse_email_notifications(entry: dict[str, Any]) -> EmailNotificationSettings:
    notifications = entry.get("notifications", {}) if isinstance(entry.get("notifications", {}), dict) else {}
    email = notifications.get("email", {}) if isinstance(notifications.get("email", {}), dict) else {}
    smtp = email.get("smtp", {}) if isinstance(email.get("smtp", {}), dict) else {}
    transport = EmailTransportSettings(
        transport=str(email.get("transport", smtp.get("transport", EmailTransportSettings().transport)) or EmailTransportSettings().transport),
        host=str(smtp.get("host", EmailTransportSettings().host) or EmailTransportSettings().host),
        port=int(smtp.get("port", EmailTransportSettings().port) or EmailTransportSettings().port),
        starttls=_as_bool(smtp.get("starttls"), default=EmailTransportSettings().starttls),
        use_ssl=_as_bool(smtp.get("use_ssl"), default=EmailTransportSettings().use_ssl),
        username_env=str(smtp.get("username_env", EmailTransportSettings().username_env) or EmailTransportSettings().username_env),
        password_env=str(smtp.get("password_env", EmailTransportSettings().password_env) or EmailTransportSettings().password_env),
        from_address=str(smtp.get("from") or smtp.get("from_address") or EmailTransportSettings().from_address),
        sendmail_path=str(smtp.get("sendmail_path", EmailTransportSettings().sendmail_path) or EmailTransportSettings().sendmail_path),
    )
    attach = _as_str_tuple(email.get("attach")) or EmailNotificationSettings().attach
    return EmailNotificationSettings(
        enabled=_as_bool(email.get("enabled"), default=EmailNotificationSettings().enabled),
        when=str(email.get("when", EmailNotificationSettings().when) or EmailNotificationSettings().when),
        to=_as_str_tuple(email.get("to")),
        attach=attach,
        subject_prefix=str(email.get("subject_prefix", EmailNotificationSettings().subject_prefix) or EmailNotificationSettings().subject_prefix),
        transport=transport,
    )


def _default_paths(processdir: Path, version: str, machine: str, config_path: Path | None) -> dict[str, Any]:
    return {
        "processdir": processdir,
        "datadir": processdir / "data",
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
    file_defaults = normalize_config_entry(payload.get("defaults"))
    selected_entry = normalize_config_entry(payload.get("configs", {}).get(chosen_machine))
    merged_entry = _deep_merge(file_defaults, selected_entry)

    if merged_entry.get("processdir"):
        defaults["processdir"] = Path(merged_entry["processdir"]).expanduser()
    processdir_path = Path(defaults["processdir"]).expanduser()
    defaults = _default_paths(processdir_path, version, chosen_machine, config_path)

    path_like_keys = {"processdir", "datadir", "public"}
    resolved = {
        **defaults,
        **{
            key: Path(value).expanduser() if key in path_like_keys else value
            for key, value in merged_entry.items()
        },
    }

    return MeopConfig(
        processdir=Path(resolved["processdir"]),
        datadir=Path(resolved["datadir"]),
        publicdir=Path(resolved["public"]),
        pdflatex=str(resolved.get("pdflatex", "pdflatex")),
        version=str(resolved.get("version", version)),
        machine=str(resolved.get("machine", chosen_machine)),
        config_path=config_path,
        diagnostics_defaults=_parse_diagnostics_defaults(merged_entry),
        batch_defaults=_parse_batch_defaults(merged_entry),
        email_notifications=_parse_email_notifications(merged_entry),
    )
