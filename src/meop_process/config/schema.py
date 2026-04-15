from __future__ import annotations

from typing import Any

REQUIRED_CONFIG_KEYS = ("datadir", "refdir", "public")


def normalize_config_entry(entry: dict[str, Any] | None) -> dict[str, Any]:
    """Return a shallow copy of a machine config entry."""

    return dict(entry or {})
