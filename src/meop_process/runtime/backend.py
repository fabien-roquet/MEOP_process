from __future__ import annotations

from typing import Any, Protocol

from ..models import MeopConfig, Selection


class Backend(Protocol):
    name: str

    def start(self, config: MeopConfig) -> None:
        ...

    def stop(self) -> None:
        ...

    def raw_command(self, command: str, verbose: bool = True) -> bool:
        ...

    def print_expression(self, expression: str) -> str:
        ...

    def init_mirounga(self, config: MeopConfig) -> Any:
        ...

    def load_info_deployment(self, config: MeopConfig, selection: Selection) -> None:
        ...

    def info_deployment_invalid(self) -> bool:
        ...

    def run_task(self, task_name: str, config: MeopConfig, selection: Selection) -> bool:
        ...
