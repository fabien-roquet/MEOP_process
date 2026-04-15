from __future__ import annotations

from ..models import MeopConfig, Selection
from ..processing.adjustments import apply_adjustments, apply_notlc
from ..processing.hr import create_hr0_python
from ..processing.ncargo import create_ncargo_python


class PythonBackend:
    """Future pure-Python backend.

    The backend now supports low-resolution creation, hr0 generation, and the no-TLC
    hr1 branch. Remaining tasks still raise ``NotImplementedError`` until their Python
    replacements are ready.
    """

    name = "python"

    def start(self, config: MeopConfig) -> None:
        _ = config

    def stop(self) -> None:
        return

    def raw_command(self, command: str, verbose: bool = True) -> bool:
        _ = verbose
        raise NotImplementedError(f"Raw Matlab command passthrough is not supported by the Python backend: {command}")

    def print_expression(self, expression: str) -> str:
        raise NotImplementedError(f"Matlab expression printing is not supported by the Python backend: {expression}")

    def init_mirounga(self, config: MeopConfig) -> dict[str, str]:
        return {"backend": "python", "processdir": str(config.processdir)}

    def load_info_deployment(self, config: MeopConfig, selection: Selection) -> None:
        _ = config
        _ = selection

    def info_deployment_invalid(self) -> bool:
        return False

    def run_task(self, task_name: str, config: MeopConfig, selection: Selection) -> bool:
        if task_name == "create_ncargo":
            return bool(create_ncargo_python(config, selection).written_files)
        if task_name == "create_hr0":
            return bool(create_hr0_python(config, selection).written_files)
        if task_name == "apply_adjustments":
            apply_adjustments(config, selection)
            return True
        if task_name == "apply_notlc":
            return bool(apply_notlc(config, selection).written_files)
        _ = config
        _ = selection
        raise NotImplementedError(f"Task {task_name!r} has not been ported to the Python backend yet.")
