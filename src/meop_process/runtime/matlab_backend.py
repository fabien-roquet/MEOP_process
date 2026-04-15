from __future__ import annotations

import io
from typing import Any

from ..models import MeopConfig, Selection


class MatlabBackend:
    """Thin Matlab-engine backend that preserves the legacy task sequence."""

    name = "matlab"

    def __init__(self) -> None:
        self._engine: Any | None = None

    def start(self, config: MeopConfig) -> None:
        if self._engine is None:
            import matlab.engine  # type: ignore

            self._engine = matlab.engine.start_matlab()
        self.raw_command(f"cd {config.processdir};", verbose=False)

    def stop(self) -> None:
        if self._engine is not None:
            self._engine.quit()
            self._engine = None

    def _require_engine(self) -> Any:
        if self._engine is None:
            raise RuntimeError("Matlab backend has not been started.")
        return self._engine

    def raw_command(self, command: str, verbose: bool = True) -> bool:
        engine = self._require_engine()
        with io.StringIO() as stdout, io.StringIO() as stderr:
            try:
                engine.eval(command, nargout=0, stdout=stdout, stderr=stderr)
            except Exception:
                if verbose:
                    message = stderr.getvalue().strip()
                    if message:
                        print(message)
                return False
            if verbose:
                message = stdout.getvalue().strip()
                if message:
                    print(message)
            return True

    def print_expression(self, expression: str) -> str:
        engine = self._require_engine()
        with io.StringIO() as output:
            engine.eval(expression, nargout=0, stdout=output, stderr=output)
            return output.getvalue()

    def init_mirounga(self, config: MeopConfig) -> Any:
        engine = self._require_engine()
        engine.addpath(str(config.processdir / "matlab"))
        conf = engine.eval("init_config();", nargout=1)
        self.raw_command("conf = init_mirounga;", verbose=False)
        return conf

    def load_info_deployment(self, config: MeopConfig, selection: Selection) -> None:
        engine = self._require_engine()
        self.init_mirounga(config)
        engine.workspace["EXP"] = selection.deployment
        engine.workspace["one_smru_name"] = selection.smru_name
        try:
            engine.eval("info_deployment=load_info_deployment(conf,EXP,one_smru_name);", nargout=0)
        except Exception:
            print("Check file list_deployment.csv for any UNKNOWN country")

    def info_deployment_invalid(self) -> bool:
        engine = self._require_engine()
        return bool(engine.eval("isfield(info_deployment,'invalid_code') && info_deployment.invalid_code", nargout=1))

    def run_task(self, task_name: str, config: MeopConfig, selection: Selection) -> bool:
        _ = config
        _ = selection
        return self.raw_command(f"{task_name}(conf,EXP,one_smru_name);")
