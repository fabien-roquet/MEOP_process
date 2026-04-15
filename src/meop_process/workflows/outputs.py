from __future__ import annotations

from ..catalog.deployments import load_info_deployment
from ..models import MeopConfig
from ..processing.hr2 import create_hr2_python


def create_hr2(config: MeopConfig, *, deployment: str = "", smru_name: str = "") -> bool:
    info = load_info_deployment(config, deployment=deployment, smru_name=smru_name)
    if info.invalid_code:
        return False
    result = create_hr2_python(config, info.selection)
    return bool(result.written_files)


def generate_doc_latex(config: MeopConfig, *, deployment: str = "", smru_name: str = "") -> bool:
    raise NotImplementedError("Legacy LaTeX report generation has been removed from the package CLI.")


def export_odv4(config: MeopConfig, *, deployment: str = "", smru_name: str = "") -> bool:
    raise NotImplementedError("ODV4 export has not yet been reimplemented in the pure-Python package.")
