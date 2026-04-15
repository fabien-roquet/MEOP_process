from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..models import MeopConfig, Selection


@dataclass(frozen=True)
class LocationAdjustmentResult:
    written_files: tuple[Path, ...] = ()
    placeholder: bool = True

    def as_dict(self) -> dict[str, object]:
        return {
            "written_files": [str(path) for path in self.written_files],
            "placeholder": self.placeholder,
        }



def apply_location_adjustment_placeholder(config: MeopConfig, selection: Selection) -> LocationAdjustmentResult:
    """Explicit no-op for the historical ``sc_adjust_locations`` step.

    The Python process boundary now makes the omission deliberate and visible: importing,
    catalog resolution, lr0, and hr0 are now pure Python, while location adjustment
    remains a future compatibility task.
    """

    _ = config
    _ = selection
    return LocationAdjustmentResult()
