from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path = [str(SRC)] + [entry for entry in sys.path if Path(entry).resolve() != SCRIPT_DIR]

from meop_process.batch.runner import main


if __name__ == "__main__":
    raise SystemExit(main())
