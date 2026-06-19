#!/usr/bin/env python3
"""Regenerate per-panel folders for model-comparison manuscript figures."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent


def run(script: str) -> None:
    subprocess.run([sys.executable, str(HERE / script)], check=True)


def main() -> None:
    for script in [
        "rebuild_figure2.py",
        "rebuild_figure3.py",
        "rebuild_figure4.py",
        "rebuild_figure5.py",
    ]:
        run(script)


if __name__ == "__main__":
    main()
