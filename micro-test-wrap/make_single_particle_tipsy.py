#!/usr/bin/env python3
"""Compatibility wrapper for the root-level IC generator.

Allows:
  python3 tools/make_single_particle_tipsy.py ...
"""

from __future__ import annotations

import runpy
from pathlib import Path


def main() -> None:
    target = Path(__file__).resolve().parents[1] / "make_single_particle_tipsy.py"
    runpy.run_path(str(target), run_name="__main__")


if __name__ == "__main__":
    main()

